"""
Background task management (django.tasks, Django 6.1+).

- Task_UpdateABaseView : enqueue the ABase compound batch upload task (modal)
- TaskResult_ListView  : paginated, filterable list of django-tasks-db results
"""
from django.contrib import messages
from django.contrib.auth.decorators import login_required, user_passes_test
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.core.exceptions import ValidationError
from django.http import HttpResponse, HttpResponseNotAllowed
from django.shortcuts import render, redirect

from django_filters import ChoiceFilter, CharFilter, DateTimeFromToRangeFilter

from django_tasks_db.models import DBTaskResult

from applib.django.base.filters import Base_Filter
from applib.django.base.views import Filtered_ListView

from dsample.tasks import task_update_ABase_Batches

#-------------------------------------------------------------------------------------------------
# Permission helper - admin rights = Django superuser or application 'Admin' permission
#-------------------------------------------------------------------------------------------------
def is_admin_user(user):
    return user.is_authenticated and (user.is_superuser or user.has_permission('Admin'))


admin_required = user_passes_test(is_admin_user, login_url='permission_not_granted')


class AdminRequiredMixin(LoginRequiredMixin, UserPassesTestMixin):
    login_url = 'permission_not_granted'

    def test_func(self):
        return is_admin_user(self.request.user)


#-------------------------------------------------------------------------------------------------
# Task status choices (values stored by django-tasks-db)
#-------------------------------------------------------------------------------------------------
TASK_STATUS_CHOICES = [
    ('READY', 'Queued'),
    ('RUNNING', 'Running'),
    ('SUCCESSFUL', 'Succeeded'),
    ('FAILED', 'Failed'),
]


#-------------------------------------------------------------------------------------------------
class TaskResult_Filter(Base_Filter):
    status = ChoiceFilter(field_name='status', choices=[], label='Status')
    task = CharFilter(field_name='task_path', lookup_expr='icontains', label='Task Name')
    queued = DateTimeFromToRangeFilter(field_name='enqueued_at', label='Queued')
    started = DateTimeFromToRangeFilter(field_name='started_at', label='Started')
    finished = DateTimeFromToRangeFilter(field_name='finished_at', label='Finished')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['status'].extra['choices'] = TASK_STATUS_CHOICES

    class Meta:
        model = DBTaskResult
        fields = ['status', 'task', 'queued', 'started', 'finished']


#-------------------------------------------------------------------------------------------------
class TaskResult_ListView(AdminRequiredMixin, Filtered_ListView):
    model = DBTaskResult
    template_name = 'apputil/task_list.html'
    filterset_class = TaskResult_Filter
    model_fields = {
        'id': 'ID',
        'task_path': 'Task',
        'status': 'Status',
        'enqueued_at': 'Queued',
        'started_at': 'Started',
        'finished_at': 'Finished',
    }
    model_name = 'TaskResult'
    app_name = 'django_tasks_db'
    ordering = ['-enqueued_at']

    def get_queryset(self):
        # Mirrors Filtered_ListView.get_queryset(), but skips the session
        # queryset caching which cannot JSON-serialize the UUID primary keys.
        queryset = super(Filtered_ListView, self).get_queryset()

        kwargs = {'deep': True} if 'applymulti' in self.request.GET else {'deep': False}
        self.filterset = self.filterset_class(
            self.request.GET,
            queryset=queryset,
            filterset_dict=self.get_filter_record(),
            **kwargs,
        )

        order = self.get_order_by()
        self.filter_count = self.filterset.qs.distinct().count()
        if order:
            order = order.replace(".", "__")
            return self.filterset.qs.distinct().order_by(order)
        return self.filterset.qs.distinct()

    def get_context_data(self, **kwargs):
        context = super(Filtered_ListView, self).get_context_data(**kwargs)

        filter_record_dict = self.get_filter_record()
        filter_record = "Selected: " + str(filter_record_dict).replace("{", "").replace("}", "") if str(filter_record_dict).replace("{", "").replace("}", "") else None

        context['filter'] = self.filterset
        context['paginate_by'] = self.get_paginate_by(self, **kwargs)
        context['fields'] = list(self.model_fields.values())
        context['filterset'] = filter_record
        context['querycount'] = self.filter_count
        return context


#-------------------------------------------------------------------------------------------------
@login_required(login_url='/')
@admin_required
def Task_UpdateABaseView(req):
    """
    GET  : render the 'Update ABase' confirmation modal
    POST : enqueue the upload task and redirect to the task monitor
    """
    if req.method == 'POST':
        test = int(req.POST.get('test') or 0)
        overwrite = req.POST.get('overwrite') == 'on'
        try:
            result = task_update_ABase_Batches.enqueue(
                test=test, upload=True, overwrite=overwrite,
            )
        except Exception as err:
            messages.error(req, f"Failed to enqueue ABase update: {err}")
            return redirect('task_current')
        messages.success(
            req,
            f"ABase update task enqueued ({result.id}). "
            f"test={test}, overwrite={overwrite}",
        )
        return redirect('task_current')

    context = {
        'backend_configured': _task_backend_configured(),
    }
    return render(req, 'modal/task_update_abase_modal.html', context)


#-------------------------------------------------------------------------------------------------
@login_required(login_url='/')
@admin_required
def Task_RemoveView(req, pk):
    """
    Remove a queued (READY) task from the database queue.

    A running (RUNNING) task cannot be cancelled with django-tasks-db, so only
    queued tasks are removed here.
    """
    if req.method != 'POST':
        return HttpResponseNotAllowed(['POST'])

    try:
        task = DBTaskResult.objects.get(id=pk)
    except (DBTaskResult.DoesNotExist, ValidationError, ValueError):
        messages.warning(req, f'Task {pk} does not exist')
        return redirect('task_current')

    if task.status == 'READY':
        task.delete()
        messages.success(req, f'Task {pk} removed from the queue')
    else:
        messages.warning(req, f'Task {pk} can no longer be removed (status: {task.status})')

    if req.headers.get('HX-Request'):
        return HttpResponse(status=200)
    return redirect('task_current')


#-------------------------------------------------------------------------------------------------
def _task_backend_configured():
    """True when a database-backed task backend is in use."""
    from django.conf import settings
    backend = settings.TASKS.get("default", {}).get("BACKEND", "")
    return "DatabaseBackend" in backend
