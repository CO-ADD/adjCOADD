"""
Views for background task management (django.tasks, Django 6.1+).

Provides:
- Task_UpdateABaseView : enqueue the ABase compound batch upload task
- Task_CurrentView     : monitor recently enqueued task results
"""
from django.apps import apps
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.decorators import login_required, user_passes_test
from django.db.utils import DatabaseError
from django.shortcuts import render, redirect

from dsample.tasks import task_update_ABase_Batches

#-------------------------------------------------------------------------------------------------
# Permission helper - admin rights = Django superuser or application 'Admin' permission
#-------------------------------------------------------------------------------------------------
def is_admin_user(user):
    return user.is_authenticated and (user.is_superuser or user.has_permission('Admin'))


admin_required = user_passes_test(is_admin_user, login_url='permission_not_granted')


#-------------------------------------------------------------------------------------------------
def _get_task_result_model():
    """Return the django-tasks-db ORM model, or None if not configured."""
    try:
        return apps.get_model("django_tasks_db", "DBTaskResult")
    except LookupError:
        return None


def _task_backend_configured():
    """True when a database-backed task backend is in use."""
    backend = settings.TASKS.get("default", {}).get("BACKEND", "")
    return "DatabaseBackend" in backend


def _recent_task_results(limit=100):
    model = _get_task_result_model()
    if model is None:
        return []
    try:
        results = list(model.objects.order_by("-enqueued_at")[:limit])
    except DatabaseError:
        # Table not migrated yet - treat as no results.
        return []
    for result in results:
        result.short_name = (
            result.task_path.rsplit(".", 1)[-1] if result.task_path else "—"
        )
    return results


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
def Task_CurrentView(req):
    """Full-page task monitor."""
    context = _task_context()
    return render(req, 'apputil/task_current.html', context)


#-------------------------------------------------------------------------------------------------
@login_required(login_url='/')
@admin_required
def Task_CurrentTableView(req):
    """Partial table (for htmx auto-refresh)."""
    return render(req, 'apputil/task_current_table.html', _task_context())


def _task_context():
    results = _recent_task_results()
    return {
        'task_results': results,
        'backend_configured': _task_backend_configured(),
    }
