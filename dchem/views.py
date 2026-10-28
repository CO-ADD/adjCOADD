from django.shortcuts import render, redirect, get_object_or_404
from django.http import HttpResponseNotFound
from django.contrib import messages
from django.contrib.auth.mixins import LoginRequiredMixin
from django.core.exceptions import FieldDoesNotExist
from django.db import transaction, IntegrityError

from rest_framework import viewsets, filters
from django_rdkit.models import MOL_FROM_SMILES
#TANIMOTO_SMILES

from applib.django.restapi import IsAppUser

from applib.django.base.views import Filtered_ListView
from apputil.models import ApplicationLog

from dchem.models import Chem_Structure, Chem_Salt, Chem_Group, Chem_Reaction, AminoAcid
from dchem.serializer import ChemStructure_Serializer
from dchem.forms import (ChemStructure_Filter,
                         ChemSalt_Filter, ChemSalt_Form, ChemSalt_UpdateForm,
                         ChemGroup_Filter, ChemGroup_Form, ChemGroup_UpdateForm,
                         ChemReaction_Filter, ChemReaction_Form, ChemReaction_UpdateForm,
                         AminoAcid_Filter, AminoAcid_Form, AminoAcid_UpdateForm)

# -----------------------------------------------------------------
class Compound_ListAPI(viewsets.ModelViewSet):
# -----------------------------------------------------------------
    queryset = Chem_Structure.objects.all()
    serializer_class = ChemStructure_Serializer
    permission_classes = [IsAppUser]

    def get_queryset(self):
        queryset = Chem_Structure.objects.all()

        # Substructure Filter
        substructure = self.request.query_params.get('substruct')
        if substructure:
            q_mol = MOL_FROM_SMILES(substructure)
            queryset = queryset.filter(molecule__contains=q_mol)

        # Similarity Filter
        similarity_smile = self.request.query_params.get('similar_to')
        threshold = float(self.request.query_params.get('threshold', 0.85))
        if similarity_smile:
            queryset = queryset.annotate(
                similarity=TANIMOTO_SMILES('molecule', similarity_smile)
            ).filter(similarity__gte=threshold).order_by('-similarity')

        return queryset


#=================================================================================================
# Chemical Structure / Salt / Group / Reaction / AminoAcid views
#=================================================================================================
CHEM_VIEWS = {
    'structure': {
        'model': Chem_Structure,
        'filter': ChemStructure_Filter,
        'create_form': None,
        'update_form': None,
        'pk': 'structure_id',
        'title': 'Chemical Structures',
        'singular': 'Chemical Structure',
        'list_name': 'chemstructure_list',
        'detail_name': 'chemstructure_detail',
        'list_template': 'dchem/chem_filter_list.html',
        'filter_first': True,
        'has_detail': True,
    },
    'salt': {
        'model': Chem_Salt,
        'filter': ChemSalt_Filter,
        'create_form': ChemSalt_Form,
        'update_form': ChemSalt_UpdateForm,
        'pk': 'salt_id',
        'title': 'Salts',
        'singular': 'Salt',
        'list_name': 'chemsalt_list',
        'create_name': 'chemsalt_create',
        'update_name': 'chemsalt_update',
        'list_template': 'dchem/chem_pag_action.html',
        'filter_first': False,
        'has_detail': False,
    },
    'group': {
        'model': Chem_Group,
        'filter': ChemGroup_Filter,
        'create_form': ChemGroup_Form,
        'update_form': ChemGroup_UpdateForm,
        'pk': 'chemgroup_id',
        'title': 'Chemical Groups',
        'singular': 'Chemical Group',
        'list_name': 'chemgroup_list',
        'create_name': 'chemgroup_create',
        'update_name': 'chemgroup_update',
        'list_template': 'dchem/chem_pag_action.html',
        'filter_first': False,
        'has_detail': False,
    },
    'reaction': {
        'model': Chem_Reaction,
        'filter': ChemReaction_Filter,
        'create_form': ChemReaction_Form,
        'update_form': ChemReaction_UpdateForm,
        'pk': 'reaction_id',
        'title': 'Chemical Reactions',
        'singular': 'Chemical Reaction',
        'list_name': 'chemreaction_list',
        'create_name': 'chemreaction_create',
        'update_name': 'chemreaction_update',
        'list_template': 'dchem/chem_pag_action.html',
        'filter_first': False,
        'has_detail': False,
    },
    'aminoacid': {
        'model': AminoAcid,
        'filter': AminoAcid_Filter,
        'create_form': AminoAcid_Form,
        'update_form': AminoAcid_UpdateForm,
        'pk': 'aa_id',
        'title': 'Amino Acids',
        'singular': 'Amino Acid',
        'list_name': 'aminoacid_list',
        'detail_name': 'aminoacid_detail',
        'create_name': 'aminoacid_create',
        'update_name': 'aminoacid_update',
        'list_template': 'dchem/chem_list.html',
        'filter_first': False,
        'has_detail': True,
    },
}


def _cfg(key):
    return CHEM_VIEWS[key]


def _display_groups(instance, group_defs):
    groups = []
    for field_names in group_defs:
        rows = []
        for name in field_names:
            try:
                field = instance._meta.get_field(name)
            except FieldDoesNotExist:
                continue
            value = getattr(instance, name, None)
            if value is None or value == '':
                value = ''
            elif isinstance(value, (list, tuple)):
                value = ', '.join(str(v) for v in value)
            else:
                value = str(value)
            rows.append({'label': field.verbose_name, 'value': value})
        groups.append(rows)
    return groups


#-------------------------------------------------------------------------------------------------
class Chem_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    app_name = 'dchem'
    key = None

    def setup(self, request, *args, **kwargs):
        super().setup(request, *args, **kwargs)
        self.cfg = _cfg(self.key)
        self.model = self.cfg['model']
        self.filterset_class = self.cfg['filter']
        self.model_fields = self.model.LIST_VIEW_FIELDS
        self.model_name = self.model.__name__
        self.ordering = ['-acreated_at']
        self.template_name = self.cfg['list_template']

    def get_queryset(self):
        # Filter-first for very large tables (e.g. Chem_Structure):
        # don't evaluate the queryset until a filter/search is supplied.
        if self.cfg.get('filter_first') and not self.get_filter_record():
            self.filterset = self.filterset_class(self.request.GET,
                                                  queryset=self.model.objects.none(),
                                                  filterset_dict={}, deep=False)
            self.filter_count = 0
            return self.model.objects.none()
        return super().get_queryset()

    def get_context_data(self, **kwargs):
        if self.cfg.get('filter_first'):
            # Replicate Filtered_ListView.get_context_data() but skip the
            # full-table COUNT() (very expensive for large tables).
            context = super(Filtered_ListView, self).get_context_data(**kwargs)
            filter_record_dict = self.get_filter_record()
            filter_record = "Selected: " + str(filter_record_dict).replace("{", "").replace("}", "") if str(filter_record_dict).replace("{", "").replace("}", "") else None
            context['filter'] = self.filterset
            context['paginate_by'] = self.request.GET.get("paginate_by", self.paginate_by)
            context['fields'] = self.model.get_fields(fields=self.model_fields)
            context['filterset'] = filter_record
            context['querycount'] = self.filter_count
            context['Count'] = None
            context['has_filters'] = bool(filter_record_dict)
            context['title'] = self.cfg['title']
            context['model_name'] = self.model_name
            context['list_url'] = self.cfg['list_name']
            return context

        context = super().get_context_data(**kwargs)
        context['title'] = self.cfg['title']
        context['model_name'] = self.model_name
        context['list_url'] = self.cfg['list_name']
        context['create_url'] = self.cfg.get('create_name')
        context['update_url'] = self.cfg.get('update_name')
        return context


#-------------------------------------------------------------------------------------------------
def Chem_DetailView(req, key, pk):
    cfg = _cfg(key)
    _object = get_object_or_404(cfg['model'], **{cfg['pk']: pk})
    context = {
        'object': _object,
        'display_groups': _display_groups(_object, cfg['model'].DETAIL_VIEW_GROUPS),
        'title': cfg['title'],
        'list_url': cfg['list_name'],
        'update_url': cfg.get('update_name'),
    }
    return render(req, 'dchem/chem_detail.html', context)


#-------------------------------------------------------------------------------------------------
def Chem_CreateView(req, key):
    cfg = _cfg(key)
    form_class = cfg.get('create_form')
    if form_class is None:
        return HttpResponseNotFound()

    form = form_class(req.POST or None)
    if req.method == 'POST':
        if form.is_valid():
            try:
                with transaction.atomic(using='dchem'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Create', str(instance.pk), 'Info', req.user, str(instance.pk), f'Create {cfg["singular"]}', 'Completed')
                if cfg.get('has_detail'):
                    return redirect(cfg['detail_name'], pk=instance.pk)
                return redirect(cfg['list_name'])
            except IntegrityError as err:
                messages.error(req, f'IntegrityError {err} happens, record may be existed!')
        else:
            messages.warning(req, f'Create failed due to {form.errors} error')
        return redirect(cfg['list_name'])

    context = {'form': form, 'modal_title': f'Create {cfg["singular"]}'}
    return render(req, 'dchem/chem_form_modal.html', context)


#-------------------------------------------------------------------------------------------------
def Chem_UpdateView(req, key, pk):
    cfg = _cfg(key)
    form_class = cfg.get('update_form')
    if form_class is None:
        return HttpResponseNotFound()

    _object = get_object_or_404(cfg['model'], **{cfg['pk']: pk})
    if req.method == 'POST':
        form = form_class(req.POST, instance=_object)
        if form.is_valid():
            try:
                with transaction.atomic(using='dchem'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Update', str(instance.pk), 'Info', req.user, str(instance.pk), f'Update {cfg["singular"]}', 'Completed')
                if cfg.get('has_detail'):
                    return redirect(cfg['detail_name'], pk=instance.pk)
                return redirect(cfg['list_name'])
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
        else:
            messages.warning(req, f'Update failed due to {form.errors} error')
        return redirect(cfg['list_name'])

    form = form_class(instance=_object)
    context = {'object': _object, 'form': form, 'modal_title': f'Update {cfg["singular"]}', 'modal_subtitle': pk}
    return render(req, 'dchem/chem_form_modal.html', context)
