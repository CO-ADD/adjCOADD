from django import forms
from applib.django.base.forms import Base_ModelForm
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from django.shortcuts import get_object_or_404
from django.forms.widgets import HiddenInput
from django.contrib.postgres.forms import SimpleArrayField
from django_filters import DateRangeFilter, CharFilter, ModelChoiceFilter, ChoiceFilter, MultipleChoiceFilter, IsoDateTimeFilter, DateFromToRangeFilter, DateFilter
from dorganism.models import  Taxonomy

from apputil.models import Dictionary, ApplicationUser, Document
from applib.django.base.filters import BaseStatus_Filter
from adjcoadd.constants import ORGANISM_CLASSES, CELL_CLASSES 
from dpeptide.models import Peptide, Peptide_Isoform, Isoform_Batch, Modification_Dictionary, Isoform_Modification, Monomer


class SequenceTextArea(forms.Textarea):
    """Multi-line, monospace textarea for amino-acid sequence notations (FASTA/BILN/HELM)."""
    def __init__(self, *args, **kwargs):
        attrs = kwargs.get('attrs', {}) or {}
        attrs.setdefault('rows', '4')
        classes = attrs.get('class', '')
        attrs['class'] = f"{classes} sequence-input font-monospace".strip()
        kwargs['attrs'] = attrs
        super().__init__(*args, **kwargs)


#=================================================================================================
# Peptide
#=================================================================================================
class Peptide_Filter(BaseStatus_Filter):
    
    peptide_id = CharFilter(field_name='peptide_id', lookup_expr='icontains')
    peptide_name = CharFilter(field_name='peptide_name', lookup_expr='icontains')
    peptide_notes = CharFilter(field_name='peptide_notes', lookup_expr='icontains')
    peptide_class = MultipleChoiceFilter(field_name='peptide_class', method='multichoices_filter', 
                                             widget=forms.CheckboxSelectMultiple(attrs={'class': 'multiselect-accord'}), choices=[])
    peptide_panel = MultipleChoiceFilter(field_name='peptide_panel', method='multichoices_filter', 
                                             widget=forms.CheckboxSelectMultiple(attrs={'class': 'multiselect-accord'}), choices=[])

    mta_status = ModelChoiceFilter(field_name='mta_status', queryset=Dictionary.objects.filter(dict_class=Peptide.DICTIONARY_FIELDS['mta_status'], astatus__gte=0))
   
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["peptide_class"].extra["choices"]=Dictionary.get_aschoices(Peptide.DICTIONARY_FIELDS['peptide_class'], showDesc = False)
        self.filters["peptide_panel"].extra["choices"]=Dictionary.get_aschoices(Peptide.DICTIONARY_FIELDS['peptide_panel'], showDesc = False)

        # Set Filter label to the Fields VerboseName or Filter Name
        for i in self.filters:
            try:
                self.filters[i].label=self.Meta.model._meta.get_field(self.filters[i].field_name).verbose_name
            except:
                self.filters[i].label=i
    class Meta:
        model=Peptide
        fields=[ 'peptide_id', 'peptide_name','peptide_notes', 'peptide_class', 'peptide_panel', 'mta_status' ]

# -----------------------------------------------------------------
class Peptide_CreateForm(Base_ModelForm):

    peptide_name= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}),required=False,)
    peptide_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}),required=False,)
    fasta = forms.CharField(widget=SequenceTextArea(), required=False)
    biln = forms.CharField(widget=SequenceTextArea(), required=False)
    helm = forms.CharField(widget=SequenceTextArea(), required=False)
    # prep_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    mta_status = forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}),required=False, queryset=Dictionary.objects.all())
    #organism_name=forms.ModelChoiceField(queryset=Taxonomy.objects.all(), widget=forms.HiddenInput(),required=False,)
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=True,)
    #collect_date = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}), required=False)
   
    def __init__(self, *args, **kwargs): 
        super(Peptide_CreateForm, self).__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name
        self.fields['peptide_class'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Peptide.DICTIONARY_FIELDS['peptide_class'])]

        self.fields['peptide_panel'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Peptide.DICTIONARY_FIELDS['peptide_panel'], showDesc=False),)
        self.fields['peptide_panel'].widget.attrs.update({'class': 'form-control', 'size':'5', 'multiple': 'true'})

        self.fields['mta_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Peptide.DICTIONARY_FIELDS['mta_status'])]
        self.create_field_groups()

        for field in self.fields.values():
            if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
                # Add the 'group-input' class to the widget attrs
                attrs = field.widget.attrs
                attrs['class'] = attrs.get('class', '') + 'input-group'
                field.widget.attrs = attrs
    
    # def clean_organism_name(self):       
    #     print(self.organism_name)
    #     data=get_object_or_404(Taxonomy, organism_name=self.organism_name)
    
    #     if data:
    #         if str(data.org_class) in CELL_CLASSES:
    #             return data
    #         else:
    #             self.add_error('org_name', 'Create failed with invalid Organism class')
                
    #     else:
    #         self.add_error('org_name', "Found No Organism")

    #     return data            


    def create_field_groups(self):
        if len(Peptide.DETAIL_VIEW_GROUPS) > 0:
            self.groups = []
            for grp in Peptide.DETAIL_VIEW_GROUPS:
                self.groups.append([self[name] for name in grp if name in self.fields])   
    
    class Meta:
        model=Peptide
        exclude=['peptide_id',  'assoc_documents'] 

# -----------------------------------------------------------------
class Peptide_UpdateForm(Peptide_CreateForm):
    peptide_id = forms.CharField(disabled=True)

    class Meta:
        model=Peptide
        exclude=['assoc_documents']
        field_order = ['peptide_id', 'peptide_name', 'peptide_notes', 'peptide_class', 'peptide_panel',
                       'fasta', 'biln', 'helm', 'mw', 'nseq', 'peptide_origin', 'reference',
                       'mta_status', 'mta_notes', 'mta_document', 'biologist']


#=================================================================================================
# Peptide Isoform
#=================================================================================================
class Isoform_Filter(BaseStatus_Filter):
    isoform_id = CharFilter(field_name='isoform_id', lookup_expr='icontains', label='Isoform ID')
    isoform_name = CharFilter(field_name='isoform_name', lookup_expr='icontains', label='Isoform Name')
    peptide = CharFilter(field_name='peptide_id__peptide_id', lookup_expr='icontains', label='Peptide ID')
    peptide_name = CharFilter(field_name='peptide_id__peptide_name', lookup_expr='icontains', label='Peptide Name')
    fasta = CharFilter(field_name='fasta', lookup_expr='icontains', label='Sequence')

    class Meta:
        model = Peptide_Isoform
        fields = ['isoform_id', 'isoform_name', 'peptide', 'peptide_name', 'fasta']


# -----------------------------------------------------------------
class Isoform_CreateForm(Base_ModelForm):
    isoform_name = forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}), required=False)
    isoform_notes = forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False)
    fasta = forms.CharField(widget=SequenceTextArea(), required=False)
    biln = forms.CharField(widget=SequenceTextArea(), required=False)
    helm = forms.CharField(widget=SequenceTextArea(), required=False)
    biologist = forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=False)

    def __init__(self, *args, **kwargs):
        super(Isoform_CreateForm, self).__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name

    class Meta:
        model = Peptide_Isoform
        exclude = ['isoform_id', 'peptide_modifications', 'assoc_documents']
        field_order = ['isoform_name', 'isoform_notes', 'peptide_id', 'fasta', 'biln', 'helm', 'mw', 'nseq', 'biologist']


# -----------------------------------------------------------------
class Isoform_UpdateForm(Isoform_CreateForm):
    isoform_id = forms.CharField(disabled=True)

    class Meta:
        model = Peptide_Isoform
        exclude = ['peptide_modifications', 'assoc_documents']
        field_order = ['isoform_id', 'isoform_name', 'isoform_notes', 'peptide_id', 'fasta', 'biln', 'helm', 'mw', 'nseq', 'biologist']


#=================================================================================================
# Isoform Batch
#=================================================================================================
class IsoBatch_CreateForm(Base_ModelForm):
    batch_notes = forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}), required=False)
    biologist = forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=False)

    def __init__(self, *args, **kwargs):
        super(IsoBatch_CreateForm, self).__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name
        self.fields['qc_status'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Isoform_Batch.DICTIONARY_FIELDS['qc_status'])]
        self.fields['source_type'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Isoform_Batch.DICTIONARY_FIELDS['source_type'])]
        self.fields['conc_unit'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Isoform_Batch.DICTIONARY_FIELDS['conc_unit'])]
        self.fields['init_amount_unit'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Isoform_Batch.DICTIONARY_FIELDS['init_amount_unit'])]

    class Meta:
        model = Isoform_Batch
        exclude = ['isoformbatch_id', 'isoform_id', 'batch_id']
        field_order = ['batch_notes', 'source', 'source_code', 'source_type', 'qc_status', 'quality_source', 'qc_record', 'stock_date', 'stock_level', 'conc', 'conc_unit', 'init_amount', 'init_amount_unit', 'biologist']


# -----------------------------------------------------------------
class IsoBatch_UpdateForm(IsoBatch_CreateForm):
    isoformbatch_id = forms.CharField(disabled=True)

    class Meta:
        model = Isoform_Batch
        exclude = ['isoform_id', 'batch_id']
        field_order = ['isoformbatch_id', 'batch_notes', 'source', 'source_code', 'source_type', 'qc_status', 'quality_source', 'qc_record', 'stock_date', 'stock_level', 'conc', 'conc_unit', 'init_amount', 'init_amount_unit', 'biologist']


#=================================================================================================
# Modification Dictionary
#=================================================================================================
class Modification_Filter(BaseStatus_Filter):
    pepmod_code = CharFilter(field_name='pepmod_code', lookup_expr='icontains', label='Code')

    class Meta:
        model = Modification_Dictionary
        fields = ['pepmod_id', 'pepmod_code']


class Modification_Form(Base_ModelForm):
    fasta = forms.CharField(widget=SequenceTextArea(), required=False)
    biln = forms.CharField(widget=SequenceTextArea(), required=False)
    helm = forms.CharField(widget=SequenceTextArea(), required=False)

    class Meta:
        model = Modification_Dictionary
        exclude = ['pepmod_id']
        field_order = ['pepmod_id', 'pepmod_code', 'fasta', 'biln', 'helm', 'mw']


class Modification_UpdateForm(Modification_Form):
    pepmod_id = forms.CharField(disabled=True)

    class Meta:
        model = Modification_Dictionary
        exclude = []


#=================================================================================================
# Isoform Modification (ordered link between Isoform and Modification_Dictionary)
#=================================================================================================
class IsoformModification_Form(Base_ModelForm):
    class Meta:
        model = Isoform_Modification
        fields = ['pepmod_id', 'mod_position', 'mod_order']

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['pepmod_id'].queryset = Modification_Dictionary.objects.filter(
            astatus__gte=0).order_by('pepmod_code')
        self.fields['mod_position'].choices = Isoform_Modification.MOD_POSITION_CHOICES


#=================================================================================================
# Monomer
#=================================================================================================
def _monomer_list_to_str(values, length=4):
    """Render an ArrayField list as comma-separated text ('H,OH,,')."""
    values = list(values or [])
    while len(values) < length:
        values.append('')
    return ','.join('' if v is None else str(v) for v in values[:length])


def _monomer_str_to_str_list(value, length=4):
    """Parse comma-separated text into a str/None list."""
    out = [None if v.strip() in ('', 'None') else v.strip()
           for v in (value or '').split(',')]
    while len(out) < length:
        out.append(None)
    return out[:length]


def _monomer_str_to_int_list(value, length=4):
    """Parse comma-separated text into an int/None list."""
    out = []
    for v in (value or '').split(','):
        v = v.strip()
        if v in ('', 'None'):
            out.append(None)
        else:
            out.append(int(v))
    while len(out) < length:
        out.append(None)
    return out[:length]


class Monomer_Filter(BaseStatus_Filter):
    symbol = CharFilter(field_name='symbol', lookup_expr='icontains', label='Symbol')
    monomer_name = CharFilter(field_name='monomer_name', lookup_expr='icontains', label='Name')
    monomer_type = CharFilter(field_name='monomer_type', lookup_expr='icontains', label='Type')
    monomer_subtype = CharFilter(field_name='monomer_subtype', lookup_expr='icontains', label='SubType')

    class Meta:
        model = Monomer
        fields = ['symbol', 'monomer_name', 'monomer_type', 'monomer_subtype']


class Monomer_Form(Base_ModelForm):
    molfile = forms.CharField(widget=SequenceTextArea(), required=False)
    rgroups = forms.CharField(required=False, help_text="Comma-separated R-group names, e.g. H,OH,,")
    rgroup_idx = forms.CharField(required=False, help_text="Comma-separated R-group atom indices, e.g. 5,6,,")
    attachment_idx = forms.CharField(required=False, help_text="Comma-separated attachment atom indices, e.g. 2,3,,")

    class Meta:
        model = Monomer
        exclude = ['monomer_id', 'mol']
        field_order = ['symbol', 'monomer_name', 'monomer_abbr', 'monomer_type',
                       'monomer_subtype', 'nat_analog', 'pdb_name',
                       'rgroups', 'rgroup_idx', 'attachment_idx', 'molfile']

    def __init__(self, *args, **kwargs):
        instance = kwargs.get('instance')
        if instance:
            initial = kwargs.setdefault('initial', {})
            initial.setdefault('rgroups', _monomer_list_to_str(instance.rgroups))
            initial.setdefault('rgroup_idx', _monomer_list_to_str(instance.rgroup_idx))
            initial.setdefault('attachment_idx', _monomer_list_to_str(instance.attachment_idx))
        super().__init__(*args, **kwargs)

    def clean_rgroups(self):
        return _monomer_str_to_str_list(self.cleaned_data.get('rgroups'))

    def clean_rgroup_idx(self):
        return _monomer_str_to_int_list(self.cleaned_data.get('rgroup_idx'))

    def clean_attachment_idx(self):
        return _monomer_str_to_int_list(self.cleaned_data.get('attachment_idx'))


class Monomer_UpdateForm(Monomer_Form):
    monomer_id = forms.CharField(disabled=True)

    class Meta:
        model = Monomer
        exclude = ['mol']
        field_order = ['monomer_id', 'symbol', 'monomer_name', 'monomer_abbr', 'monomer_type',
                       'monomer_subtype', 'nat_analog', 'pdb_name',
                       'rgroups', 'rgroup_idx', 'attachment_idx', 'molfile']

