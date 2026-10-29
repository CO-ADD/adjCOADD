from django import forms
from applib.django.base.forms import Base_ModelForm
from django_filters import CharFilter, ChoiceFilter

from apputil.models import Dictionary
from applib.django.base.filters import BaseStatus_Filter

from dchem.models import Chem_Structure, Chem_Salt, Chem_Group, Chem_Reaction, AminoAcid

#=================================================================================================
# Chem Structure
#=================================================================================================
class ChemStructure_Filter(BaseStatus_Filter):
    structure_id = CharFilter(field_name='structure_id', lookup_expr='icontains', label='Structure ID')
    structure_name = CharFilter(field_name='structure_name', lookup_expr='icontains', label='Structure Name')
    structure_type = ChoiceFilter(field_name='structure_type', choices=[], label='Type')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['structure_type'].extra['choices'] = [(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Chem_Structure.DICTIONARY_FIELDS['structure_type'])]

    class Meta:
        model = Chem_Structure
        fields = ['structure_id', 'structure_name', 'structure_type']


#=================================================================================================
# Chem Salt
#=================================================================================================
class ChemSalt_Filter(BaseStatus_Filter):
    salt_id = CharFilter(field_name='salt_id', lookup_expr='icontains', label='Salt ID')
    salt_name = CharFilter(field_name='salt_name', lookup_expr='icontains', label='Salt Name')
    salt_type = ChoiceFilter(field_name='salt_type', choices=[], label='Type')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['salt_type'].extra['choices'] = [(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Chem_Salt.DICTIONARY_FIELDS['salt_type'])]

    class Meta:
        model = Chem_Salt
        fields = ['salt_id', 'salt_name', 'salt_type']


class ChemSalt_Form(Base_ModelForm):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name
        self.fields['salt_type'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Chem_Salt.DICTIONARY_FIELDS['salt_type'])]

    class Meta:
        model = Chem_Salt
        exclude = ['salt_id']
        field_order = ['salt_name', 'salt_type', 'smiles', 'mf', 'mw', 'natoms', 'charge', 'h_equiv']


class ChemSalt_UpdateForm(ChemSalt_Form):
    salt_id = forms.CharField(disabled=True)

    class Meta:
        model = Chem_Salt
        exclude = []


#=================================================================================================
# Chem Group
#=================================================================================================
class ChemGroup_Filter(BaseStatus_Filter):
    chemgroup_code = CharFilter(field_name='chemgroup_code', lookup_expr='icontains', label='Code')
    chemgroup_name = CharFilter(field_name='chemgroup_name', lookup_expr='icontains', label='Name')
    chemgroup_set = CharFilter(field_name='chemgroup_set', lookup_expr='icontains', label='Set')
    chemgroup_type = ChoiceFilter(field_name='chemgroup_type', choices=[], label='Type')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['chemgroup_type'].extra['choices'] = [(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Chem_Group.DICTIONARY_FIELDS['chemgroup_type'])]

    class Meta:
        model = Chem_Group
        fields = ['chemgroup_code', 'chemgroup_name', 'chemgroup_set', 'chemgroup_type']


class ChemGroup_Form(Base_ModelForm):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name
        self.fields['chemgroup_type'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Chem_Group.DICTIONARY_FIELDS['chemgroup_type'])]

    class Meta:
        model = Chem_Group
        exclude = ['chemgroup_id']
        field_order = ['chemgroup_code', 'chemgroup_name', 'chemgroup_set', 'chemgroup_type', 'smarts', 'allowed_min', 'allowed_max']


class ChemGroup_UpdateForm(ChemGroup_Form):
    chemgroup_id = forms.CharField(disabled=True)

    class Meta:
        model = Chem_Group
        exclude = []


#=================================================================================================
# Chem Reaction
#=================================================================================================
class ChemReaction_Filter(BaseStatus_Filter):
    reaction_code = CharFilter(field_name='reaction_code', lookup_expr='icontains', label='Code')
    reaction_name = CharFilter(field_name='reaction_name', lookup_expr='icontains', label='Name')
    reaction_type = ChoiceFilter(field_name='reaction_type', choices=[], label='Type')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['reaction_type'].extra['choices'] = [(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Chem_Reaction.DICTIONARY_FIELDS['alert_type'])]

    class Meta:
        model = Chem_Reaction
        fields = ['reaction_code', 'reaction_name', 'reaction_type']


class ChemReaction_Form(Base_ModelForm):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name
        self.fields['reaction_type'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Chem_Reaction.DICTIONARY_FIELDS['alert_type'])]

    class Meta:
        model = Chem_Reaction
        exclude = ['reaction_id']
        field_order = ['reaction_code', 'reaction_name', 'reaction_type', 'smarts', 'run_status']


class ChemReaction_UpdateForm(ChemReaction_Form):
    reaction_id = forms.CharField(disabled=True)

    class Meta:
        model = Chem_Reaction
        exclude = []


#=================================================================================================
# Amino Acid
#=================================================================================================
class AminoAcid_Filter(BaseStatus_Filter):
    aa_id = CharFilter(field_name='aa_id', lookup_expr='icontains', label='AA ID')
    aa_code = CharFilter(field_name='aa_code', lookup_expr='icontains', label='Code')
    aa_letter = CharFilter(field_name='aa_letter', lookup_expr='icontains', label='Letter')

    class Meta:
        model = AminoAcid
        fields = ['aa_id', 'aa_code', 'aa_letter']


class AminoAcid_Form(Base_ModelForm):
    class Meta:
        model = AminoAcid
        exclude = ['aa_id']
        field_order = ['aa_code', 'aa_letter', 'smiles', 'mf', 'mw', 'link_cterm', 'link_nterm', 'link_ester', 'link_disulphide']


class AminoAcid_UpdateForm(AminoAcid_Form):
    aa_id = forms.CharField(disabled=True)

    class Meta:
        model = AminoAcid
        exclude = []
