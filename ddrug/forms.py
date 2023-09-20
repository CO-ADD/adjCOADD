import django_filters
from rdkit import Chem
from django_rdkit.models import *

from django import forms
from django.core.exceptions import ValidationError
from django.contrib.postgres.forms import SimpleArrayField
from django.shortcuts import get_object_or_404

from django.contrib.postgres.search import TrigramSimilarity
from django.db.models.functions import Greatest

from apputil.models import Dictionary, ApplicationUser
from apputil.utils.filters_base import Filterbase
from .models import Drug, VITEK_Card, VITEK_AST, VITEK_ID, MIC_COADD, MIC_Pub, Breakpoint
from adjcoadd.constants import *


#========================================Drug Form================================================================
class Drug_form(forms.ModelForm):
    drug_type = forms.ModelChoiceField(widget=forms.Select(attrs={'class':'form-select'}), required=False, queryset=Dictionary.objects.all())
    max_phase = forms.ChoiceField(widget=forms.Select(attrs={'class':'form-select'}), required=False, choices= [], )
    drug_codes= SimpleArrayField(forms.CharField(), required=False)
    drug_othernames = SimpleArrayField(forms.CharField(), required=False)
    drug_note= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    approval_note=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    drug_id=forms.CharField(widget=forms.HiddenInput(), required=False)
    # smol=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '12'}),)
   


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name
        self.fields['drug_panel'].widget = forms.CheckboxSelectMultiple(choices= [])# Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_panel'], showDesc=False),)
        self.fields['drug_panel'].widget.attrs.update({'class': 'form-select', 'size':'5', 'multiple': 'true'})
        self.fields['drug_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Drug.Choice_Dictionary['drug_type'])]
        self.fields['max_phase'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Drug.Choice_Dictionary['max_phase'])]
        self.create_field_groups()
        for field in self.fields.values():
            if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
                # Add the 'group-input' class to the widget attrs
                attrs = field.widget.attrs
                attrs['class'] = attrs.get('class', '') + 'input-group'
                field.widget.attrs = attrs
    
    def create_field_groups(self):
        self.group1 = [self[name] for name in ("drug_othernames", "drug_codes", "drug_type", 'n_compounds',"drug_panel",'mw','mf',"drug_note", )]
        self.group2 = [self[name] for name in ("drug_class", "drug_subclass", "drug_target", "drug_subtarget", 'moa', 'antimicro', 'antimicro_class','max_phase','approval_note','admin_routes','application',)]
        self.group3 = [self[name] for name in ('chembl', 'drugbank', 'cas', 'pubchem', 'chemspider','unii', 'kegg', 'comptox', 'echa', 'chebi', 'uq_imb', 'vendor', 'vendor_catno')]

    class Meta:
        model =Drug
        fields='__all__'
        exclude=['ffp2', 'torsionbv', 'mfp2', 'smol' ]
       
    

    def clean_smol(self):
        data=self.cleaned_data['smol']

        if data:
            data=Chem.MolFromMolBlock('M- \n'+data)
        else:
            self.add_error('smol', 'Provide smol value, currently is None')
        
        return data

    # def clean_mfp2(self):
    #     data=self.cleaned_data['mfp2']
    #     data=Chem.AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(self.instance.smiles),radius=2, bitInfo={})
    #     print(data)
    #     return data


# -------------fitlerset Forms---------------------------------------------------------------

class Drug_filter(Filterbase):
    Drug_Name = django_filters.CharFilter(field_name='drug_name', lookup_expr='icontains')
    Drug_Type=django_filters.ChoiceFilter(field_name='drug_type',widget=forms.RadioSelect, choices=[], empty_label=None)
    Target=django_filters.ChoiceFilter(field_name='drug_target', choices=[])
    Drug_Class=django_filters.ChoiceFilter(field_name='drug_class', choices=[])
    Antimicro=django_filters.ChoiceFilter(field_name='antimicro', choices=[])
    Other_Name = django_filters.CharFilter(field_name='drug_othernames', lookup_expr='icontains')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["Drug_Type"].extra['choices']=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Drug.Choice_Dictionary['drug_type'])]
        self.filters['Drug_Name'].label='Drug Name'
        self.filters['Drug_Type'].label='Drug Type'
        self.filters['Target'].label='Drug Target'
        self.filters['Drug_Class'].label='Drug Class'
        self.filters['Antimicro'].label='Antimicro'
        self.filters['Target'].extra["choices"] = self.Meta.model.get_field_choices(field_name='drug_target')
        self.filters['Drug_Class'].extra["choices"] = self.Meta.model.get_field_choices(field_name='drug_class')
        self.filters['Antimicro'].extra["choices"] = self.Meta.model.get_field_choices(field_name='antimicro')
    
    class Meta:
        model=Drug
        fields=['Drug_Name', 'Drug_Type', 'Target', 'Drug_Class', 'Antimicro', 'Other_Name']
        # fields=list(model.HEADER_FIELDS.keys())


#=================================================================================================
# Vitek Data
#=================================================================================================

# -----------------------------------------------------------------
# VitekCard
# -----------------------------------------------------------------
class VitekCard_Filter(Filterbase):
    f_OrgID = django_filters.CharFilter(field_name='card_barcode__orgbatch_id__organism_id__organism_id', lookup_expr='icontains',label="Organism ID")
    #card_barcode = django_filters.CharFilter(lookup_expr='icontains')
    card_code = django_filters.ChoiceFilter(field_name='card_code', choices=[], label="Card Code")
    card_type = django_filters.ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=VITEK_Card.Choice_Dictionary['card_type']))
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['card_code'].extra["choices"] = self.Meta.model.get_field_choices(field_name='card_code')

    class Meta:
        model=VITEK_Card
        fields=['f_OrgID']
        fields += list(model.HEADER_FIELDS.keys())
        exclude = ['orgbatch_id.organism_id.organism_id',
                   'orgbatch_id.batch_id',
                   ]


# -----------------------------------------------------------------
# Vitek AST
# -----------------------------------------------------------------
class VitekAST_Filter(Filterbase):
# -----------------------------------------------------------------
    f_OrgID = django_filters.CharFilter(field_name='card_barcode__orgbatch_id__organism_id__organism_id', lookup_expr='icontains',label="Organism ID")
    f_OrgName = django_filters.CharFilter(field_name='card_barcode__orgbatch_id__organism_id__organism_name', lookup_expr='icontains',label='Organism Name')
    f_OrgBatchID = django_filters.CharFilter(field_name='card_barcode__orgbatch_id__batch_id', lookup_expr='icontains',label='Batch')
    f_DrugName = django_filters.CharFilter(field_name='drug_id__drug_name', lookup_expr='icontains',label='Drug Name')

    # bp_profile = django_filters.ChoiceFilter(field_name = 'bp_profile', choices=[], label = 'BP')
    bp_source = django_filters.ChoiceFilter(field_name = 'bp_source', choices=[], label = 'Source')

    codes = django_filters.CharFilter(field_name='drug_id__drug_codes', lookup_expr='icontains',label="Drug Code")
    #codes = django_filters.ChoiceFilter(field_name = 'drug_id__drug_codes', choices=[], label = 'Code')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self.filters['bp_profile'].extra["choices"] = self.Meta.model.get_field_choices(field_name='bp_profile')
        self.filters['bp_source'].extra["choices"] = self.Meta.model.get_field_choices(field_name='bp_source')
        
        # code is Foreighkey field
        #choice_query = Drug.objects.order_by().values_list('drug_codes').distinct()
        #choices = [(str(i[0]), str(i[0])) for i in choice_query]
        #self.filters['codes'].extra["choices"] = choices

    class Meta:
        model=VITEK_AST
        fields=['f_OrgID','f_OrgBatchID','f_OrgName','f_DrugName', 'codes']
        fields +=list(model.HEADER_FIELDS.keys())
        exclude = ['card_barcode.orgbatch_id.organism_id.organism_id',
                   'card_barcode.orgbatch_id.batch_id',
                   'card_barcode.orgbatch_id.organism_id.organism_name',
                   'drug_id.drug_name',
                   'drug_id.drug_codes',
                   ]

# -----------------------------------------------------------------
# Vitek ID
# -----------------------------------------------------------------
class VitekID_Filter(Filterbase):
# -----------------------------------------------------------------
    fOrg_ID = django_filters.CharFilter(field_name='card_barcode__orgbatch_id__organism_id__organism_id', lookup_expr='icontains',label="Organism ID")
    fBatch_ID = django_filters.CharFilter(field_name='card_barcode__orgbatch_id__batch_id', lookup_expr='icontains',label="Batch")
    fOrg_Name = django_filters.CharFilter(field_name='card_barcode__orgbatch_id__organism_id__organism_name', lookup_expr='icontains',label="Organism")
    id_confidence = django_filters.ChoiceFilter(field_name = 'id_confidence', choices=[], label = 'ID Confidence')
    process = django_filters.ChoiceFilter(field_name = 'process', choices=[], label = 'Vitek Process')
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['id_confidence'].extra["choices"] = self.Meta.model.get_field_choices(field_name='id_confidence')
        self.filters['process'].extra["choices"] = self.Meta.model.get_field_choices(field_name='process')

    class Meta:
        model=VITEK_ID
        fields = ['fOrg_ID','fBatch_ID','fOrg_Name']
        fields +=list(model.HEADER_FIELDS.keys())
        exclude = ['card_barcode.orgbatch_id.organism_id.organism_id',
                   'card_barcode.orgbatch_id.batch_id',
                   'card_barcode.orgbatch_id.organism_id.organism_name',
                   ]

#=================================================================================================
# AntiBiogram
#=================================================================================================

# -----------------------------------------------------------------
class MIC_COADDfilter(Filterbase):
    fOrgBatch_ID = django_filters.CharFilter(field_name='orgbatch_id__orgbatch_id', lookup_expr='icontains',label="OrgBatch ID")
    fBatch_ID = django_filters.CharFilter(field_name='orgbatch_id__batch_id', lookup_expr='icontains',label="Batch")
    fOrg_Name = django_filters.CharFilter(field_name='orgbatch_id__organism_id__organism_name', lookup_expr='icontains',label="Organism")
    fDrug_Name = django_filters.CharFilter(field_name="drug_id__drug_name", lookup_expr='icontains',label="Drug Name")
    mic = django_filters.CharFilter(lookup_expr='icontains', label="MIC")
    bp_profile = django_filters.ChoiceFilter(field_name = 'bp_profile', choices=[], label = 'Break Point')

    class Meta:
        model=MIC_COADD
        fields = ['fOrgBatch_ID','fBatch_ID','fOrg_Name','fDrug_Name']
        fields +=list(model.HEADER_FIELDS.keys())
        exclude = ['orgbatch_id.organism_id.organism_id',
                   'orgbatch_id.batch_id',
                   'orgbatch_id.organism_id.organism_name',
                   'drug_id.drug_name',
                   ]
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['bp_profile'].extra["choices"] = self.Meta.model.get_field_choices(field_name='bp_profile')
    
    def filter_all_fields(self, queryset, name, value):
        if value:
            # print(f"filtr fields is {self._meta.model._meta.fields})")
            fields=[f.name for f in self._meta.model._meta.fields]
            value=value         
            similarity = Greatest(
                TrigramSimilarity('drug_id__drug_name', value),
                TrigramSimilarity('mic', value),
                TrigramSimilarity('orgbatch_id__organism_id__organism_name', value),
                TrigramSimilarity('mic_unit', value),
                TrigramSimilarity('mic_type__dict_value', value),
                TrigramSimilarity('bp_profile', value),
                TrigramSimilarity('bp_source', value),
                TrigramSimilarity('run_id', value),
                TrigramSimilarity('testplate_id', value),
                TrigramSimilarity('testwell_id', value),
                TrigramSimilarity('plate_size__dict_value', value),
                TrigramSimilarity('plate_material__dict_value', value),
                TrigramSimilarity('media', value),
                TrigramSimilarity('dye', value),
               )
            queryset=queryset.annotate(similarity=similarity)#(similarity=TrigramSimilarity('drug_id__drug_name', value),)

            return queryset.filter(similarity__gt=0.1)#(q_object)
        return queryset
    
    def filter_all_fields_deep(self, queryset, name, value):
        queryset=self.filter_all_fields(queryset, name, value)
        return queryset

    
    
# -----------------------------------------------------------------
class MIC_Pubfilter(Filterbase):
# -----------------------------------------------------------------
    fOrg_ID = django_filters.CharFilter(field_name='orgbatch_id', lookup_expr='icontains',label="Org ID")
    fOrg_Name = django_filters.CharFilter(field_name='organism_id__organism_name', lookup_expr='icontains',label="Organism")
    fDrug_Name = django_filters.CharFilter(field_name="drug_id__drug_name", lookup_expr='icontains',label="Drug Name")
    mic_type = django_filters.ChoiceFilter(field_name = 'mic_type', choices=[], label = 'Type', empty_label=None)
    bp_profile = django_filters.ChoiceFilter(field_name = 'bp_profile', choices=[], label = 'BP')
    source = django_filters.ChoiceFilter(field_name = 'source', choices=[], label = 'Source')

    def filter_all_fields(self, queryset, name, value):
        if value:
            # print(f"filtr fields is {self._meta.model._meta.fields})")
            fields=[f.name for f in self._meta.model._meta.fields]
            value=value         
            similarity = Greatest(
                TrigramSimilarity('organism_id__organism_name', value),
                TrigramSimilarity('drug_id__drug_name', value),
                TrigramSimilarity('mic', value),
                TrigramSimilarity('mic_unit', value),
                TrigramSimilarity('zone_diameter', value),
                TrigramSimilarity('mic_type__dict_value', value),
                TrigramSimilarity('source', value),
                TrigramSimilarity('bp_profile', value),
                TrigramSimilarity('bp_source', value),
               )
            queryset=queryset.annotate(similarity=similarity)#(similarity=TrigramSimilarity('drug_id__drug_name', value),)

            return queryset.filter(similarity__gt=0.1)#(q_object)
        return queryset
    
    def filter_all_fields_deep(self, queryset, name, value):
        queryset=self.filter_all_fields(queryset, name, value)
        return queryset

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['mic'].label='MIC'
        self.filters["mic_type"].extra['choices']=[('', ''),] + [(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(MIC_Pub.Choice_Dictionary['mic_type'])]
        self.filters['bp_profile'].extra["choices"] = self.Meta.model.get_field_choices(field_name='bp_profile')
        self.filters['source'].extra["choices"] = self.Meta.model.get_field_choices(field_name='source')

    
    class Meta:
        model=MIC_Pub
        fields=["fOrg_ID", "fOrg_Name", "fDrug_Name"]
        fields +=list(model.HEADER_FIELDS.keys())
        exclude = ['organism_id.organism_id',
                   'organism_id.organism_name',
                   'drug_id.drug_name',
                   ]


# -----------------------------------------------------------------
class Breakpointfilter(Filterbase):
# -----------------------------------------------------------------
    drug_name = django_filters.CharFilter(field_name='drug_id__drug_name', lookup_expr='icontains', label="Drug")
    bp_type=django_filters.ChoiceFilter(field_name='bp_type', choices=[], empty_label=None)
    notorg_rank=django_filters.ChoiceFilter(field_name='notorg_rank', choices=[], empty_label=None)
    org_rank=django_filters.ChoiceFilter(field_name='org_rank', choices=[], empty_label=None)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["bp_type"].extra['choices']=[('', ''),] + [(obj.dict_value, obj) for obj in Dictionary.get_filterobj(Breakpoint.Choice_Dictionary['bp_type'])]
        self.filters["notorg_rank"].extra['choices']=[('', ''),] + [(obj.dict_value, obj) for obj in Dictionary.get_filterobj(Breakpoint.Choice_Dictionary['notorg_rank'])]
        self.filters["org_rank"].extra['choices']=[('', ''),] + [(obj.dict_value, obj) for obj in Dictionary.get_filterobj(Breakpoint.Choice_Dictionary['org_rank'])]

    class Meta:
        model=Breakpoint
        fields = ['drug_name',]
        fields += list(model.HEADER_FIELDS.keys()) 
        exclude=['drug_id.drug_name']
 
