# import re
# import csv
import os
# from os.path import exists
# from django.contrib.admin.views.decorators import staff_member_required
# from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.core.paginator import Paginator
# from django.core import management
from django.shortcuts import HttpResponse, render, redirect
# from django_rdkit.models import * 
from aa_chem.models import Mytest, Organisms, Taxonomy #Drugbank,,Organisms,
# import aa_chem.models
from app.models import Dictionaries, Mytest2
from app.utils import molecule_to_svg, clearIMGfolder
from rdkit import Chem

from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic import ListView
from django.views.generic.detail import DetailView
from django.urls import reverse_lazy
from .forms import CreateNewOrgForm, UpdateNewOrgForm
# from model_utils import Choices
# from django.forms import modelform_factory
from django.shortcuts import get_object_or_404
from django.http import JsonResponse




# # =======================================Taxo View=============================================================================#

# # Create your views here.
# @login_required
def home(req): 
    clearIMGfolder()
    # search function
    if req.method=='POST':
        search =req.POST.get('search')
        field=req.POST.get('field')
        if field=='Organism_Name':
            result=Taxonomy.objects.filter(Organism_Name__contains=search)
    else:
        result=Taxonomy.objects.all()
        
    objects_all=result
    p=Paginator(objects_all, 24)
    page_number = req.GET.get('page')
    page_obj=p.get_page(page_number)
    for object_ in page_obj:
        m=Chem.MolFromSmiles('Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1')
        molecule_to_svg(m, object_.Organism_Name)
    
    context={
        'page_obj':page_obj,
        'chose':Taxonomy.Choice_Dictionaries   
       
    }
  
    return render(req, 'aa_chem/chem.html', context)




class TaxoListView(ListView):
    model=Taxonomy
    paginate_by = 24
    fields='__all__'
    template_name = 'aa_chem/taxoListview.html'

    def get_context_data(self, **kwargs):
        for filename in os.listdir("static/images"):
                file_path=os.path.join("static/images", filename)
                try:
                    os.unlink(file_path)
                    print("removed!")
                except Exception as err:
                    print(err)

        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        objects_all=[object_ for object_ in context["objects"]]
        p=Paginator(objects_all, 24)
        
        for object_ in p.get_page(self.request.GET.get('page')):

            m=Chem.MolFromSmiles('Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1')
            molecule_to_svg(m, object_.Organism_Name)
            
        return context



class TaxoCreateView(CreateView):
    model=Taxonomy
    fields='__all__'
    template_name = 'aa_chem/taxoCreate.html'
    success_url = reverse_lazy('compounds')

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
    
        context["objects"]=self.model.objects.all()
        return context

class TaxoUpdateView(UpdateView):
    model=Taxonomy
    fields='__all__'
    template_name = 'aa_chem/taxoUpdate.html'
    success_url = reverse_lazy('compounds')

# # ===============================================================OrgView==============================================#
# # class OrgCreateView(CreateView):
# #     model=Organisms
# #     form_class=CreateNewOrgForm
# #     template_name = 'aa_chem/orgCreate2.html'
# #     success_url = reverse_lazy('compounds')


# ============================Create new Organism======================================================#
# =============================1. Ajax Call search Taxo=============================================#

def searchTaxo(req):
    if req.headers.get('x-requested-with') == 'XMLHttpRequest':
        res=None
        taxo=req.POST.get('inputtext')
        # print(taxo)
        qs=Taxonomy.objects.filter(Organism_Name__icontains=taxo)
        if len(qs)>0 and len(taxo)>0:
            data=[]
            for i in qs:
                item={
                    'name':i.Organism_Name,
                }
                data.append(item)
            res=data
        else:
            res='No games found...'
        
        return JsonResponse({'data':res})
    return JsonResponse({})

# =============================2. Create New Organism based Taxonomy====================================================================================================#
    Choice_Dictionaries = {
        'Risk_Group':'Risk_Group',
        'Pathogen_Group':'Pathogen_Group',
        'Bio_Approval':'Bio_Approval',
        'Oxygen_Pref':'Oxygen_Preference',
        'MTA_Status':'License_Status',
        'Strain_Type':'Strain_Type',
    }
# =======================================================================================================================================================================

def querysetToChoiseList_Dictionaries(model_name, field_name):
    options=model_name.objects.filter(Dictionary_ID=field_name).values('Dict_Value', 'Dict_Desc')
    choices=[tuple(d.values()) for d in options]
    return choices

Strain_Type_choices=querysetToChoiseList_Dictionaries(Dictionaries, 'Strain_Type')
Oxygen_Pref_choices=querysetToChoiseList_Dictionaries(Dictionaries, 'Oxygen_Preference') #[tuple(d.values()) for d in Oxygen_Pref_options]
Risk_Group_choices=querysetToChoiseList_Dictionaries(Dictionaries, 'Risk_Group')

Pathogen_Group_choices=querysetToChoiseList_Dictionaries(Dictionaries, 'Pathogen_Group')

# ================================================================================================================================================================
def createOrgnisms(req):
    '''
    Function View Create new Organism table row with foreignkey: Taxonomy and Dictionary. 
    '''
    # Strain_Type=Dictionaries.objects.filter(Dictionary_ID='Strain_Type') #===multi choice
    Strain_Type_choices=querysetToChoiseList_Dictionaries(Dictionaries, 'Strain_Type')
    Oxygen_Pref_choices=querysetToChoiseList_Dictionaries(Dictionaries, 'Oxygen_Preference') #[tuple(d.values()) for d in Oxygen_Pref_options]
    print(Oxygen_Pref_choices)
    Risk_Group_choices=querysetToChoiseList_Dictionaries(Dictionaries, 'Risk_Group')
    # MTA_Status_choices=querysetToChoiseList_Dictionaries(Dictionaries, '')
    # Biological_Approval_choices=querysetToChoiseList_Dictionaries(Dictionaries, '')
    Pathogen_Group_choices=querysetToChoiseList_Dictionaries(Dictionaries, 'Pathogen_Group')
    print(Pathogen_Group_choices)

    # Retreive Values for each column========================
    if req.method=='POST':
        form=CreateNewOrgForm(Strain_Type_choices, Oxygen_Pref_choices, Risk_Group_choices, Pathogen_Group_choices, req.POST)
        Organism_Name=req.POST.get('Organism_Name')
        Organism_Name_fk=get_object_or_404(Taxonomy, Organism_Name=Organism_Name)
        Strain_Type_list=req.POST.getlist('Strain_Type')
       
        try:
            if form.is_valid():
                instance=form.save(commit=False)
                instance.Organism_Name=Organism_Name_fk
                instance.save()
                print("saved")
                return redirect("/")
        
        except Exception as err:
            print(err)
    else:
        form=CreateNewOrgForm(Strain_Type_choices, Oxygen_Pref_choices, Risk_Group_choices, Pathogen_Group_choices)
 
    return render(req, 'aa_chem/createForm/Organism.html', { 'form':form})


#=======================================================================================================================================================

def organismDetail(req, Organism_ID):
    obj=get_object_or_404(Organisms, Organism_ID=Organism_ID)
    context={'Organism': obj}
    return render(req, "aa_chem/readForm/Organism_detail.html", context)



def updateOrganism(req, Organism_ID):
    context={}
    obj=get_object_or_404(Organisms, Organism_ID=Organism_ID)
    form=CreateNewOrgForm(Strain_Type_choices, Oxygen_Pref_choices, Risk_Group_choices, Pathogen_Group_choices,instance=obj)
    if req.method=='POST':
        Organism_Name=req.POST.get('Organism_Name')
        Organism_Name_fk=get_object_or_404(Taxonomy, Organism_Name=Organism_Name)
        obj.delete()
        form=CreateNewOrgForm(Strain_Type_choices, Oxygen_Pref_choices, Risk_Group_choices, Pathogen_Group_choices, req.POST, instance=obj)    
        if form.is_valid():
            instance=form.save(commit=False)
            instance.Organism_Name=Organism_Name_fk
            instance.save()
            print('save updated')
            return redirect("/")
    
    context["form"]=form
    # context['Organisms_Name']=obj.Organism_Name
    return render(req, "aa_chem/updateForm/Organism.html", context)

def deleteOrganism(req, Organism_ID):
    
    obj=get_object_or_404(Organisms, Organism_ID=Organism_ID)
    # if req.method=='POST':
    #     obj.delete()
    #     print("deleted!")
    #     return redirect("/")
    # context={'obj':obj}
    obj.delete()
    return redirect("/")
    # return render(req, "aa_chem/deleteForm/Organism_del.html", context)

# ==============================List View ===============================================================
class OrgListView(ListView):
    model=Organisms
    paginate_by = 24
    fields='__all__'
    template_name = 'aa_chem/readForm/Organism_list.html'
    
    def get_context_data(self, **kwargs):
        for filename in os.listdir("static/images"):
                file_path=os.path.join("static/images", filename)
                try:
                    os.unlink(file_path)
                    print("removed!")
                except Exception as err:
                    print(err)

        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        objects_all=[object_ for object_ in context["objects"]]
        p=Paginator(objects_all, 24)
        
        for object_ in p.get_page(self.request.GET.get('page')):
       
            # if exists(f"static/images/{object_.id}.png"):
            #     print("file exists")
            #     return context
            # else:

            
            m=Chem.MolFromSmiles('Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1')
        
            molecule_to_svg(m, object_.Organism_ID)
            
        return context

class OrgCardView(OrgListView):
  
    template_name = 'aa_chem/readForm/Organism_card.html'

   



# # ==============================Import Excel files===========================================================#
import pandas as pd
from django.conf import settings
from django.core.files.storage import FileSystemStorage

def import_excel_taxo(req):
    print('importing....')
    try:
        if req.method=='POST' and req.FILES['myfile']:
            myfile=req.FILES['myfile']
            fs=FileSystemStorage()
            filename=fs.save(myfile.name, myfile)
            uploaded_file_url=fs.url(filename)
            excel_file=uploaded_file_url
            print(excel_file)
            exmpexceldata=pd.read_csv("."+excel_file, encoding='utf-8')
            print(type(exmpexceldata))
            dbframe=exmpexceldata
            for dbframe in dbframe.itertuples():
                class_fkey=Dictionaries.objects.filter(Dict_Desc=dbframe.ORGANISM_CLASS)
                print(class_fkey)
                division_fkey=Dictionaries.objects.filter(Dict_Desc=dbframe.DIVISION)
                linea=dbframe.LINEAGE.split(",")
                print(division_fkey)
                # fromdata_time_obj=dt.datetime.strptime(dbframe.DOB, '%d-%m-%Y')
                obj, created=Taxonomy.objects.get_or_create(Organism_Name=dbframe.ORGANISM_NAME, Other_Names=dbframe.ORGANISM_NAME_OTHER, Code=dbframe.ORGANISM_CODE, 
                    Class=class_fkey[0], Tax_ID=dbframe.TAX_ID,Parent_Tax_ID=dbframe.PARENT_TAX_ID, 
                    Tax_Rank=dbframe.TAX_RANK, Division=division_fkey[0], Lineage=linea
                    )
                print(type(obj))
                # obj.save()
            
            return render(req, 'aa_chem/importexcel_taxo.html', {'uploaded_file_url': uploaded_file_url})
    except Exception as err:
        print(err)
    return render(req, 'aa_chem/importexcel_taxo.html', {})

def import_excel_dict(req):
    print('importing....')
    try:
        if req.method=='POST' and req.FILES['myfile']:
            myfile=req.FILES['myfile']
            fs=FileSystemStorage()
            filename=fs.save(myfile.name, myfile)
            uploaded_file_url=fs.url(filename)
            excel_file=uploaded_file_url
            print(excel_file)
            exmpexceldata=pd.read_csv("."+excel_file, encoding='utf-8')
            print(type(exmpexceldata))
            dbframe=exmpexceldata
            for dbframe in dbframe.itertuples():
                print("iam here")             
                # int_order=int(0+dbframe.Dict_View_Order)
                # fromdata_time_obj=dt.datetime.strptime(dbframe.DOB, '%d-%m-%Y')
                obj, created=Dictionaries.objects.get_or_create(Dictionary_ID=dbframe.Dictionary_ID, Dictionary_Class=dbframe.Dictionary_Class, Dict_Value=dbframe.Dict_Value, Dict_Desc =dbframe.Dict_Desc, Dict_Value_Type =dbframe.Dict_Value_Type, Dict_View_Order =dbframe.Dict_View_Order
                    )
                print(type(obj))
                # obj.save()
            
            return render(req, 'aa_chem/importexcel_dict.html', {'uploaded_file_url': uploaded_file_url})
    except Exception as err:
        print(err)
    return render(req, 'aa_chem/importexcel_dict.html', {})

#======================================================Export Data Views Function==================================================#  
# # @user_passes_test(lambda u: u.is_admin)
# @permission_required('importdata')
def exportCSV(req):
    response=HttpResponse(content_type='text/csv')
   
    writer=csv.writer(response)
    writer.writerow(['S', 'O','C', "Cl", "NC", "NP", "T","D", "Di", "Lineaage"])#['id', 'drug_name','drug_mol']
    query=Taxonomy.objects.all()
    comp_list=[comp for comp in query]

    for comp in comp_list:
        writer.writerow(['S', 'O','C', "Cl", "NC", "NP", "T","D", "Di", "Lineaage"]) #[comp.id, comp.drug_name, comp.smiles]

    response['Content-Disposition']='attachment; filename="Taxo_export.csv"'
    return response
   