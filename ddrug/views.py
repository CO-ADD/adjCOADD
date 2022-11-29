from django.shortcuts import render

# Create your views here.
#===================== A Example =================================================================================#
# def home(req): 
#     clearIMGfolder()
#     # search function
#     if req.method=='POST':
#         search =req.POST.get('search')
#         field=req.POST.get('field')
#         if field=='Organism_Name':
#             result=Taxonomy.objects.filter(astatus=1, Organism_Name__contains=search)
#     else:
#         result=Taxonomy.objects.filter(astatus=1)
        
#     objects_all=result
#     p=Paginator(objects_all, 24)
#     page_number = req.GET.get('page')
#     page_obj=p.get_page(page_number)
#     for object_ in page_obj:
#         m=Chem.MolFromSmiles('Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1')
#         molecule_to_svg(m, object_.Organism_Name)
    
#     context={
#         'page_obj':page_obj,
#         'chose':Taxonomy.Choice_Dictionary   
       
#     }
  
#     return render(req, 'organism/chem.html', context)
