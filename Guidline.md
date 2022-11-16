# Coadd Coding Guidline

Coding guideline is to ensure the code is consistent, minimize the erros and help maintain easlily, which provides guidance for building a project and functionalities structure, naming rules and coding reuseable code blocks in different section. 
Also provids certain constructions explaination in the last section.

## Prerequisites
### Documents:
- Django Documents for name convention
- GetStartedWithCoadd: Computer setting, Environment setting
- Test Documents with IssueLog

### Project Overview:


## Section 1. Project Directories 
The Directory includs project level and application level, displayed as the following:
### 1 Project level
```
Django Project Root
├── aa_django-core(application level refer below...)
├── django-app-utilities
├── django-app1
...
├── static
├── templates/base.html
├── .gitignore
├── manage.py
├── readme.md
└── requirements.txt
```

### 2 Application level
```
Django Project Root
├── django-core
│   ├── __init__.py
│   ├── asgi.py
│   ├── settings.py
│   ├── urls.py
│   ├── wsgi.py
├── django-app-utilities
│   ├── migrations
│   ├── templates
│   ├── __init__.py
│   ├── admin.py
│   ├── app.py
│   ├── forms.py
│   ├── models.py
│   ├── signals.py
│   ├── views.py
...(maybe more)
├── django-app01
│   ├── migrations
│   ├── templates
│   ├── __init__.py
│   ├── admin.py
│   ├── app.py
│   ├── forms.py
│   ├── dbRouter.py
│   ├── signals.py(maybe)
│   ├── models.py
│   ├── utils.py(maybe)
│   ├── views.py
...(maybe more)
```

### 3 In database Postgres
```
Database name(projectname)
├── app schema (for django migration project level info(extension, session, contenttype...); django-app-utilities)
├── application schema 1 (per one application)
....
├── application schema N
```
in Django Settings.py the code is presented like the following:
```

    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'OPTIONS':{'options': '-c search_path=app'},
        'NAME': 'orgdb',
        # 'USER': 'postgres', #os.environ.get('db_user'),
        # 'PASSWORD':os.environ.get('db_password', 'password'),
        'HOST': 'Localhost',
        'PORT': '5432',
    },
    'drugs_db': {
        ...
        'OPTIONS':{'options': '-c search_path=aa_chem,app'},
        ...
    },
      'test_sch': {
        ...
        'OPTIONS':{'options': '-c search_path=test_sch,app'},
        ...
    }
```

## Section 2. Naming Rules

### 1 Variable Name
- General Variables:
    - context: variable collection will be send to html template in render
    - i or item  : item in iterable collection
    - instance: model instance used in function based create, update views
    - objects : a queryset objects from a model, used to put query results to html template.
    - object_:  a single model query result
    - qs  : queryset
    - req : request
    - _fk: indicate Foreign Key type
    - _list, _tul, _dict...: variable is list, tuple and dictionary... type
    
- Special variables will be commented in Functions and Classes

### 2 Project level naming
A project name start with a_, therefore the core file name will be with "a_..." positioned on the first line under project root. For applications' name following the below: 
- utility application (for admin management, functions sharing...purpose)starts with "a"
- django applications start with "d"
other django files and files following django and python name convention: templates, static, uploads...

### 3 Application level naming
Applications' module name following django and python name convention: admin.py, models.py, utils.py ....
In coding process following the below rules:
- Model naming: 
    - class names in models.py are in Singular
    - Field Names can contain a underscore starting with a capital letter
    - CharField with max_length: 10, 15, 25, 50, 125, 150, 200, 220, 500 ...
    - Model variable is in camelcase, starts with "semantic word(purpose) in small letters", follows with referred model field name if exists, ends with dataType. 
    - ForeignKey setting: <b>need to discussing with DO_Nothing or Set_Null??</b>
- Forms name is in CamelCase with (form purpose if needs to do: Create, Update...)(Modelname)_form

- View Functions variable names
    - Views in Class will be defined with (Modelname)(Functions: List/Card/Update/Create...)(View); Views in Function will be defined with (Modelname in small letters)(Functions: List/Card/Update/Create...)
    - For django build-in function variables following django convention e.g. for pagination, using pag_num ...
    - An objects collection or object of model instances, which need to pass to template, will be defined with <b>objects / object_(model name in small letter in template)</b>
    - For other variables and bridge variables passed to Templates/Model instance/Forms or received from Templates input in camelCase. with belonged model name at the beginning all letters in small case with following semantic word (capital letter at the beginning) and ending with data type e.g. ```organismNameFkChar```
- Urls names:
   - Pass primary key by using ```.../<int:pk>``` or ```.../<str:pk>```. refer to django document [link](https://docs.djangoproject.com/en/4.0/ref/models/instances/#the-pk-property).
   - path defined with ```/(projectname)/(modelclassname)/(overview_list/detail/delete/update) or (ID_name)``` except homepage with ```/(projectname)/home```
- Js files name and variables name...

<br>

## Section 3. Programming and Code Blocks guidance
This provides guidance for programming each module.
1. In python module file whenever need to import other modules, importing codes are following the order: 
<b>Other packages: e.g. os, rdkit... -> Django buildin -> project local modules: e.g. from app.models import comebined with an alphabetical order.</b>

2. Environment variable in settings.py:
In develop mode: set environment variable in Django shell is enough. (previously exp... set variables in file .venv/bin/activate)
In Product refer to [link](https://stackoverflow.com/questions/44693485/where-do-i-set-environment-variables-for-django):

3. Models.py for model Classes <span style="color: yellow">...not completed</span>

4. Forms.py.
- Using ModelForm to inherite models fields. 
- Avoid making programs in Views Function and move these to Form class methods.
e.g.:
```
class CreateOrganism_form(ModelForm):  
    
   
    Organism_Desc=forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    ...
    Bio_Approval = forms.ChoiceField(choices=(('',''),('','')),  widget=forms.Select(attrs={'class':'form-select'}))
    
    Organism_Name=forms.ModelChoiceField(queryset=Taxonomy.objects.all(), widget=forms.HiddenInput(),required=False,)

   # programs for get value from Models and Views
   
    def __init__(self, Strain_Type_choices,  *args, **kwargs):      
        super(CreateOrganism_form, self).__init__(*args, **kwargs)
        self.strainTypeChoices= Strain_Type_choices
        self.fields['Strain_Type'].widget = forms.CheckboxSelectMultiple(choices=self.strainTypeChoices)
        self.fields['Oxygen_Pref'].choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Oxygen_Pref'])
        ...
    
    def get_object(self, Organism_Name):
        Organism = get_object_or_404(Taxonomy, Organism_Name=Organism_Name)
        print(Organism.Class.Dict_Value)
        
        self.fields['Organism_Name'].queryset=Organism
        print(f'Thisis from def get_object: {self.fields["Organism_Name"]}')
       

    # inheriting model fields       
    class Meta:
        model=Organisms
        exclude = ['Organism_ID']
```

5. In views.py. 
For simple model Class (only defined with modelfields ) using GenericView for list/card View. Detail, Create, Update and Delete view using function based view(at this stage). For complexed model Class list/card View may use Function based view. Keep Views functions as minimal.
- function based view functions for each model class:
    - example(create new entry):
```
def createOrgnisms(req):  # function name refer to [namerules](#application-level-naming)
    '''
    Function View Create new Organism table row with foreignkey: Taxonomy and Dictionary.   
    '''
    #define function variables directly here if possible
    Strain_Type_choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Strain_Type']) # 
    kwargs={}
    kwargs['user']=req.user 
    ...

    if req.method=='POST':
        form=CreateOrganism_form(Strain_Type_choices,  req.POST,)
        Strain_Type_list=req.POST.getlist('Strain_Type')
        
        try:                                                 #code block
            if form.is_valid():
                print("form is valid")  
                Organism_Name=req.POST.get('Organism_Name')
                form.get_object(Organism_Name) 
                instance=form.save(commit=False)
                instance.save(**kwargs)
                print("saved")
                return redirect("org_list")
            else:
                print(f'something wrong...{form.errors}')     #debug print
                return redirect(req.META['HTTP_REFERER'])  
        except Exception as err:
            print(f'error is {form.errors} with {err}')
            return redirect(req.META['HTTP_REFERER'])  

    else:
        form=CreateOrganism_form(Strain_Type_choices,)
 
    return render(req, 'aa_chem/createForm/Organism.html', { 'form':form, 'Strain_Type':Strain_Type_choices})
```
  <br> 

- "Debug prints" will be used in development stage, some will be deleted in production stage.
- If only to find a queryset has at least one recorder, using querySet.exists() instead of querySet. e.g. replace ```if querySet``` to ```if querySet.exists()``` [link](https://docs.djangoproject.com/en/4.1/ref/models/querysets/#django.db.models.query.QuerySet.exists) 


<br>

## Section 4. Functional testing:
<span style="color:yellow">Not completed</span>

## Section 5. Specific construction
1. To solve multiple users access data entries problem (concurrency), <span style="color: yellow">using django build_in ```with transaction... select_and_update....```.</span>

2. option for Choicesfield in models [not implemented]:
 choice fields can be placed in settings.py and modified via python environment variables.

3. Ldap login information with signals.py

4. multi_schema set up in postgresql:
```
testschema=# create schema app;
CREATE SCHEMA
testschema=# create schema aa_chem;
CREATE SCHEMA
testschema=# create schema test_sch;
CREATE SCHEMA
testschema=# grant usage on schema app to public;
GRANT
testschema=# grant execute on all functions in schema app to public;
GRANT
testschema=# alter default privileges in schema app
testschema-# grant execute on functions to public;
ALTER DEFAULT PRIVILEGES
testschema=# alter default privileges in schema app
testschema-# grant usage on types to public;
ALTER DEFAULT PRIVILEGES
testschema=# create extension if not exists rdkit schema app;
CREATE EXTENSION
```
