# Coadd Guidline

## pre-requisites

## Structure
Django project structure includs project level and application level.
### Project level
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
### Application level
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
### In database Postgres
```
Database name(projectname)
├── public schema (for django migration project level info(extension, session, contenttype...); django-app-utilities)
├── application schema 1 (per one application)
....
├── application schema N
```
## Name rules
### project level names
A project name start with a_, therefore the core file name will be with "a_..." positioned on the first line under project root. For applications' name following the below: 
- utility application (for admin management, functions sharing...purpose)starts with "a"
- django applications start with "d"
other django files and files following django and python name convention: templates, static, uploads...

### application level names
Applications' module name following django and python name convention: admin.py, models.py, utils.py ....
In coding process following the below rules:
- Model naming: 
    - class names in models.py are in Singular
    - Field Names can contain a underscore starting with a capital letter
    - CharField with max_length: 10, 50, 150, ...
- Forms name is in CamelCase with (Modelname)(form purpose: Create, Update...)_form

- View Functions variable names
    - For django build-in variables refer to convention e.g. pagination pag_num ...
    - For varibles and bridge variables passed to Templates/Model instance/Forms or received from Templates input in camelCase. with belonged model name at the beginning all letters in small case with following semantic word (capital letter at the beginning) and ending with data type e.g. ```organismNameFkChar```
- Urls names:
   - Pass primary key by using ```.../<int:pk>``` or ```.../<str:pk>```. refer to django document [link](https://docs.djangoproject.com/en/4.0/ref/models/instances/#the-pk-property).
   - path defined with ```/(projectname)/(modelclassname)/(overview_list/detail/delete/update) or (ID_name)``` except homepage with ```/(projectname)/home```
- Js files name and variables name...

## 
