### Development Issues and Solutions:

1. (i) migration error about contenttype:
   (s) delete public schema, migrate, add public, migrate again.

2. (i) File "/home/jzhong/.local/lib/python3.10/site-packages/rdkit/Chem/Draw/**init**.py", line 14, in <module>
   from rdkit.Chem.Draw import rdMolDraw2D
   ImportError: libXrender.so.1: cannot open shared object file: No such file or directory

   (s)fixed it by installing libXrender:

   sudo apt-get install libxrender1
   
   or install Cairo / GTK (or both)

3. (i)Recommanded steps for multi schemas migrate:

   (s)
   step 1 run `manage.py migrate <yourapp> --database <wanted database name>`
   step 2 run `manage.py migrate` (can be deleted, after switch app to public schema!!!)

4. (i)Error: migrations.exceptions.InvalidBasesError...This can happen if you are inheriting ....

   (s)check if wrongly deleted directory "migrations", if it is, solve it by creating a migrations directory at the root of your app and adding an empty **init**.py file might resolve your issue.

5. (i) When save a new record in OrgDB model with foreignkey Taxonomy, DB ERROR: column strain_organisms.Organism_Name_id does not exist at character 310, in Django :return self.cursor.execute(sql, params)
   psycopg2.errors.UndefinedColumn: column strain_organisms.Organism_Name_id does not exist
   LINE 1: ...deleted_by_id", "strain_organisms"."Organism_ID", "strain_or...

   (s)This happens when wrongly deleted some user and records. To solve it, - have to clear this table records from DB, - run migrate app(contains this table) zero, - then again, run migrate app.


6. (i) implement signal function to class Audit  and let each model inherit this function

   (i) clear up html files, .py files...

   (i) creating new org takes long time (200k Taxo Foreign Key)

   (i) test remote server DB
   (i) templates in the individual apps
   (i) choice: value with explain text upon choosing, after choosing only pich up value. refer OrgDB input. 

7. Error: (fields.E304) Reverse accessor'Group.user_set' for 'app.ApplicationUser.groups' clashes Reverse accessor clashes with reverse accessor for 'app.User.groups' 
Hint: add or change a related_name argument to the definition for 'app.ApplicationUser.groups' or 'app.User.groups'.
(Error happens during inherite class and with foreignkey)


8. (i) Cannot DELETE or UPDATE with on_delete=models.DO_NOTHING
   (s) Change to PROTECT

#### on_delete for Foreignkey change to on_delete=models.PROTECT.
CASCADE
Cascade emulates the SQL constraint of ON DELETE CASCADE. Whenever the referenced object (post) is deleted, the objects referencing it (comments) are deleted as well. 

PROTECT argument of the ForeignKey on_delete option prevents the referenced object from being deleted if it already has an object referencing it in the database. Put simply, Django will prevent a post from deletion if it already has comments. we tried deleting this POST that already has a comment, it will raise PROTECTEDERROR, remove the comment firstly then remove the post.

SET_NULL argument of the ForeignKey on_delete option is only available to you when you have set the null option on the ForeignKey field to True. When you use this argument, and, in our case, delete a post, it is going to leave the comments in the database without deleting it.

SET_DEFAULT
This argument on the ForeignKey on_delete option requires you to set a default value when defining the relationship. When you delete a post that has comments, the comments are automatically assigned to a default post you had set when creating the model.

DO_NOTHING
As the name implies, it does nothing when a referenced object is deleted. This is essentially discouraged because it defeats the purpose of an RDBMS. <b> HERE I CANNOT DELETE OR UPDATE any recorder!</b>

9. (i)Difficulty in AuditModel to retrieve logged or request user information
(s) in the view funtion: create, update and delete to set the value, e.g. instance.acreate_by=request.user 

10. (i) exception: python manage.py runserver exception in thread django-main-thread typeerror issubclass() arg 1 must be a class app_config_class
     (s) settings.py - installed_Apps=[...'app',...] change to [...'app.apps.AppConfig',...]


