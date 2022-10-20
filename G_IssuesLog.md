### Development Issues and Solutions:

1. (i) migration error about contenttype:
   (s) delete public schema, migrate, add public, migrate again.

2. (i) File "/home/jzhong/.local/lib/python3.10/site-packages/rdkit/Chem/Draw/**init**.py", line 14, in <module>
   from rdkit.Chem.Draw import rdMolDraw2D
   ImportError: libXrender.so.1: cannot open shared object file: No such file or directory

   (s)fixed it by installing libXrender:

   sudo apt-get install libxrender1

3. (i)Recommanded steps for multi schemas migrate:

   (s)
   step 1 run `manage.py migrate <yourapp> --database <wanted database name>`
   step 2 run `manage.py migrate` (can be delete, after switch app to public schema!!!)

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


   (i) audit model including Foreign key User in aa_chem app result in migrating User models to the same schema
       plan: using django entry_log. create a new app for audit purpose

   (i) sequence short name 

   (i) double check indexes setting, in django db_index

7. Error: (fields.E304) Reverse accessor'Group.user_set' for 'app.ApplicationUser.groups' clashes Reverse accessor clashes with reverse accessor for 'app.User.groups' 
Hint: add or change a related_name argument to the definition for 'app.ApplicationUser.groups' or 'app.User.groups'.
(Error happens during inherite class and with foreignkey)