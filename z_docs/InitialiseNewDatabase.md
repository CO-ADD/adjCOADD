#
#
#

1. Drop (Cascade) all tables from apputil, dorganism, ddrug, dgene, ... 
2. Remove all migration folders from applications (apputil, dorganism, ddrug, dgene, ...)
3. Makemigrations and Migrate each application separately, starting with apputil, and finish with general migration

    manage.py makemigrations apputil
    manage.py migrate apputil
    ...
    for <app> in [dorganism,ddrug,dgene,..]:
        manage.py makemigrations <app>
        manage.py migrate <app> --database <app>
    ...

    manage.py migrate

4. Upload apputil data
    Check data in impdata/Data/ApplicationData_vnn.xlsx [User] and [Dictionary]
    Upload data with
        ->impdata
        python upload_OrgDB_Data.py -t User --upload
        python upload_OrgDB_Data.py -t Dictionary --upload

5. Upload Application data from Oracle:OrgDB 
    Upload data with
        ->impdata
        python upload_OrgDB_Data.py -t Taxonomy --upload
            ~ 45min for 207665 entries
        python upload_OrgDB_Data.py -t Organism --upload
            ~ 
        python upload_OrgDB_Data.py -t OrganismBatch --upload
            ~

6. Upload Drug data
    Check data in impdata/Data/DrugData_vnn.xlsx [Drug]
    Upload data with
        ->impdata
        python upload_OrgDB_Data.py -t Drug --upload

7. Upload Antibiogram data from Oracle:CastDB and Xls for LMIC strains
    Check data in impdata/Data/LMIC_Data_vnn.xlsx [Drug]
    Upload data with
        ->impdata
        python upload_OrgDB_Data.py -t MICCOADD --runid PMC045_R01  --upload
        python upload_OrgDB_Data.py -t MICCOADD --runid PMC036_R01  --upload
        python upload_OrgDB_Data.py -t MICCOADD --runid PMC045_FDB1 --upload
        python upload_OrgDB_Data.py -t MICCOADD --runid AntiBio     --upload

        python upload_OrgDB_Data.py -t MICPub --upload

        python upload_OrgDB_Data.py -t MICCollab --upload

8. Upload Vitek data from PDF