#
# Application Constants/Settings 
#

# -dOrganism Settings ---------------------------------------------------
ORGANISM_CLASSES = ['GN','GP','MB','FG','MA']
ORGANSIM_SEP = "_"
ORGBATCH_SEP = "_"
# column name can be edited here 
# make a dictioinary  with Key and value, if value is none choose verbose name else choose the dictionary name.
# TAXONOMY_FIELDs = ['Specie',  'Code',  'Class',  'NCBI Parent Tax ID','Taxonomy Rank','Division', ]
TAXONOMY_FIELDs = {
    'organism_name':'Specie_test',  
    'code':'Code_test',  
    'org_class':'Class_test',
    'tax_rank':'Taxonomy Rank',
    'division':'Division', 
    }

ORGANISM_FIELDs = {
    'organism_id':'Organism ID', 
    'organism_name':'Organism Name',  
    'risk_group':'Risk Group',  
    'pathogen_group':'Pathogen', 
    'strain_code':'Strain Code'
    }

ORGANISM_STOCK_FIELDs={
    "orgbatch_id":"OrgBatch ID",
    "stock_id":"Stock ID",
    "stock_note":"Stock Note",
    "stock_type":"Stock Type",
    "stock_date":"Stock Date",
    "biologist":"Biologist"
}

ORGANISM_BATCH_FIELDs = {
    "orgbatch_id":"OrgBatch ID",
    "supplier":"Supplier",
    "stock_date":"Stock Date",
    "stock_level":"Stock Levels",
    "biologist":"Biologist"
}

# -dDrug Settings ---------------------------------------------------
