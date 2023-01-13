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
    'organism_name':'Organism Name',  
    'code':'Code', 
    'lineage':'Lineage', 
    'tax_rank':'Rank',
    'division':'Division', 
    'org_class':'Class',
    }

ORGANISM_FIELDs = {
    'organism_name':'Organism Name',
    'strain_ids':'Strain IDs',
    'strain_type':'Strain Type',
    'strain_panel':'Panel',
    'res_property':'Phenotype',  
    'gen_property':'Genotype', 
    'biologist':'Biologist',
    'strain_origin':'Origin',
    'organism_id':'Organism ID', 
    }

ORGANISM_BATCH_FIELDs = {
    "orgbatch_id":"OrgBatch ID",
    "supplier":"Supplier",
    "supplier_code":"Supplier Code",
    "stock_date":"Stock Date",
    "stock_level":"Stock Levels",
    "qc_status":"QC_Status",
    "batch_notes":"Batch Notes",
    "biologist":"Biologist"
}
ORGANISM_STOCK_FIELDs={
    "orgbatch_id":"OrgBatch ID",
    "stock_id":"Stock ID",
    "stock_note":"Stock Note",
    "stock_type":"Stock Type",
    "stock_date":"Stock Date",
    "biologist":"Biologist"
}


ORGANISM_CULTR_FIELDs = {
    
    "organism_id":"Organism ID",
    "culture_type":"Culture Type",
    "media_use":"Media Use",
    "atmosphere":"Atmosphere",
    "temperature":"Temperature",
    "labware":"Labware",
    "notes":"Media",
    "biologist":"Biologist"
}

# -dDrug Settings ---------------------------------------------------
DRUG_FIELDs = {
    
    "drug_id":"Drug ID",
    "drug_name":"Drug Name",
    "drug_othernames":"Other Name",
    "drug_codes":"Drug Code",
    "drug_type":"Drug Type",
    "drug_note":"Drug Note"
   
}