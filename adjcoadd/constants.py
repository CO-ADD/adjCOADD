#
# Application Constants/Settings 
#

# -dOrganism Settings ---------------------------------------------------
ORGANISM_CLASSES = ['GN','GP','MB','FG','MA']
ORGANSIM_SEP = "_"
ORGBATCH_SEP = "_"
# column name can be edited here 
# make a dictioinary  with Key and value, if value is none choose verbose name else choose the dictionary name.
DICTIONARY_FIELDs = {
    'dict_value':'Value', 
    'dict_class':'Class',  
    'dict_desc':'Description', 
    
    }


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

VITEKCARD_FIELDs = {
    
    "orgbatch_id_pk":"orgbatch_id",
    "card_barcode":"Barcode",
    "card_type":"Card Type",
    "card_code":"Card Code",
    "expiry_date":"expiry_date",
    "instrument":"instrument",
    "proc_date":"proc_date",
    "analysis_time":"Analysis with",
    "acreated_at":"acreated_at"
   
}

VITEKID_FIELDs = {
    
    "card_barcode":"Barcode",
    "process":"Process",
    "id_organism":"ID organism",
    "id_probability":"ID Probability",
    "id_confidence":"ID Confidence",
    "id_source":"Source",
    "filename":"PDF Name",
    "page_no":"PDF PageNo"
   
}

VITEKAST_FIELDs = {
    
    "card_barcode":"Barcode",
    "drug_id":"Drug",
    "mic":"MIC",
    "process":"Vitek Process",
    "bp_profile":"Break Point",
    "bp_comment":"Comment",
    "bp_source":"Source",
    "selection":"Selection",
    "organism":"Organism",
    "filename":"PDF Filename",
    "page_no":"PDF pageNo"
   
}