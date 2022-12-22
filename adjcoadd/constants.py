#
# Application Constants/Settings 
#

# -dOrganism Settings ---------------------------------------------------
ORGANISM_CLASSES = ['GN','GP','MB','FG','MA']
ORGANSIM_SEP = "_"
ORGBATCH_SEP = "_"

TAXONOMY_FIELDs = ['Specie', 'Other Names',  'Code',  'Class',  'NCBI Parent Tax ID','Taxonomy Rank','Division', ]
ORGANISM_FIELDs = ['Organism ID', 'Organism Name',  'Risk Group',  'Pathogen', 'Strain Code']
ORGANISM_BATCH_FIELDs = ["OrgBatch ID", "Batch No","Batch Notes","Supplier","Biologist"]
ORGANISM_BATCH_modelFIELDs=["orgbatch_id","batch_no","batch_notes","supplier","biologist"]
ORGANISM_STOCK_modelFIELDs=["orgbatch_id","stock_id","stock_note","stock_type","stock_date","biologist"]
# ORGANISM_BATCH_FIELDs = ["OrgBatch ID","Organism ID","Batch No","Batch Notes","QC Notes","QC Records","Supplier","Supplier Code","Supplier PO","Stock Date","Stock Levels","Biologist"]
# ORGANISM_BATCH_modelFIELDs = ["orgbatch_id","organism_id","batch_no","batch_notes","qc_status","qc_record","supplier","supplier_code","supplier_po","stock_date","biologist"]

# -dDrug Settings ---------------------------------------------------
