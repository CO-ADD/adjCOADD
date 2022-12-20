#
# Application Constants/Settings 
#

# -dOrganism Settings ---------------------------------------------------
ORGANISM_CLASSES = ['GN','GP','MB','FG','MA']
ORGANSIM_SEP = "_"
ORGBATCH_SEP = "_"

TAXONOMY_FIELDs = ['Specie', 'Other Names',  'Code',  'Class', 'NCBI Tax ID', 'NCBI Parent Tax ID','Taxonomy Rank','Division', 'lineage']
ORGANISM_FIELDs = ['Organism ID', 'Organism Name',  'Risk Group',  'Pathogen', 'Lab Restriction', 'Bio Approval']
ORGANISM_BATCH_FIELDs = ["OrgBatch ID","Organism ID","Batch No","Batch Notes",]
ORGANISM_BATCH_modelFIELDs=["orgbatch_id","organism_id","batch_no","batch_notes",]
# ORGANISM_BATCH_FIELDs = ["OrgBatch ID","Organism ID","Batch No","Batch Notes","QC Notes","QC Records","Supplier","Supplier Code","Supplier PO","Stock Date","Stock Levels","Biologist"]
# ORGANISM_BATCH_modelFIELDs = ["orgbatch_id","organism_id","batch_no","batch_notes","qc_status","qc_record","supplier","supplier_code","supplier_po","stock_date","biologist"]

# -dDrug Settings ---------------------------------------------------
