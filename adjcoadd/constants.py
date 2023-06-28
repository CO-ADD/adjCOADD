#
# Application Constants/Settings 
#

# -dOrganism Settings ---------------------------------------------------
ORGANISM_CLASSES = ['GN','GP','MB','FG','MA']
ORGANSIM_SEP = "_"
ORGBATCH_SEP = "_"

COMPOUND_SEP = '|'
# column name can be edited here 
# make a dictioinary  with Key and value, if value is none choose verbose name else choose the dictionary name.

# Links:
LinkList={
"urlname": '/dorganism/taxonomy/',
'organism_id': '/dorganism/organism/',
'gene_id':'/dgene/gene/',
'tax_id': 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=',
'drug_id': "/ddrug/drug/",
'chembl':'https://www.ebi.ac.uk/chembl/compound_report_card/',
'drugbank':'https://go.drugbank.com/drugs/',

}

# App and Model
# app_model_list=['apputil_ApplicationUser_', ]

