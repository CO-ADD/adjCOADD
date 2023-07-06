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
    'taxonomny':    '/dorganism/taxonomy/',
    'tax_id':       'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=',
    'organism_id':  '/dorganism/organism/',
    'drug_id':      '/ddrug/drug/',
    'cas':          'https://commonchemistry.cas.org/detail?cas_rn=',
    'pubchem':      'https://pubchem.ncbi.nlm.nih.gov/compound/',
    'drugbank' :    'https://www.drugbank.ca/drugs/',
    'chemspider':   'https://www.chemspider.com/Chemical-Structure.{VALUE}.html',
    'unii':         'https://precision.fda.gov/uniisearch/srs/unii/',
    'kegg':         'https://www.kegg.jp/entry/',
    'chebi':        'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=',
    'chembl':       'https://www.ebi.ac.uk/chembldb/index.php/compound/inspect/',
    'comptox':      'https://comptox.epa.gov/dashboard/chemical/details/',
    'echa':         'https://echa.europa.eu/substance-information/-/substanceinfo/',
    'gene_id':      '/dgene/gene/',

}

# App and Model
# app_model_list=['apputil_ApplicationUser_', ]

