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
    'taxonomny':    '/dorganism/taxonomy/{VALUE}',
    'tax_id':       'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={VALUE}',
    'organism_id':  '/dorganism/organism/{VALUE}',
    'drug_id':      '/ddrug/drug/{VALUE}',
    'cas':          'https://commonchemistry.cas.org/detail?cas_rn={VALUE}',
    'pubchem':      'https://pubchem.ncbi.nlm.nih.gov/compound/{VALUE}',
    'drugbank' :    'https://www.drugbank.ca/drugs/{VALUE}',
    'chemspider':   'https://www.chemspider.com/Chemical-Structure.{VALUE}.html',
    'unii':         'https://precision.fda.gov/uniisearch/srs/unii/{VALUE}',
    'kegg':         'https://www.kegg.jp/entry/{VALUE}',
    'chebi':        'https://www.ebi.ac.uk/chebi/searchId.do?chebiId={VALUE}',
    'chembl':       'https://www.ebi.ac.uk/chembldb/index.php/compound/inspect/{VALUE}',
    'comptox':      'https://comptox.epa.gov/dashboard/chemical/details/{VALUE}',
    'echa':         'https://echa.europa.eu/substance-information/-/substanceinfo/{VALUE}',
    'gene_id':      '/dgene/gene/',
    'nctc':         'https://www.culturecollections.org.uk/products/bacteria/detail.jsp?refId=NCTC+{NUMVALUE}&collection=nctc',
    'atcc':         'https://www.atcc.org/products/{NUMVALUE}',
    'cdc':          'https://wwwn.cdc.gov/ARIsolateBank/Panel/IsolateDetail?IsolateID={VALUE1}&PanelID={VALUE2}}',
    'ncbi_project': 'https://www.ncbi.nlm.nih.gov/bioproject/?term={VALUE}',
    'ncbi_assembly':'https://www.ncbi.nlm.nih.gov/assembly/{VALUE}',
    'ncbi_nuccore': 'https://www.ncbi.nlm.nih.gov/nuccore/{VALUE}',
}

# App and Model
# app_model_list=['apputil_ApplicationUser_', ]

