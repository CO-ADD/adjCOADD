#
# Application Constants/Settings 
#

# -dOrganism Settings ---------------------------------------------------
ORGANISM_CLASSES = ['GN','GP','MB','FG']
ORGANSIM_SEP = "_"
ORGBATCH_SEP = "_"

# -dCell Settings 
CELL_CLASSES = ['MA']
CELL_SEP = "_"
CELLBATCH_SEP = "_"

COMPOUND_SEP = '|'
# column name can be edited here 
# make a dictioinary  with Key and value, if value is none choose verbose name else choose the dictionary name.

# Links:
LinkList={
    'taxonomny':    '/dorganism/taxonomy/{VALUE1}',
    "urlname":      '/dorganism/taxonomy/{VALUE1}',
    'tax_id':       'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={VALUE1}',
    'organism_id':  '/dorganism/organism/{VALUE1}',
    'cell_id':  '/dcell/cell/{VALUE1}',
    'drug_id':      '/ddrug/drug/{VALUE1}',
    'cas':          'https://commonchemistry.cas.org/detail?cas_rn={VALUE1}',
    'pubchem':      'https://pubchem.ncbi.nlm.nih.gov/compound/{VALUE1}',
    'drugbank' :    'https://www.drugbank.ca/drugs/{VALUE1}',
    'chemspider':   'https://www.chemspider.com/Chemical-Structure.{VALUE1}.html',
    'unii':         'https://precision.fda.gov/uniisearch/srs/unii/{VALUE1}',
    'kegg':         'https://www.kegg.jp/entry/{VALUE1}',
    'chebi':        'https://www.ebi.ac.uk/chebi/searchId.do?chebiId={VALUE1}',
    'chembl':       'https://www.ebi.ac.uk/chembldb/index.php/compound/inspect/{VALUE1}',
    'comptox':      'https://comptox.epa.gov/dashboard/chemical/details/{VALUE1}',
    'echa':         'https://echa.europa.eu/substance-information/-/substanceinfo/{VALUE1}',
    'gene_id':      '/dgene/gene/',
    'seq_id':       '/dgene/seq/',
    'nctc':         'https://www.culturecollections.org.uk/products/bacteria/detail.jsp?refId=NCTC+{NUMVALUE}&collection=nctc',
    'atcc':         'https://www.atcc.org/products/{NUMVALUE}',
    'cdc':          'https://wwwn.cdc.gov/ARIsolateBank/Panel/IsolateDetail?IsolateID={VALUE1}&PanelID={VALUE2}}',
    'ncbi_project': 'https://www.ncbi.nlm.nih.gov/bioproject/?term={VALUE1}',
    'ncbi_assembly':'https://www.ncbi.nlm.nih.gov/assembly/{VALUE1}',
    'ncbi_nuccore': 'https://www.ncbi.nlm.nih.gov/nuccore/{VALUE1}',
    'seq_id':       '',
}

# Cache Name
# •	Antiogram . Break Point
# •	MIC Pub . BP
# •	MIC Pub . Type
# •	MIC Pub . Source 
# •	Vitek ID . Vitek Process
# •	Vitek ID .  ID Confidence
# •	Vitek Cards . Card Type
# •	Vitek Cards . Card Code 
# •	Vitek AST . BP 
# •	Vitek AST . Source 
# •	Vitek AST . Codes (maybe) Foreignkey field!
# •	Drug . Drug Target
# •	Drug . Drug Class
# •	Drug . Antimicro
# CharToChoice_filterList=["log_code","bp_profile","mic_type", "process", "id_confidence",  "card_type", 
# "card_code", "bp_source", "drug_target", "drug_class", " antimicro"
# ]
CharToChoice_filterList = ["bp_profile",]

