from django.core.validators import RegexValidator

#
# Application Constants/Settings 
#

#
# Use models.DecimalField(max_digits=x, decimal_places=x) for Floats
#
# Concentrations:  models.DecimalField(max_digits=12, decimal_places=4)
#

#
# Use models.DecimalField(max_digits=x, decimal_places=x) for Floats
#
# Concentrations:  models.DecimalField(max_digits=12, decimal_places=4)
#



# -dOrganism Settings ---------------------------------------------------
ORGANISM_CLASSES = ['GN','GP','MB','FG']
ORGANSIM_SEP = "_"
ORGBATCH_SEP = "_"

# -dCell Settings 
CELL_CLASSES = ['MA']
CELL_SEP = "_"
CELLBATCH_SEP = "_"

# -dPeptide Settings 
PEPTIDE_CLASSES = ['MA']
PEPTIDE_SEP = "_"
PEPBATCH_SEP = "_"

# -dChem Settings ---------------------------------------------------
SAMPLE_SEP = "_"
SAMPLEBATCH_SEP = "_"
COMPOUND_SEP = '|'

# -dScreen Settings ---------------------------------------------------
RUN_CLASSES = ['PSR','HCR','QCR']
RUN_SEP = ""

# -dSample Settings ---------------------------------------------------
PROJECT_COMPOUND_STATUS = ['MissingStructureData','NoStructureData','ToImportStructureData',
                           'NoCompoundData','NoCompoundReceived',]   
PROJECT_SCREEN_STATUS   = ['PS','HC','HV']
PROJECT_DATA_STATUS     = ['NoScreenData','DoNotYet_makePublic','DoNot_makePublic','Limited','Confidential']
PROJECT_REPORT_STATUS   = ['PS','HC','HV']

AlphaNumeric = RegexValidator(r'^[0-9a-zA-Z]*$', 'Only alphanumeric characters are allowed.')

# column name can be edited here 
# make a dictioinary  with Key and value, if value is none choose verbose name else choose the dictionary name.

# Links:
LinkList={
    'taxonomny':        '/dorganism/taxonomy/{VALUE1}',
    "urlname":          '/dorganism/taxonomy/{VALUE1}',
    'organism_id':      '/dorganism/organism/{VALUE1}',
    
    'cell_id':          '/dcell/cell/{VALUE1}',
    
    'drug_id':          '/ddrug/drug/{VALUE1}',
    
    'screenrun_id':     '/dscreen/screenrun/{VALUE1}',
    
    'project_id':       '/dsample/project/{VALUE1}',
    
    'testplate_id':     '/dplate/testplate/{VALUE1}',
    'masterplate_id':   '/dplate/masterplate/{VALUE1}',
    
    'peptide_id':       '/dpeptide/peptide/{VALUE1}',
    
    'organisation_id':  '/dcollab/organisation/{VALUE1}',
    'user_id':          '/dcollab/collaborator/{VALUE1}',
    'group_id':         '/dcollab/organisation/{VALUE1}',
    
    'gene_id':          '/dgene/gene/',
    'seq_id':           '/dgene/seq/',
    
    'tax_id':           'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={VALUE1}',
    'cas':              'https://commonchemistry.cas.org/detail?cas_rn={VALUE1}',
    'pubchem':          'https://pubchem.ncbi.nlm.nih.gov/compound/{VALUE1}',
    'drugbank' :        'https://www.drugbank.ca/drugs/{VALUE1}',
    'chemspider':       'https://www.chemspider.com/Chemical-Structure.{VALUE1}.html',
    'unii':             'https://precision.fda.gov/uniisearch/srs/unii/{VALUE1}',
    'kegg':             'https://www.kegg.jp/entry/{VALUE1}',
    'chebi':            'https://www.ebi.ac.uk/chebi/searchId.do?chebiId={VALUE1}',
    'chembl':           'https://www.ebi.ac.uk/chembldb/index.php/compound/inspect/{VALUE1}',
    'comptox':          'https://comptox.epa.gov/dashboard/chemical/details/{VALUE1}',
    'echa':             'https://echa.europa.eu/substance-information/-/substanceinfo/{VALUE1}',
    'nctc':             'https://www.culturecollections.org.uk/products/bacteria/detail.jsp?refId=NCTC+{NUMVALUE}&collection=nctc',
    'atcc':             'https://www.atcc.org/products/{NUMVALUE}',
    'cdc':              'https://wwwn.cdc.gov/ARIsolateBank/Panel/IsolateDetail?IsolateID={VALUE1}&PanelID={VALUE2}}',
    'ncbi_project':     'https://www.ncbi.nlm.nih.gov/bioproject/?term={VALUE1}',
    'ncbi_assembly':    'https://www.ncbi.nlm.nih.gov/assembly/{VALUE1}',
    'ncbi_nuccore':     'https://www.ncbi.nlm.nih.gov/nuccore/{VALUE1}',
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

