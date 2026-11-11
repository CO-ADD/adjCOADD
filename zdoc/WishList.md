
# Wish list App


Add Documents for most Items : Organism, Cell, COADD_Compound, ABase_Compound, Biologic_Isoform 
 with document_type 

Split app_dictionary into app-based dictionary, dorganism, dbiologics, 
Add access permission for different app-groups : OrgDB, CastDB, MoleculeDB, BiologicDB 

On generating a new Organism, automatically generate a [Std] Assay entry

Calculated Fields in Extra Status Line

Masterplate card view in ScreenRun_DetailView

# To DO Data


# To Fix

Reduce Process_Progress init to (self, request, progress_cb=None, progress_key=None, progress_timeout=None)

def __init__(self, request, dirname=None, file_list=None, appuser=None,
                 progress_cb=None, progress_key=None, progress_timeout=None):

Create Screen Run revisit the create model, and on finish could the List be refreshed

Project detail view - display nC, nS, nB 
Compound list - display nB 