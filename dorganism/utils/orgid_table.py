import pandas as pd
from django.apps import apps
from ddrug.models import VITEK_Card, VITEK_ID 
from dgene.models import WGS_CheckM, ID_Sequence 

def get_org_identification_summary(OrgID):

    # Source            Vitek          WGS                 WGS
    # Method           <card_code>     Kraken             CheckM
    # Identification    <id_organism>   <kraken_organism>   

    orgID = []
    nIdx = 1

    vitekID = VITEK_ID.objects.filter(card_barcode__orgbatch_id__organism_id=OrgID, astatus__gte=0)
    for v in vitekID:
        aDict = {}
        aDict['Batch ID'] = v.card_barcode.orgbatch_id.batch_id
        aDict['Identification'] = f"{v.id_organism} ({v.id_probability})"
        aDict['Method'] = f"Vitek {v.card_barcode.card_code}"
        orgID.append(aDict)
        nIdx += 1

    seqID = ID_Sequence.objects.filter(seq_id__orgbatch_id__organism_id=OrgID, astatus__gte=0)
    for s in seqID:
        aDict = {}
        aDict['Batch ID'] = s.seq_id.orgbatch_id.batch_id
        aDict['Identification'] = '; '.join(s.kraken_organisms)
        aDict['Method'] = "WGS Kraken"
        orgID.append(aDict)
        aDict = {}
        aDict['Batch ID'] = s.seq_id.orgbatch_id.batch_id
        aDict['Identification'] = f"{s.gtdbtk_class} - {s.gtdbtk_fastani}"
        aDict['Method'] = "WGS GTDB-TK"
        orgID.append(aDict)

        if s.mlst_scheme != '-':
            aDict = {}
            aDict['Batch ID'] = s.seq_id.orgbatch_id.batch_id
            aDict['Identification'] = f"{s.mlst_scheme} (MLST: {s.mlst_seqtype})"
            aDict['Method'] = "WGS MLST"
            orgID.append(aDict)
            nIdx += 1

    checkID = WGS_CheckM.objects.filter(seq_id__orgbatch_id__organism_id=OrgID, astatus__gte=0)
    for c in checkID:
        aDict = {}
        aDict['Batch ID'] = c.seq_id.orgbatch_id.batch_id
        aDict['Identification'] = f"{c.marker_lineage} (completness: {c.completeness}% contamination: {c.contamination}%)"
        aDict['Method'] = f"WGS CheckM {c.assembly}"
        orgID.append(aDict)


    return pd.DataFrame(orgID)