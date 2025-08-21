import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import configargparse
from functools import reduce
from pathlib import Path

from tqdm import tqdm
# from zUtils import zData

import matplotlib.ticker as tic
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import django

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Sum_CmpBatch_DoseResp"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------

def main(prgArgs,djDir):

    django.setup()

    from dplate.models import Labware, TestPlate, TestWell
    from dsample.models import COADD_Compound, Compound_Batch
    from ddrug.models import Drug
    from dsummary.utils.summary_data import sum_cmpbatch_dr
    from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run, Assay
    from dsummary.models import Summary_CmpBatch, Summary_CmpBatch_Doseresp
    from adjcoadd.constants import COMPOUND_SEP
    #from django_pandas.io import read_frame

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir['djPrj']}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")
    

    

   # AssayData MIC -------------------------------------------------------------

    def plotDR(DrData,PlotTitle,PlotFile):
        pltColors = [mcolors.TABLEAU_COLORS[c] for c in mcolors.TABLEAU_COLORS]

        n_data = len(DrData)

        TW_VALUES = ['plate_id','well_id','cmpbatch_id','conc_lst','conc_unit_lst','readouts', 'inhibition']

        #subTitle = f"{self.run_id} - {self.tp_id}.{_nset:02d}"
        fig, ax = plt.subplots(figsize=(12,8))
        fig.text(0.05,0.91,PlotTitle, fontsize=19, ha = 'left')
        #fig.text(0.97,0.91,subTitle,fontsize=12, color = 'darkgrey', ha = 'right')

        nplot=0
        for _DR in DrData:
            _color = pltColors[nplot]
            qryTW = TestWell.objects.filter(cmpbatch_id__cmpbatch_id__contains = _DR['cmpbatch_id'],
                                        plate_id__assay_id__assay_id__contains= _DR['assay_id'],
                                        plate_id =_DR['testplate_id']).values(*TW_VALUES)
            _data = []
            for _tw in list(qryTW):
                _data.append({'cmpbatch_id':_tw['cmpbatch_id'], 
                                     'conc':_tw['conc_lst'][0], 'conc_unit':_tw['conc_unit_lst'][0],
                                     'inhibition':_tw['inhibition'],
                                     'readout':_tw['readouts'][0]}
                                     )
            _df = pd.DataFrame(_data)

            labText = f"{_DR['cmpbatch_id']} {_DR['assay_id']}\n {_DR['run_id']} : {_DR['testplate_id']}\n MIC: {_DR['mic']} {_DR['mic_unit']} ({_DR['inhibit_max']:.1f}%) \n IC50: {_DR['ic50']}\n QC: {_DR['data_quality']} {_DR['testplate_id__plate_quality']}\n"
            sns.lineplot(data=_df, ax=ax, x='conc',y='inhibition',label=labText,marker='o',linestyle='dotted',color=_color,markersize=10)
            nplot += 1
        
        ax.set(xscale="log")
        ax.set_ylim(-20, 120)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.80, box.height])
        ax.xaxis.set_major_formatter(tic.FormatStrFormatter('%g'))
        ax.yaxis.set_major_formatter(tic.FormatStrFormatter('%g'))
        plt.legend(bbox_to_anchor=(1.01, 1),loc=2, borderaxespad=0)
        plt.xlabel(f"Concentration ", fontsize= 12)
        plt.ylabel(f"Inhibition [%]", fontsize= 12)
        print(f" [DR-Plot] {PlotFile}.jpeg")
        fig.savefig(f"{PlotFile}.jpeg")


    compound_lst = ['MCC_000636','MCC_000094']

    if prgArgs.csvfile:
        SEQ = pd.read_csv(prgArgs.csvfile)
        SEQ.columns = [x.upper() for x in SEQ.columns]
        for idx,row in SEQ.iterrows():
            _orgid = "_".join(row['ORGBATCH_ID'].split('_')[0:2])
            
            if Assay.objects.filter(organism_id=_orgid).exists():
                djAssay = Assay.objects.get(organism_id=_orgid)
                if djAssay:
                    _pub_id = djAssay.organism_id.pub_id
                else:
                    _pub_id = '-'

                # Get MIC Data
                MIC_VALUES = ['testplate_id','testwell_id','assay_id','run_id','cmpbatch_id','mic','mic_unit', 'inhibit_max','data_quality','ic50','ic50_quality','testplate_id__plate_quality']
                lstMIC = []
                for _cmpid in compound_lst:

                    djDrug = Drug.objects.get(uq_imb=_cmpid)
                    if djDrug:
                        _drug_name = djDrug.drug_name
                    else:
                        _drug_name = '-'

                    qryMIC = AssayData_MIC.objects.filter(cmpbatch_id__cmpbatch_id__contains = _cmpid,
                                                    assay_id__assay_id__contains= _orgid).values(*MIC_VALUES)
                    if qryMIC.count() > 0:
                        _lst = list(qryMIC)
                        for _l in _lst:
                            _l['drug_name'] = _drug_name
                        lstMIC += _lst

                if len(lstMIC) > 0:
                    plotDR(lstMIC,f"{prgArgs.orgid} ({_pub_id}) - COL (094) PmxB (636) ",f"{_pub_id}_PmxB_Col")
                else:
                    print(f" [DR-Plot] {_orgid} NO MIC Data") 
            else:
                print(f" [DR-Plot] {_orgid} NO Assay Data") 

    elif prgArgs.orgid:

        djAssay = Assay.objects.get(organism_id=prgArgs.orgid)
        if djAssay:
            _pub_id = djAssay.organism_id.pub_id
        else:
            _pub_id = '-'

        # Get MIC Data
        MIC_VALUES = ['testplate_id','testwell_id','assay_id','run_id','cmpbatch_id','mic','mic_unit', 'inhibit_max','data_quality','ic50','ic50_quality','testplate_id__plate_quality']
        lstMIC = []
        for _cmpid in compound_lst:

            djDrug = Drug.objects.get(uq_imb=_cmpid)
            if djDrug:
                _drug_name = djDrug.drug_name
            else:
                _drug_name = '-'

            qryMIC = AssayData_MIC.objects.filter(cmpbatch_id__cmpbatch_id__contains = _cmpid,
                                            assay_id__assay_id__contains= prgArgs.orgid).values(*MIC_VALUES)
            if qryMIC.count() > 0:
                _lst = list(qryMIC)
                for _l in _lst:
                    _l['drug_name'] = _drug_name
                lstMIC += _lst

        if len(lstMIC) > 0:
            plotDR(lstMIC,f"{prgArgs.orgid} ({_pub_id}) - COL (094) PmxB (636) ",f"{_pub_id}_PmxB_Col")
        

#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='plot_DR', 
                                description="Plot Doseresponse data")
    prgParser.add_argument("-t",default='Plot',required=False, dest="table", action='store', help="Table to upload [Plot]")
    # prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    # prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    # prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    # prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
    # prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

    # prgParser.add_argument("-c","--cid",default=None,required=False, dest="compoundid", action='store', help="CompoundID")
    # prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="RunID")
    prgParser.add_argument("-o","--orgid",default=None,required=False, dest="orgid", action='store', help="OrganismID")
    prgParser.add_argument("-c","--csvfile", default=None,required=False, dest="csvfile", action='store', help="CSVFile")
    prgParser.add_argument("--plate",default=None,required=False, dest="plateid", action='store', help="TestPlate")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")

    prgParser.add_argument("--django",default='Local',required=False, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    #prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    prgArgs = prgParser.parse_args()

    from zDjango.djUtils import init_django_dir

    # Django -------------------------------------------------------------
    djDir = init_django_dir(prgArgs,"adjCOADD")
    if djDir:
        print(djDir)
        main(prgArgs,djDir)

#==============================================================================