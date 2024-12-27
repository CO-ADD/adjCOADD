#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import configargparse
from pathlib import Path

from tqdm import tqdm
# from zUtils import zData

from zSql import zSqlConnector
from zUtils import zData
#from zUtils.zPlates import Plate, conv_well_multi_to_list, conv_well_list_to_multi
from zBio.zBioUtils import CmpdSep,SampleSep,BatchSep

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "getPublic"
#logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)
#-----------------------------------------------------------------------------


# ------------------------------------------------------------------------------
def openCoaddDB(User='coadd', Passwd='MtMaroon23',DataBase="coadd",verbose=1):
    dbPG = zSqlConnector.PostgreSQL()
    dbPG.open(User,Passwd,"imb-coadd-db.imb.uq.edu.au",DataBase,verbose=verbose)
    return(dbPG)

def openChemDB(User='chemdb', Passwd='chemdb',DataBase="chemdb",verbose=1):
    dbPG = zSqlConnector.PostgreSQL()
    dbPG.open(User,Passwd,"imb-coadd-db.imb.uq.edu.au",DataBase,verbose=verbose)
    return(dbPG)

def openCastDB(User='castdb', Passwd='coaddAdmin',verbose=1):
    db = zSqlConnector.Oracle()
    db.open(User,Passwd,"imb-coadd-db.imb.uq.edu.au","1521","coadb",verbose=verbose)
    return(db)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
def get_Projects(dataset = 'Public', test=0):

    if dataset == 'Public':
        prjSQL = """
            Select project_id, project_type, pub_date, pub_status, pub_name,
                compound_status, data_status, group_id
            From  dsample.project 
            Where project_type in ('CO-ADD', 'WADI') and pub_status = 'Public'
            """

    if test>0:
        prjSQL += f" Fetch First {test} Rows Only "

    pgDB = openCoaddDB(verbose=0)
    _DF = pd.DataFrame(pgDB.get_dict_list(prjSQL))
    logger.info(f"[Projects] {len(_DF)} ")
    pgDB.close()

    return(_DF)

# ------------------------------------------------------------------------------
def get_Compounds(prj_df, dataset = 'Public', test=0):

    cmpSQL = """
        Select c.compound_id, c.compound_code, c.compound_name, c.compound_type, c.compound_source,
            c.ora_compound_id,
            c.std_status, c.std_nfrag, c.std_smiles, c.std_mw, c.std_mf, c.std_issues,
            c.reg_smiles, c.reg_mw, c.reg_mf,
            b.structure_id, b.structure_type, b.full_mw, b.full_mf,
            c.pub_status
        From  dsample.coadd_compound c
         Left Join dsample.cmpbatch b on c.cmpbatch_id = b.cmpbatch_id
        """

    prjid_lst = prj_df['project_id'].to_list()

    pgDB = openCoaddDB(verbose=0)
    cmp_lst = []
    for prjid in tqdm(prjid_lst,desc="Cmpds: "):
        _sql = f"{cmpSQL} Where project_id = '{prjid}'"
        if test>0:
            _sql += f" Fetch First {test} Rows Only "
        cmp_lst = cmp_lst + (pgDB.get_dict_list(_sql))

    _DF = pd.DataFrame(cmp_lst)
    logger.info(f"[Compounds] {len(_DF)} ")
    pgDB.close()

    return(_DF)


# ------------------------------------------------------------------------------
def get_Inhibition_ora(cmp_df, dataset = 'Public', test=0):

    

    inhibSQL = """
        Select tw.compound_id, 
            tp.Readout_ID, tp.Plate_ID TestPlate_ID, tw.Well_ID TestWell_ID,
            tp.AssayType_ID, a.AssayType_Code, a.AssayType_Class, a.Organism, a.Strain,
            tp.Run_ID, tp.Plate_Quality,
            tw.Conc, tw.Conc_Unit,
            tw.Inhibition, tw.MScore, tw.Active, 
            tw.Data_Status
        From  castdb.TestWell tw
            Left Join castdb.TestPlate tp on tw.Plate_ID = tp.Plate_ID
            Left Join castdb.AssayType a on tp.AssayType_ID = a.AssayType_ID
        Where tw.n_compounds = 1
         And (tp.Status > 0  And tw.Status > 0)
         And ((tw.isSample > 0 And tw.isControl = 0)
         And (tw.isNegControl = 0 And tw.isPosControl = 0)
	     And (tw.isSkip is NULL OR tw.isSkip = 0))
         And (tw.Compound_ID <> 'DMSO' And tw.Compound_ID is not NULL)
         And tp.Result_Type = 'Inhibition'
         And tw.Inhibition is not Null
         And a.AssayType_Class in ('Standard','Mutant_Membrane')
        """

    # inhibSQL = """
    #     Select compound1_id compound_id, 
    #         Readout_ID,
    #         AssayType_ID, AssayType_Code, AssayType_Class, Organism, Strain,
    #         Run_ID, Plate_Quality,
    #         Conc1 Conc, Conc1_Unit Conc_unit,
    #         Inhibition, MScore, Active, 
    #         Data_Status
    #     From  castdb.vAssayInhibition
    #     Where n_compounds = 1
    #      And AssayType_Class in ('Standard','Mutant_Membrane')
    #     """

    oraDB = openCastDB(verbose=0)

    cmpid_lst = cmp_df['ora_compound_id'].to_list()
    inhib_lst = []
    for cmpid in tqdm(cmpid_lst,desc="Inhib: "):
        _sql = f"{inhibSQL} And tw.compound_id = '{cmpid}'"
        if test>0:
            _sql += f" Fetch First {test} Rows Only "
        inhib_lst = inhib_lst + (oraDB.get_dict_list(_sql))

    _DF = pd.DataFrame(inhib_lst)
    logger.info(f"[Inhibitions] {len(_DF)} ")

    oraDB.close()
    return(_DF)


# ------------------------------------------------------------------------------
def get_DoseResponse_ora(cmp_df, dataset = 'Public', test=0):

    micSQL = """
        Select dr.compound1_id compound_id,
            tp.Readout_ID, dr.TestPlate_ID, dr.TestWell_ID,
            tp.AssayType_ID, a.AssayType_Code, a.AssayType_Class, a.Organism, a.Strain,
            tp.Run_ID, tp.Plate_Quality,
            'MIC' Result_Type, 
            dr.MIC DR, dr.MIC_Unit DR_Unit,
            dr.Active, dr.ActScore, dr.pScore,
            dr.Data_Quality
        From  castdb.AssayData_MIC dr
            Left Join castdb.TestPlate tp on dr.TestPlate_ID = tp.Plate_ID
            Left Join castdb.AssayType a on tp.AssayType_ID = a.AssayType_ID
        Where dr.n_compounds = 1
         And a.AssayType_Class in ('Standard','Mutant_Membrane')
         And dr.Valid > 0
         And (dr.Data_Quality like 'Valid%' or dr.Data_Quality like 'Retest%')
        """

    cc50SQL = """
        Select dr.compound1_id compound_id,
            tp.Readout_ID, dr.TestPlate_ID, dr.TestWell_ID,
            tp.AssayType_ID, a.AssayType_Code, a.AssayType_Class, a.Organism, a.Strain,
            tp.Run_ID, tp.Plate_Quality,
            'CC50' Result_Type, 
            dr.CC50 DR, dr.CC50_Unit DR_Unit,
            dr.Active, dr.ActScore, dr.pScore,
            dr.Data_Quality
        From  castdb.AssayData_CC50 dr
            Left Join castdb.TestPlate tp on dr.TestPlate_ID = tp.Plate_ID
            Left Join castdb.AssayType a on tp.AssayType_ID = a.AssayType_ID
        Where dr.n_compounds = 1
         And a.AssayType_Class in ('Standard','Mutant_Membrane')
         And dr.Valid > 0
         And (dr.Data_Quality like 'Valid%' or dr.Data_Quality like 'Retest%')
        """

    hc50SQL = """
        Select dr.compound1_id compound_id,
            tp.Readout_ID, dr.TestPlate_ID, dr.TestWell_ID,
            tp.AssayType_ID, a.AssayType_Code, a.AssayType_Class, a.Organism, a.Strain,
            tp.Run_ID, tp.Plate_Quality,
            'HC50' Result_Type, 
            dr.HC50 DR, dr.HC50_Unit DR_Unit,
            dr.Active, dr.ActScore, dr.pScore,
            dr.Data_Quality
        From  castdb.AssayData_HC50 dr
            Left Join castdb.TestPlate tp on dr.TestPlate_ID = tp.Plate_ID
            Left Join castdb.AssayType a on tp.AssayType_ID = a.AssayType_ID
        Where dr.n_compounds = 1
         And a.AssayType_Class in ('Standard','Mutant_Membrane')
         And dr.Valid > 0
         And (dr.Data_Quality like 'Valid%' or dr.Data_Quality like 'Retest%')
        """

    oraDB = openCastDB(verbose=0)

    cmpid_lst = cmp_df['ora_compound_id'].to_list()
    dr_lst = []
    for cmpid in tqdm(cmpid_lst,desc="DoseResp: "):
        # MIC
        _sql = f"{micSQL} And compound1_id = '{cmpid}'"
        if test>0:
            _sql += f" Fetch First {test} Rows Only "
        dr_lst = dr_lst + (oraDB.get_dict_list(_sql))

        # CC50
        _sql = f"{cc50SQL} And compound1_id = '{cmpid}'"
        if test>0:
            _sql += f" Fetch First {test} Rows Only "
        dr_lst = dr_lst + (oraDB.get_dict_list(_sql))

        # HC50
        _sql = f"{hc50SQL} And compound1_id = '{cmpid}'"
        if test>0:
            _sql += f" Fetch First {test} Rows Only "
        dr_lst = dr_lst + (oraDB.get_dict_list(_sql))

    _DF = pd.DataFrame(dr_lst)
    logger.info(f"[DoseResponse] {len(_DF)} ")

    oraDB.close()
    return(_DF)


# ===============================================================================
def main(prgArgs):
    
    OutPut = ['CSV']

    if prgArgs.directory is None:
        PubDir = 'C:/Data/COADD/DataSet'
    else:
        PubDir = prgArgs.directory
    if not os.path.exists(PubDir):
        os.makedirs(PubDir)

    # Projects
    PrjDF = get_Projects(dataset = 'Public',test=int(prgArgs.test))

    if 'Excel' in OutPut or prgArgs.to_excel:
        xlFile = "COADD_Public_Projects.xlsx"
        logger.info(f"[Projects] {len(PrjDF)} -> {xlFile}")
        PrjDF.to_excel(os.path.join(PubDir,xlFile),sheet_name='Projects')
    if 'CSV' in OutPut:
        csvFile = "COADD_Public_Projects.csv.gz"
        logger.info(f"[Projects] {len(PrjDF)} -> {csvFile}")
        PrjDF.to_csv(os.path.join(PubDir,csvFile),index=False,compression='gzip')

 
    # Compounds
    CmpDF = get_Compounds(PrjDF,dataset = 'Public',test=int(prgArgs.test))

    if 'Excel' in OutPut or prgArgs.to_excel:
        xlFile = "COADD_Public_Compounds.xlsx"
        logger.info(f"[Compounds] {len(CmpDF)} -> {xlFile}")
        CmpDF.to_excel(os.path.join(PubDir,xlFile),sheet_name='Compounds')
    if 'CSV' in OutPut:
        csvFile = "COADD_Public_Compounds.csv.gz"
        logger.info(f"[Compounds] {len(CmpDF)} -> {csvFile}")
        CmpDF.to_csv(os.path.join(PubDir,csvFile),index=False,compression='gzip')


     # Activity - DoseResponse
    DRDF = get_DoseResponse_ora(CmpDF,dataset = 'Public',test=int(prgArgs.test))

    if 'Excel' in OutPut or prgArgs.to_excel:
        xlFile = "COADD_Public_DoseResponse.xlsx"
        logger.info(f"[DoseResponse] {len(CmpDF)} -> {xlFile}")
        DRDF.to_excel(os.path.join(PubDir,xlFile),sheet_name='DoseResponse')
    if 'CSV' in OutPut:
        csvFile = "COADD_Public_DoseResponse.csv.gz"
        logger.info(f"[DoseResponse] {len(CmpDF)} -> {csvFile}")
        DRDF.to_csv(os.path.join(PubDir,csvFile),index=False,compression='gzip')

   # Activity - Inhibition
    InhibDF = get_Inhibition_ora(CmpDF,dataset = 'Public',test=int(prgArgs.test))

    if 'Excel' in OutPut or prgArgs.to_excel:
        xlFile = "COADD_Public_Inhibitions.xlsx"
        logger.info(f"[Inhibition] {len(CmpDF)} -> {xlFile}")
        InhibDF.to_excel(os.path.join(PubDir,xlFile),sheet_name='Inhibition')
    if 'CSV' in OutPut:
        csvFile = "COADD_Public_Inhibitions.csv.gz"
        logger.info(f"[Inhibition] {len(CmpDF)} -> {csvFile}")
        InhibDF.to_csv(os.path.join(PubDir,csvFile),index=False,compression='gzip')



# ------------------------------------------------------------------------------
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")

    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='get Public Dataset', 
                                description="Get Public CO-ADD Dataset")
#    prgParser.add_argument("-t","--table",default=None,required=True, dest="table", action='store', help="Table to upload [ProjectID]")
    # prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    # prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
#    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")

    prgParser.add_argument("--excel",default=False,required=False, dest="to_excel", action='store_true', help="Save to Excel")
    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
#    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")

    # prgParser.add_argument("--django",default='Local',required=False, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    # prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    flRun = True
    try:
        prgArgs = prgParser.parse_args()
    except:
        prgParser.print_help()
        flRun = False

    if flRun:
        main(prgArgs)
        print("-------------------------------------------------------------------")
 