#

import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse

from PIL import Image
# from zUtils import zData

#import django
#from djCOADD import djOrgDB
# from oraCastDB import oraCastDB
#-----------------------------------------------------------------------------

import logging
#-----------------------------------------------------------------------------

# Logger ----------------------------------------------------------------
logTime= datetime.datetime.now()
logName = "Process_BMP"
#logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------------
def reformat_OrganismID(OrgID):
#-----------------------------------------------------------------------------------
    xStr = OrgID.split("_")
    return(f"{xStr[0]}_{int(xStr[1]):04d}")

#-----------------------------------------------------------------------------------
def conv_BMP_JPEG(BmpFolder=None,upload=False,uploaduser=None):
#-----------------------------------------------------------------------------------

    # Get PDF Files
    nBMP = 0
    if BmpFolder:
        DirFiles = os.listdir(BmpFolder)
        BmpFiles = [f for f in DirFiles if f.endswith(".bmp")]
        nBMP = len(BmpFiles)

    if nBMP>0:
        vImages = []
        JpegFolder = os.path.join(BmpFolder,"jpeg")
        if not os.path.exists(JpegFolder):
            os.makedirs(JpegFolder)

        logger.info("-------------------------------------------------------------------------")
        for i in range(nBMP):
            xOrg = reformat_OrganismID(BmpFiles[i].replace(".bmp",""))

            inBmp = os.path.join(BmpFolder,BmpFiles[i])
            outJpeg = os.path.join(JpegFolder,f"{xOrg}.jpeg")

            logger.info(f"[Image    ] {i+1:3d}/{nBMP:3d} - {xOrg} {BmpFiles[i]} {inBmp}")
            im = Image.open(inBmp)
            rgb_im = im.convert("RGB")
            rgb_im.save(outJpeg)
            # update_VitekPDF(PdfFile=PdfFiles[i],VitekFolder=VitekFolder,ProcessedFolder=ProcessedFolder,OrgBatchID=None,
            #                 upload=upload,appuser=appuser)
    else:
        logger.info(f"[Image    ] NO PDF to process in {BmpFolder}  ")


#-----------------------------------------------------------------------------
def main():

    # ArgParser -------------------------------------------------------------
    prgParser = argparse.ArgumentParser(prog='upload_OrgDB_Data', 
                                description="Uploading data to adjCOADD from Oracle or Excel")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [User]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("-f","--bmpfolder",default=None,required=False, dest="bmpfolder", action='store', help="BMP Folder to parse")
    prgArgs = prgParser.parse_args()

    if prgArgs.table == 'Image':
        conv_BMP_JPEG(prgArgs.bmpfolder)

#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")
    main()
    print("...................................................................")

#==============================================================================
