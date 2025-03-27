import sys, os
import datetime
import numpy as np
import pandas as pd
import logging

#-----------------------------------------------------------------------------
# Scoring Functions
#-----------------------------------------------------------------------------
ActScore_Cutoff = {
    'Invalid':  {'Score':-1,'Code':'-','Desc':'wrong Data'},
    'Inactive': {'Score':0, 'Code':'I','Desc':'>32/20 & DMax<50%'},
    'Partial':  {'Score':1, 'Code':'P','Desc':'>32/20 & DMax>=50%'},
    'LowActive':{'Score':2, 'Code':'L','Desc':'=32/20'},
    'Active':   {'Score':3, 'Code':'A','Desc':'<32/20', 'CutOff':{'uM':20,'ug/mL':32,'pct':50}},
    'Hit':      {'Score':4, 'Code':'H','Desc':'<16/10', 'CutOff':{'uM':10,'ug/mL':16,'pct':20}},
    'SuperHit': {'Score':5, 'Code':'S','Desc':'<2/1',   'CutOff':{'uM':1,'ug/mL':2,'pct':10}}
}


