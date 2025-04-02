#
import pandas as pd


def resort_pivtable(df,nLevel=0):
    _code = list(set([c[nLevel] for c in df.columns]))
    _order = df.columns.reindex(_code, level=0)
    return(df.reindex(columns=_order[nLevel]))
