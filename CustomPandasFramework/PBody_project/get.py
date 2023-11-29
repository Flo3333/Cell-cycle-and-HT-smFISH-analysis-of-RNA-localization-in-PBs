"""
Submodule used to get index from DF.
"""
import pandas as pd
import re
import numpy as np


def get_GenesWithlessThanXCells(Df : pd.DataFrame, CellNumber: int, CellKey = 'CellId') :
    """
    Returns list of genes (str) which has less than 'CellNumber' cells.
    """

    if not isinstance(CellNumber, int) : raise TypeError("CellNumber should be an int, it is a {0}".format(type(CellNumber)))

    if 'rna name' in Df.columns :
        count = Df.groupby('rna name')[CellKey].count()
        

    elif 'rna name' in Df.index.names and CellKey in Df.index.names :
        gene_level = Df.index.names.index('rna name')
        group = Df.groupby("rna name", axis= 0, level= gene_level)
        count = group.count()
        if isinstance(count, pd.DataFrame) :
            count = count.iloc[:,0]
        else :
            assert isinstance(count, pd.Series)

    else : raise KeyError("Could find an axis to group on 'rna name'.")
    res = list(count[count < CellNumber].index)
    return res


def get_pbody_spot_distance_parameters(Pbody : pd.DataFrame)-> list :
    """
    Return list of different distance parameters used to compute spot number in nucleus.
    """
    res = [int(re.findall('(\d+)nm', str(col))[0]) for col in Pbody.filter(regex= ".*\d*nm.*")]
    res = list(np.unique(res))
    return res

def get_minCell_columns() :
    return [
        'id', 'AcquisitionId','nb_rna_out_nuc',
       'nb_rna_in_nuc','cell_area', 'nuc_area','nucleus_mip_mean_signal', 'nucleus_mip_max_signal',
       'nucleus_mip_min_signal', 'nucleus_mip_median_signal',
       'nucleus_mean_mean_signal', 'nucleus_mean_max_signal',
       'nucleus_mean_min_signal', 'nucleus_mean_median_signal',
       'nucleus area (px)', 'nucleus area (nm^2)', 'malat1 spots in nucleus',
       'malat1 spots in cytoplasm','pbody number', 'count pbody in nucleus', 'count pbody in cytoplasm','rna count', 'malat1 count'
    ]

def get_minSpots_columns():
    return [
        'id', 'AcquisitionId','spots_type','InNucleus', 'PbodyId', 'CellId'
    ]

def get_minPbody_columns():
    return [
        'id', 'AcquisitionId','InNucleus', 'rna 0nm count', 'malat1 0nm count',
       'rna 100nm count', 'malat1 100nm count', 'rna 200nm count',
       'malat1 200nm count', 'rna 400nm count', 'malat1 400nm count',
       'rna 600nm count', 'malat1 600nm count', 'rna 800nm count',
       'malat1 800nm count', 'rna 1000nm count', 'malat1 1000nm count',
       'rna 1500nm count', 'malat1 1500nm count', 'rna 2000nm count',
       'malat1 2000nm count', 'CellId'
    ]