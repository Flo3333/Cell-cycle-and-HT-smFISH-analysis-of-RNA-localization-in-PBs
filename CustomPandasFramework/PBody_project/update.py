import pandas as pd
import warnings
import os
import numpy as np
from pbwrap.integrity import check_parameter
from ..operations import foreign_key
from ..integrity import check_expectedcolumns, MissingColumnError
from ..utils import is_contained, hist_maximum
from ..utils import is_contained

"""module aiming to update values saved in pbody project dataframes"""

def JoinCellAcquisition(Acquisition: pd.DataFrame, Cell: pd.DataFrame, Cell_columns: list = None, Acquisition_columns: list = None) :
    """
    Perform left join Cell - Acquisition. If 'Cell_columns' or 'Acquisiton_columns' is None all columns from corresponding DataFrame are returned. 
    """
    
    if Cell_columns == None : Cell_columns = list(Cell.columns)
    else : check_expectedcolumns(Cell, Cell_columns)
    if "AcquisitionId" in Cell_columns : Cell_columns.remove("AcquisitionId")
    if Acquisition_columns == None : Acquisition_columns = list(Acquisition.columns)
    else : check_expectedcolumns(Acquisition, Acquisition_columns)
    if "id" in Acquisition_columns : Acquisition_columns.remove("id")

    JoinFrame = pd.merge(left= Cell.loc[:, ["AcquisitionId"] + Cell_columns], right= Acquisition.loc[:,["id"] + Acquisition_columns], how= 'inner', left_on= 'AcquisitionId', right_on= 'id')
    JoinFrame = JoinFrame.drop(axis=1, labels= ["id_y"]).rename(columns={"id_x": "id"})
    return JoinFrame

def AddSpotsCount(Cell: pd.DataFrame, Spots: pd.DataFrame) :
    """
    """
    malat_idx = Spots.query("spots_type == 'malat1'").index
    rna_idx = Spots.query("spots_type == 'rna'").index
    malat_group = Spots.loc[malat_idx,["CellId","spots_coords"]].groupby(["CellId"]).count().reset_index().rename(columns= {"spots_coords" : "malat_count"})
    rna_group = Spots.loc[rna_idx,["CellId","spots_coords"]].groupby(["CellId"]).count().reset_index().rename(columns= {"spots_coords" : "rna_count"})
    
    spot_count = pd.merge(rna_group, malat_group, 'outer', on= "CellId")
    JoinFrame = pd.merge(Cell, spot_count, 'left', left_on="id", right_on="CellId")

    return JoinFrame


def remove_acquisitions(Acquisition: pd.DataFrame, Cell: pd.DataFrame, AcquisitionIds: 'list[int]') -> 'list[pd.DataFrame]' :
    """Remove acquisition from Acquisition and Cell tables.

    Parameters
    ----------
        Acquisition : pd.DataFrame
        Cell : pd.DataFrame
        AcquisitionId : list[int] or pd.Index
            ids of acquisition that are to be removed.
            
    Returns
    -------
        new_frames : list[pd.DataFrame]
            [new_Acquisition, new_Cell] : DataFrame with acquisition removed
    """

    Acquisition_drop_index = Acquisition.query("id in {0}".format(list(AcquisitionIds))).index
    Cell_drop_index = Cell.query("AcquisitionId in {0}".format(list(AcquisitionIds))).index

    new_Acquisition = Acquisition.drop(axis= 0, index= Acquisition_drop_index)
    new_Cell = Cell.drop(axis= 0, index= Cell_drop_index)
    res = [new_Acquisition, new_Cell]

    return res




def from_detectionthreshold_remove_acquisition(Acquisition: pd.DataFrame, Cell: pd.DataFrame, limit: float, keep= 'greater') -> 'list[pd.DataFrame]' :
    """
    Parameters
    ----------
        limit : float
            Will only keep acquisition with (strictly) greater or lower detection threshold than limit.
        keep : str
            'greater' or 'lower'
    """
    if keep == 'greater' : comparison = '<=' #comparison is reversed on purpose since it used to fetch acquisitions NOT to keep.
    elif keep == 'lower' : comparison = '>='
    else : raise ValueError("keep should be either 'greater' or 'lower'.")

    drop_id = Acquisition.query('`RNA spot threshold` {0} {1}'.format(comparison, limit))["id"]

    new_Acquisition, new_Cell = remove_acquisitions(Acquisition, Cell, AcquisitionIds= drop_id)

    return new_Acquisition, new_Cell




def from_Input_remove_damagedacquisitions(Input: pd.DataFrame, path= 'default', print_deletion_number= False, oligopool = 10):
    """
    Remove from Input DataFrame acquisitions which are reported as damaged in the 'Oligopool_annotation.xlsx' file.
    Make sure the computer has access to this location, in other words that the server is mounted.
    
    File location :
    ---------------
        Oligopool 8 
            "/run/user/1001/gvfs/smb-share:server=archive.igh.internal,share=bertrand/Commun-Bertrand/Oriane/Floric images pbodies/Oligopool8_annotation.xlsx"
        Oligopool 8.1 
            "/run/user/1001/gvfs/smb-share:server=archive.igh.internal,share=bertrand/Commun-Bertrand/Oriane/Floric images pbodies/Images new dapi/Oligopool8_newdapi_annotation.xlsx"
        Oligopool 10
            "/run/user/1001/gvfs/smb-share:server=archive.igh.internal,share=bertrand/Commun-Bertrand/Oriane/Oligpools/Oligopool 10/Oligopool10_annotation.xlsx"

    """

    if path == 'default' : 
        if oligopool == 8 : path = "/run/user/1001/gvfs/smb-share:server=archive.igh.internal,share=bertrand/Commun-Bertrand/Oriane/Floric images pbodies/Oligopool8_annotation.xlsx"
        elif oligopool == 8.1 : path =  "/run/user/1001/gvfs/smb-share:server=archive.igh.internal,share=bertrand/Commun-Bertrand/Oriane/Floric images pbodies/Images new dapi/Oligopool8_newdapi_annotation.xlsx"
        elif oligopool == 10 : path =  "/run/user/1001/gvfs/smb-share:server=archive.igh.internal,share=bertrand/Commun-Bertrand/Oriane/Oligpools/Oligopool 10/Oligopool10_annotation.xlsx"
    if not os.path.isfile(path) : raise FileNotFoundError("Could not find the annotation file. Please make sure the sure the server is mounted or that the given path is correct.")

    annotation = pd.read_excel(path)
    
    #Detection
    remarques_banned_words = ['damaged', 'no cells', 'pb dapi']
    signal_banned_words = ['no signal']
    Gene_to_ban = []

    for word in remarques_banned_words :
        query = "remarques.str.contains('{0}')".format(word)
        Gene_to_ban.extend(annotation.dropna(axis=0, subset= "remarques").query(query).loc[:,"GeneName"])
    for word in signal_banned_words :
        Gene_to_ban.extend(annotation.dropna(axis=0, subset= "signal quality").query('`signal quality`.str.contains("{0}")'.format(word)).loc[:,"GeneName"])
    
    
    #Drop
    drop_index = Input.query("`rna name` in {0}".format(Gene_to_ban)).index
    new_Input = Input.drop(index= drop_index, axis=0)
    if print_deletion_number : print("{0} files were removed from pipeline according to the annotation file report.".format(len(Input)- len(new_Input)))
    return new_Input




def from_pbodynum_remove_cell(Cell: pd.DataFrame, limit = 0, keep= 'greater') -> pd.DataFrame:
    """keep : 'greater', 'lower', 'equal' or 'different'."""
    if keep == 'greater' : comparison = '<=' #comparison is reversed on purpose since it used to fetch acquisitions NOT to keep.
    elif keep == 'lower' : comparison = '>='
    elif keep == 'equal' : comparison = '!='
    elif keep == 'different' : comparison = '=='
    else : raise ValueError("keep should be one of the following : 'greater', 'lower', 'equal' or 'different'.")

    drop_index = Cell.query("`pbody number` {0} {1}".format(comparison, limit)).index
    new_Cell = Cell.drop(axis=0, index= drop_index)

    return new_Cell



def from_malat_remove_acquisition(Acquisition: pd.DataFrame, Cell: pd.DataFrame, limit= 10) :

    temp_Cell = Cell.copy()
    temp_Cell["total malat spots"] = temp_Cell["malat1 spots in nucleus"] + temp_Cell['malat1 spots in cytoplasm']
    group_by_total = temp_Cell.loc[:,["AcquisitionId", "total malat spots"]].groupby(["AcquisitionId"]).mean()
    drop_index_cell = list(group_by_total[group_by_total["total malat spots"] <= limit].index)

    new_Acquisition, new_Cell = remove_acquisitions(Acquisition, Cell, AcquisitionIds= drop_index_cell)
    return new_Acquisition, new_Cell





def from_spike_value_compute_CellullarCycleGroup(Cell: pd.DataFrame, spikeg1: float= 2.2e7)-> pd.DataFrame :
    """
    Add a new column to Cell DataFrame : "Cellular Cycle Group" which sorts cells between the g1 phase and g2 phase.
    """
    measurement = "nucleus_mean_mean_signal"
    new_Cell = Cell.copy()
    measure = Cell.loc[:, measurement]/spikeg1
    new_Cell["Mean signal / spike"] = measure
    new_Cell.loc[:,new_Cell["Mean signal / spike"] >= 1.5] = 'g2'
    new_Cell.loc[:,new_Cell["Mean signal / spike"] < 1.5] = 'g1'

    return new_Cell



def from_nucleus_malat_proportion_compute_CellullarCycleGroup(Cell: pd.DataFrame, threshold = 0.5)-> pd.DataFrame :

    """
    creates column : 'Cellular_cycle (malat proportion)' ; 'malat proportion in nucleus'
    """
    proportion = Cell.loc[:, "malat1 spots in nucleus"] / (Cell.loc[:, "malat1 spots in nucleus"] + Cell.loc[:, "malat1 spots in cytoplasm"])
    new_Cell = Cell.copy()
    new_Cell["malat proportion in nucleus"] = proportion
    new_Cell["Cellular_cycle (malat proportion)"] = 'g1'
    index = new_Cell.query("`malat proportion in nucleus` > {0}".format(threshold)).index
    new_Cell.loc[index,["Cellular_cycle (malat proportion)"]] = 'g2'

    return new_Cell

def from_IntegratedSignal_ranking_compute_CellularCycleGroup(Cell_raw: pd.DataFrame, g1 = 0.502, s = 0.225, g2 = 0.165, column_name= 'cellular_cycle') :
    """
    Default value are for Hek Cells : g1 = 0.53, s = 0.22, g2 = 0.25
    New key : column_name (default = 'cellular_cycle')
    """
    if column_name in Cell_raw.index : 
        warnings.warn("{0} column already in dataframe, returning dataframe unchanged.".format(column_name))
        return Cell_raw
    
    if g1 >1 or g2>1 or s > 1 : raise ValueError("g1, g2 and s must be defined between 0 and 1")
    if g1 <0 or g2<0 or s < 0 : raise ValueError("g1, g2 and s must be defined between 0 and 1")
    if g1 + g2 + s > 1 : raise ValueError("g1+g2+s should equal 1 or be inferior to 1.")

    Cell = Cell_raw.copy()
    if 'rna name' in Cell.index.names :
        level = Cell.index.names.index("rna name")
        rna_list = Cell.index.get_level_values(level).unique()
    elif 'rna name' in Cell.columns : 
        rna_list = Cell.value_counts(subset= 'rna name').index
    else :    
        raise MissingColumnError("rna name column wasn't found in Cell DataFrame consider using update.AddRnaName")
    
    if "IntegratedSignal" not in Cell.columns :
        Cell = compute_IntegratedSignal(Cell)
    
    for rna in rna_list :
        
        rna_idx = Cell.query("`rna name` == '{0}'".format(rna)).index

        if g1 + g2 + s == 1 :
            Cell.loc[rna_idx,column_name] = "s"
        else : 
            Cell.loc[rna_idx,column_name] = np.NaN
            if s != 0 :
                limit_inf = g1 + (1 - g1 - g2 - s)/2
                limit_sup = 1 - g2 - (1 - g1 - g2 - s)/2
                assert limit_inf < 1 and limit_sup < 1
                assert limit_inf >= 0 and limit_sup >= 0
                quantile_inf = float(Cell.loc[rna_idx,["IntegratedSignal"]].quantile(limit_inf))
                quantile_sup = float(Cell.loc[rna_idx,["IntegratedSignal"]].quantile(limit_sup))
                s_idx = Cell.loc[rna_idx,:].query("IntegratedSignal >= {0} and IntegratedSignal < {1}".format(quantile_inf,quantile_sup)).index
                Cell.loc[s_idx, column_name] = "s"

        g1_quantile = float(Cell.loc[rna_idx,["IntegratedSignal"]].quantile(g1))
        g1_idx = Cell.loc[rna_idx,:].query("IntegratedSignal < {0}".format(g1_quantile)).index
        g2_quantile = float(Cell.loc[rna_idx,["IntegratedSignal"]].quantile(1-g2))
        g2_idx = Cell.loc[rna_idx,:].query("IntegratedSignal >= {0}".format(g2_quantile)).index
        Cell.loc[g1_idx, [column_name]] = "g1"
        Cell.loc[g2_idx, [column_name]] = "g2"

    return Cell


def from_IntegratedSignal_spike_compute_CellularCycleGroup(Cell : pd.DataFrame, column_name= 'cellular_cycle', shift = -0.12e12, auto_shift= True, multiplicator = 2, bins= 100,) :

    """
    Add new column to Cell DF : column_name (default = 'cellular_cycle')
        Values can be either one of ['g1', 'g2' or np.nan]
    """
    if column_name in Cell.index : 
        warnings.warn("{0} column already in dataframe, returning dataframe unchanged.".format(column_name))
        return Cell

    Cell[column_name] = np.NaN
    if "IntegratedSignal" not in Cell.columns :
        Cell.loc[:,["IntegratedSignal"]] = Cell['nucleus_mean_mean_signal'] * Cell["nucleus area (nm^2)"]
    distribution = Cell["IntegratedSignal"]
    median_value = np.median(Cell['IntegratedSignal'])
    maximum_value = hist_maximum(np.histogram(distribution, bins=bins))
    box_size = abs(median_value - maximum_value)

    if auto_shift : shift = (-np.percentile(distribution, 0.5)/2)

    g2_spike = (maximum_value + shift)*multiplicator


    indexg1 = Cell[Cell["IntegratedSignal"] < maximum_value + box_size].index
    indexg2 = Cell[Cell["IntegratedSignal"] > g2_spike].index

    Cell.loc[indexg1,column_name] = 'g1'
    Cell.loc[indexg2,column_name] = 'g2'


    return Cell

def from_IntegratedSignal_distribution_compute_Cellular_CycleGroup(Cell: pd.DataFrame, g1_interval = [7e7,9.5e7], g2_interval= [1.5e8,2e8], drop_IntegratedSignal = False):

    check_parameter(Cell = (pd.DataFrame), g1_interval = (list, tuple), g2_interval = (list, tuple))
    if len(g1_interval) != 2 : raise ValueError("g1_interval : expected length 2 : g1min, g1max")
    if len(g2_interval) != 2 : raise ValueError("g2_interval : expected length 2 : g2min, g2max")

    if "IntegratedSignal" not in Cell.columns : Cell["IntegratedSignal"] = Cell["nucleus_mean_mean_signal"] * Cell["nucleus area (nm^2)"]

    g1min, g1max = g1_interval
    g2min, g2max = g2_interval

    Cell["cellular_cycle"] = np.NaN
    Cell[np.logical_and(Cell["IntegratedSignal"] >= g1min, Cell["IntegratedSignal"] <= g1max)]["cellular_cycle"] = "g1"
    Cell[np.logical_and(Cell["IntegratedSignal"] >= g2min, Cell["IntegratedSignal"] <= g2max)]["cellular_cycle"] = "g2"
    Cell[np.logical_and(Cell["IntegratedSignal"] > g1max, Cell["IntegratedSignal"] < g2min)]["cellular_cycle"] = "?"

    if drop_IntegratedSignal : Cell = Cell.drop(columns= "IntegratedSignal", axis= 1)
    return Cell

def compute_corrected_dapi_signal_measures(Cell, measurement = None) :
    """
    During analysis absolute dapi signal is measured in nucleus as well as in background.
    This function compute new columns for each measure as relative_'measure' = measure(nucleus) - measure(background)
    """

    dapi_signal_measurements = ["nucleus_mean_mean_signal"]
    
    if isinstance(measurement, list) :
        if is_contained(measurement, dapi_signal_measurements) :
            dapi_signal_measurements = measurement
        else : raise ValueError("Currently only {0} mean signal are supported".format(dapi_signal_measurements))
    
    for measurement in dapi_signal_measurements :
        new_col = "relative_" + measurement
        Cell[new_col] = Cell[measurement] - Cell["background dapi signal mean"]
    
    return Cell


def Pbody_AddCellFK(Pbody: pd.DataFrame, Cell: pd.DataFrame, drop_cell_label= True, how= 'inner') :
    """
    Add FK to Pbody_DataFrame referencing Cell_DataFrame.

    This is achieved using a left join operation using AcquisitionId and cell_label columns.
    """
    check_expectedcolumns(Cell, ['AcquisitionId', 'label'])
    check_expectedcolumns(Pbody, ['AcquisitionId', 'cell_label', "id"])

    if 'label' in Pbody.columns :
        Pbody = Pbody.rename(columns= {'label' : 'own_label'})

    res = foreign_key(Pbody, "CellId", Cell, left_key= ['AcquisitionId','cell_label'], right_key= ['AcquisitionId', 'label'], how= how)
    # res = res.reset_index(drop= True)
    if drop_cell_label : res = res.drop(['cell_label', 'label'], axis = 1)
    else : res = res.drop(['label'], axis= 1)
    if 'own_label' in res.columns : res = res.rename(columns= {'own_label' : 'label'})
    return res



def Spots_AddFK(Spots: pd.DataFrame, Cell: pd.DataFrame, Pbody: pd.DataFrame, drop_cell_label= False, drop_Pbody_label= False) :
    """
    Add FK to Cell and Pbody table. Spots in Pbodies will get same cellid as the Pbody, else a join is performed on (AcquisitionId, Cell_label).
    """
    
    Spots = Spots_AddPbodyFK(Spots, Pbody, drop_Pbody_label)
    
    if "CellId" in Spots.columns : Spots = Spots.drop("CellId", axis=1)
    Spots = pd.merge(Spots, Pbody.loc[:,["id","CellId"]], 'left', left_on= "PbodyId", right_on= 'id').drop(["id_y"], axis= 1).rename(columns={"id_x" : "id"})
    ext_idx = Spots.query("CellId.isnull()").index
    int_idx = Spots.query("not CellId.isnull()").index
    Spots_int = Spots.loc[int_idx,:]
    Spots_ext = Spots_AddCellFK(Spots.loc[ext_idx, :], Cell, drop_cell_label)



    Spots_res  = pd.concat([Spots_int, Spots_ext], axis= 0).sort_index().reset_index(drop=True)
    drop_idx = Spots_res.query("PbodyId.isna() and CellId.isna()").index
    Spots_res = Spots_res.drop(drop_idx, axis= 0).reset_index(drop=True)

    return Spots_res




def Spots_AddPbodyFK(Spots: pd.DataFrame, Pbody: pd.DataFrame, drop_pbody_label= True) :
    """
    Add FK to Pbody_DataFrame referencing Cell_DataFrame.

    This is achieved using a left join operation using AcquisitionId and cell_label columns.
    """
    check_expectedcolumns(Spots, ['AcquisitionId', 'Pbody_label'])
    check_expectedcolumns(Pbody, ['AcquisitionId', 'label', "id"])

    if 'label' in Spots.columns :
        Spots = Spots.rename(columns= {'label' : 'own_label'})

    res = foreign_key(Spots, "PbodyId", Pbody, left_key= ['AcquisitionId','Pbody_label'], right_key= ['AcquisitionId', 'label'], validate='many_to_one')
    if drop_pbody_label : res = res.drop(['Pbody_label', 'label'], axis = 1)
    else : res = res.drop(['label'], axis= 1)
    if 'own_label' in res.columns : res = res.rename(columns= {'own_label' : 'label'})
    return res


def Spots_AddCellFK(Spots: pd.DataFrame, Cell: pd.DataFrame, drop_cell_label= True) :
    """
    Add FK to Pbody_DataFrame referencing Cell_DataFrame.

    This is achieved using a left join operation using AcquisitionId and cell_label columns.
    """
    check_expectedcolumns(Spots, ['AcquisitionId', 'cell_label'])
    check_expectedcolumns(Cell, ['AcquisitionId', 'label', "id"])
    if 'label' in Spots.columns :
        Spots = Spots.rename(columns= {'label' : 'own_label'})

    res = foreign_key(Spots, "CellId", Cell, left_key= ['AcquisitionId','cell_label'], right_key= ['AcquisitionId', 'label'], validate=('many_to_one'))
    if drop_cell_label : res = res.drop(['cell_label', 'label'], axis = 1)
    else : res = res.drop(['label'], axis= 1)
    if 'own_label' in res.columns : res = res.rename(columns= {'own_label' : 'label'})
    return res


def AddRnaName(Df,Acquisition) :

    if "AcquisitionId" not in Df.columns : raise MissingColumnError("AcquisitionId was not found in DataFrame, cannot join with Acquisition.")
    if 'rna name' in Df.columns :
        return Df

    check_length = len(Df)
    Df = pd.merge(Df, Acquisition.loc[:,["id", "rna name"]], "inner", left_on='AcquisitionId', right_on= 'id')
    Df = Df.drop("id_y", axis= 1).rename(columns= {"id_x" : "id"})

    assert len(Df) == check_length, "Reference from DataFrame to Acquisition was not complete or not unique, joined frame resulted in a different number of lines."
    
    return Df

def AddCellularCycle(Df,Cell) :
    """
    Add to Df 'cellular_cycle' key from Cell.
    If cellular_cycle is already in Df, Df is not modified.
    """
    if 'cellular_cycle' in Df.columns : return Df
    if not 'cellular_cycle' in Cell.columns : raise KeyError("Cellular_cycle key wasn't computed in Cell DataFrame. Please Compute cellular_cycle in Cell DataFrame first.")
    if 'CellId' not in Df.columns : raise KeyError("Not FK to Cell DataFrame found : missing CellId key in Df argument.")

    check_length = len(Df)
    Df = pd.merge(Df,Cell.loc[:,['id','cellular_cycle']], 'inner', left_on= 'CellId', right_on= 'id', validate= 'many_to_one').drop("id_y", axis= 1).rename(columns={'id_x' : 'id'})
    assert len(Df) <= check_length, "length mismatch, before : {0} ; after : {1}".format(check_length, len(Df))

    return Df



def removeCellsWithoutPbodies(Cell: pd.DataFrame, Pbody: pd.DataFrame) :
    """
    """
    CellIdsWithPbody = list(Pbody.loc[:,'CellId'].unique())
    CellWithPbodyIds = Cell.query("id in {0}".format(CellIdsWithPbody)).index
    new_Cell = Cell.loc[CellWithPbodyIds,:]

    return new_Cell


def removeCellsWithoutSpotsInPodies(Cell: pd.DataFrame, Spots: pd.DataFrame) :
    Trueidx = Spots.query("not PbodyId.isna()").index
    CellIdsTrue = list(Spots.loc[Trueidx,'CellId'].unique())
    CellTrue = Cell.query("id in {0}".format(CellIdsTrue)).index

    new_Cell = Cell.loc[CellTrue,:]
    return new_Cell

def compute_IntegratedSignal(Cell: pd.DataFrame) :
    Cell_copy = Cell.copy()
    Cell_copy.loc[:,["IntegratedSignal"]] = Cell['nucleus_mean_mean_signal'] * Cell["nucleus area (nm^2)"]
    return Cell_copy