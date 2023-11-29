"""
This submodule is meant to load results data, in particular from different plate and put it together. 
"""

import pandas as pd
import numpy as np
import CustomPandasFramework.PBody_project.update as update
import CustomPandasFramework.integrity as integrity
from .get import get_minCell_columns, get_minPbody_columns, get_minSpots_columns


def load_results(input_path_list:list ,plate_name_list: list, load_Acquisitions= True, load_Cells= True, load_Pbodys = True, load_Spots = True, display_steps = False, load_min_columns= True)-> dict :
    """
    Load results tables data from different plates. Add a new key to each table : 'plate_name'.
    Every couple of (plate_name; id) is unique. On the contrary beware that usual 'id' are not unique across different plates.
    /!\ load_min_columns (default : True) reduces number of informations loaded to have better memory performance.

    Returns
    -------
        res : dict
            dictionary with concatenated DataFrames, keys are 'Acquisition', 'Cell', 'Spots', 'Pbody'.
    """
    if len(input_path_list) != len(plate_name_list) : raise ValueError("Length of input_path_list and plate_name_list parameters must match.")
    
    Acquisition = pd.DataFrame()
    Cell = pd.DataFrame()
    Pbody = pd.DataFrame()
    Spots = pd.DataFrame()
    
    for plate_name, path in zip(plate_name_list, input_path_list) :

        if display_steps : print("Loading plate {0}".format(plate_name))
        
        if load_Acquisitions :
            loaded_Acquisition = pd.read_feather(path + 'Acquisition')
            loaded_Acquisition.loc[:,["plate_name"]] = plate_name
            loaded_Acquisition = make_InterPlates_id(loaded_Acquisition, ['id'])
            Acquisition = pd.concat([Acquisition, loaded_Acquisition], axis=0)

        if load_Cells :
            if load_min_columns : loaded_Cell = pd.read_feather(path + 'Cell', columns= get_minCell_columns())
            else : loaded_Cell = pd.read_feather(path + 'Cell')
            loaded_Cell.loc[:,["plate_name"]] = plate_name
            loaded_Cell = make_InterPlates_id(loaded_Cell, ['id','AcquisitionId'])
            Cell = pd.concat([Cell, loaded_Cell], axis=0)
            
        if load_Pbodys :
            if load_min_columns : loaded_Pbody = pd.read_feather(path + 'Pbody', columns= get_minPbody_columns())
            else : loaded_Pbody = pd.read_feather(path + 'Pbody')
            loaded_Pbody.loc[:,["plate_name"]] = plate_name
            loaded_Pbody = make_InterPlates_id(loaded_Pbody, ['id', 'AcquisitionId', 'CellId'])
            Pbody = pd.concat([Pbody, loaded_Pbody], axis=0)

        if load_Spots :
            if load_min_columns : loaded_Spots = pd.read_feather(path + 'Spots', columns= get_minSpots_columns())
            else : loaded_Spots = pd.read_feather(path + 'Spots')
            loaded_Spots.loc[:,["plate_name"]] = plate_name
            loaded_Spots = make_InterPlates_id(loaded_Spots, ['id', 'AcquisitionId', 'CellId', 'PbodyId'])
            Spots = pd.concat([Spots, loaded_Spots], axis=0)

    res = {
        'Acquisition' : Acquisition,
        'Cell' : Cell,
        'Pbody' : Pbody,
        'Spots' : Spots
    }

    return res



def make_InterPlates_id(Dataframe: pd.DataFrame, id_list:list = ['id']) :
    """
    if an id is NA then the new id is NA as well.
    """
    if 'plate_name' not in Dataframe.columns : raise KeyError("plate_name key wasn't found in DataFrame, can't make interplate id.")

    df = Dataframe.copy()

    for id_name in id_list :
        if id_name not in df.columns : raise KeyError('{0} column from id_list wasnt found in dataframe.'.format(id_name))
        na_idx = df.query("{0}.isna()".format(id_name)).index
        df[id_name] = list(zip(df['plate_name'], df[id_name]))
        df.loc[na_idx,id_name] = np.NaN
    
    return df



def add_rna_names(DataFrames:dict) :

    """
    Add 'rna name' key to all Frames inside DataFrames, must contain Acquisition DataFrame.
    """

    if DataFrames['Acquisition'].empty : raise KeyError("Acquisition table is empty : cannot add rna names.")

    for table in ['Cell', 'Pbody', 'Spots'] :
        DataFrame = DataFrames[table].copy()
        if DataFrame.empty : continue
        if 'rna name' in DataFrame : continue

        DataFrames[table] = update.AddRnaName(DataFrame, DataFrames['Acquisition'])

    return DataFrames