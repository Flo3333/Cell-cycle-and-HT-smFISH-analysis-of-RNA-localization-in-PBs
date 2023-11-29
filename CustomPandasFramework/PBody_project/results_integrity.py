"""This submodule contains function to test coherence of Data between results Tables

    So far 4 results tables are computed


                                    Acquisition (Field of View)
                                              
                                                |
                                                | (1,N)
                                                V

                                              Cell

                                                |
                                                | (1,N)
                                                V

                                   Pbody <--- (1,N) ---> Spots
"""


import pandas as pd
from ..integrity.datashape import is_empty, is_primarykey, check_expectedcolumns
from ..utils import check_parameter
from .DataFrames import newframe_Spots, newframe_Acquisitions, newframe_Cell, newframe_Pbody



def Run_Results_Integrity_checks(Acquisition_result: pd.DataFrame= None, Cell_result : pd.DataFrame=None, Pbody_result : pd.DataFrame=None, Spots_result : pd.DataFrame=None, Print_failed_checks= False) :
    """
    Check data coherence between all results tables and saves a log with results at 'log_path'.
    """
    
    if type(Acquisition_result) == type(None) : Acquisition_result = newframe_Acquisitions()
    if type(Cell_result) == type(None) : Cell_result = newframe_Cell()
    if type(Pbody_result) == type(None) : Pbody_result = newframe_Pbody()
    if type(Spots_result) == type(None) : Spots_result = newframe_Spots()
    
    dict_report = {}
    Cell_result["rna count"] = Cell_result['nb_rna_in_nuc'] + Cell_result['nb_rna_out_nuc']
    Cell_result["malat1 count"] = Cell_result['malat1 spots in nucleus'] + Cell_result['malat1 spots in cytoplasm']

    #Emptiness
    dict_report["Acquisition is empty"] = is_empty(Acquisition_result)
    dict_report["Cell is empty"] = is_empty(Cell_result)
    dict_report["Pbody is empty"] = is_empty(Pbody_result)
    dict_report["Spots is empty"] = is_empty(Spots_result)

    #Id is primary key
    dict_report["Acquisition id is valid"] = is_primarykey(Acquisition_result, "id")
    dict_report["Cell id is valid"] = is_primarykey(Cell_result, "id")
    dict_report["Pbody id is valid"] = is_primarykey(Pbody_result, "id")
    dict_report["Spots id is valid"] = is_primarykey(Spots_result, "id")

    #Referencement relation
    dict_report["Acquisition defines (1,N) relation with Cell"] = get_referencement_relation(Acquisition_result, "id", Cell_result, "AcquisitionId") == ('1','N')
    dict_report["Acquisition defines (1,N) relation with Pbody"] = get_referencement_relation(Acquisition_result, "id", Pbody_result, "AcquisitionId") == ('1','N')
    dict_report["Acquisition defines (1,N) relation with Spots"] = get_referencement_relation(Acquisition_result, "id", Spots_result, "AcquisitionId") == ('1','N')
    dict_report["Cell defines (1,N) relation with Pbody"] = get_referencement_relation(Cell_result, "id", Pbody_result, "CellId") == ('1','N')
    dict_report["Cell defines (1,N) relation with Spots"] = get_referencement_relation(Cell_result, "id", Spots_result, "CellId") == ('1','N')
    dict_report["Pbody defines (1,N) relation with Spots"] = get_referencement_relation(Pbody_result, "id", Spots_result, "PbodyId") == ('1','N')

    #Coherence between counts
    dict_report["Acquisition cell number matches Cell table"] = bool_object_count(Cell_result, 'AcquisitionId', Acquisition_result, "cell number",Print_failed_checks=Print_failed_checks)
    dict_report["Cell Pbody number matches Pbody Table"] = bool_object_count(Pbody_result, 'CellId', Cell_result, "pbody number", Print_failed_checks=Print_failed_checks)
    dict_report["Pbody rna count matches Spots Table"] = bool_object_count(Spots_result.query("spots_type == 'rna'"), 'PbodyId', Pbody_result, "rna 0nm count", Print_failed_checks=Print_failed_checks)
    dict_report["Pbody malat1 count matches Spots Table"] = bool_object_count(Spots_result.query("spots_type == 'malat1'"), 'PbodyId', Pbody_result, "malat1 0nm count", Print_failed_checks=Print_failed_checks)

    return dict_report


def bool_object_count(DataFrame: pd.DataFrame, FK: str, check_array: pd.DataFrame, check_key: pd.Series, Print_failed_checks=False) :
    """
    Group DataFrame by 'objectname' and returns true if counts equal values in check_key.

    Parameters
    ----------
        FK : int
            Foreign key to check array

    """
    if DataFrame.empty : return 0
    # check_parameter(DataFrame = pd.DataFrame, groupkey = str, check_series= (pd.Series,pd.DataFrame))
    if 'label' in check_array.columns :
        Truth_table = check_array.loc[:,["id","label", check_key]].sort_values("id").astype(int).rename(columns={"id" : FK, check_key : "count"})
    else : Truth_table = check_array.loc[:,["id", check_key]].sort_values("id").astype(int).rename(columns={"id" : FK, check_key : "count"})
    sum_df = DataFrame.groupby([FK])["id"].count().reset_index(drop= False).sort_values(FK).astype(int).rename(columns={"id" : "count"})
    
    join_frame = pd.merge(sum_df, Truth_table, how='outer', on= FK, indicator=True).query('count_x != 0 and count_y != 0 and count_y != count_x').rename(columns={'count_x' : 'Table aggregate', 'count_y' : 'count from measure'})
    res = join_frame.empty
    if not res and Print_failed_checks : print("Bool object count failed : \n", join_frame)
    return join_frame.empty


def get_referencement_relation(localDataFrame,local_key, ForeignFrame, foreign_key ):
    """Checks that keycolumn 1 from DataFrame1 forms a (order1, order2) referencement relationship with key column2 from DataFrame2.
    Raise Exception if it fails.
    
    Parameters
    ----------
        DataFrame1 : pd.DataFrame
        keycolumn1 : str
        order1 : str
            Either "1", "N" or "n"
        DataFrame2 : pd.DataFrame
        keycolumn2 : str
        order2 ; str
            Either "1", "N" or "n"
    """
    
    #Integrity
    check_parameter(localDataFrame= (pd.DataFrame), local_key= (str), ForeignFrame= (pd.DataFrame), foreign_key= (str))
    check_expectedcolumns(localDataFrame, [local_key])
    check_expectedcolumns(ForeignFrame, [foreign_key])
    if localDataFrame.empty : return 'empty local df'
    if ForeignFrame.empty : return 'empty foreign frame'

    #removing NULL rows
    tempForeignFrame = ForeignFrame.copy()
    tempForeignFrame = tempForeignFrame.dropna(subset= [foreign_key]).loc[:,[foreign_key]]

    tempLocalFrame = localDataFrame.copy()
    tempLocalFrame = localDataFrame.dropna(subset = [local_key]).loc[:, [local_key]]

    #group by keys
    group_fk = tempForeignFrame.value_counts(subset= foreign_key).reset_index().drop(0, axis= 1)
    group_lk = tempLocalFrame.value_counts(subset= local_key).reset_index().drop(0, axis= 1)

    #join tables
    left_join = pd.merge(left= group_lk, right= tempForeignFrame, how= 'left', left_on= local_key, right_on= foreign_key)
    right_join = pd.merge(left= tempLocalFrame, right= group_fk, how= 'right', left_on= local_key, right_on= foreign_key)


    #match testing
    if len(left_join) == len(group_lk) : foreign_order = "1"
    else : foreign_order = "N"

    if len(right_join) == len(group_fk) : local_order = "1"
    else : local_order = "N"

    if any(pd.isna(left_join[local_key])) : raise Warning("All foreign keys don't have a match in the local DataFrame.")
    if any(pd.isna(right_join[foreign_key])) : raise Warning("All local keys don't have a match in the foreign DataFrame.")

    referencement_relation = (local_order, foreign_order)
    return referencement_relation