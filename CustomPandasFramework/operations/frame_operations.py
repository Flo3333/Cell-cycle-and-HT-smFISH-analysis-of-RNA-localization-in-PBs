import pandas as pd
import numpy as np
from bigfish.stack import check_parameter
from ..integrity import check_id, has_samedatashape, has_id, is_primary, is_primarykey, has_column, check_expectedcolumns, MergeError

def get_missing_column(dataframe: pd.DataFrame, reference: pd.DataFrame) -> 'list[str]':
    """Returns missing columns from dataframe given reference columns as model."""

    missing_column = []
    if has_samedatashape(DataFrame1= dataframe, DataFrame2= reference) : pass
    else :
        for column in reference.columns :
            if not column in dataframe.columns : missing_column += [column]
    
    return missing_column



def set_missingcolumns_toNA(dataframe: pd.DataFrame, missing_columns : 'list[str]') -> pd.DataFrame :
    """Create missing columns into 'dataframe' and fill values with Null."""

    res_frame = dataframe.copy()
    if type(missing_columns) == pd.DataFrame : missing_columns = get_missing_column(dataframe, missing_columns)

    for column in missing_columns:
        if column in dataframe.columns : raise Exception("{0} already exists in dataframe. Consider using get_missing_column".format(column))
        res_frame[column] = np.NaN
    return res_frame



def add_data(DataFrame, newDataFrame, increament_id = True) :
    """Add new data to an existing (pandas) DataFrame automatically updating the Ids of the added data.
    
    Parameters
    ----------
        DataFrame : pd.DataFrame
            Existing data frame.
        newDataFrame : pd.DataFrame
            Data to add to DataFrame.
            
    Returns
    -------
        res : pd.DataFrame
            DataFrame with same frame as DataFrame with new data added."""
    
    check_parameter(DataFrame = (pd.DataFrame), newDataFrame = (pd.DataFrame))
    check_id(DataFrame)

    #Creating temp id in newdataframe if not found.
    if not has_id(newDataFrame) :
        newDataFrame["id"] = np.arange(0,len(newDataFrame))

    if not has_samedatashape(DataFrame, newDataFrame) : raise Exception("Cannot add 2 dataframes with different columns")

    #Separating past and new instances in case of trying to add dataframe to self.
    if id(DataFrame) == id(newDataFrame) : newDataFrame = DataFrame.copy()
    else : newDataFrame = newDataFrame.copy()

    if increament_id :
        if np.isnan(DataFrame["id"].max()) : last_id = 0
        else : 
            last_id = int(DataFrame["id"].max() + 1)
        newDataFrame["id"] = np.arange(start= last_id, stop= last_id + len(newDataFrame))
    
    res = pd.concat([DataFrame, newDataFrame])
    res = res.reset_index(drop= True)
    return(res)




def resetID (DataFrame, drop_id = False) :
    """
    Reset the column id conserving the current id order (ascending). By default the old id is kept.
    WARNING : this should not be applied without handling properly table referencing.
    To do so it is recommended to keep the old id and use foreign_key() to create a new reference using old ids. 
    """

    res = DataFrame.copy()
    if not drop_id : res["old_id"] = res["id"]
    res["id"] = np.arange(len(res))
    return(res)




def foreign_key(dataframe, fk_name, foreign_dataframe, left_key, right_key, how= 'left', validate = 'many_to_many') :
    #A finir : tester que la table right_key fait bien référence à une clef primaire (càd relation 1,1 avec l'id).
    """
    Update or create a foreign key column in {dataframe} referencing the primary key of {foreing_dataframe}.
    The functions uses a join with dataframe as the left table and foreign_datafram as the right table on dataframe.left_key = foreign_dataframe.right_key.
    This function should not be used to define (n,n) relationship but only (n,1) relation ship with (n : dataframe; 1 foreign_dataframe).
    
    Parameters
    ----------
        dataframe : pd.Dataframe
            Dataframe in which the foreign key is to be created.
        fk_name : str
            Name of the column to update/create in dataframe.
        foreign_dataframe : pd.Dataframe
            Dataframe which is referenced by the foreign key.
        left_key : str
            Column from dataframe on which to perform the table joining.
        right_key : str
            Column from foreign_dataframe on which to perform the table joining.
        validate : str, optional

            If specified, checks if merge is of specified type.
            "one_to_one" or "1:1": check if merge keys are unique in both left and right datasets.
            "one_to_many" or "1:m": check if merge keys are unique in left dataset.
            "many_to_one" or "m:1": check if merge keys are unique in right dataset.
            "many_to_many" or "m:m": allowed, but does not result in checks.
            

    Returns
    -------
        dataframe : pd.Dataframe
            Dataframe with new/updated fk column.

    """
    dataframe_copy: pd.DataFrame = dataframe.copy()

    #Integrity checks
    if not has_id(foreign_dataframe) : raise Exception("Foreign table has no id column")
    if not is_primary(foreign_dataframe, right_key) : raise Exception("{0} is not primary which makes unique referencing impossible".format(right_key))
    if not is_primarykey(foreign_dataframe, "id") : raise Exception("Foreign dataframe id is not a primary key, proceeding would result in non unique referencement.")

    #computing foreign key
    if has_column(dataframe_copy, fk_name) :
        dataframe_copy = dataframe_copy.drop(fk_name, axis= 1)
    columns = ["id"]
    columns.extend(right_key)
    right_index = foreign_dataframe.loc[:,columns]
    dataframe_copy = pd.merge(dataframe_copy,    right_index,     how= how,     left_on= left_key,     right_on= right_key, validate=validate)
    dataframe_copy = dataframe_copy.rename(columns={"id_x" : "id", "id_y" : fk_name})

    operation_was_valid = get_referencement_relation(dataframe_copy, fk_name, foreign_dataframe, 'id') in [('1','1'), ('N','1')]
    if not operation_was_valid : raise ReferenceError("Foreign key computation resulted in a {0} relation.".format(get_referencement_relation(dataframe_copy, fk_name, foreign_dataframe, 'id')))
    if len(dataframe) != len(dataframe_copy)  and how == 'left': raise MergeError("Merge operation lead to line duplication during Foreign key computing. Length before FK addition : {0}; length after : {1}".format(len(dataframe), len(dataframe_copy)))
    return(dataframe_copy)


def keep_columns(Dataframe, columns):
    """Drop all columns from Dataframe except {columns}.
    
    Parameters
    ----------
        Dataframe : pd.DataFrame
        columns : List[str]
            columns to keep.
            
    Returns
    -------
        res = pd.DataFrame
    """
    check_parameter(Dataframe = (pd.DataFrame), columns = list)
    for elmt in columns : check_parameter(elmt = (str))
    col2drop = list(Dataframe.keys())

    for column2keep in columns :
        if column2keep not in col2drop : continue
        else: col2drop.remove(column2keep)
    
    res = Dataframe.drop(columns= col2drop, axis= 1)

    return res




def get_elmtindex(elmt,List) :
    """Returns index (position) of elmt in list
    
    Parameters
    ----------
        elmt : any
        list : list
    Returns
    -------
        res : int
    """

    check_parameter(List = (list))

    for idx in range(0,len(List)) :
        if List[idx] == elmt : return idx
    
    raise Exception("Could not find elmt in List")




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

    #removing NULL rows
    tempForeignFrame = ForeignFrame.copy()
    tempForeignFrame = tempForeignFrame.dropna(subset= [foreign_key])
    nonnull_fk_number = len(tempForeignFrame)
    tempForeignFrame = keep_columns(tempForeignFrame, [foreign_key])

    tempLocalFrame = localDataFrame.copy()
    tempLocalFrame = localDataFrame.dropna(subset = [local_key])
    nonnull_lk_number = len(tempLocalFrame)
    tempLocalFrame = keep_columns(tempLocalFrame, [local_key])

    #group by keys
    group_fk = group_by(tempForeignFrame, foreign_key)
    group_lk = group_by(tempLocalFrame, local_key)

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





def remove_fromlist(elmt, List) :
    """Removes an element from a list

    Parameters
    ----------
        elmt : any
        list : list
    Returns
    -------
        res : list
    """

    check_parameter(List = (list))

    idx = get_elmtindex(elmt,List)
    res = List[0:idx] + List[idx+1:]

    return res


def group_by(Dataframe, columns):
    """Usual SQL Group-by function"""

    check_parameter(Dataframe= (pd.DataFrame), columns = (str, list))
    if type(columns) == str : 
        columns = [columns]
    else : 
        for clumn in columns : 
            check_parameter(clumn = (str))

    summarize = Dataframe.value_counts(subset= columns)
    summarize = summarize.reset_index().drop(0, axis= 1)

    return summarize
