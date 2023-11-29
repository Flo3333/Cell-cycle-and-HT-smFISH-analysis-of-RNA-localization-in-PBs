"""This modules contains function for boolean test and datashape checks"""


import pandas as pd
from bigfish.stack import check_parameter
from ..errors import *



def has_samedatashape(DataFrame1, DataFrame2):
    """Returns True if DataFrame1 has the same datashape than DataFrame2"""

    res = list_match(list(DataFrame1.columns), list(DataFrame2.columns))
    return res



def has_id(Dataframe) :
    """Returns true if the dataframe has a column named id"""
    return "id" in Dataframe.columns




def has_column(Dataframe, column) :
    """Returns true if the dataframe has a column named 'column''"""
    if column in list(Dataframe.columns) : res = True
    else : res = False
    
    return(res)


def is_primary(Dataframe:pd.DataFrame, columns : 'list[str]') :
    """Returns true if all non null values of column are unique in the Dataframe
        
    Parameters
    ----------
        DataFrame : pd.DataFrame
        columns : str or list[str]
            if list is given the group is tested as a primary element.
        
    Returns
    -------
        res : bool
    """
    if type(columns) == list : pass
    elif type(columns) == str : columns = [columns]
    else : raise TypeError("columns argument should be of type list['str'] or str. It is {0}".format(type(columns)))
    
    if len(Dataframe.dropna(subset=columns)) == 0 : return True
    groupby = Dataframe.groupby(columns, dropna= True)
    return len(Dataframe.dropna(subset= columns)) == len(groupby.count())




def is_primarykey(Dataframe: pd.DataFrame, column) :
    """Returns true if all values of columns are non null and unique. If Df is empty will return True
    
    Parameters
    ----------
        DataFrame : pd.DataFrame
        column : str
        
    Returns
    -------
        res : bool
    """
    if Dataframe.empty : return True
    return is_primary(Dataframe, column) and Dataframe.count().at[column] == Dataframe.shape[0]


def is_contained(list1, list2) : 
    """ Returns True if all list1 elements are in list2
    
    Parameters
    ----------
        list1 : list
        list2 : list
        
    Returns
    -------
        res : bool
        
    """

    check_parameter(list1 = (list), list2 = (list))
    truth = []

    for elmt in list1 : truth += [elmt in list2]
    res = all(truth)

    return res


def is_empty(DataFrame):
    res = len(DataFrame) == 0
    return res




def list_match(list1, list2) :
    """Returns true if both lists have the same elements without taking order into accounts
    
    Parameters
    ----------
        list1 : str
        list 2: str
    
    Return
    ------
        res : bool
        
    """

    check_parameter(list1 = (list), list2 = (list))

    res = is_contained(list1, list2) and is_contained(list2, list1)

    return res
    

def check_id(*DataFrames) : 
    """Check if DataFrame has valid id column, raises exception otherwise.

    Parameter
    ---------
        DataFrames : pd.DataFrame
    """

    for DataFrame in DataFrames :
        if not has_id(DataFrame) : raise NoIdError("No column nammed 'id' found in DataFrame")
        if not is_primarykey(DataFrame, "id") : 
            raise IdIsNotPrimaryError("'id' column is not a valid primary key column.")


def check_expectedcolumns(DataFrame, expectedcolumns) :
    """Raises exception if expected columns are not found in DataFrame
    
    Parameters
    ----------
        DataFrame : pd.DataFrame
        expectedcolumns : List[str]
        
    """

    found_columns = list(DataFrame.keys())
    if not is_contained(expectedcolumns, found_columns) : raise MissingColumnsError("Expected columns were not found in DataFrame")


def check_samedatashape(DataFrame1, DataFrame2):
    """Raise Exception if DataFrame1 and DataFrame2 haven't got equals datashapes.

    Parameters
    ----------
        DataFrame1 : pd.DataFrame
        DataFrame2 : pd.DataFrame
    """

    check_parameter(DataFrame1 = (pd.DataFrame), DataFrame2 = (pd.DataFrame))
    if len(DataFrame1.columns) != len(DataFrame2.columns): raise MissingColumnsError("Dataframes have not the same number of columns.")

    if not has_samedatashape(DataFrame1, DataFrame2): raise DataShapeError("DataFrames have not the same datashape.")




def check_isnotempty(*DataFrames) :
    
    truth_table = []
    for DataFrame in DataFrames:
        truth_table += [is_empty(DataFrame)]

    if any(truth_table) : raise EmptyFrameError("At least one DataFrame is empty.")