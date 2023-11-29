import pandas as pd
import numpy as np
import CustomPandasFramework.PBody_project.update as update
from CustomPandasFramework.integrity.Errors import MissingColumnError

def SpotsCell_view(Spots: pd.DataFrame) :

    groupby_obj = Spots.groupby(["spots_type","CellId" ])
    count = groupby_obj["id"].count()
    count_in_nuc = groupby_obj["InNucleus"].sum()
    count_in_cyto = count - count_in_nuc
    count_in_Pbody = groupby_obj["PbodyId"].count()

    view = pd.DataFrame({
        "count" : count,
        "count_in_nuc" : count_in_nuc,
        "count_in_cyto" : count_in_cyto,

        "count_in_Pbody": count_in_Pbody
    })
    

    return view

def SpotsRna_view(Spots: pd.DataFrame) :
    if "rna name" not in Spots.columns : raise  MissingColumnError("Missing `rna name` column consider update.AddRnaName")

    groupby_obj = Spots.groupby(["spots_type","rna name",])
    count = groupby_obj["id"].count()
    count_in_nuc = groupby_obj["InNucleus"].sum()
    count_in_cyto = count - count_in_nuc
    count_in_Pbody = groupby_obj["PbodyId"].count()
    
    view = pd.DataFrame({
        "count" : count,
        "count_in_nuc" : count_in_nuc,
        "count_in_cyto" : count_in_cyto,
        "count_in_Pbody": count_in_Pbody
    })

    return view


def Spots_view(Spots: pd.DataFrame) :
    if "rna name" not in Spots.columns : raise  MissingColumnError("Missing `rna name` column consider update.AddRnaName")

    groupby_obj = Spots.groupby(["spots_type","rna name", "CellId"])
    count = groupby_obj["id"].count()
    count_in_nuc = groupby_obj["InNucleus"].sum()
    count_in_cyto = count - count_in_nuc
    count_in_Pbody = groupby_obj["PbodyId"].count()
    
    view = pd.DataFrame({
        "count" : count,
        "count_in_nuc" : count_in_nuc,
        "count_in_cyto" : count_in_cyto,
        "count_in_Pbody": count_in_Pbody
    })

    return view


def Cell_view(Cell:pd.DataFrame, col_list= None) :
    if "rna name" not in Cell.columns : raise  MissingColumnError("Missing `rna name` column consider update.JoinCellAcquisition")
    if type(col_list) == type(None) : col_list = ["id", "rna name","AcquisitionId", "nuc_area", "nucleus_mean_mean_signal"] 
    if "id" not in col_list : col_list += ["id"]
    if "rna name" not in col_list : col_list += ["rna name"]

    view = Cell.loc[:, col_list].rename(columns= {"id" : "CellId"}).set_index(["rna name","CellId"])
    return view



def detection_view(Cell, Spots, Pbody: pd.DataFrame) : 
    """
    View containing data needed to check detection.
    """
    cell_rna_view = _detection_rna_view(Cell)
    cell_malat_view = _detection_malat_view(Cell)

    cell_view = pd.concat([cell_rna_view, cell_malat_view], axis= 0)
    spot_view = SpotsCell_view(Spots)

    Pbody_view = _detection_pbody_view(Pbody)

    detection_view = spot_view.join(cell_view)
    detection_view = detection_view.join(Pbody_view).rename(columns={"id" : "Pbody number"}).fillna(0).reorder_levels(["spots_type","rna name", "CellId"])
    detection_view["BigFish_total"] = detection_view["BigFish_count_in_nuc"] +  detection_view["BigFish_count_in_cyto"]
    detection_view["difference"] = detection_view["count"] - detection_view["BigFish_total"]

    return detection_view.sort_index()


def _detection_malat_view(Cell) : 
    
    if "rna name" not in Cell.columns : raise  MissingColumnError("Missing `rna name` column consider update.JoinCellAcquisition")
    
    cell_malat_view = Cell.loc[:,["AcquisitionId","id","rna name", "malat1 spots in nucleus", "malat1 spots in cytoplasm"]].rename(columns={"id" : "CellId","malat1 spots in cytoplasm" : "BigFish_count_in_cyto", "malat1 spots in nucleus" : "BigFish_count_in_nuc"})
    cell_malat_view["spots_type"] = "malat1"
    cell_malat_view = cell_malat_view.set_index(["spots_type","rna name", "CellId"])

    return cell_malat_view

def _detection_rna_view(Cell) : 
    
    if "rna name" not in Cell.columns : raise  MissingColumnError("Missing `rna name` column consider update.JoinCellAcquisition")
    
    cell_rna_view = Cell.loc[:,["AcquisitionId","id","rna name","nb_rna_in_nuc", "nb_rna_out_nuc"]].rename(columns={"id" : "CellId", "nb_rna_in_nuc" : "BigFish_count_in_nuc", "nb_rna_out_nuc" : "BigFish_count_in_cyto"})
    cell_rna_view["spots_type"] = "rna"
    cell_rna_view = cell_rna_view.set_index(["spots_type", "rna name", "CellId"])

    return cell_rna_view

def _detection_pbody_view(Pbody: pd.DataFrame) :
    
    if "rna name" not in Pbody.columns : raise  MissingColumnError("Missing `rna name` column consider using update.AddRnaName")

    view = Pbody.loc[:,["id", "CellId", "rna name"]]
    view = view.groupby(["rna name", "CellId"])["id"].count()
    return view


def CellularCycle_view(Cell, Spots= None, spots_view = None) : 
    """
    """

    if type(Spots) == type(None) == type(spots_view) : raise ValueError("At least one the two arguments 'Spots', 'Spots_view' must be passed.")

    if type(spots_view) == type(None) :
        spots_view = Spots_view(Spots)

    spots_view = spots_view.loc["malat1"]

    cell_view = Cell_view(Cell)
    cell_view["IntegratedSignal"] = cell_view["nuc_area"] * cell_view["nucleus_mean_mean_signal"]
    cellularcycle_view = cell_view.join(spots_view).sort_index()
    cellularcycle_view = update.from_IntegratedSignal_ranking_compute_CellularCycleGroup(cellularcycle_view)

    return cellularcycle_view

def CellularCycle_distribution_view(CellularCycle_view: pd.DataFrame) :
    """
    Computes view showing cellular_cycle classification repartition.
    """
    if 'CellId' not in CellularCycle_view.columns : X = 'id'
    else : X = 'CellId'
    Df = CellularCycle_view.loc[:,['rna name','cellular_cycle',X]].fillna('na')
    Df = Df.reset_index(drop=False).groupby(['rna name', 'cellular_cycle'])[X].count().reset_index(level=1, drop= False)
    total_DF = Df.groupby(axis= 0, level= 0)[X].sum().rename("total")
    Df = pd.merge(Df,total_DF, how= 'left', left_index=True, right_index=True).rename(columns= {X : 'count'})
    Df['proportion'] = Df['count'] / Df['total'] * 100
    Df = Df.set_index('cellular_cycle', append=True)
    
    return Df

def CellularCycle_distribution_overallmeasure(CellularCycle_distribution_view: pd.DataFrame) :
    """
    Computes weighted average (weigth = cell number) of cellular cycle distribution.

    Returns
    -------
        res : dict
            keys : 'g1_mean', 'g1_std','g2_mean', 'g2_std','S_mean', 'S_std'
    """
    res = {}
    for step in ['g1', 'g2', 's'] :

        if step not in CellularCycle_distribution_view.index.get_level_values(1) :
            res[step + '_mean'] = 0
            res[step + '_std'] = 0
        else :
            Df = pd.DataFrame()
            Df["normal"] = CellularCycle_distribution_view.loc[:,step,:]["count"] * CellularCycle_distribution_view.loc[:,step,:]['proportion']
            Df["squared"] = CellularCycle_distribution_view.loc[:,step,:]["count"] * CellularCycle_distribution_view.loc[:,step,:]['proportion']**2
            total_weigth = CellularCycle_distribution_view.loc[:,step,:]["count"].sum()
            mean = Df["normal"].sum() / total_weigth
            squared_mean = Df["squared"].sum() / total_weigth
            std = np.sqrt(np.abs(squared_mean - mean**2))
            res[step + '_mean'] = round(mean)
            res[step + '_std'] = round(std)

    return res



def Cell_detection_G1G2(Cell: pd.DataFrame, Spots: pd.DataFrame, classifier= update.from_IntegratedSignal_spike_compute_CellularCycleGroup, **classifier_kargs) :

    """"
    Creates view with following keys : RNA name and cellular_cycle informations are in index and the counts are in columns.
    index : ['plate_name', 'rna name', 'cellular_cycle', 'CellId'] or ['rna name', 'cellular_cycle', 'CellId']
    columns : ['mean number of rna in cell', 'mean number of rna in pbodies', 'mean fraction of rna in pbodies', 'mean number of rna in cytoplasm', 'mean number of malat1 in cell', 'mean number of malat1 in pbodies', 'mean fraction of malat1 in pbodies', 'mean number of malat1 in cytoplasm']
    
    -RNA name

    -mean number of RNA in cell (total), G1

    -mean number of RNA in cell (total), G2

    -mean number of RNA in cytoplasm, G1

    -mean number of RNA in cytoplasm, G2

    -mean number of RNA in PB,G1

    -mean number of RNA in PB,G2

    -mean fraction of RNA in PB, G1

    -mean fraction of RNA in PB, G2

    -same for MALAT1

    """

    if 'rna name' not in Cell.columns : raise KeyError('rna_name was not found in Cell df please use update.AddRnaName.')
    if 'rna name' not in Spots.columns : raise KeyError('rna_name was not found in Spots df please use update.AddRnaName.')
    if 'cellular_cycle' in Cell.columns : Cell = Cell.drop('cellular_cycle', axis= 0)
    if 'cellular_cycle' in Spots.columns : Spots = Spots.drop('cellular_cycle', axis= 0)
    Cell = classifier(Cell, **classifier_kargs)
    Spots = update.AddCellularCycle(Spots, Cell)


    if 'plate_name' in Cell.columns and 'plate_name' in Spots.columns :
        IsMultiplate = True
    else :
        IsMultiplate = False 

    rna_idx = Spots.query('spots_type == "rna"').index
    malat1_idx = Spots.query('spots_type == "malat1"').index
    rna_cyto_idx = Spots.loc[rna_idx,:].query('not InNucleus').index
    malat1_cyto_idx = Spots.loc[malat1_idx,:].query('not InNucleus').index

    if IsMultiplate : 
        grouping = ['plate_name', 'rna name', 'cellular_cycle', 'CellId']
        cellular_cycle_level = 3
    else : 
        grouping = ['rna name', 'cellular_cycle', 'CellId']
        cellular_cycle_level = 2
    
    grouper_rna, grouper_malat1, grouper_rna_cyto, grouper_malat1_cyto = [Spots.loc[idx,:].groupby(grouping) for idx in [rna_idx, malat1_idx, rna_cyto_idx, malat1_cyto_idx]]

    rna_number_cell= grouper_rna["id"].count()#Spots total par cellule
    rna_pbody_number_cell= grouper_rna["PbodyId"].count() #Spots dans pbody par cellule
    rna_in_pbody_fraction = (rna_pbody_number_cell / rna_number_cell)
    assert all(rna_in_pbody_fraction <= 1)
    malat1_number_cell= grouper_malat1["id"].count()#Spots total par cellule
    malat1_pbody_number_cell= grouper_malat1["PbodyId"].count()#Spots dans pbody par cellule
    malat1_in_pbody_fraction = malat1_pbody_number_cell / malat1_number_cell
    assert all(malat1_in_pbody_fraction <= 1)
    levels = [level for level in range(0,cellular_cycle_level)]

    mean_rna_in_pbody_fraction = rna_in_pbody_fraction.groupby(level=levels).mean().rename('mean fraction of rna in pbodies')
    mean_rna_number_cell= rna_number_cell.groupby(level=levels).mean().rename('mean number of rna in cell')
    mean_rna_number_cyto= grouper_rna_cyto["id"].count().groupby(level=levels).mean().rename('mean number of rna in cytoplasm')
    mean_rna_number_pbody = grouper_rna["PbodyId"].count().groupby(level=levels).mean().rename('mean number of rna in pbodies')
    
    mean_malat1_in_pbody_fraction = malat1_in_pbody_fraction.groupby(level=levels).mean().rename('mean fraction of malat1 in pbodies')
    mean_malat1_number_cell= malat1_number_cell.groupby(level=levels).mean().rename('mean number of malat1 in cell')
    mean_malat1_number_cyto= grouper_malat1_cyto["id"].count().groupby(level=levels).mean().rename('mean number of malat1 in cytoplasm')
    mean_malat1_number_pbody = grouper_malat1["PbodyId"].count().groupby(level=levels).mean().rename('mean number of malat1 in pbodies')

    view = pd.DataFrame(mean_rna_number_cell).join([mean_rna_number_cyto, mean_rna_number_pbody, mean_rna_in_pbody_fraction, mean_malat1_number_cell, mean_malat1_number_cyto, mean_malat1_number_pbody, mean_malat1_in_pbody_fraction], how= 'outer', on= grouping.remove('CellId'))
    return view


def Number_G1G2Cells_computed(Cell: pd.DataFrame, classifier, **classifier_kargs) :
    
    
    if 'rna name' not in Cell.columns : raise KeyError('rna_name was not found in Cell df please use update.AddRnaName.')

    if 'cellular_cycle' in Cell.columns : Cell = Cell.drop('cellular_cycle', axis= 1)
    if 'IntegratedSignal' not in Cell.columns : Cell = update.compute_IntegratedSignal(Cell)

    if len(Cell.at[0,'id']) == 2 :
        IsMultiplate = True
        columns = ['id','plate_name', 'rna name', 'IntegratedSignal', 'cellular_cycle']
    else :
        IsMultiplate = False 
        columns = ['id', 'rna name', 'IntegratedSignal', 'cellular_cycle']

    columns += ['IntegratedSignal']

    Cell = classifier(Cell, **classifier_kargs)
    Cell_df = Cell.loc[:,columns]

    if IsMultiplate : 
        grouping = ['plate_name', 'rna name', 'cellular_cycle']
    else : 
        grouping = ['rna name', 'cellular_cycle']

    CellCount = Cell_df.groupby(grouping)['id'].count().rename('CellNumber').reset_index(level=len(grouping)-1, drop= False)
    total = Cell_df.groupby(grouping[:len(grouping) - 1])['id'].count().rename('Total')
    CellCount = pd.merge(CellCount, total, how= 'left', on= grouping[:len(grouping) - 1])#TODO : this has high chances not to work --> check it
    CellCount = CellCount.set_index('cellular_cycle',append=True, verify_integrity=True).sort_index()
    CellCount['fraction'] = CellCount['CellNumber'] / CellCount['Total']
    assert all(CellCount["fraction"] <= 1)
    
    return CellCount



def cell_level_view(Cell:pd.DataFrame, Spots: pd.DataFrame) :
    if 'rna name' not in Cell : raise KeyError("rna name wasn't found in Cell dataframe.")
    if 'plate_name' in Cell.columns : grouping = ['plate_name', 'rna name', 'id']
    else : grouping = ['rna name', 'id']
    Cell_df = update.from_IntegratedSignal_spike_compute_CellularCycleGroup(Cell=Cell, column_name= 'cellular_cycle_spike_classification')
    Cell_df = update.from_IntegratedSignal_ranking_compute_CellularCycleGroup(Cell_raw=Cell_df, column_name= 'cellular_cycle_ranking_classification')
    columns = grouping + ['pbody number', "count pbody in nucleus", "count pbody in cytoplasm", "cellular_cycle_spike_classification", "cellular_cycle_ranking_classification"]
    view = Cell_df.loc[:,columns].rename(columns={'id' : 'CellId'})
    grouping = grouping[:len(grouping) - 1] + ['CellId'] #changing 'id' to 'CellId'
    view = view.set_index(grouping).sort_index()
    
    Cytospots_grouper = Spots[(Spots['spots_type'] == 'rna') & (~Spots['InNucleus'])].groupby(grouping)
    spots_grouper = Spots[(Spots['spots_type'] == 'rna')].groupby(grouping)
    cytospots_count_inpbody = Cytospots_grouper['PbodyId'].count().rename('nb_rna_cyto_in_pbody')
    cytospots_count = Cytospots_grouper['id'].count().rename('nb_rna_in_cyto')
    spots_count = spots_grouper['id'].count().rename('nb_rna_total')
    spots_count_inpbody = spots_grouper['PbodyId'].count().rename('nb_rna_in_pbody')
    view= view.join([spots_count, cytospots_count, cytospots_count_inpbody, spots_count_inpbody], how= 'left')
    view['nb_rna_in_nucleus'] = view['nb_rna_total'] - view["nb_rna_in_cyto"]
    view['rna_in_cytoplasm_fraction'] = view["nb_rna_in_cyto"] / view['nb_rna_total']
    view['cyto_rna_in_pbody_fraction'] = view['nb_rna_cyto_in_pbody'] / view['nb_rna_in_cyto']
    view['rna_in_pbody_fraction'] = view['nb_rna_in_pbody'] / view['nb_rna_total']


    assert all(view["rna_in_cytoplasm_fraction"].dropna() <= 1)
    assert all(view["cyto_rna_in_pbody_fraction"].dropna() <= 1)
    assert all(view["rna_in_pbody_fraction"].dropna() <= 1)

    return view