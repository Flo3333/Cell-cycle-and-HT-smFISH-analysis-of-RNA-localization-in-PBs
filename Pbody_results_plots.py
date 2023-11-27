import pandas as pd
import numpy as np
import os
import CustomPandasFramework.PBody_project.update as update
import CustomPandasFramework.PBody_project.views as views
import pbwrap.plot.histogram as histogram
import pbwrap.plot.bar as bar
import pbwrap.plot.scatter as scatter
import pbwrap.plot.box as box
import matplotlib.pyplot as plt
from pbwrap.plot import get_colors_list
from CustomPandasFramework.PBody_project import from_IntegratedSignal_spike_compute_CellularCycleGroup


################################################################################
### SETTINGS
################################################################################

do_AddFK = False #Should be False except when working on running analysis
do_G1G2_quantification = True
do_dapi_plots = True
do_rna_plots = True
do_malat1_plots = True
do_pbody_plots = True
banned_genes = ['C-', 'DGKH','ELOA', 'SYNE1', 'C_FLAPXY_FLAPCD', 'MALAT1_FLAPXY'] #Oligopool 8
# banned_genes = ['AC005288_1', 'AC006329_1', 'AC006504_5', 'AC016708_1', 'AC021092_1', 'AC093010_3', 'AL021068_1', 'AL590428_1', 'ARHGAP11A', 'AURKA', 'BOLA3-AS1'] #Oligopool 5

################################################################################
### PIPELINE
################################################################################


print("Starting results plots pipeline")
# /media/floricslimani/SSD 4To/SSD_floricslimani/1_P_body/O8/stack_O8_p19p22/output/20231002 14-59-33/result_tables/
# /media/floricslimani/SSD 4To/SSD_floricslimani/1_P_body/O8/stack_O8_p20/output/20231003 10-08-44/result_tables/
# /media/floricslimani/SSD 4To/SSD_floricslimani/1_P_body/O8/stack_O8_p21/output/20231004 14-19-59/result_tables/

#Input
in_path = "/media/floricslimani/SSD 4To/SSD_floricslimani/1_P_body/O5/output/results O5p1/result_tables"
output_path = "/media/floricslimani/SSD 4To/SSD_floricslimani/1_P_body/O5/output/results O5p1/result_plots/"
os.makedirs(output_path, exist_ok=True)
if not in_path.endswith('/') : in_path += '/'
if not output_path.endswith('/') : output_path += '/'

#Data loading
print("Loading and preprocessing data")
Acquisition = pd.read_feather(in_path + 'Acquisition')
Cell = pd.read_feather(in_path + 'Cell')
Pbody = pd.read_feather(in_path + 'Pbody')
Spots = pd.read_feather(in_path + "Spots")


print("Banned genes list is currently set to {0}".format(banned_genes))
with open(output_path + 'banned_genes.txt', 'w') as txt :
    txt.write("These genes are removed from plots because they are damaged or not confluent enough.\n")
    txt.write(str(banned_genes))
drop_index = Acquisition.query("`rna name` in {0}".format(banned_genes)).index
drop_ids = list(Acquisition.loc[drop_index, "id"])

Acquisition = Acquisition.drop(drop_index, axis=0).sort_values("rna name").reset_index(drop= True)
cell_drop_id = Cell.query("AcquisitionId in {0}".format(drop_ids)).index
Cell = Cell.drop(cell_drop_id, axis=0)
pbody_drop_id = Pbody.query("AcquisitionId in {0}".format(drop_ids)).index
Pbody = Pbody.drop(pbody_drop_id, axis=0)
spots_drop_id = Spots.query("AcquisitionId in {0}".format(drop_ids)).index
Spots = Spots.drop(spots_drop_id, axis=0)

if do_AddFK :
    Pbody = update.Pbody_AddCellFK(Pbody,Cell)
    Spots = update.Spots_AddFK(Spots, Cell, Pbody)

Pbody = update.AddRnaName(Pbody, Acquisition)
Spots = update.AddRnaName(Spots, Acquisition)
Cell = update.AddRnaName(Cell,Acquisition)

#Decorators
g1g2RawData = histogram.RawData
dapi_signal = histogram.dapi_signal
dapi_density = histogram.dapi_density

#############
### PLOTS ###
#############

print("\n##################\n# Starting plots #\n##################\n")

bar.threshold(Acquisition, path_output= output_path + "threshold plot", show= False)
bar.total_cell_number(Acquisition, path_output= output_path + "total cell number computed", show= False)
box.cell_number(Acquisition= Acquisition, path_output= output_path + "BoxPlot_CellNumber", title= "Cell number", show= False)

########################
# G1G2 quantification # 
########################
if do_G1G2_quantification :
    print("plotting G1G2 quantification...")
    #A cellular_cycle classifier is a function taking Cell dataframe as parameter and returns Cell daframe with new key : 'cellular_cycle'
    classifier_dict = {'G1G2_spike_classification' : update.from_IntegratedSignal_spike_compute_CellularCycleGroup, 
                       'Ranking_classification' : update.from_IntegratedSignal_ranking_compute_CellularCycleGroup
                       }

    for name, classifier in  classifier_dict.items() :
        print("\n{0} classifier".format(name))
        g1g2_path = output_path + "G1G2_quantifications/{0}/".format(name)
        os.makedirs(g1g2_path, exist_ok= True)
        if 'cellular_cycle' in Cell.columns : Cell = Cell.drop('cellular_cycle', axis= 1)
        Cell = classifier(Cell)

        dapi_path = output_path + "DapiPlots/"
        os.makedirs(dapi_path, exist_ok= True)
        rnas = Cell.value_counts(subset="rna name").index
        plot_number = len(rnas)

        root = np.sqrt(plot_number)
        if root - np.floor(root) == 0 :
            n_lin = n_col = root
        elif root - np.floor(root) < 0.5 :
            n_lin = np.floor(root)
            n_col = np.floor(root) + 1
        else :
            n_lin = np.floor(root) + 1
            n_col = np.floor(root) + 1

        plt.close()
        fig, ax = plt.subplots(nrows= int(n_lin), ncols= int(n_col), figsize= (12*n_col, 12*n_lin))

        line = 0,
        col = 0
        plot_idx=1
        color_list = iter(get_colors_list(len(rnas)))
        for rna in rnas :
            color = next(color_list)
            gene_index = Cell.query("`rna name` == '{0}'".format(rna)).index
            gene_Cell = from_IntegratedSignal_spike_compute_CellularCycleGroup(Cell.loc[gene_index,:], auto_shift= True)
            plt.subplot(n_lin, n_col, plot_idx)
            #Plotting Hists
            histogram.dapi_signal(gene_Cell.query("cellular_cycle == 'g1'"), projtype= 'mean', summarize_type= 'mean', title= "{0} distribution histogram".format(rna), show= False, bins= 100, auto_bins= False, close = False, reset= False, label= "g1", color= 'blue', edgecolor = 'white', linewidth= 2)
            histogram.dapi_signal(gene_Cell.query("cellular_cycle == 'g2'"), projtype= 'mean', summarize_type= 'mean',  title= "{0} Histogram Dapi MeanSignal MeanProj".format(rna), show= False, bins= 100, auto_bins= False, close = False, reset= False, label= "g2", color= 'red', edgecolor = 'white', linewidth= 2)
            distribution = dapi_signal(gene_Cell, projtype= 'mean', summarize_type= 'mean',  title= "{0}".format(rna), xlabel= "", show= False, bins= 100, close = False, reset= False, label= rna, color= color, edgecolor = 'black', linewidth= 2)
            xmin,xmax,ymin,ymax = plt.axis()

            #Plotting Dean-Jett Model cellular cycle classification
            g1_expected_limitsup = np.quantile(distribution, 0.502)
            g2_expected_limitinf = np.quantile(distribution, 1-0.165)
            Xg1 = [g1_expected_limitsup] * 2
            Xg2 = [g2_expected_limitinf] * 2
            Y = [ymin,ymax]
            plt.plot(Xg1,Y, color= 'blue')
            plt.plot(Xg2,Y, color= 'red')

            plt.legend(loc= "upper right")
            plot_idx +=1

        fig.text(0.5, 0.04, 'Dapi integrated signal (area * pixel intensity)', ha='center', fontsize= 65)
        plt.suptitle("Dapi Individual Distributions", fontsize= 80, fontweight= 'bold')
        plt.savefig(g1g2_path + "individual_distributions_mosaic.png")
        plt.close()

        #rna plots
        if any(Cell['cellular_cycle'] == 's') :
            cellular_cycle_combination = [('g1', 'g2'), ('g1','s'), ('s','g2')]
        else :
            cellular_cycle_combination = [('g1','g2')]
        for cellular_cycle in cellular_cycle_combination :
            print("cellular_cycle combination : ", cellular_cycle)
            cellular_cycle_x, cellular_cycle_y = cellular_cycle

            scatter.G1G2_spot_area_Quantif(Cell, Spots, spots_type= 'rna', cellular_cycle_x=cellular_cycle_x, cellular_cycle_y=cellular_cycle_y,
                                           title= 'Rna spots per μm² quantification', show= False, reset= True, path_output= g1g2_path + 'rna_Spots_Area_{0}_vs_{1}'.format(cellular_cycle_x, cellular_cycle_y) , edgecolors= 'black')
            
            scatter.G1G2_Spots_Quantif(Cell, Spots, spots_type= 'rna', cellular_cycle_x=cellular_cycle_x, cellular_cycle_y=cellular_cycle_y,
                                       title= 'Rna spots gene quantification', show= False, reset= True, path_output= g1g2_path + 'rna_Spots_{0}_vs_{1}'.format(cellular_cycle_x, cellular_cycle_y) , edgecolors= 'black')
            
            scatter.G1G2_CytoSpotsInPbody_Quantif(Cell, Spots, Pbody=Pbody, spots_type= 'rna',  cellular_cycle_x=cellular_cycle_x, cellular_cycle_y=cellular_cycle_y,
                                       title= 'Cytoplasmic rna spots gene quantification', show= False, reset= True, path_output= g1g2_path + 'rna_CytoSpotsInPbody_{0}_vs_{1}'.format(cellular_cycle_x, cellular_cycle_y) , edgecolors= 'black')
            
            for distance in [0,100,200,400,600,800,1000,1500,2000] :
                print('G1G2 rna in pbody quantif distance = {0} '.format(distance), end='')    
                
                scatter.G1G2_SpotsInPbody_Quantif(Cell=Cell, Spots= Spots, Pbody= Pbody, spots_type= 'rna', distance_SpotPbody= distance, cellular_cycle_x=cellular_cycle_x, cellular_cycle_y=cellular_cycle_y,
                                                   title= "Rnas in P-bodies quantification per genes\nmaximal distance = {0}nm".format(distance), show= False, reset= True,path_output= g1g2_path + "rna_SpotsInPbody_{1}vs{2}_{0}".format(distance, cellular_cycle_x, cellular_cycle_y) , edgecolors= 'black')    
                print("ok")

            #malat1 plots
            scatter.G1G2_spot_area_Quantif(Cell, Spots, spots_type= 'malat1',  cellular_cycle_x=cellular_cycle_x, cellular_cycle_y=cellular_cycle_y,
                                        title= 'Malat1 spots per μm² quantification', show= False, reset= True, path_output= g1g2_path + 'malat1_Spots_Area_{0}_vs_{1}'.format(cellular_cycle_x, cellular_cycle_y) , edgecolors= 'black')
            
            scatter.G1G2_Spots_Quantif(Cell, Spots, spots_type= 'malat1',  cellular_cycle_x=cellular_cycle_x, cellular_cycle_y=cellular_cycle_y,
                                       title= 'Malat1 spots gene quantification', show= False, reset= True, path_output= g1g2_path + 'malat1_Spots_{0}_vs_{1}'.format(cellular_cycle_x, cellular_cycle_y), edgecolors= 'black')
            
            scatter.G1G2_CytoSpotsInPbody_Quantif(Cell, Spots, Pbody=Pbody, spots_type= 'malat1',  cellular_cycle_x=cellular_cycle_x, cellular_cycle_y=cellular_cycle_y,
                                       title= 'Cytoplasmic malat1 spots gene quantification', show= False, reset= True, path_output= g1g2_path + 'malat1_CytoSpotsInPbody_{0}_vs_{1}'.format(cellular_cycle_x, cellular_cycle_y), edgecolors= 'black')
            
            for distance in [0,100,200,400,600,800,1000,1500,2000] :    
                print('G1G2 malat1 in pbody quantif distance = {0} '.format(distance), end='')    
                
                scatter.G1G2_SpotsInPbody_Quantif(Cell=Cell, Spots= Spots, Pbody= Pbody, spots_type= 'malat1', distance_SpotPbody= distance,  cellular_cycle_x=cellular_cycle_x, cellular_cycle_y=cellular_cycle_y,
                                                title= "Malat1 in P-bodies quantification per genes\nmaximal distance = {0}nm".format(distance), show= False, reset= True, path_output= g1g2_path + 'malat1_SpotsInPbody_{1}vs{2}_{0}'.format(distance, cellular_cycle_x, cellular_cycle_y), edgecolors= 'black')
                print('ok')

        #classifier plots
        bar.CellularCycle_classification([Cell], [''], classifier,title= 'Cellular cycle repartition' ,show= False, path_output= g1g2_path + 'cellular_cycle_repartition')
    print("done")

########
# DAPI # 
########
if do_dapi_plots :
    print("plotting DAPI plots...", end= '')
    dapi_path = output_path + "DapiPlots/"
    os.makedirs(dapi_path, exist_ok= True)

    scatter.DapiSignal_vs_CellNumber(Cell, show= False, reset= True, path_output= dapi_path + "DapiSignal_VS_CellNumber", title= "Mean Dapi Signal VS Cell Number")
    histogram.RawData(Cell, variable_name= "nucleus area (nm^2)", show= False, reset= True, path_output= dapi_path + "Histogram_nucleus_area", title= "Nucleus area")
    histogram.RawData(Cell, variable_name= 'nucleus_mean_mean_signal', color= 'red', show= False, reset= True, path_output= dapi_path + "Histogram_mean_signal", title= "Mean dapi signal")
    histogram.dapi_signal(Cell, projtype= 'mean', summarize_type= 'mean', path_output= dapi_path + "Histogram_IntegratedSignal", title= "Integrated dapi signal \nzstack mean projection", show= False)
    box.dapi_signal(Cell= Cell,integrated_signal= True, show= False, path_output= dapi_path + "BoxPlot_Dapi_IntegratedSignal", title= "Integrated dapi signal \nzstack mean projection")
    scatter.dapi_signal(Cell= Cell, integrated_signal= True, show= False, path_output= dapi_path + "ScatterPlot_Dapi_Integrated_MeanSignal_MeanProj", title= "Dapi Integrated signal (mean projection-mean signal)")

    print(" done")

###############
# MALAT PLOTS #
###############
if do_malat1_plots :
    print("plotting malat1 plots...", end= "")
    cellcycle_view = views.CellularCycle_view(Cell, Spots)
    malat_path = output_path + "MalatPlots/"
    os.makedirs(malat_path, exist_ok= True)


    histogram.in_nuc_malat_proportion(cellcycle_view, path_output = malat_path + "MalatNucleus_proportion_hist", show= False, title= "Proportion of malat spots detected inside nucleus")
    histogram.malat_count(cellcycle_view, location= 'nucleus', path_output = malat_path + "Histogram_Malat_Nucleus", show= False, title= "Malat spots count in nucleus", bins= 250)
    histogram.malat_count(cellcycle_view, location= 'cytoplasm', path_output = malat_path + "Histogram_Malat_Cytoplasm", show= False, title= "Malat spots count in cytoplasm", bins= 250)
    box.count_Malat_per_Cell(CellCellular_cycle= cellcycle_view, path_output = malat_path + "BoxPlot_MalatSpotsPerCell", show= False, title= "Malat spots detected per cell")

    scatter.Malat_inNuc_asDapiIntensity(Cell,Spots, path_output= malat_path + "MalatSpots_DapiIntensity_MeanSignal_MeanProj", show= False, plot_linear_regression= True)
    scatter.count_Malat_per_Cell(CellCellular_cycle= cellcycle_view, path_output = malat_path + "ScatterPlot_MalatSpotsPerCell", show= False, title= "Malat spots detected per cell")
    print(" done")

#############
# RNA PLOTS #
#############
if do_rna_plots :
    print("plotting rna plots...", end= '')
    rna_path = output_path + "RnaPlots/"
    os.makedirs(rna_path, exist_ok= True)

    detection_view = views.detection_view(Cell, Spots, Pbody)
    histogram.rna_in_pbodies(Pbody, path_output= rna_path + "mean RNA in pbody", show= False)
    bar.cytoRNA_proportion_in_pbody(detection_view, path_output= rna_path + "cytoplasmic RNA proportion in pbody", title= "cytoplasmic RNA proportion in pbody", show= False)
    bar.RNA_proportion_in_pbody(detection_view, path_output= rna_path + "RNA proportion in pbody", title= "RNA proportion in pbody", show= False)
    bar.RNApercentage_in_nucleus(detection_view, path_output= rna_path + "RNA percentage inside nucleus", show= False, plot_in_and_out_bars= False)
    box.count_rna_per_Cell(detection_view= detection_view, path_output=None, show= False, title= "Rna spots detected per cell", close= False)
    scatter.count_rna_per_Cell(detection_view= detection_view, path_output= rna_path + "ScatterPlot_RnaSpotsPerCell", show= False, title= "Rna spots detected per cell")
    print(" done")

################
# P-Body Plots #
################
if do_pbody_plots :
    print("plotting P-bodies plots...", end= '')

    pbody_path = output_path + "PbodyPlots/"
    os.makedirs(pbody_path, exist_ok= True)
    bar.P_body_detect_inside_nucleus(Cell, show= False, path_output= pbody_path + "Pbodies_detected_inside_nucleus")
    box.count_pbody_per_Cell(Pbody= Pbody, show= False, path_output= pbody_path + "BoxPlot_PbodyPerCell", title= "Pbody detected number per cell")

    print(" done")
print("End of pipeline. Plots were saved at {0}".format(output_path))