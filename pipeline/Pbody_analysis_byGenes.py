import os, time, warnings, signal
import pbwrap.data as data
import pbwrap.preprocessing as prepro
import pbwrap.detection as detection
import pbwrap.quantification as quant
import CustomPandasFramework.PBody_project as framework
import CustomPandasFramework.PBody_project.DataFrames as DataFrame
import bigfish.stack as stack
import bigfish.multistack as multistack
import bigfish.plot as plot
import numpy as np
from random import sample 
from pbwrap.segmentation import Nucleus_segmentation, Cytoplasm_segmentation, pbody_segmentation
from pbwrap.plot import plot_detection_steps, plot_cell, output_spot_tiffvisual, nucleus_signal_control, plot_labels
from pbwrap.integrity import detectiontimeout_handler
from pbwrap.errors import DetectionError, CellnumberError, TooManySpotsError, SegmentationError, NoSpotError
from CustomPandasFramework.operations import add_data
from CustomPandasFramework.PBody_project import Run_Results_Integrity_checks
from bigfish.detection import detect_clusters    

####################
#### Parameters ####
####################

#Path
path_in = "/media/floricslimani/SSD 4To/SSD_floricslimani/1_P_body/O5/input/"
path_out = "/media/floricslimani/SSD 4To/SSD_floricslimani/1_P_body/O5/output/"
# path_in = "/home/floricslimani/Documents/Projets/1_P_body/Workshop/input/"
# path_out = "/home/floricslimani/Documents/Projets/1_P_body/Workshop/output/"
datetime = data.get_datetime()
result_path = "{0}{1}/".format(path_out, datetime)
pandas_path = "{0}/result_tables/".format(result_path)
os.makedirs(pandas_path, exist_ok= True)

#Images
channels_list = ["DAPI", "EGFP", "Cy3", "Alexa 647"]
acquisition_start = 0
z_min = 5
z_max = None
crop_zstack = [9, 15] #Cropping used ONLY during threshold computation. WARNING : if croping with z_min/z_max parameters, crop_zstack will do an additional crop for detection threshold computation

#Segmentation
nucleus_diameter = 80 #in pixels 
cytoplasm_diameter = 105 #in pixels
allow_only_one_nucleus_in_a_cell = True
allow_cell_without_nucleus = False
pbody_kernel_size = 2.25
min_cell_number = 1
model_cytosegmentation= 'cyto2'

#Detection
RNA_THRESHOLD= None
MALAT_THRESHOLD = None
rna_penalty = 1.4 #(malus) multiplicative factor applied to bigfish autothreshold.
malat1_penalty = 0.2
pbody_penalty = 1
max_rna_number = 1000 #1000
min_rna_number = 5 #5
voxelsz = (300,103,103)
spot = 0.75
rna_spot_radius = (int(round(350*spot)),int(round(150*spot)),int(round(150*spot)))
malat_spot_radius = (int(round(350*spot)),int(round(150*spot)),int(round(150*spot)))
alpha_rna = 0.5  # alpha impacts the number of spots per candidate region for dense region decomposition
beta_rna = 1
alpha_malat1 = 0.1
beta_malat1 = 1  # beta impacts the number of candidate regions to decomposeget_images_as_gen
rna_cluster_radius = 350 #nm
malat1_cluster_radius = 350 #nm
ndim= 3

#Measurements
#P-body
spot_distance = [0,100,200,400,600,800,1000,1500,2000]

############
# Settings #
############
#Cleaning
do_clean_damaged_acquisition = False #If true please make sure the server "Bertrand-Commun" is mounted on the computer
oligopool = 8.1
#Plots
do_plot_cells = False
do_steps_plot = True
#Visuals
do_rna_create_tif_spotvisual = True
do_malat1_create_tif_spotvisual = True
do_dapi_signal_control = True


###################################
########### Main script ###########
###################################

### Init ###
signal.signal(signal.SIGALRM, detectiontimeout_handler) #Initiating timeout handling
if not path_in.endswith('/') : path_in += '/'
if not path_out.endswith('/') : path_out += '/'

#### Input ####
Input = framework.get_Input(path_in, channels_list)


if do_clean_damaged_acquisition : Input = framework.from_Input_remove_damagedacquisitions(Input, print_deletion_number= True, oligopool= oligopool)
acquisition_num = data.get_acquisition_num(Input)
Acquisition_result = DataFrame.newframe_Acquisitions()
Cell_result = DataFrame.newframe_Cell()
Spots_result = DataFrame.newframe_Spots()
Pbody_result = DataFrame.newframe_Pbody(distance=spot_distance)

#Preparing dirs for cell/steps ploting
if do_plot_cells :
    cell_path = result_path + "cell_plots/"
    os.makedirs(cell_path)
if do_steps_plot :
    steps_path = result_path + "steps_plots/"
    os.makedirs(steps_path)

#Preparing logsmalat_threshold
log_path = path_out + datetime + "/"
os.makedirs(log_path, exist_ok= True)
parameters_log = data.parameter_log('parameters', log_path)
parameters_log.add_parameters(log_path,rna_penalty,malat1_penalty, RNA_THRESHOLD, MALAT_THRESHOLD, rna_spot_radius,malat_spot_radius,alpha_rna,beta_rna,alpha_malat1,beta_malat1,
                              rna_cluster_radius,malat1_cluster_radius,pbody_kernel_size,ndim, z_max, z_min, min_cell_number, max_rna_number, min_rna_number)
parameters_log.write()
error_log = data.error_log('errors_log', log_path)
log = data.run_log('log', log_path)


#########################
### Analysis pipeline ###
#########################
print("{0} acquisitions (fov) were found. Begining analysis pipeline...".format(acquisition_num))

#### Gene loop ####
rna_computed = []
for rna in Input.value_counts(subset= "rna name").index :

    #Global threshold computation and spot detection
    print("Starting RNA spot detection for all {0} gene's fovs.".format(rna))
    acquisition_list = list(Input.query("`rna name` == '{0}'".format(rna)).value_counts(subset= "acquisition index").index)
    fov_list = data.get_images_as_gen(path_in, Input, acquisition_list, "Cy3", z_min= z_min, z_max= z_max)
    spots, rna_threshold = detection.detect_spots(fov_list, ndim= 3, threshold= RNA_THRESHOLD, threshold_penalty= rna_penalty, crop_zstack = crop_zstack, return_threshold= True, voxel_size = voxelsz, spot_radius = rna_spot_radius)
    spots_gen = iter(spots)
    del fov_list,spots
    print("RNA threshold gene {0} : {1}".format(rna, rna_threshold))

    #malat detection
    print("Starting malat1 spot detection for all {0} gene's fovs.".format(rna))
    malat_list = data.get_images_as_gen(path_input= path_in, Input= Input, acquisition_list= acquisition_list, channels_list= "Alexa 647", z_min= z_min, z_max= z_max)
    malat_spots, malat_threshold = detection.detect_spots(malat_list, ndim= 3, threshold= MALAT_THRESHOLD, threshold_penalty= malat1_penalty, crop_zstack= [2,5], return_threshold= True, voxel_size = voxelsz, spot_radius = malat_spot_radius)
    malat_spots_gen = iter(malat_spots)
    del malat_spots, malat_list
    print("Malat1 threshold gene {0} : {1}".format(rna, malat_threshold))

    #field of view loop
    images = data.get_images_as_gen(path_input= path_in, Input= Input, acquisition_list= acquisition_list, channels_list= channels_list, z_min= z_min, z_max= z_max)
    for acquisition in acquisition_list :
        new_threshold = None
        #Preparing DataFrames
        Acquisition_Cell_frame = DataFrame.newframe_Cell()
        Acquisition_Spots_frame = DataFrame.newframe_Spots()
        Acquisition_id = acquisition
        rootfilename = data.get_rootfilename(acquisition, Input)
        rna_name = Input.at[Input.query("`acquisition index` == {0}".format(acquisition)).index[0], "rna name"]

        #Updating run log
        log.update(filename= rootfilename, rna_computed= rna_computed)
        print("Processing acquisition number {0} : {1}".format(acquisition, rootfilename))  

        #Get images
        try:
            alexa647 = next(images)
            cy3 = next(images)
            dapi = next(images)
            egfp = next(images)
            rna_spots = next(spots_gen)
            malat1_spots = next(malat_spots_gen)
            print("rna spots len : ",len(rna_spots))

            #2D segmentation image projection
            dapi_proj = stack.mean_projection(dapi)
            cy3_seg = stack.mean_projection(cy3)
            cy3_maxproj = stack.maximum_projection(cy3)
            egfp_proj = stack.mean_projection(egfp)
            alexa647_proj = stack.mean_projection(alexa647)
            image_contrasted = stack.rescale(cy3_maxproj, channel_to_stretch= 0)

            # Cell/Nucleus segmentation #
            nucleus_label = Nucleus_segmentation(dapi_proj, nucleus_diameter, min_cell_number= 10, model_type= 'Nuc_proj1.0',
                                                  use_gpu= True)
            nucleus_mask = nucleus_label.astype(bool)
            cytoplasm_label = Cytoplasm_segmentation(cy3_seg, model_type= model_cytosegmentation,
                                                         dapi= dapi_proj, diameter= cytoplasm_diameter, min_cell_number=10, use_gpu= True)
            
            nucleus_label, cytoplasm_label = multistack.match_nuc_cell(nucleus_label, cytoplasm_label, 
                                                                       allow_only_one_nucleus_in_a_cell, allow_cell_without_nucleus) 
            
            # Segmentation/Detection failure raising
            if cytoplasm_label.max() < min_cell_number :
                raise CellnumberError("{0} cells were segmented which is below threshdold ({1}).".format(cytoplasm_label.max(),min_cell_number))
            if max_rna_number != None :
                rna_percell = len(rna_spots)/cytoplasm_label.max()
                if rna_percell > max_rna_number : raise TooManySpotsError("{0} rna spots have been detected for this fov with threshold {1}. Thresholding or Fish targeting have probably failed.".format(rna_percell, rna_threshold)) 
            if min_rna_number != None :
                rna_percell = len(rna_spots)/cytoplasm_label.max()
                if rna_percell < min_rna_number : 
                    print("{0} rna spots have been detected for this fov with threshold {1}. Thresholding or Fish targeting have probably failed, trying individual detection".format(rna_percell, rna_threshold))
                    rna_spots, new_threshold = detection.detect_spots(cy3, ndim= 3,threshold_penalty= rna_penalty, crop_zstack = crop_zstack, return_threshold= True, voxel_size = voxelsz, spot_radius = rna_spot_radius)
                    rna_percell = len(rna_spots)/cytoplasm_label.max()
                    if rna_percell < min_rna_number : raise NoSpotError("{0} rna spots have been detected for this fov with threshold {1}. Thresholding or Fish targeting have probably failed, trying individual detection".format(rna_percell, new_threshold))
                    elif rna_percell > max_rna_number : raise TooManySpotsError("{0} rna spots have been detected for this fov with threshold {1}. Thresholding or Fish targeting have probably failed.".format(rna_percell, new_threshold))
                    else : print("Sucess : {0} RNAs per cell detected, proceeding with analysis.".format(rna_percell))

            # Pbody segmentation :
            egfp_rescale = stack.rescale(egfp_proj)
            pbody_label = pbody_segmentation(egfp_rescale, sigma= pbody_kernel_size , threshold= None, small_object_size= 10, fill_holes= True)
            # Spot detection 
            cy3_rmvbackgrnd = prepro.remove_mean_gaussian_background(cy3, 5)
            alexa647_rmvbackgrnd = prepro.remove_mean_gaussian_background(alexa647, 5)
            signal.alarm(300) #Due to high processing time in higly noised image we're setting a timeout after 5 mins

            # Dense region decomposition
            try : rna_postdecomp = detection.spot_decomposition_nobckgrndrmv(cy3_rmvbackgrnd, 
                                                                       spots= rna_spots, spot_radius= rna_spot_radius, voxel_size_nm= voxelsz, alpha= alpha_rna, beta= beta_rna)
            except (ValueError, RuntimeError) as error : # Value Error may occured when too few spots are detected. BigFish will raise an error while trying to decompose dense regions where no dots have been detected. Even tough this should be extremely rare regarding rna detection compared to malat1 detection, see explanation in comment below.
                warnings.warn("Error was caugth during rna spot decomposition, process will resume with spots as they were before decomposition. \nError caught : {0}".format(error))
                rna_postdecomp = rna_spots
            try : malat1_postdecomp = detection.spot_decomposition_nobckgrndrmv(alexa647_rmvbackgrnd, 
                                                                          spots= malat1_spots, spot_radius= rna_spot_radius, voxel_size_nm= voxelsz, alpha= alpha_malat1, beta= beta_malat1)
            except (ValueError, RuntimeError) as error : # Value Error may occured when too few spots are detected. BigFish will raise an error while trying to decompose dense regions where no dots have been detected. Such cases occur when Malat1 Fish has failed 'dense region' might appear because  fish failure leads to the median spot having low intensity. Even tough due to global thresholding, this region will not have any spots.
                warnings.warn("Error was caugth during malat spot decomposition, process will resume with spots as they were before decomposition. \nError caught : {0}".format(error))
                malat1_postdecomp = malat1_spots

            Acquisition_Pbody_frame = quant.compute_Pbody(Acquisition_id, Pbody_label= pbody_label, cell_label= cytoplasm_label, nucleus_mask= nucleus_mask, rna_coords= rna_postdecomp, malat1_coords= malat1_postdecomp, distance= spot_distance)
            Acquisition_Spots_frame = quant.compute_Spots(AcquisitionId= Acquisition_id, Cell_label= cytoplasm_label, Nucleus_mask=nucleus_mask, spots_dictionary= {'rna' : rna_postdecomp, 'malat1' : malat1_postdecomp}, Pbody_label= pbody_label)

            # Cluster detection
            cy3_spots_postclustering, cy3_clusters = detect_clusters(spots= rna_postdecomp, radius= rna_cluster_radius,  voxel_size= voxelsz, nb_min_spots=4)
            alexa647_spots_postclustering, alexa647_clusters = detect_clusters(spots= malat1_postdecomp, radius = malat1_cluster_radius, voxel_size= voxelsz)
        
            #Transcription sites / Focis
            rna_spots_no_ts, foci, transcription_sites = multistack.remove_transcription_site(cy3_spots_postclustering, 
                                                                                          clusters= cy3_clusters, nuc_mask= nucleus_mask, ndim= ndim)
            signal.alarm(0)

            #Steps ploting
            try :
                if do_steps_plot:
                    rnasteps_path = steps_path + "{0}/".format(rna_name)
                    os.makedirs(rnasteps_path, exist_ok= True)
                    #Cell segmentation
                    plot.plot_segmentation_boundary(cy3_seg, 
                                                cell_label= cytoplasm_label, nuc_label= nucleus_label, boundary_size= 3, contrast=True, title= "{0} : Cell Segmentation".format(rna_name), path_output= rnasteps_path + rootfilename + "_Cell_segmentation",  show= False),
                    #Pbody segmentation
                    plot.plot_segmentation(egfp_proj, pbody_label, 
                                                contrast= True, title= "{0} : Pbody_segmentation".format(rna_name),  show= False ,path_output= rnasteps_path + rootfilename + "_Pbody_segmentation")
                    plot_labels([cytoplasm_label,pbody_label],arrows=[False, True], show= False, path_output= rnasteps_path + rootfilename + "_cells_labels")

                #Spots visual ploting
                if do_rna_create_tif_spotvisual and len(rna_postdecomp) != 0:
                    visual_path = result_path + "visuals/{0}/".format(rna_name)
                    os.makedirs(visual_path + "detection/",exist_ok= True)

                    #Retrieving spots in pbodies
                    Z,Y,X = zip(*rna_postdecomp)
                    spots_inpbody = pbody_label.copy()
                    spots_inpbody[Y,X] *= -1
                    rna_in_pbodies = np.array(list(zip(*np.nonzero(spots_inpbody < 0))))
                    output_spot_tiffvisual(cy3_maxproj, [rna_spots, rna_postdecomp, rna_in_pbodies],
                                            dot_size= 2, path_output= visual_path + "detection/rna_{0}.tif".format(rootfilename))

                if do_malat1_create_tif_spotvisual and len(malat1_postdecomp) != 0:
                    visual_path = result_path + "visuals/{0}/".format(rna_name)
                    os.makedirs(visual_path + "detection/",exist_ok= True)

                    #Retrieving spots in pbodies
                    Z,Y,X = zip(*malat1_postdecomp)
                    spots_inpbody = pbody_label.copy()
                    spots_inpbody[Y,X] *= -1
                    malat1_in_pbodies = np.array(list(zip(*np.nonzero(spots_inpbody < 0))))

                    output_spot_tiffvisual(stack.maximum_projection(alexa647), [malat1_spots, malat1_postdecomp, malat1_in_pbodies], 
                                           dot_size= 2, path_output= visual_path + "detection/malat1_{0}.tif".format(rootfilename))
                    


            except Exception as error:
                raise error
                warnings.warn("An error occured during plots plotting : {0}".format(error))
                error_log.add_error(filename = rootfilename, error= 'plot', msg= "An error occured during plots plotting.")
                warnings.warn("Filename was added to error log. Process will try to perform cell analysis nonetheless.")
                


            
        ################
        # FOV Features #
        ################

            ####################
            # Cell computation #
            ####################
            other_coords = {"foci" : foci, "transcription_site" : transcription_sites, "malat1_coord" : malat1_postdecomp}
            fov_results = multistack.extract_cell(cell_label= cytoplasm_label, 
                                                  ndim=ndim, nuc_label= nucleus_label, rna_coord= cy3_spots_postclustering, others_coord= other_coords, image= image_contrasted, others_image= {"dapi" : dapi_proj, "smfish" : cy3_seg})
            cell_number = len(fov_results)
            print("number of cells identified: {0}".format(cell_number))


            if do_plot_cells :
                try :
                    rnacell_path = cell_path + "{0}/".format(rna_name)
                    os.makedirs(rnacell_path, exist_ok= True)

                except:
                    warnings.warn("An error occured during cell plots folders creation")

            #Cell analysis loop
            warnings.simplefilter("ignore")
            for cell_id, cell in enumerate(fov_results) :
                new_Cell = quant.compute_Cell(Acquisition_id, cell, cell_label= cytoplasm_label, voxel_size= voxelsz, dapi= dapi, Pbody_Acquisition=Acquisition_Pbody_frame)
                new_Cell["plot index"] = cell_id
                if not Cell_result.empty : cell_id += Cell_result["id"].max()
                Acquisition_Cell_frame = add_data(Acquisition_Cell_frame, new_Cell)
                if do_plot_cells :
                    try :
                        cell_title = rootfilename + "{0}_CellPlot_".format(rna_name) + str(cell_id)
                        plot_cell(cell, title= cell_title, show= False, path_output= rnacell_path + cell_title)

                    except:
                        warnings.warn("An error occured during cell plots plotting")
                        error_log.add_error(filename = rootfilename, error= 'plot', msg= "An error occured during cell plots plotting.")
                        warnings.warn("Filename was added to error log. Process will try to perform cell analysis nonetheless.")
            
            warnings.simplefilter("default")
            
            ###############
            #Signal_visual#
            ###############
            if do_dapi_signal_control and len(malat1_postdecomp) != 0:
                visual_path = result_path + "visuals/{0}/".format(rna_name)
                os.makedirs(visual_path,exist_ok= True)
                dapi_mip = stack.maximum_projection(dapi)
                centroids_coordinates = list(Acquisition_Cell_frame.loc[:,"cell_coordinates"])
                #Integrated signal visual
                measurement = "Integrated_Signal"
                g1spike_value = 9.45e11 # 2nd spike is around 1.6e08 
                os.makedirs(visual_path + measurement +'/', exist_ok=True)
                measures = list(Acquisition_Cell_frame.loc[:,"nucleus area (nm^2)"] * Acquisition_Cell_frame.loc[:,"nucleus_mean_mean_signal"])
                visu_title = "{0}/{1}_Dapi_Signal_Control".format(measurement,rootfilename)
                nucleus_signal_control(dapi_mip, nucleus_label, measures, centroids_coordinates, malat1_postdecomp, value_multiplicator= 1/g1spike_value, use_scientific_notation= False, show= False, title= measurement, path_output= visual_path + visu_title)



        #####################
        ### Error Catcher ###
        #####################


        except SegmentationError as error :
            warnings.warn(str(error))
            error_log.add_error(filename= rootfilename, error= 'cell segmentation', msg= str(error))
            warnings.warn("Filename was added to error log. Process will now compute next acquisition.")
            log.failure()
            continue

        except DetectionError as error:
            warnings.warn(str(error))
            error_log.add_error(filename= rootfilename, error= 'detection', msg= str(error))
            warnings.warn("Filename was added to error log. Process will now compute next acquisition.")
            log.failure()
            continue

        except AssertionError as error:
            raise error

        except Exception as error : 
            warnings.warn("Unexpected Exception was raised within analysis loop.")
            warnings.warn(str(error))
            error_log.add_error(filename= rootfilename, error = 'CRITICAL', msg= str(error))
            log.failure()
            pandas_path = "{0}/result_tables/".format(result_path)
            os.makedirs(pandas_path, exist_ok= True)
            if not Acquisition_result.empty : Acquisition_result.to_feather(pandas_path + "Acquisition")
            if not Cell_result.empty : Cell_result.to_feather(pandas_path + "Cell")
            raise error
            continue
            

        else :
            log.sucess()
            if type(new_threshold) != type(None) : 
                new_Acquisition = quant.compute_fov(Acquisition_id, rootfilename, rna_name, cell_number, new_threshold, malat_threshold, dapi= dapi, nucleus_mask = nucleus_mask)
            else : 
                new_Acquisition = quant.compute_fov(Acquisition_id, rootfilename, rna_name, cell_number, rna_threshold, malat_threshold, dapi= dapi, nucleus_mask = nucleus_mask)
            Acquisition_result = add_data(Acquisition_result, new_Acquisition, increament_id= False)
            Acquisition_Cell_frame["AcquisitionId"] = Acquisition_id
            Cell_result = add_data(Cell_result, Acquisition_Cell_frame)
            Spots_result = add_data(Spots_result, Acquisition_Spots_frame)
            Pbody_result = add_data(Pbody_result, Acquisition_Pbody_frame)

        finally :
            signal.alarm(0)
            warnings.simplefilter("default")

    rna_computed += [rna]
    if not Acquisition_result.empty : Acquisition_result.to_feather(pandas_path + "Acquisition")
    if not Cell_result.empty : Cell_result.to_feather(pandas_path + "Cell")
    if not Spots_result.empty : Spots_result.to_feather(pandas_path + "Spots")
    if not Pbody_result.empty : Pbody_result.to_feather(pandas_path + "Pbody")


##############
### OUTPUT ###
##############

Pbody_result = framework.Pbody_AddCellFK(Pbody_result, Cell_result, drop_cell_label=False, how= 'inner')
Spots_result = framework.Spots_AddFK(Spots_result, Cell_result, Pbody_result)

# Spots_result = framework.Spots_AddCellFK(Spots_result, Cell_result, drop_cell_label=False)
# Spots_result = framework.Spots_AddPbodyFK(Spots_result, Pbody_result, drop_pbody_label=False)

error_log.output_errors() #saving all filenames which encountered errors to be processed indivudually.
#Folders making and saving parameters
print("Analysis done.")

process_time = time.process_time()

#Runnig full integrity checks on result tables
print("Running integrity checks...")
log_report = Run_Results_Integrity_checks(Acquisition_result, Cell_result, Pbody_result= Pbody_result, Spots_result= Spots_result, Print_failed_checks=True)
log_report["run time"] = process_time
log_report["cell number"] = len(Cell_result)

#Saving tables
print("Saving results...")
Acquisition_result.reset_index(drop=True).to_feather(pandas_path + "Acquisition")
Cell_result.reset_index(drop=True).to_feather(pandas_path + "Cell")
Pbody_result.reset_index(drop=True).to_feather(pandas_path + "Pbody")
Spots_result.reset_index(drop=True).to_feather(pandas_path + "Spots")
Input.reset_index(drop=True).to_feather(pandas_path + "Input")
print("Pipeline finished after {0}s (process time).".format(time.process_time()))

log.endrun(log_report)
print("Log was saved at path : {0}".format(log_path))