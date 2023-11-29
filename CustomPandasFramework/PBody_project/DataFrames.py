""" Making DataFrames structures tuned for our P-body project."""

import pandas as pd
import numpy as np
import os, re
from bigfish.stack import check_parameter
from bigfish.classification import get_features_name

#Temporary frame making
def get_Input(input_path, channels_list, check_completeness= True) :
    """ Returns a panda dataFrame filled with input folder info. To learn about expected data informations/shape see newFrame_Input.   
     
     Parameters
    ----------
    input_path : string.
        full path to the input folder as a string.
    channels_list : List[str].
        List of strings indicating the different channels of an acquisition. Thoses strings MUST be contained in the filenames otherwise they will be filtered.

    Returns
    -------
    Input : pd.Dataframe
        Dataframe containing informations on files within the input folder. See newFrame_Input for more informations.
        
        """     
     
    #Integrity checks
    check_parameter(input_path = (str), channels_list = (list))
    for string in channels_list : check_parameter(string = (str))
    
    #filename
    Input = newFrame_Input()
    filenames = os.listdir(input_path)
    Input["filename"] = filenames

    #Channel labelling
    for channel in channels_list :
        Input.loc[
        Input.filename.str.match(".*({0}).*".format(channel)), ["channel"]
        ]= channel
    Input = Input.dropna(subset=["channel"])

    #extension
    extensions = []
    ext_regex = "\.\w*$"
    for file in Input.loc[:,"filename"] :
        extensions += re.findall(ext_regex, file)
    Input.loc[:, "file extension"] = extensions

    #root filename
    root_filenames = []
    for line in Input.loc[:,["filename","channel"]].itertuples(index= False) :
        filename, channel = line
        regex = "({0})|(\.\w*$)".format(channel)
        root_filenames += [re.sub(r'\.(\d)',r'_\1',re.sub(regex, "", filename))]
    Input.loc[:, "root filename"] = root_filenames
    ##Acquisition index
    #Integrity
    channel_num = len(channels_list)
    Input_groupbyroot = Input.sort_values(["filename"]).value_counts(subset= "root filename")
    if not all(Input_groupbyroot == channel_num) and check_completeness :
        print(Input_groupbyroot[Input_groupbyroot != channel_num].index) 
        raise Exception("Some acquisition are missing at least one channel. Please check the completeness of files placed in input.")

    #Computing acquisition index
    Input_acquisitionindex = Input_groupbyroot.reset_index(drop = False)
    Input_acquisitionindex = Input_acquisitionindex.sort_values("root filename").reset_index(drop= True)
    Input_acquisitionindex = Input_acquisitionindex.reset_index(drop = False)
    Input_acquisitionindex = Input_acquisitionindex.rename(columns={"index" : "acquisition index"})
    Input_acquisitionindex = Input_acquisitionindex.drop(0, axis= 1)
    Input = Input.drop(["acquisition index", "rna name"], axis= 1)
    
    #Computing rna name
    regex = "-([\.\w-]*)f\d{2}--"
    rna_names = []
    for index in Input_acquisitionindex.index :
        rootfilename = Input_acquisitionindex.at[index, "root filename"]
        rna_name = re.findall(regex,rootfilename)[0]
        rna_names += [rna_name]
    Input_acquisitionindex.loc[:,"rna name"] = rna_names
    Input = pd.merge(Input, Input_acquisitionindex,
            how= "left", left_on= "root filename", right_on= "root filename")
    Input = Input.sort_values(["filename", "rna name","channel"]).reset_index(drop= True)

    return Input




#Empty frames making
def newFrame_Input(filename= None, channel= None, rootfilename= None, file_extension= None, acquisition_index=None, rna_name= None, **kwargs) :
    """Returns an empty pandas DataFrame with expected input data shape. 
    This frame is used for navigating through files put in the input folder during the segmentation process and isn't meant to be stored, therefor it does not contains any key.
        
    Returns
    -------
    
        new_Input : pd.Dataframe
            filename : full name of files within the input folder.
            channel : name of the channel such as dapi or egfp.
            root filename : filename without channel and extension, refering to the acquisition. In our project should be like : RNAspecies-wellID--sk1fk1fl1. It does NOT contains channel and extension informations.
            file extension : extension of the file.
            acquisition : Int. Computed, should be one per root filename, it aims at grouping different channels into the same acquisition.
            
    """

    if type(filename) == type(None) : filename = []
    if type(channel) == type(None) : channel = np.NaN
    if type(rootfilename) == type(None) : rootfilename = np.NaN
    if type(file_extension) == type(None) : file_extension = np.NaN
    if type(acquisition_index) == type(None) : acquisition_index = np.NaN
    if type(rna_name) == type(None) : rna_name = np.NaN

    dic = {
        "filename" : filename,
        "channel" : channel,
        "root filename" : rootfilename,
        "file extension" : file_extension,
        "acquisition index" : acquisition_index,
        "rna name" : rna_name,
        }

    Input_dic = dict(dic,**kwargs)
    new_Input = pd.DataFrame(Input_dic)
    
    return(new_Input)




def newframe_Acquisitions(pk= None, rootfilename= None, rna_name=None, cell_number=None, RNA_spot_threshold= None, malat1_spot_threshold= None, background_dapi_signal_mean=None, background_dapi_signal_std= None, **kwargs) :
    """Returns an empty pandas DataFrame with expected Acquisition data shape
        id (int, primary key) : Unique identifier to each acquisition
        filenameroot (str) : root of the file name meaning the part shared by all chanels filenames. 
    """

    if type(pk) == type(None) : pk= []
    if type(rootfilename) == type(None) : rootfilename= np.NaN
    if type(rna_name) == type(None) : rna_name= np.NaN
    if type(cell_number) == type(None) : cell_number= np.NaN
    if type(RNA_spot_threshold) == type(None) : RNA_spot_threshold= np.NaN
    if type(malat1_spot_threshold) == type(None) : malat1_spot_threshold= np.NaN
    if type(background_dapi_signal_mean) == type(None) : background_dapi_signal_mean= np.NaN
    if type(background_dapi_signal_std) == type(None) : background_dapi_signal_std = np.NaN


    dic = {
        "id" : pk, 
        "rootfilename" : rootfilename,
        "rna name" : rna_name,
        "cell number" : cell_number,
        "RNA spot threshold" : RNA_spot_threshold,
        "malat1 spot threshold" : malat1_spot_threshold,
        "background dapi signal mean" : background_dapi_signal_mean,
        "background dapi signal std" : background_dapi_signal_std
        }

    dic_Acquisition = dict(dic, **kwargs)
    res = pd.DataFrame(dic_Acquisition)    
    return(res)




def newframe_Cell(names_features_distance= True,
                    names_features_area= True,
                    names_features_centrosome= True,
                    names_features_dispersion= True,
                    names_features_foci= True,
                    names_features_intranuclear= True,
                    names_features_protrusion= True,
                    names_features_topography= True,
                    names_features_signal= True,
                    names_features_malat1= True,
                    names_features_pbody= True,
                    plot_index= True) :
    
    """Returns an empty pandas DataFrame with expected Cell data shape. """
    
    #Datamodel
    header = ['id', 'AcquisitionId']
    #Bigfish built-in features
    header += get_features_name(names_features_distance= names_features_distance,
                                names_features_area= names_features_area,
                                names_features_centrosome= names_features_centrosome,
                                names_features_dispersion= names_features_dispersion,
                                names_features_foci= names_features_foci,
                                names_features_intranuclear= names_features_intranuclear,
                                names_features_protrusion= names_features_protrusion,
                                names_features_topography= names_features_topography)
    #Custom features

    #Custom signal features
    if names_features_signal :
        header += [ "cell_coordinates", "label","label_bis", "bbox", 
                "nucleus_mip_mean_signal","nucleus_mip_max_signal","nucleus_mip_min_signal","nucleus_mip_median_signal",
                "nucleus_mean_mean_signal","nucleus_mean_max_signal","nucleus_mean_min_signal","nucleus_mean_median_signal",'nucleus area (px)','nucleus area (nm^2)']
    
    #Other Custom features
    if names_features_malat1 :
        header += ['malat1 spots in nucleus', 'malat1 spots in cytoplasm', 'cluster number']
               
    
    if names_features_pbody :
        header += ["pbody coordinates", 'pbody number', "count pbody in nucleus", "count pbody in cytoplasm"]

    #Plot index
    if plot_index : header += ['plot index'] # computed on its own during Cell Computation
    res = pd.DataFrame(columns= header)
    return(res)



def newframe_Pbody(distance = [0]):
    """
    Returns an empty pandas DataFrame with expected Pbody data shape
    """
    if isinstance(distance, (int, float)) : distance = [distance]
    elif not isinstance(distance, list): raise TypeError('distance parameter should be a list, it is a {0}'.format(type(distance)))

    header = [
        "id",
        "AcquisitionId",
        "CellId",
        "centroid_coordinates",
        "area",
        "volume",
        "label",
        "cell_label",
        "InNucleus"
        ]
    for dist in distance :
        header += ['rna {0}nm count'.format(dist)]
        header += ['malat1 {0}nm count'.format(dist)]
    res = pd.DataFrame(columns= header)
    
    return res




def newframe_Spots() :
    """
    Returns an empty pandas DataFrame with expected Spots data shape
    
    Columns
    -------
        'id' : int
            unique indentifier for spots table (PK)
        'AcquisitionId' : int
            FK to Acquisition table. (fov data level)
        'CellId' : int
            FK to Cell table (cell level result) : this means each rna with one CellId belong the this one cell.
        'spots coordinates' : tupple(int,int,int)
            coordinates of computed cell (z,y,x)
        'type' : str
            refers to the biological object wich is detected. Such as 'rna' or 'malat1'.
    
    """
    res = pd.DataFrame(
        {"id" : [],
        "AcquisitionId" : [],
        "CellId" : [],
        "PbodyId" : [],
        "spots_coords" : [],
        "spots_type" : [],
        "cell_label" : [],
        "Pbody_label" : [],
        "InNucleus" : []

    })
    res.loc[:,'InNucleus'] = res.loc[:,'InNucleus'].astype(bool)

    return res


def newframe_Errors() -> pd.DataFrame :
    res = pd.DataFrame({
        "rootfilename" : [],
        "error" : [],
        "message" : []
    })