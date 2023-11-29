import datetime as dt
import numpy as np
import bigfish.stack as stack
from bigfish.stack import check_parameter
from scipy.ndimage import binary_dilation

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

def hist_maximum(hist:tuple) :
    highest_count = np.array(hist[0]).max()
    bins_num = list(hist[0])
    index = bins_num.index(highest_count)
    if index < len(bins_num) : index +=1
    return hist[1][index]



def get_datetime():
    return dt.datetime.now().strftime("%Y%m%d %H-%M-%S")


def create_tiff_detection_check(cy3: np.ndarray, spots, path_output, dot_size = 2):
    """
    Creates 3D tiff image with cy3 channel and spot detected.
    """

    if cy3.ndim == 3 : 
        channel = stack.maximum_projection(cy3)

    z,y,x = list(zip(*spots))
    spots_mask = np.zeros_like(channel)
    spots_mask[y, x] = 1

    
    #enlarging dots
    if dot_size > 1 : spots_mask = binary_dilation(spots_mask, iterations= dot_size-1)


    spots_mask = stack.rescale(np.array(spots_mask, dtype = channel.dtype))
    
    im = np.zeros([2] + list(channel.shape))
    im[0,:,:] = channel
    im[1,:,:] = spots_mask

    stack.rescale(channel, channel_to_stretch= 0)
    stack.save_image(im, path_output, extension= 'tif')