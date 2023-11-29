# -*- coding: utf-8 -*-
# Author: Floric Slimani <floric.slimani@live.fr>


"""
This subpackages handles all Custom pandas framework related to the P-body project.
"""

from .DataFrames import newframe_Errors,newframe_Acquisitions,newframe_Cell,newFrame_Input,newframe_Pbody,newframe_Spots
from .DataFrames import get_features_name, get_Input

from .update import from_spike_value_compute_CellullarCycleGroup, from_nucleus_malat_proportion_compute_CellullarCycleGroup, from_IntegratedSignal_spike_compute_CellularCycleGroup
from .update import Pbody_AddCellFK, Spots_AddPbodyFK, Spots_AddCellFK, Spots_AddFK
from .update import from_detectionthreshold_remove_acquisition, from_Input_remove_damagedacquisitions, from_malat_remove_acquisition, from_pbodynum_remove_cell
from .update import remove_acquisitions, JoinCellAcquisition
from .results_integrity import Run_Results_Integrity_checks
