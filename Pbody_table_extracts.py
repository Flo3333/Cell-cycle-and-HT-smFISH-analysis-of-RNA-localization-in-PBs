"""
This script is meant to print specific extractions to Excel files.
"""

import os
from pbwrap.data import get_datetime
import CustomPandasFramework.PBody_project.load as load
import CustomPandasFramework.PBody_project.views as views
import CustomPandasFramework.PBody_project.update as update

#Settings
do_cell_level_extract = True
do_cell_detection_G1G2 =  True
do_Number_G1G2Cells_computed = True

classifiers = {'ranking_classifier' : update.from_IntegratedSignal_ranking_compute_CellularCycleGroup, 'spike_classifier' : update.from_IntegratedSignal_spike_compute_CellularCycleGroup}

p19p22 = '/media/floricslimani/SSD 4To/SSD_floricslimani/1_P_body/O8/stack_O8_p19p22/output/20231002 14-59-33/result_tables/'
p20 = '/media/floricslimani/SSD 4To/SSD_floricslimani/1_P_body/O8/stack_O8_p20/output/20231003 10-08-44/result_tables/'
p21 = '/media/floricslimani/SSD 4To/SSD_floricslimani/1_P_body/O8/stack_O8_p21/output/20231004 14-19-59/result_tables/'


#preparing dirs
date = get_datetime()
output_path = '/media/floricslimani/SSD 4To/SSD_floricslimani/1_P_body/O8/results/extractions/'
if not output_path.endswith('/') : output_path += '/'
output_path += date + '/'
os.makedirs(output_path)


#Pipeline

#Loading data
input_path_list = [p19p22,p20,p21]
plate_list = ['p19p22', 'p20', 'p21']

print("Loading data from input result tables...")
dataframes = load.load_results(input_path_list= input_path_list, plate_name_list=plate_list, display_steps = True)
dataframes['Acquisition'] = dataframes['Acquisition']
print('adding rna names...')
dataframes = load.add_rna_names(dataframes)

Acquisition = dataframes['Acquisition']
Cell = dataframes['Cell']
Pbody = dataframes['Pbody']
Spots = dataframes['Spots']

del dataframes

print("Extracts will be saved at {0}".format(output_path))
print("Extracting data...")
if do_cell_level_extract:
    print("Cell_level extract ",end='')
    cell_level_extract = views.cell_level_view(Cell, Spots)
    if len(cell_level_extract) > 900000 :
        print('warning extract is too long for excel file')
        pass
    else :
        cell_level_extract.reset_index(drop=False).to_excel(output_path + 'cell_level_extract_{0}.xlsx'.format(date))

    print('done')

if do_cell_detection_G1G2 :
    for name, classifier in classifiers.items() :
        print("classifier : {0} ; cell detection G1G2 extract ".format(name),end='')
        cell_detection_g1g2_extract = views.Cell_detection_G1G2(Cell, Spots, classifier)
        if len(cell_detection_g1g2_extract) > 900000 :
            print('warning extract is too long for excel file')
        else :
            cell_detection_g1g2_extract.reset_index(drop=False).to_excel(output_path + 'cell_detection_g1g2_extract_{1}_{0}.xlsx'.format(date, name))

        print('done')

if do_Number_G1G2Cells_computed :
    print("Number_G1G2Cells extract ",end='')
    for name, classifier in classifiers.items() :
        print("classifier : {0} ; cell detection G1G2 extract ".format(name),end='')
        Number_G1G2Cells = views.Number_G1G2Cells_computed(Cell,classifier)
        if len(Number_G1G2Cells) > 900000 :
            print('warning extract is too long for excel file')
        else :
            Number_G1G2Cells.reset_index(drop=False).to_excel(output_path + 'Number_G1G2Cells{1}_{0}.xlsx'.format(date, name))
        print('done')
