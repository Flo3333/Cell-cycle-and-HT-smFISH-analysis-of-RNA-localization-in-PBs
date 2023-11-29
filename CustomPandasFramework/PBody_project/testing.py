import pandas as pd
import numpy as np
import os
import CustomPandasFramework.PBody_project as framework
import CustomPandasFramework.operations as dataop
import CustomPandasFramework.PBody_project.update as update

Rooms = pd.DataFrame({
    "id" : [16, 18],
    "spots" : [9,9]
})

Workers = pd.DataFrame({
    "id" : list(np.arange(1,15)),
    "Roomid" : [18]*7 + [16]*7,
    "name" : ["Vera", "Davide", "Janna", "Hussein", "Rachel", "Ed", "MCsol", "Heloïse", "Flavia", "Oriane", "Céline", "Séverine", "Sylvain", "Soha"],
    "language" : [('english','czech','french'),('english','italian','french'),('english','arabic','french'),('english','arabic'),('english','indian','french'),('english','french'),('english','french'),('english','quebecois','french'),('english','italian','french'),('english','french'),('english','french'),('english','french'),('english','french'),('english','arabic','french')]
})

Birthdays = pd.DataFrame({
    "id" : list(np.arange(1,15)),
    "Workersid" : list(np.arange(1,15)),
    "date" : ["30/10", "04/12", "18/07", "20/08", "12/08", "28/05", "19/06",  "06/07", "27/02", "27/04", "18/09", "15/07", "02/05", "17/02"]
})

Names = pd.DataFrame({
    'id' : np.arange(len(Workers["name"])),
    'name' : Workers["name"]
})



Murkers = pd.DataFrame({
    "id" : list(np.arange(1,15)),
    "Roomid" : [18]*7 + [16]*7,
    "name" : ["Vera", "Davide", "Janna", "Hussein", "Rachel", "Ed", "MCsol", "Heloïse", "Flavia", "Orianne", "Céline", "Séverine", "Sylvain", "Soha"],
    "language" : [('english','czech','french'),('english','italian','french'),('english','arabic','french'),('english','arabic'),('english','indian','french'),('english','french'),('english','french'),('english','quebecois','french'),('english','italian','french'),('english','french'),('english','french'),('english','french'),('english','french'),('english','arabic','french')]
})



path = "/media/floricslimani/SSD 4To/SSD_floricslimani/1_P_body/O8/stack_O8_p19p22 part 1 new dapi stack/output/20230830 11-57-24/result_tables/"
Acquisition = pd.read_feather(path + "Acquisition")

Acquisition.to_excel(path + "Acquisition.xlsx")