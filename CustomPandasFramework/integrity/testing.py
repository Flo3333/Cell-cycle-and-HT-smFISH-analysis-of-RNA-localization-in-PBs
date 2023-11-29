import pandas as pd
import numpy as np
from CustomPandasFramework.integrity import is_primary, is_primarykey
# path = "/home/floricslimani/Documents/Projets/1_P_body/stack_O8_p21/output/20230531 17-01-21/result_tables/"
# Acquisition = pd.read_feather(path + "Acquisition")
# Cell = pd.read_feather(path + "Cell")

Rooms = pd.DataFrame({
    "id" : [16, 18],
    "spots" : [9,9]
})

Workers = pd.DataFrame({
    "id" : list(np.arange(1,15)),
    "Roomid" : [18]*7 + [16]*7,
    "name" : ["Vera", "Davide", "Janna", "Hussein", "Rachel", "Ed", "MCsol", "Heloïse", "Flavia", "Oriane", "Céline", "Séverine", "Sylvain", "Soha"]
})

Birthdays = pd.DataFrame({
    "id" : list(np.arange(1,15)),
    "Workersid" : list(np.arange(1,15)),
    "date" : ["30/10", "04/12", "18/07", "20/08", "12/08", "28/05", "19/06",  "06/07", "27/02", "27/04", "18/09", "15/07", "02/05", "17/02"]
})

fruits = ['mela', 'pomme', 'apple', 'orange', 'abricot', 'arancia', 'peperone', 'poivron', 'chou','cavolo', 'orange']
country = ['Italia', 'France','England', 'France', 'France', 'Italia', 'Italia', 'France', 'France', 'Italia', 'England']

Fruits = pd.DataFrame({
    'id' : np.arange(len(fruits)),
    'name' : fruits,
    'country' : country
})

Void = pd.DataFrame({'id' : []})

print(Fruits.value_counts(subset="name"))