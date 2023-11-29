### testing operations subpackage
import pandas as pd
import numpy as np
from CustomPandasFramework.operations import get_referencement_relation


Rooms = pd.DataFrame({
    "id" : [16, 18],
    "spots" : [9,9]
})

Workers = pd.DataFrame({
    "id" : list(np.arange(1,16)),
    "Roomid" : [18]*7 + [16]*7+ [np.nan],
    "name" : ["Vera", "Davide", "Janna", "Hussein", "Rachel", "Ed", "MCsol", "Heloïse", "Flavia", "Oriane", "Céline", "Séverine", "Sylvain", "Soha",""]
})

Birthdays = pd.DataFrame({
    "id" : list(np.arange(1,15)),
    "Workersid" : list(np.arange(2,16)),
    "date" : ["30/10", "04/12", "18/07", "20/08", "12/08", "28/05", "19/06",  "06/07", "27/02", "27/04", "18/09", "15/07", "02/05", "17/02"]
})


print(Workers)
idx = Workers.query("Roomid.isnull()").index
new_df = Workers.loc[idx, :]
Workers.loc[idx,"id"] = new_df["id"] + 10
print(Workers)