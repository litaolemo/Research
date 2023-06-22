#!/usr/bin/env python
# coding: utf-8

# In[1]:
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import numpy
import scipy
from scipy import stats
import warnings
warnings.filterwarnings("ignore")


# In[4]:


import statsmodels.api as sm
import itertools
import pandas as pd
import warnings
import matplotlib.pyplot as plt
# hr2 = pd.read_csv("file_resource/ecg-id-database-1.0.0/Person_03",names=["hr"])
# print(hr2.head())
# create new datatime index to hr2, every 0.5s

import wfdb
import matplotlib.pyplot as plt

# Specify the path to your downloaded data
path_to_data = 'file_resource/ecg-id-database-1.0.0/Person_03/'

# The record name is the filename without the extension
record_name = 'rec_1'

# Use the 'rdrecord' function to read the ECG data
record = wfdb.rdrecord(f'{path_to_data}/{record_name}')
# Plot the ECG data
plt.figure(figsize=(10, 4))
plt.plot(record.p_signal[:,1])
plt.title('ECG Signal')
plt.xlabel('Time (samples)')
plt.ylabel('Amplitude')
#plt.show()
hr2 = pd.DataFrame(record.p_signal[:,1],columns=["hr"])
hr2.to_csv("heart_rate.csv")

datetime_series = pd.Series(
    pd.date_range("2000-01-01", periods=10000, freq="2ms")
)
#
hr2.set_index(datetime_series,inplace=True)

hr2_train = hr2[:8000]
hr2_test = hr2[8000:]

hr2[0:800].plot()
plt.figure(figsize=(10, 4))
hr2[0:800].rolling(15,min_periods=1).mean().plot()

from statsmodels.tools.sm_exceptions import ValueWarning
warnings.simplefilter('ignore', ValueWarning)



# In[11]:
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score,mean_absolute_percentage_error


# # Compute the mean square error



def cal_sarima(q,d,p):
    mod = sm.tsa.SARIMAX(hr2_train[0:1000],
                                                order=(q,d,p),
                                                seasonal_order=(0, 1, 0, 420),
                                                enforce_stationarity=True,
                                                enforce_invertibility=True)

    results = mod.fit(maxiter=20)
    # results.plot_diagnostics(figsize=(15, 12))
    print("aic: "+ str(results.aic))
    print("mse: " + str(results.mse))
    pred = results.get_prediction(start=1000, end=2999, dynamic=False)


    y_forecasted = pred.predicted_mean

    y_truth = hr2_train[1000:3000]
    print(len(y_forecasted),len(y_truth))
    mse = mean_squared_error(y_truth, y_forecasted)
    # rmse = mean_squared_error(y_truth, y_forecasted, squared=False)
    mae = mean_absolute_error(y_truth, y_forecasted)
    mape = mean_absolute_percentage_error(y_truth, y_forecasted)
    r2=r2_score(y_truth, y_forecasted)
    with open("sarima_result.txt", "a") as f:
        f.write("sarima: "+str(q)+" "+str(d)+" "+str(p)+"\n")
        f.write("train aic: "+str(results.aic)+"\n")
        f.write("train mse: "+str(results.mse)+"\n")
        f.write("mse: "+str(mse)+"\n")
        f.write("mape: "+str(mape)+"\n")
        f.write("mae: "+str(mae)+"\n")
        f.write("r2: "+str(r2)+"\n")
        f.write("\n")
    # print('The Mean Squared Error of our forecasts is {}'.format(round(mse, 8)))
    del mod
    del results
    del pred




# In[ ]:

para_lsit = [
    (0,0,0),
    (0,0,1),
    (0,0,2),
    (0,0,3),
    (0,0,4),
    (0,0,5),
    (1, 0, 0),
    (1, 0, 1),
    (1, 0, 2),
    (1, 0, 3),
    (1, 0, 4),
    (1, 0, 5),
(2,0,0),
    (2,0,1),
    (2,0,2),
    (2,0,3),
    (2,0,4),
    (2,0,5),
(3,0,0),
    (3,0,1),
    (3,0,2),
    (3,0,3),
    (3,0,4),
    (3,0,5),
(4,0,0),
    (4,0,1),
    (4,0,2),
    (4,0,3),
    (4,0,4),
    (4,0,5),
    (5,0,0),
    (5,0,1),
    (5,0,2),
    (5,0,3),
    (5,0,4),
    (5,0,5),
]



# In[ ]:

for p,d,q in para_lsit:
    res = cal_sarima(p,d,q)
    del res



# In[ ]:


