# Stress-Predict Dataset

This study aims to develop a **stress-predict dataset** and perform descriptive, statistical and classification analysis of biophysiological data collected from healthy individuals, who underwent various induced emotional states, to assess the relative sensitivity and specificity of common biophysiological indicators of stress and provide a stepping-stone towards the development of an accurate stress monitoring device. In this study, 35 healthy volunteers performed three different stress inducing tasks (i.e., Stroop color word test, interview session and hyperventilation session) with baseline/relax period in-between each task, for 60 minutes. Blood volume pulse (BVP), Inter-beat-intervals, and heart rate were being continuously recorded using Empatica watches while respiratory rate was estimated using PPG-based respiratory rate estimation algorithm.

## Dataset Files

#### |-- Processed_data 

            |-----Improved_All_Combined_hr_rsp_binary.csv (contain information of heart rates and respiratory rates of all the participants along with timestamps and labels (for nonstress/baseline and 1 for stress task duration))

            |-----Questionnaires_scores.xlsx (contains information about the PSS and STAI questionnaire scores of each participant)

            |-----Time_logs.xlsx (contain date and start/end time of each task for each participant, Irish standard time)

#### |-- Raw_data

            |-----SX Folder (folders with raw files from Empatica E4. Where X is participant number)

                  |-----ACC.csv (contains accelerometer data (x, y, z axis))

                  |-----BVP.csv (contains raw BVP data)

                  |-----EDA.csv (contains EDA data (skin conductance))

                  |-----HR.csv (contains heart rate data)

                  |-----IBI.csv (contains inter-beat-interval data)

                  |-----info.txt (contains information about all the csv file and sampling rate)
     
                  |-----tags_SX.csv (contains timestamp tags. start-end time of each task)

                  |-----TEMP.csv (contains skin temperature data)

## Libraries

Following libraries were used for analysis:

- Statistical Analysis (Rstudio):
  - tidyverse
  - gtsummary
  - sjPlot
  - knitr
  - nlme
  - tidyr
- Descriptive Analysis (python):
  - pandas
  - numpy
  - seaborn
  - matplotlib
- Classification Analysis (python):
  - Numpy
  - Tensorflow
  - Pandas
  - Scikitlearn
  - Scipy
  - Pickle
  - Matplotlib
  - Keras

## References

**When using this dataset, please cite the following:**

1. Talha Iqbal, Andrew J. Simpkin, Davood Roshan, Nicola Glynn, John Killilea, Jane Walsh, Gerard Molloy, Sandra Ganly, Hannah Ryman, Eileen Coen, Adnan Elahi, William Wijns, and Atif Shahzad. 2022. _"Stress Monitoring Using Wearable Sensors: A Pilot Study and Stress-Predict Dataset"_, Sensors 22, no. 21: 8135. https://doi.org/10.3390/s22218135
2. Talha Iqbal, Adnan Elahi, Sandra Ganly, William Wijns, and Atif Shahzad. _"Photoplethysmography-Based Respiratory Rate Estimation Algorithm for Health Monitoring Applications."_ Journal of medical and biological engineering 42, no. 2 (2022): 242-252. https://doi.org/10.1007/s40846-022-00700-z
