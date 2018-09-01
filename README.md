# DataQualityMonitor
Code to control the quality of data

DQM_Coincidence.py Instructions
Run Dat2Root first
Use the root file from Dat2Root, run DQM_Coincidence.py
Example of how to run code
python DQM_Coincidence.py S13360-3025_10450_57V_anotherSiPM_57V_Na22_Coincidence_version2.root -S /home/lab-cpt03/DataQualityMonitor-master -N 0 0
-S is the saving location
-N is the number of runs, we only need to run once, put 2 numbers that are the same
Default config file: Na22_config.txt
need to add what graphs to draw in the first row, example: IL_2
