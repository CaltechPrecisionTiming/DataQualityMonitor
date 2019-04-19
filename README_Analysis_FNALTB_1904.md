# Instructions to analyze runs from lxplus

## Set the environment
First of all gather all the code and create all the directories.
```
mkdir FNALTB_1904
cd FNALTB_1904

ln -s /eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Apr2019 data
git clone https://github.com/CaltechPrecisionTiming/DataQualityMonitor.git

mkdir results
```

Optional: what is possible to do it to replace the folder results with a symbolic link to the eos www folder. A possible example is:
```
rmdir results
ln -s ~/cernbox/www/FNALTB_1904
```

## Run the analysis

the analysis runs on couple of channels and produce plots bars-orinted.
It requires as input the creation of a special config file which states the details of the analysis. Example of the config file can be found [here](https://github.com/CaltechPrecisionTiming/DataQualityMonitor/blob/master/config/FNAL_TB_1904/Analysis_SiPM_Bar_v3.txt).

When the config file is ready, run the tile analysis with:
```
python Analysis_SiPM_Bar.py -C config.txt -S results_dir -i file_or_template -N runs_interval
```
where ``<vA>`` is the analysis version. The standard is to use ``<vA> = <v><NA>``, where ``<NA>`` is an incremental number identifying the subversion of the analysis, if any (e.g. ``<vA> = vr1`` or ``<vA> = vr12``).
It is save to run different analysis versions on the same reco version since different directories will be created.

A practical example:
```
python Analysis_SiPM_Bar.py -C config/FNAL_TB_1904/Analysis_SiPM_Bar_v3.txt -S ../results/v3/ -i ../data/VME/RecoData/RecoWithTracks/v3/RawDataSaver0CMSVMETiming_RunXXX_0_Raw.root -N 6369 6401
```
