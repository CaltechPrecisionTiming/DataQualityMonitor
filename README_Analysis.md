# Instructions to analyze runs from FNAL TB

## Set the environment
First of all gather all the code and create all the directories.
```
mkdir TB_dir
cd TB_dir
mkdir -p data/VME/RAW
mkdir data/VME/RECO
mkdir data/Tracks
mkdir data/NimPlus

git clone https://github.com/CaltechPrecisionTiming/TimingDAQ.git
git clone https://github.com/CaltechPrecisionTiming/DataQualityMonitor.git

cd TimingDAQ
make -j8
cd -
```

Copy in the proper data folder all the necessary collected data.

## Decode the raw data and create root files

Create config files to run on the desired channels. A new version of config files is supposed to be created for each different set of analysis (i.e. whenever the analyzed channel content is not consistent).
Two config files need to be created and placed in the proper directories. One for ``TimingDAQ`` with the name ``VME_<v>.config`` (like this one [here](https://github.com/CaltechPrecisionTiming/TimingDAQ/blob/master/config/FNAL_TestBeam_1811/VME_vr3.config)) and on for the DQM with the name ``VME_<v>.txt`` (like this one [here](https://github.com/CaltechPrecisionTiming/DataQualityMonitor/blob/master/config/FNAL_TB_1811/VME_vr3.txt)), where ``<v>`` is the version name. The standard is to use ``v<N>`` for TB versions and ``vr<N>`` for post TB versions, where ``<N>`` is an increasing number.
When creating the ``TimingDAQ`` config make sure that it is placed in the same directory stated [here](https://github.com/CaltechPrecisionTiming/TimingDAQ/blob/master/automation/DecodeData.py#L14).
