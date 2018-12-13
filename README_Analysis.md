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

mkdir results

git clone https://github.com/CaltechPrecisionTiming/TimingDAQ.git
git clone https://github.com/CaltechPrecisionTiming/DataQualityMonitor.git

cd TimingDAQ
make -j8
cd -
```

Copy in the proper data folder all the necessary collected data.

Optional: what is possible to do it to replace the folder results with a symbolic link to the eos ww folder. A possible example is:
```
rmdir results
ln -s ~/cernbox/www/TB results
```

## Decode the raw data and create root files

Create config files to run on the desired channels. A new version of config files is supposed to be created for each different set of analysis (i.e. whenever the analyzed channel content is not consistent). We will work in the assumption that for the set of considered runs the hardware has not been modified (OV, cable length, position w.t.r. to the tracker, etc.).
Two config files need to be created and placed in the proper directories. One for ``TimingDAQ`` with the name ``VME_<v>.config`` (like this one [here](https://github.com/CaltechPrecisionTiming/TimingDAQ/blob/master/config/FNAL_TestBeam_1811/VME_vr3.config)) and on for the DQM with the name ``VME_<v>.txt`` (like this one [here](https://github.com/CaltechPrecisionTiming/DataQualityMonitor/blob/master/config/FNAL_TB_1811/VME_vr3.txt)), where ``<v>`` is the version name. The standard is to use ``<v> = v<N>`` for TB versions and ``<v> = vr<N>`` for post TB versions, where ``<N>`` is an increasing number.
When creating the ``TimingDAQ`` config make sure that it is in a directory consistent with  [this](https://github.com/CaltechPrecisionTiming/TimingDAQ/blob/master/automation/DecodeData.py#L14).

### Configurations validation

Before running the decoding on the full set of runs it is necessary to be sure that the two config are properly set. Particular attention needs to be given to the baseline computation interval.
For this purpose, decode one single condition-safe run (let's assume it has a run number ``<RN>``), version you are willing to use.
```
cd TimingDAQ
python automation/DecodeData.py --vVME <v> -f -R <RN>
```
In the DQM config, make sure to add the ``WaveColor`` plot in the print statement. Then run:
```
mkdir results/<v>
cd DataQualityMonitor
python DQM_SiPM.py -C config/<TB_dir>/VME_<v>.txt -S ../results/<v>/ -i ../data/VME/RECO/<v>/DataVMETiming_RunXXX.root -N <RN> --No_list
```
Check the DQM output plot to validate the config file. In particular check form ``WaveColor_ch*.png`` taht the baseline is properly computed.

When the config is finalized, remove ``WaveColor`` from the print statement.

### Decode the whole set of runs

To create the root output file for all the runs in the interval ``<N_st>`` - ``<N_end>``:
```
cd TimingDAQ
python automation/run_REdecode.py --opt_DecodeData "xxNO_save_meas" --v_fast <v> -R <N_st> <N_end> --run_DQM
```
It might be usefull to use add ``xxforce`` at the string passed to ``--opt_DecodeData`` in order to force the re-processing.
Inside the ``results/<v>`` folder two files conteinning the good runs an bad runs will be automatically created if ``List`` is present in the ``TracksConsistency`` statement of the DQM config.

### Run the analysis

Two different analysis scripts exists. The first one is meant to run on tiles: it runs the analysis on single channels and produce plots tile-orinted. The second one is meant to run on bars: it runs the analysis on couple of channels and produce plots bars-orinted.
Both of them require the creation of a special config file which states the details of the analysis. Example of the config file can be found [here](https://github.com/CaltechPrecisionTiming/DataQualityMonitor/blob/master/config/FNAL_TB_1811/Analysis_SiPM_Tile_vr21.txt) for tiles analysis and [here](https://github.com/CaltechPrecisionTiming/DataQualityMonitor/blob/master/config/FNAL_TB_1811/Analysis_SiPM_Bar_vr4.txt) for bars analysis.

When the config file is ready, run the tile analysis with:
```
python Analysis_SiPM_Tile.py -C config/<TB_dir>/Analysis_SiPM_Tile_<vA>.txt -S ../results/<v>/ -i ../data/VME/RECO/<v>/DataVMETiming_RunXXX.root -N ../results/<v>/TracksConsistency_Good.txt
```
where ``<vA>`` is the analysis version. The standard is to use ``<vA> = <v><NA>``, where ``<NA>`` is an incremental number identifying the subversion of the analysis, if any (e.g. ``<vA> = vr1`` or ``<vA> = vr12``).
It is save to run different analysis versions on the same reco version since different directories will be created.
Similarly, to run the bar analysis:
```
python Analysis_SiPM_Bar.py -C config/<TB_dir>/Analysis_SiPM_Bar_<vA>.txt -S ../results/<v>/ -i ../data/VME/RECO/<v>/DataVMETiming_RunXXX.root -N ../results/<v>/TracksConsistency_Good.txt
```
