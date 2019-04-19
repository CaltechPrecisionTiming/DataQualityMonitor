if [ ! -d "../results/online_dqm" ]; then
  mkdir -p ../results/online_dqm
fi

python DQM_SiPM.py -C config/FNAL_TB_1904/VME_v1.txt -S ../results/online_dqm/ -i ../data/VME/RecoData/RecoWithTracks/v3/RawDataSaver0CMSVMETiming_RunXXX_0_Raw.root -N $1
