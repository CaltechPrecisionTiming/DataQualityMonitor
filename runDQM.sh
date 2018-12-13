source ~/.bash_profile
cd $3/DataQualityMonitor

if [ ! -d "../results/$2" ]; then
  mkdir -p ../results/$2
fi

python DQM_SiPM.py -C config/FNAL_TB_1811/VME_$2.txt -S ../results/$2/ -i ../data/VME/RECO/$2/DataVMETiming_RunXXX.root -N $1
