source ~/.bash_profile
cd $3/DataQualityMonitor
python DQM_SiPM.py -C config/FNAL_TB_1811/VME_$2.txt -S ~/cernbox/www/FNAL_TB_1811/$2/ -i ../data/VME/RECO/$2/DataVMETiming_RunXXX.root -N $1
