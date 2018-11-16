source ~/.bash_profile
cd /Users/fqnet1/FNAL_TB/18_11/DataQualityMonitor
python DQM_SiPM.py -C config/FNAL_TB_1811/VME_$2.txt -S ~/cernbox/cpena/www/FNAL_TB_1811/ -i ../data/VME/RECO/$2/DataVMETiming_RunXXX.root -N $1
