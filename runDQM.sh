source ~/.bash_profile
cd /Users/fqnet1/FNAL_TB/18_11/DataQualityMonitor
python DQM_SiPM.py -C config/FNAL_TB_1811/VME_v3.txt -S ~/cernbox/ocerri/www/FNAL_TB_1811/ -i ../data/VME/RECO/v3/DataVMETiming_RunXXX.root -N $1
