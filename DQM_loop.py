import argparse
import os, re
from glob import glob
import numpy as np
import time

def parsing():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file_template", type=str, default="../data/VME/RecoData/RecoWithTracks/v3/RawDataSaver0CMSVMETiming_RunXXX_0_Raw.root")
    parser.add_argument("-C", "--config", type=str, default='config/FNAL_TB_1904/VME_v1.txt', help="Config file")
    parser.add_argument("-S", "--save_loc", type=str, default='../results/online_dqm/', help="Saving location")
    parser.add_argument("-g", "--group", type=int, default=5, help="Number of runs to be grouped")
    args = parser.parse_args()
    return args

def getDirsRunN(s):
    out = re.search('[0-9]+_DQM', s)
    if hasattr(out, 'group'):
        return int(out.group(0)[:-4])
    else:
        return np.nan

def getFileRunN(s):
    out = re.search('Run[0-9]+', s)
    if hasattr(out, 'group'):
        return int(out.group(0)[3:])
    else:
        return np.nan

args = parsing()
if args.save_loc[-1] != '/':
    args.save_loc += '/'

while(True):
    dir_list = glob(args.save_loc + '*')
    last_run_processed = int(np.nanmax(map(getDirsRunN, dir_list)))
    file_list = glob(args.input_file_template.replace('XXX', '*'))
    newest_run = int(np.nanmax(map(getFileRunN, file_list)))

    if newest_run - last_run_processed > args.group+1:
        cmd = 'python DQM_SiPM.py -C ' + args.config
        cmd += ' -S ' + args.save_loc
        cmd += ' -i ' + args.input_file_template
        cmd += ' -N {} {}'.format(last_run_processed + 1, newest_run - 1)
        os.system(cmd)
    else:
        print time.time(), 'Nothing to be done'
        time.sleep(5)
