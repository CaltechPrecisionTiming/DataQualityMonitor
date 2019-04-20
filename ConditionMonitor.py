#!/usr/bin/env python

import numpy as np
from glob import glob
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--condition_location', type=str, default='../data/RP/', help='Name of output file')
args = parser.parse_args()

loc_condition = args.condition_location
flist = glob(loc_condition+'*.txt')

fields2import = ['RunN', 'BoxVoltage', 'BarVoltage', 'BoxHum', 'BoxTemp', 'x_stage', 'y_stage', 'z_rotation']
data = []
for fname in flist:
    f = open(fname)

    rn = int(re.search('n[0-9]+_', fname).group(0)[1:-1])
    ddic = {}
    for l in f.readlines():
        idx = l.find(':')
        ddic[l[:idx]] = l[idx+2:]

    aux = []
    for n in fields2import:
        if n == 'RunN':
            aux.append(rn)
        elif not n in ddic.keys():
            aux.append(np.nan)
        else:
            try:
                a = float(ddic[n])
            except ValueError:
                a = np.nan
            aux.append(a)
    f.close()
    data.append(tuple(aux))

dt = [(n, '<i4'if n=='RunN'else'<f4') for n in fields2import]
data = np.array(data, dtype=dt)


matplotlib.rcParams.update({'font.size': 16})
for i,n in enumerate(fields2import):
    if n == 'RunN': continue
    plt.figure(i, figsize=(8,6))
    plt.plot(data['RunN'], data[n], '.', lw=0)
    plt.xlim(6200, 50+np.max(data['RunN']))
    units = ''
    if 'Voltage' in n:
        units=' [V]'
    elif 'stage' in n:
        units=' [mm]'
    elif 'Hum' in n:
        units=' [%]'
    elif 'Temp' in n:
        units=' [C]'
    elif 'rota' in n:
        units=' [deg]'

    plt.ylabel(n+units)
    plt.xlabel('Run Number')
    plt.grid()
    plt.savefig('../../../www/FNALTB_1904/StageCondition/{}.png'.format(n), bbox_inches='tight')
