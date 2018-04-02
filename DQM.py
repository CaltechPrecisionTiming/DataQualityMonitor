import numpy as np
import os, re, shutil
import argparse
import ROOT as rt
from root_numpy import root2array
from lib.histo_utilities import create_TH1D, create_TH2D
from lib.cebefo_style import cebefo_style


# input_file = '../data/RawDataSaver0CMSVMETiming_Run81_0_Raw.root'

def parsing():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str, help="input root file")
    # parser.add_argument("-a", "--add_to_tree", default=False, action='store_true', help="Add output to input tree")
    parser.add_argument("-c", "--config", type=str, default='config/VME_test.txt', help="Config file")
    parser.add_argument("-s", "--save_loc", type=str, default='./', help="Saving location")

    # parser.add_argument("-T", "--train", type=str, help="variable to train", nargs='+')
    # parser.add_argument("-P", "--predict", type=str, help="variable to predict", nargs='+')
    # parser.add_argument("--addCPX", default=False, action='store_true', help="Add CPX variables")


    args = parser.parse_args()
    return args


class Config:
    def __init__(self, infile):
        self.raw_conf = open('config/VME_test.txt', 'r').readlines()
        self.channel = {}
        self.ch_ordered = []
        self.plots = []
        self.labels = []

        for i, l in enumerate(self.raw_conf):
            if l[0] == '#':
                continue
            l = l[0:-1].split(' ')
            if '-->Pri' in l[0]:
                self.plots = l[1:]
            elif len(self.labels) == 0:
                self.labels = l[1:]
            else:
                n = int(l[0])
                self.ch_ordered.append(n)
                l = l[1:]
                self.channel[n] = {}
                for k, val in zip(self.labels, l):
                    if 'idx' in k:
                        self.channel[n][k] = int(val)
                    else:
                        self.channel[n][k] = val


if __name__ == '__main__':
    cebefo_style()

    args = parsing()

    configurations = Config(args.config)

    input_file = args.input_file
    aux = re.search(r'Run[0-9]+', input_file)
    run_number = int(aux.group(0)[3:])

    save_loc = args.save_loc
    if os.path.isdir(save_loc):
        if save_loc[-1] != '/':
            save_loc += '/'
        out_dir = save_loc + 'Run' + str(run_number) + 'plots'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.mkdir(out_dir)
    else:
        print 'Save location not existing:', save_loc
        raise


    canvas = {}
    canvas['amp'] = {}
    canvas['int'] = {}
    canvas['wave'] = {}
    canvas['pos'] = {}
    canvas['w_pos'] = {}
    canvas['t_res_raw'] = {}

    branches = ['amp', 'channel', 'integral', 'time', 'x_dut', 'y_dut']
    for conf in configurations.channel.values():
        aux = conf['var_ref']
        if not aux in branches:
            branches.append(aux)

    root_file = rt.TFile.Open( input_file,"READ");
    if not root_file:
        print "[ERROR]: Input file not found"
        raise
    tree_name = root_file.GetListOfKeys().At(0).GetName()
    data = root2array(input_file, tree_name, branches=branches)

    for k in configurations.ch_ordered:
        conf = configurations.channel[k]


        '''=========================== Amplitude ==========================='''
        name = 'h_amp_'+str(k)
        title = 'Amplitude channel '+str(k)
        amp = (data['amp'].T)[k]
        h = create_TH1D(amp, name, title,
                        binning = [100, 0, 550],
                        axis_title = ['Peak amplitude [mV]', 'Events / 5.5 mV'])

        h.GetXaxis().SetRange(int(40/5.5)+1,int(450/5.5)+1)
        i_max = h.GetMaximumBin()
        h.GetXaxis().SetRange(1,100)
        peak = h.GetBinCenter(i_max)
        conf['amp_range'] = [max(20,peak*0.6), min(450, peak*1.6)]
        res = h.Fit('landau','LQSR', '', conf['amp_range'][0], conf['amp_range'][1])
        #     h.SetOptStat
        lowLimit = rt.TLine(conf['amp_range'][0], 0, conf['amp_range'][0],99) 
        highLimit = rt.TLine(conf['amp_range'][1], 0, conf['amp_range'][1],99) 
        lowLimit.SetLineWidth(3)
        highLimit.SetLineWidth(3)
        lowLimit.SetLineStyle(9)
        highLimit.SetLineStyle(9)


        if 'Amp' in configurations.plots:
            canvas['amp'][k] = rt.TCanvas('c_amp_'+str(k), 'c_amp_'+str(k), 800, 600)
            if(h.GetMaximum() - h.GetMinimum() > 5000):
                canvas['amp'][k].SetLogy()
            h.DrawCopy('E1')
            lowLimit.Draw("same")
            highLimit.Draw("same")
            canvas['amp'][k].donotdelete = [lowLimit, highLimit]
            canvas['amp'][k].Update()
            canvas['amp'][k].SaveAs(out_dir + '/Amp_ch'+str(k)+'.png')

        # Compute the selection and save it to the
        conf['amp_sel'] = np.logical_and(np.greater(amp,conf['amp_range'][0]), np.less(amp, conf['amp_range'][1]))


        '''=========================== Integral ==========================='''
        if 'Int' in configurations.plots:
            canvas['int'][k] = rt.TCanvas('c_int_'+str(k), 'c_int_'+str(k), 800, 600)

            name = 'h_int_'+str(k)
            title = 'Integral channel '+str(k)
            integral = -(data['integral'].T)[k]
            h = create_TH1D(integral, name, title,
                            binning = [100, np.min(integral), np.max(integral)],
                            axis_title = ['Integral [pC]', 'Events'])

            if(h.GetMaximum() - h.GetMinimum() > 5000):
                canvas['int'][k].SetLogy()
            h.DrawCopy('E1')
            canvas['int'][k].Update()
            canvas['int'][k].SaveAs(out_dir + '/Int_ch'+str(k)+'.png')


        '''=========================== Waveform color chart ==========================='''
        if 'WaveColor' in configurations.plots:
            name = 'h_wave_'+str(k)
            title = 'Waveform color chart channel '+str(k)

            ch = data['channel'][:,k].flatten()
            t = data['time'][:,conf['idx_time']].flatten()

            h = create_TH2D(np.column_stack((t,ch)), name, title,
            binning = [250, 0, np.max(t), 250, np.min(ch), np.max(ch)],
            axis_title = ['Time [ns]', 'Voltage [mV]']
            )

            h.SetStats(0)
            canvas['wave'][k] = rt.TCanvas('c_wave_'+str(k), 'c_wave_'+str(k), 800, 600)
            if(h.GetMaximum() - h.GetMinimum() > 5000):
                canvas['wave'][k].SetLogz()
            h.DrawCopy('colz')
            canvas['wave'][k].Update()
            canvas['wave'][k].SaveAs(out_dir + '/Waveform_ch'+str(k)+'.png')


        '''=========================== Track position ==========================='''
        if conf['idx_dut'] >= 0:
            name = 'h_pos_'+str(k)
            title = 'Track position at z_DUT[' + str(conf['idx_dut']) + '], for channel amp['+str(k)+'] in {'
            title += str(conf['amp_range'][0]) + ', ' + str(conf['amp_range'][1]) + '} mV'


            x = (data['x_dut'].T)[conf['idx_dut']]
            y = (data['y_dut'].T)[conf['idx_dut']]

            sel = np.logical_and(conf["amp_sel"], np.greater(x, -998))
            sel = np.logical_and(sel, np.greater(y, -998))

            x = x[sel]
            y = y[sel]

            h = create_TH2D(np.column_stack((x,y)), name, title,
                            binning = [100, np.min(x), np.max(x)+0.3*np.std(x), 100, np.min(y), np.max(y)+0.3*np.std(x)],
                            axis_title = ['x [mm]', 'y [mm]']
                           )

            canvas['pos'][k] = rt.TCanvas('c_pos_'+str(k), 'c_pos_'+str(k), 800, 600)
            h.DrawCopy('colz')

            if 'PosRaw' in configurations.plots:
                canvas['pos'][k].Update()
                canvas['pos'][k].SaveAs(out_dir + '/PositionXY_raw_ch'+str(k)+'.png')

            name = 'h_weight_pos_'+str(k)
            title = 'Track position at z_DUT[' + str(conf['idx_dut']) + '] weighted with channel amp['+str(k)+'] in {'
            title += str(conf['amp_range'][0]) + ', ' + str(conf['amp_range'][1]) + '} mV'
            h = create_TH2D(np.column_stack((x,y)), name, title, weights=amp[sel],
                            binning = [250, np.min(x), np.max(x)+0.3*np.std(x), 250, np.min(y), np.max(y)+0.3*np.std(x)],
                            axis_title = ['x [mm]', 'y [mm]']
                           )


            if 'PosWeight' in configurations.plots:
                canvas['w_pos'][k] = rt.TCanvas('c_w_pos_'+str(k), 'c_w_pos_'+str(k), 800, 600)
                h.DrawCopy('colz')
                canvas['w_pos'][k].Update()
                canvas['w_pos'][k].SaveAs(out_dir + '/PositionXY_amp_weight_ch'+str(k)+'.png')


        '''=========================== Raw time resolution ==========================='''
        if conf['idx_ref'] >= 0:
            time_stamp = (data[conf['var_ref']].T)[k]

            var_chref = configurations.channel[conf['idx_ref']]['var_ref']
            time_ref = (data[var_chref].T)[conf['idx_ref']]

            delta_t_amp_selections = np.logical_and(conf['amp_sel'], configurations.channel[conf['idx_ref']]['amp_sel'])

            delta_t = time_stamp - time_ref
            delta_t = delta_t[delta_t_amp_selections]

            name = 'h_delta_t_raw_'+str(k)
            title = 'Time resolution for channel '+str(k)

            x_axis_title = conf['var_ref'] + '[{}]'.format(k) + ' - ' + var_chref + '[{}]'.format(conf['idx_ref']) + '  [ns]'

            median = np.percentile(delta_t, 50)
            width = np.abs(np.percentile(delta_t, 20) - np.percentile(delta_t, 80))

            h = create_TH1D(delta_t, name, title,
                                binning = [ None, median-2*width, median+2*width],
                                axis_title = [x_axis_title, 'Events'])

            res = h.Fit('gaus', 'LQSR', '', np.percentile(delta_t, 20), np.percentile(delta_t, 80))

            if 'TimeResRaw' in configurations.plots:
                canvas['t_res_raw'][k] = rt.TCanvas('c_t_res_raw_'+str(k), 'c_t_res_raw_'+str(k), 800, 600)
                h.DrawCopy('E1')
                canvas['t_res_raw'][k].Update()
                canvas['t_res_raw'][k].SaveAs(out_dir + '/TimeResolution_raw_ch'+str(k)+'.png')

        '''=========================== Time resolution vs impact pointgit ==========================='''
