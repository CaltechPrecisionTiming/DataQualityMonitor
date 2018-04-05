import numpy as np
import os, re, shutil
import argparse
import ROOT as rt
from root_numpy import root2array, tree2array
from lib.histo_utilities import create_TH1D, create_TH2D
from lib.cebefo_style import cebefo_style



def parsing():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str, help="input root file", nargs='+')
    parser.add_argument("-C", "--config", type=str, default='config/VME_test.txt', help="Config file")
    parser.add_argument("-S", "--save_loc", type=str, default='./', help="Saving location")

    parser.add_argument("-B", "--batch", default=False, action='store_true', help="Root batch mode")

    parser.add_argument("-N", "--runs_interval", default=None, type=int, help="Runs to run", nargs='+')

    args = parser.parse_args()
    return args


class Config:
    def __init__(self, infile):
        self.raw_conf = open(infile, 'r').readlines()
        self.channel = {}
        self.ch_ordered = []
        self.plots = []
        self.labels = []
        self.xy_center = [0,0, 15]

        for i, l in enumerate(self.raw_conf):
            if l[0] == '#':
                continue
            l = l[0:-1].split(' ')
            if '-->Pri' in l[0]:
                self.plots = l[1:]
            elif '-->XYcenter' in l[0]:
                self.xy_center = [float(l[1]), float(l[2]), float(l[3])]
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

    if args.batch:
        rt.gROOT.SetBatch()

    configurations = Config(args.config)

    if args.runs_interval==None:
        aux = re.search(r'Run[0-9]+', args.input_file[0])
        flag = aux.group(0)[3:]
        if len(args.input_file) > 1:
            aux = re.search(r'Run[0-9]+', args.input_file[-1])
            flag += '_'+aux.group(0)[3:]
    elif len(args.runs_interval)==2:
        deduced_file_list = []
        for i in range(args.runs_interval[0], args.runs_interval[1]+1):
            aux = args.input_file[0].replace('XXX', str(i))
            deduced_file_list.append(aux)

        args.input_file = deduced_file_list
        flag = str(args.runs_interval[0]) + '_' + str(args.runs_interval[1])
    elif len(args.runs_interval) > 2:
        deduced_file_list = []
        for i in args.runs_interval:
            aux = args.input_file[0].replace('XXX', str(i))
            deduced_file_list.append(aux)
        args.input_file = deduced_file_list
        flag = str(args.runs_interval[0]) + '-' + str(args.runs_interval[-1])
    else:
        raise

        print 'Directory flag: ', flag
    save_loc = args.save_loc
    if os.path.isdir(save_loc):
        if save_loc[-1] != '/':
            save_loc += '/'
        out_dir = save_loc + 'Run' + str(flag) + '_plots'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.mkdir(out_dir)
        shutil.copyfile('./lib/index.php', out_dir+'/index.php')
    else:
        print 'Save location not existing:', save_loc
        raise



    chain = rt.TChain('pulse')
    chain.tree_name = tree_name = 'pulse'
    for i, f in enumerate(args.input_file):
        root_file = rt.TFile.Open( f,"READ");
        if not root_file:
            print "[ERROR]: Input file not found:", f
            continue
        tree_name = root_file.GetListOfKeys().At(0).GetName()
        if i == 0:
            if tree_name != 'pulse':
                print 'Change'
                chain = rt.TChain(tree_name)
                chain.tree_name = tree_name
        elif tree_name != chain.tree_name:
            print 'Skipping', f, ' for tree name incompatibility'
            continue
        chain.Add(f)

    canvas = {}
    canvas['amp'] = {}
    canvas['int'] = {}
    canvas['risetime'] = {}
    canvas['wave'] = {}
    canvas['pos'] = {}
    canvas['w_pos'] = {}
    canvas['t_res_raw'] = {}

    rt.gStyle.SetStatY(0.98);
    rt.gStyle.SetStatX(0.999);
    rt.gStyle.SetStatW(0.15);
    rt.gStyle.SetStatH(0.1);

    for k in configurations.ch_ordered:
        conf = configurations.channel[k]

        '''=========================== Amplitude ==========================='''
        name = 'h_amp_'+str(k)
        title = 'Amplitude channel '+str(k)

        h = rt.TH1D(name, title, 100, 3, 550)
        h.SetXTitle('Peak amplitude [mV]')
        h.SetYTitle('Events / '+str(h.GetBinWidth(1))+' mV')
        chain.Project(name, 'amp['+str(k)+']')

        h.GetXaxis().SetRange(int(40/5.5)+1,int(450/5.5)+1)
        i_max = h.GetMaximumBin()
        h.GetXaxis().SetRange(1,100)
        peak = h.GetBinCenter(i_max)
        conf['amp_range'] = [max(40,peak*0.6), min(450, peak*1.7)]
        conf['amp_sel'] = '(amp['+str(k)+'] < ' + str(conf['amp_range'][1])
        conf['amp_sel'] += ' && '
        conf['amp_sel'] += 'amp['+str(k)+'] > ' + str(conf['amp_range'][0]) + ')'

        res = h.Fit('landau','LQSR', '', conf['amp_range'][0], conf['amp_range'][1])

        if 'Amp' in configurations.plots:
            canvas['amp'][k] = rt.TCanvas('c_amp_'+str(k), 'c_amp_'+str(k), 800, 600)
            if(h.GetMaximum() - h.GetMinimum() > 5000):
                canvas['amp'][k].SetLogy()
            h.DrawCopy('E1')

            line = rt.TLine()
            line.SetLineColor(6)
            line.SetLineWidth(2)
            line.SetLineStyle(10)
            line.DrawLine(conf['amp_range'][0], 0, conf['amp_range'][0], h.GetBinContent(i_max))
            line.DrawLine(conf['amp_range'][1], 0, conf['amp_range'][1], h.GetBinContent(i_max))

            canvas['amp'][k].Update()
            canvas['amp'][k].SaveAs(out_dir + '/Amp_ch'+str(k)+'.png')


        '''=========================== Integral ==========================='''
        if 'Int' in configurations.plots:
            canvas['int'][k] = rt.TCanvas('c_int_'+str(k), 'c_int_'+str(k), 800, 600)

            name = 'h_int_'+str(k)
            title = 'Integral channel '+str(k)
            h = rt.TH1D(name, title, 100, 0, 550)
            h.SetXTitle('Integral [pC]')
            h.SetYTitle('Events / '+str(h.GetBinWidth(1))+' mV')
            chain.Project(name, '-integral['+str(k)+']')

            if(h.GetMaximum() - h.GetMinimum() > 5000):
                canvas['int'][k].SetLogy()
            h.DrawCopy('E1')
            canvas['int'][k].Update()
            canvas['int'][k].SaveAs(out_dir + '/Int_ch'+str(k)+'.png')

        '''=========================== Risetime ======================================='''
        if 'Risetime' in configurations.plots:
            canvas['risetime'][k] = rt.TCanvas('c_risetime_'+str(k), 'c_int_'+str(k), 800, 600)

            t_low  = .0
            t_high = 10;#20 ns range
            nbins  = 100
            name   = 'h_risetime_' + str(k)
            title  = 'Risetime ' + str(k)
            h      = rt.TH1D(name, title, nbins, t_low, t_high)
            h.SetXTitle('risetime [ns]')
            h.SetYTitle('events /' + str(h.GetBinWidth(1)) + 'ns')
            if conf['idx_ref'] == -1:#do not apply cut again on reference channels
                cut = conf['amp_sel']
            else:
                cut = conf['amp_sel'] + '&&' + configurations.channel[conf['idx_ref']]['amp_sel']

            chain.Project(name, 'risetime['+str(k)+']', cut)

            if (h.GetMaximum() - h.GetMinimum() > 50000):
                canvas.SetLogy()

            ##ploting histogram
            h.DrawCopy('E1')
            canvas['risetime'][k].Update()
            canvas['risetime'][k].SaveAs(out_dir+'/risetime_ch'+str(k)+'.png')


        '''=========================== Waveform color chart ==========================='''
        if ('WaveColor' in configurations.plots):
            name = 'h_wave_'+str(k)
            title = 'Waveform color chart channel '+str(k)

            h = rt.TH2D(name, title, 250, 0, 210, 250, -600, 600)
            h.SetXTitle('Time [ns]')
            h.SetYTitle('Voltage [mV]')

            chain.Project(name, 'channel[{}]:time[{}]'.format(k,conf['idx_time']))

            h.SetStats(0)
            canvas['wave'][k] = rt.TCanvas('c_wave_'+str(k), 'c_wave_'+str(k), 800, 600)
            if(h.GetMaximum() - h.GetMinimum() > 5000):
                canvas['wave'][k].SetLogz()
            h.DrawCopy('colz')
            canvas['wave'][k].Update()
            canvas['wave'][k].SaveAs(out_dir + '/Waveform_ch'+str(k)+'.png')


        '''=========================== Track position ==========================='''
        if conf['idx_dut'] >= 0:
            var = 'y_dut[{}]:x_dut[{}]'.format(conf['idx_dut'], conf['idx_dut'])
            dy = configurations.xy_center[1]
            dx = configurations.xy_center[0]
            width = configurations.xy_center[2]
            if 'PosRaw' in configurations.plots:
                name = 'h_pos_'+str(k)
                title = 'Track position at z_DUT[' + str(conf['idx_dut']) + '], for channel amp['+str(k)+'] in {'
                title += str(conf['amp_range'][0]) + ', ' + str(conf['amp_range'][1]) + '} mV'

                h = rt.TH2D(name, title, 100, -width+dx, width+dx, 100, -width+dy, width+dy)
                h.SetXTitle('x [mm]')
                h.SetYTitle('y [mm]')

                chain.Project(name, var, conf['amp_sel'])
                # chain.Project(name, var, 'amp[{}]>200'.format(k))


                canvas['pos'][k] = rt.TCanvas('c_pos_'+str(k), 'c_pos_'+str(k), 800, 600)
                h.DrawCopy('colz')

                canvas['pos'][k].Update()
                canvas['pos'][k].SaveAs(out_dir + '/PositionXY_raw_ch'+str(k)+'.png')

            if 'PosWeight' in configurations.plots:
                name = 'h_weight_pos_'+str(k)
                title = 'Track position at z_DUT[' + str(conf['idx_dut']) + '] weighted with channel amp['+str(k)+'] in {'
                title += str(conf['amp_range'][0]) + ', ' + str(conf['amp_range'][1]) + '} mV'
                h_w = rt.TH2D(name, title, 100, -width+dx, width+dx, 100, -width+dy, width+dy)
                h_w.SetXTitle('x [mm]')
                h_w.SetYTitle('y [mm]')

                weights = '('+ conf['amp_sel'] +') * amp[' + str(k) + ']'
                chain.Project(name, var, weights)

                if 'PosRaw' in configurations.plots:
                    h_w.Divide(h)

                canvas['w_pos'][k] = rt.TCanvas('c_w_pos_'+str(k), 'c_w_pos_'+str(k), 800, 600)
                h_w.DrawCopy('colz')
                canvas['w_pos'][k].Update()
                canvas['w_pos'][k].SaveAs(out_dir + '/PositionXY_amp_weight_ch'+str(k)+'.png')


        '''=========================== Raw time resolution ==========================='''
        if conf['idx_ref'] >= 0:
            time_var_chref = configurations.channel[conf['idx_ref']]['var_ref']+'[{}]'.format(conf['idx_ref'])
            time_var = conf['var_ref']+'[{}]'.format(k)
            var = time_var + ' - ' + time_var_chref

            selection = '('+ conf['amp_sel'] +') * (' + configurations.channel[conf['idx_ref']]['amp_sel'] +')'

            delta_t = (tree2array(chain, var, selection).view(np.recarray).T)[0]
            # print data

            name = 'h_delta_t_raw_'+str(k)
            title = 'Time resolution for channel '+str(k)

            x_axis_title = var + '  [ns]'

            if ( len(delta_t) ==0):
                print 'Empty delta'
                continue
            median = np.percentile(delta_t, 50)
            width = np.abs(np.percentile(delta_t, 10) - np.percentile(delta_t, 90))

            h = create_TH1D(delta_t, name, title,
                                binning = [ None, median-3*width, median+3*width],
                                axis_title = [x_axis_title, 'Events'])

            low_edge = min(h.GetBinCenter(h.GetMaximumBin()-5), np.percentile(delta_t, 10))
            upper_edge = min(h.GetBinCenter(h.GetMaximumBin()+8), np.percentile(delta_t, 95))

            res = h.Fit('gaus', 'LQSR', '', low_edge, upper_edge)

            if 'TimeResRaw' in configurations.plots:
                canvas['t_res_raw'][k] = rt.TCanvas('c_t_res_raw_'+str(k), 'c_t_res_raw_'+str(k), 800, 600)
                h.DrawCopy('E1')
                canvas['t_res_raw'][k].Update()
                canvas['t_res_raw'][k].SaveAs(out_dir + '/TimeResolution_raw_ch'+str(k)+'.png')

            # '''=========================== Time resolution vs impact point ==========================='''
            # if conf['idx_dut'] >= 0 and 'ImpactCorrection' in configurations.plots:
            #     i_s = conf['idx_dut']
            #
            #     selection += selection+' && ntracks>0'
            #     bounds = {}
            #     h = {}
            #
            #     for c in ['x','y']:
            #         cc = (tree2array(chain, c+'_dut['+str(i_s)+']', selection).view(np.recarray).T)[0]
            #         bounds[c] = [np.percentile(cc, 5), np.percentile(cc, 95)]
            #         var = ''
            #         h[c] = TH2D('h_'+c+'_space_dip', 'h_'+c+'_space_dip', 50, bounds[c][0], bounds[c][1], 50, low_edge, upper_edge)




    # Compile the php index.php file
    current_dir = os.getcwd()
    os.chdir(out_dir)
    os.system('php ./index.php > web.html')
    os.chdir(current_dir)
