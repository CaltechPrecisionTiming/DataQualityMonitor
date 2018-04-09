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
    canvas['space_corr'] = {}
    canvas['t_res_space'] = {}
    canvas['dt_vs_amp'] = {}
    canvas['dtcorr_vs_amp'] = {}
    canvas['dt_corr'] = {}

    rt.gStyle.SetStatY(0.98);
    rt.gStyle.SetStatX(1.);
    rt.gStyle.SetStatW(0.2);
    rt.gStyle.SetStatH(0.1);
    rt.gStyle.SetHistLineWidth(2);


    for k in configurations.ch_ordered:
        conf = configurations.channel[k]

        '''=========================== Amplitude ==========================='''
        name = 'h_amp_'+str(k)
        title = 'Amplitude channel '+str(k)

        amp_aux = (tree2array(chain, 'amp['+str(k)+']').view(np.recarray).T)[0]
        h = rt.TH1D(name, title, 80, min(float(conf['amp_min']), np.percentile(amp_aux, 1)), max(np.percentile(amp_aux, 99), float(conf['amp_max'])))
        h.SetXTitle('Peak amplitude [mV]')
        h.SetYTitle('Events / {:.1f} mV'.format(h.GetBinWidth(1)))
        chain.Project(name, 'amp['+str(k)+']')

        h.GetXaxis().SetRange(int(float(conf['amp_min'])/h.GetBinWidth(1))+1,int(float(conf['amp_max'])/h.GetBinWidth(1))+1)
        i_max = h.GetMaximumBin()
        h.GetXaxis().SetRange(1,200)
        peak = h.GetBinCenter(i_max)
        # conf['amp_range'] = [max(float(conf['amp_min']),peak*0.6), min(float(conf['amp_max']), peak*1.7)+20]
        conf['amp_range'] = [float(conf['amp_min']), float(conf['amp_max'])]
        conf['amp_sel'] = '(amp['+str(k)+'] < ' + str(conf['amp_range'][1])
        conf['amp_sel'] += ' && '
        conf['amp_sel'] += 'amp['+str(k)+'] > ' + str(conf['amp_range'][0]) + ')'

        res = h.Fit('landau','LQSR', '', conf['amp_range'][0], conf['amp_range'][1])

        if 'Amp' in configurations.plots:
            canvas['amp'][k] = rt.TCanvas('c_amp_'+str(k), 'c_amp_'+str(k), 800, 600)
            # if(h.GetMaximum() - h.GetMinimum() > 5000):
            canvas['amp'][k].SetLogy()
            h.DrawCopy('LE')

            line = rt.TLine()
            line.SetLineColor(6)
            line.SetLineWidth(2)
            line.SetLineStyle(10)
            line.DrawLine(conf['amp_range'][0], 0, conf['amp_range'][0], h.GetMaximum())
            line.DrawLine(conf['amp_range'][1], 0, conf['amp_range'][1], h.GetMaximum())

            canvas['amp'][k].Update()
            canvas['amp'][k].SaveAs(out_dir + '/Amp_ch'+str(k)+'.png')


        '''=========================== Integral ==========================='''
        if 'Int' in configurations.plots:
            canvas['int'][k] = rt.TCanvas('c_int_'+str(k), 'c_int_'+str(k), 800, 600)

            name = 'h_int_'+str(k)
            title = 'Integral channel '+str(k)
            int_aux = (tree2array(chain, '-integral['+str(k)+']').view(np.recarray).T)[0]
            h = rt.TH1D(name, title, 100, np.percentile(int_aux, 5), np.percentile(int_aux, 99.9))
            h.SetXTitle('Integral [pC]')
            h.SetYTitle('Events / {:.1f} pC'.format(h.GetBinWidth(1)))
            chain.Project(name, '-integral['+str(k)+']')

            if(h.GetMaximum() - h.GetMinimum() > 5000):
                canvas['int'][k].SetLogy()
            h.DrawCopy('LE')
            canvas['int'][k].Update()
            canvas['int'][k].SaveAs(out_dir + '/Int_ch'+str(k)+'.png')

        '''=========================== Risetime ======================================='''
        if 'Risetime' in configurations.plots:
            canvas['risetime'][k] = rt.TCanvas('c_risetime_'+str(k), 'c_int_'+str(k), 800, 600)

            nbins  = 100
            name   = 'h_risetime_' + str(k)
            title  = 'Risetime ch' + str(k)
            rise_aux = (tree2array(chain, 'risetime['+str(k)+']').view(np.recarray).T)[0]
            bin_w = 200.0/1024
            lower_bound = np.percentile(rise_aux, 5) - bin_w*0.5
            upper_bound = np.percentile(rise_aux, 99) + bin_w*0.5
            N_bins = (upper_bound - lower_bound) / bin_w
            h      = rt.TH1D(name, title, int(N_bins), lower_bound, upper_bound)
            h.SetXTitle('Risetime [ns]')
            h.SetYTitle('Events / {:.1f} ns'.format(h.GetBinWidth(1)))
            if conf['idx_ref'] == -1:#do not apply cut again on reference channels
                cut = conf['amp_sel']
            else:
                cut = conf['amp_sel'] + '&&' + configurations.channel[conf['idx_ref']]['amp_sel']

            chain.Project(name, 'risetime['+str(k)+']', cut)

            if (h.GetMaximum() - h.GetMinimum() > 50000):
                canvas['risetime'][k].SetLogy()

            ##ploting histogram
            h.DrawCopy('lE')
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
            var_dT = time_var + ' - ' + time_var_chref

            selection = conf['amp_sel'] +' && ' + configurations.channel[conf['idx_ref']]['amp_sel']

            delta_t = (tree2array(chain, var_dT, selection).view(np.recarray).T)[0]
            if ( len(delta_t) ==0):
                print 'Empty delta'
                continue
            median = np.percentile(delta_t, 50)
            width = np.abs(np.percentile(delta_t, 10) - np.percentile(delta_t, 90))

            name = 'h_delta_t_raw_'+str(k)
            title = 'Time resolution for channel '+str(k)+', width 10-90 = {:.2f} ns'.format(width)
            h = create_TH1D(delta_t, name, title,
                                binning = [ None, median-2*width, median+2*width],
                                axis_title = [var_dT + ' [ns]', 'Events'])

            if 'TimeResRaw' in configurations.plots:
                canvas['t_res_raw'][k] = rt.TCanvas('c_t_res_raw_'+str(k), 'c_t_res_raw_'+str(k), 800, 600)
                h.Fit('gaus', 'LQR','', median-width, median+width)
                h = h.DrawCopy('LE')
                line = rt.TLine()
                line.SetLineColor(6)
                line.SetLineWidth(2)
                line.SetLineStyle(7)
                line.DrawLine(median-width, 0, median-width, h.GetMaximum())
                line.DrawLine(median+width, 0, median+width, h.GetMaximum())
                canvas['t_res_raw'][k].Update()
                canvas['t_res_raw'][k].SaveAs(out_dir + '/TimeResolution_raw_ch'+str(k)+'.png')

            selection += '&& ({}>{} && {}<{})'.format(var_dT, median-width, var_dT, median+width)
            delta_t = (tree2array(chain, var_dT, selection).view(np.recarray).T)[0]
            arr = {}
            '''=========================== Time resolution vs impact point ==========================='''
            if conf['idx_dut'] >= 0:
                i_s = conf['idx_dut']

                name = 'c_space_corr'+str(k)
                canvas['space_corr'][k] = rt.TCanvas(name, name, 1000, 600)
                canvas['space_corr'][k].Divide(2)
                add_sel = ''

                selection += ' && ntracks > 0'
                delta_t = (tree2array(chain, var_dT, selection).view(np.recarray).T)[0]

                for i, c in enumerate(['x', 'y']):
                    conf[c] = {}
                    pos = (tree2array(chain, c+'_dut[0]', selection).view(np.recarray).T)[0]
                    conf[c]['pl'] = np.percentile(pos, 10)
                    conf[c]['ph'] = np.percentile(pos, 90)
                    add_sel += ' && '+c+'_dut['+str(i_s)+'] > ' + str(conf[c]['pl']) + ' && '+c+'_dut['+str(i_s)+'] < ' + str(conf[c]['ph'])

                    canvas['space_corr'][k].cd(i+1)
                    h = create_TH2D(np.column_stack((pos, delta_t)), name=c+name, title='Time resolution '+c+' dependence',
                                    binning=[50, conf[c]['pl']-4, conf[c]['ph']+4, 50, median-2*width, median+2*width],
                                    axis_title=[c+' [mm]', '#DeltaT [ns]']
                                   )

                    h.DrawCopy('COLZ')

                    prof = h.ProfileX('prof_'+c )
                    prof.SetLineColor(2)
                    prof.SetLineWidth(2)

                    f = rt.TF1(c+'_fit','[0]+[1]*x+[2]*x^2',conf[c]['pl'], conf[c]['ph'])

                    aux_t = delta_t[np.logical_and(pos>conf[c]['pl'], pos<conf[c]['ph'])]
                    pos = pos[np.logical_and(pos>conf[c]['pl'], pos<conf[c]['ph'])]
                    coeff, r, rank, s = np.linalg.lstsq(np.column_stack((0*pos+1, pos, pos**2)), aux_t)
                    for j,a in enumerate(coeff):
                        f.SetParameter(j, a)
                    conf[c]['coeff'] = np.flip(coeff, 0)
                    f.SetLineColor(8)
                    f.DrawCopy('SAMEL')

                    prof.DrawCopy('SAMEE')

                canvas['space_corr'][k].Update()
                canvas['space_corr'][k].SaveAs(out_dir + '/TimeResolution_Position_dependece_ch'+str(k)+'.png')


                line = rt.TLine()
                line.SetLineColor(6)
                line.SetLineStyle(7)
                line.SetLineWidth(3)
                for can in [canvas['pos'][k], canvas['w_pos'][k]]:
                    can.cd()
                    line.DrawLine(conf['x']['pl'], conf['y']['pl'], conf['x']['pl'], conf['y']['ph'])
                    line.DrawLine(conf['x']['ph'], conf['y']['pl'], conf['x']['ph'], conf['y']['ph'])
                    line.DrawLine(conf['x']['pl'], conf['y']['pl'], conf['x']['ph'], conf['y']['pl'])
                    line.DrawLine(conf['x']['pl'], conf['y']['ph'], conf['x']['ph'], conf['y']['ph'])
                canvas['w_pos'][k].SaveAs(out_dir + '/PositionXY_amp_weight_ch'+str(k)+'.png')
                canvas['pos'][k].SaveAs(out_dir + '/PositionXY_raw_ch'+str(k)+'.png')

                selection += add_sel
                arr['x'] = (tree2array(chain, 'x_dut[0]', selection).view(np.recarray).T)[0]
                arr['y'] = (tree2array(chain, 'y_dut[0]', selection).view(np.recarray).T)[0]
                delta_t = (tree2array(chain, var_dT, selection).view(np.recarray).T)[0]

                dt_space_corrected = np.copy(delta_t)
                for c in ['x', 'y']:
                    dt_space_corrected -= np.polyval(conf[c]['coeff'], arr[c])

                h = create_TH1D(dt_space_corrected, 'h_delta_space_corr'+str(k), 'Time resolution space corrected ch '+str(k),
                                binning = [ None, np.min(dt_space_corrected), np.max(dt_space_corrected)],
                                axis_title = [var_dT+' [ns]', 'Events'])

                canvas['t_res_space'][k] = rt.TCanvas('c_t_res_space'+str(k), 'c_t_res_raw'+str(k), 700, 500)
                h.DrawCopy('E1')
                canvas['t_res_space'][k].Update()
                canvas['t_res_space'][k].SaveAs(out_dir + '/TimeResolution_space_ch'+str(k)+'.png')

            '''=========================== Time resolution vs amplitude ==========================='''
            conf['amp'] = {}
            arr['amp'] = (tree2array(chain, 'amp['+str(k)+']', selection).view(np.recarray).T)[0]
            canvas['dt_vs_amp'][k] = rt.TCanvas('dt_vs_amp'+str(k), 'dt_vs_amp'+str(k), 1200, 600)
            canvas['dt_vs_amp'][k].Divide(2)

            h = create_TH2D(np.column_stack((arr['amp'], delta_t)), name='h_amp_dip', title='h_amp_dip',
                                binning=[50, conf['amp_range'][0], conf['amp_range'][1], 50, np.min(delta_t), np.max(delta_t)],
                                axis_title=['Amp [mV]', '#DeltaT [ns]']
                               )
            canvas['dt_vs_amp'][k].cd(1)
            h.DrawCopy('colz')
            prof = h.ProfileX('prof_amp')
            prof.SetLineColor(6)
            prof.SetLineWidth(2)
            prof.DrawCopy('SAMEE1')

            f = rt.TF1('amp_fit'+str(k),'[0]+[1]*x+[2]*x^2', conf['amp_range'][0], conf['amp_range'][1])
            f.DrawCopy('SAMEL')

            coeff, r, rank, s = np.linalg.lstsq(np.column_stack((0*arr['amp']+1, arr['amp'], arr['amp']**2)), delta_t)
            for j,a in enumerate(coeff):
                f.SetParameter(j, a)
            conf['amp']['coeff'] = np.flip(coeff, 0)
            f.SetLineColor(8)
            f.DrawCopy('SAMEL')


            dt_amp_corrected = np.copy(delta_t) - np.polyval(conf['amp']['coeff'], arr['amp'])

            h = create_TH1D(dt_amp_corrected, 'h_delta_amp_corr'+str(k), 'Time resolution amp corrected',
                            binning = [ None, np.min(dt_amp_corrected), np.max(dt_amp_corrected)],
                            axis_title = [var_dT+' [ns]', 'Events'])
            canvas['dt_vs_amp'][k].cd(2)
            h.Fit('gaus', 'LQR','', np.percentile(dt_amp_corrected, 1), np.percentile(dt_amp_corrected, 99))
            h.DrawCopy('E1')

            canvas['dt_vs_amp'][k].Update()
            canvas['dt_vs_amp'][k].SaveAs(out_dir + '/TimeResolution_amp_ch'+str(k)+'.png')

            '''=========================== Time resolution vs amplitude and space ==========================='''
            if conf['idx_dut'] >= 0:
                canvas['dtcorr_vs_amp'][k] = rt.TCanvas('dtcorr_vs_amp'+str(k), 'dtcorr_vs_amp'+str(k), 1000, 500)
                canvas['dtcorr_vs_amp'][k].Divide(2)

                h = create_TH2D(np.column_stack((arr['amp'], dt_space_corrected)), name='h_space_amp_dip'+str(k), title='h_space_amp_dip'+str(k),
                                    binning=[50, conf['amp_range'][0], conf['amp_range'][1], 50, np.min(dt_space_corrected), np.max(dt_space_corrected)],
                                    axis_title=['Amp [mV]', '#DeltaT [ns]']
                                   )
                canvas['dtcorr_vs_amp'][k].cd(1)
                h.DrawCopy('colz')
                prof = h.ProfileX('prof_corr_amp')
                prof.SetLineColor(6)
                prof.SetLineWidth(2)
                prof.DrawCopy('SAMEE1')

                f = rt.TF1('f_tcorr_amp_fit'+str(k),'[0]+[1]*x+[2]*x^2', conf['amp_range'][0], conf['amp_range'][1])
                coeff, r, rank, s = np.linalg.lstsq(np.column_stack((0*arr['amp']+1, arr['amp'], arr['amp']**2)), delta_t)
                for j,a in enumerate(coeff):
                    f.SetParameter(j, a)
                conf['amp']['coeff_space_corr'] = np.flip(coeff, 0)
                f.SetLineColor(8)
                f.DrawCopy('SAMEL')


                dt_amp_corrected = np.copy(dt_space_corrected) - np.polyval(conf['amp']['coeff_space_corr'], arr['amp'])

                h = create_TH1D(dt_amp_corrected, 'h_delta_amp_corr', 'Time resolution amp corrected',
                                binning = [ None, np.min(dt_amp_corrected), np.max(dt_amp_corrected)],
                                axis_title = [var_dT+' [ns]', 'Events'])
                f = rt.TF1('f_all_corr'+str(k),'gaus', np.percentile(dt_amp_corrected, 5), np.percentile(dt_amp_corrected, 95))
                canvas['dtcorr_vs_amp'][k].cd(2)
                h.Fit('f_all_corr'+str(k), 'LQR+')
                h.DrawCopy('E1')
                f.SetLineColor(2)
                f.DrawCopy('SAME L')

                canvas['dtcorr_vs_amp'][k].Update()
                canvas['dtcorr_vs_amp'][k].SaveAs(out_dir + '/TimeResolution_SpaceAmp_ch'+str(k)+'.png')

                '''=========================== Time resolution w/ one-shot corrections ==========================='''
                def create_regression_input(x, y, amp):
                    out = (np.ones_like(x), x, y, amp, x**2, y**2, amp**2, x*y, amp*x, amp*y)
                    return np.column_stack(out)

                inputs = create_regression_input(arr['x'], arr['y'], arr['amp'])
                coeff, r, rank, s = np.linalg.lstsq(inputs, delta_t)

                dt_corr = delta_t - np.dot(inputs, coeff)

                h = create_TH1D(dt_corr, 'h_dt_corr'+str(k), 'Time resolution one-shot correction',
                                binning = [ None, np.min(dt_corr), np.max(dt_corr)],
                                axis_title = [var_dT+' [ns]', 'Events'])

                f = rt.TF1('f_corr'+str(k),'gaus', np.percentile(dt_corr, 3), np.percentile(dt_corr, 97))
                h.Fit('f_corr'+str(k), 'LQR+')
                canvas['dt_corr'][k] = rt.TCanvas('c_dt_corr'+str(k), 'c_dt_corr'+str(k), 800, 600)
                h.DrawCopy('E1')
                f.SetLineColor(2)
                f.DrawCopy('SAMEL')

                canvas['dt_corr'][k].Update()
                canvas['dt_corr'][k].SaveAs(out_dir + '/TimeResolution_OneShot_ch'+str(k)+'.png')




    # Compile the php index.php file
    current_dir = os.getcwd()
    os.chdir(out_dir)
    os.system('php ./index.php > web.html')
    os.chdir(current_dir)
