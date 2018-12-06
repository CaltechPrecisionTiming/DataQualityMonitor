import numpy as np
import os, re, shutil
import argparse
import itertools
from prettytable import PrettyTable

import ROOT as rt
from root_numpy import tree2array, tree2rec
from lib.histo_utilities import create_TH1D, create_TH2D, quantile
from lib.cebefo_style import cebefo_style, Set_2D_colz_graphics

donotdelete = []



def parsing():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str, help="input root file, if -N is given XX is replaced with runNumber", nargs='+')
    parser.add_argument("-C", "--config", type=str, default='config/FNAL_TB_1811/VME_SiPM.txt', help="Config file")
    parser.add_argument("-S", "--save_loc", type=str, default='./out_plots/', help="Saving location")

    parser.add_argument("-B", "--batch", default=True, action='store_false', help="Root batch mode")

    parser.add_argument("-N", "--runs_interval", default=None, help="Runs to run", nargs='+')

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
            l = l[0:-1].split()
            if '-->Print' in l[0]:
                self.plots = l[1:]
            elif '-->XYcenter' in l[0]:
                self.xy_center = [float(l[1]), float(l[2]), float(l[3])]
            elif '-->TracksCleaning' in l[0]:
                self.TracksCleaning = {}
                cuts = ''
                for aux in l[1:]:
                    if aux[:7] == 'CutOnCh':
                        self.TracksCleaning['ch'] = int(aux[7:])
                        print 'Cuts from ch', int(aux[7:]), 'will be used for cleaning tracks'
                    else:
                        cuts += aux + ' && '
                if len(cuts) > 3:
                    self.TracksCleaning['cuts'] = cuts[:-4]
                    print 'Cleaning tracks with:', cuts[:-4]
                else:
                    self.TracksCleaning['cuts'] = ''
            elif '-->TracksConsistency' in l[0]:
                self.TracksConsistency = {}
                self.TracksConsistency['Bad'] = 0
                self.TracksConsistency['Checked'] = 0
                for aux in l[1:]:
                    if aux[:5] == 'RefCh':
                        self.TracksConsistency['RefCh'] = int(aux[5:])
                        # print 'Tracks consistency: ch', int(aux[5:]), 'will be checked'

                    if aux == 'List':
                        self.TracksConsistency['List'] = 1
            elif '-->Var' in l[0]:
                if not hasattr(self, 'Var'):
                    self.Var = {}
                self.Var[l[0][7:-1]] = l[1:]
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

class Bauble:
    def __init__(self):
        pass

def define_range_around_peak(h, perc = [0.4, 0.4], Range=[0.0, 99999.0]):
    h_aux = h.Clone('h_aux')
    for i in range(1, h_aux.GetNbinsX()+1):
        if h_aux.GetBinCenter(i) < Range[0] or h_aux.GetBinCenter(i) > Range[1]:
            h_aux.SetBinContent(i, 0.)

    SS = rt.TSpectrum()
    n_pks = SS.Search(h_aux, 2, "", 0.2)

    i_max = 0
    max = -1
    y = SS.GetPositionY()
    x = SS.GetPositionX()
    found = False
    for i in range(n_pks):
        if y[i] > max and x > Range[0] and x < Range[1]:
            max = y[i]
            i_max = i
            found = True

    if found:
        n_pk = h.FindBin(SS.GetPositionX()[i_max])
    else:
        n_pk = h_aux.GetMaximumBin()

    thr = perc[0] * h.GetBinContent(n_pk)
    n_low = n_pk
    while h.GetBinContent(n_low) > thr and n_low > 1:
        n_low -= 1
    x_low = h.GetBinLowEdge(n_low)

    thr = perc[1] * h.GetBinContent(n_pk)
    n_up = n_pk
    while h.GetBinContent(n_up) > thr and n_up < h.GetNbinsX():
        n_up += 1
    x_up = h.GetBinLowEdge(n_up) + h.GetBinWidth(n_up)
    return x_low, x_up, n_pk

def rootTH2_to_np(h, cut = None, Norm = False):
    nx = h.GetNbinsX()
    ny = h.GetNbinsY()

    arr = np.zeros((ny, nx))
    pos = np.zeros((ny, nx, 2))

    for ix in range(nx):
        for iy in range(ny):
            x = h.GetXaxis().GetBinCenter( ix+1 );
            y = h.GetYaxis().GetBinCenter( iy+1 );
            z = h.GetBinContent(h.GetBin(ix+1, iy+1))

            if cut == None:
                arr[iy, ix] = z
            else:
                arr[iy, ix] = z if z > cut else 0
            pos[iy, ix] = [x,y]
    return arr, pos

def circle_filter(h_arr, cx, cy, rx, ry):
    p = 0
    for i,j in itertools.product(np.arange(h_arr.shape[0]), np.arange(h_arr.shape[1])):
        if (float(cx-j)/rx)**2 + (float(cy-i)/ry)**2 < 1:
            p += h_arr[i-1, j-1]

    return p/(np.pi*rx*ry)

def fill_TimeResHisto(k, dt, h2D, out_list, tagin, title_tag, canvas):
    q_up, e_up = quantile(1000*dt, 0.15)
    q_dwn, e_dwn = quantile(1000*dt, 0.85)
    disp_est = 0.5*np.abs(q_up - q_dwn)
    disp_unc = 0.5*np.hypot(e_up, e_dwn)

    out_list.append([disp_est, disp_unc])
    ix = h2D.GetXaxis().FindBin(0.5*(xu+xd))
    iy = h2D.GetYaxis().FindBin(0.5*(yu+yd))
    ig = h2D.GetBin(ix, iy)
    h2D.SetBinContent(ig, disp_est)
    h2D.SetBinError(ig, disp_unc)

    median = np.percentile(dt, 50)
    width = np.abs(np.percentile(dt, 10) - np.percentile(dt, 90))
    tag = tagin+'_'+str(k)+'_'+str(ib)
    h = create_TH1D(dt, 'h_'+tag,
                    title_tag + ' time resolution - x #in [{:.1f}, {:.1f}], y #in [{:.1f}, {:.1f}]'.format(xd, xu, yd, yu),
                    binning = [ None, median-2*width, median+2*width],
                    axis_title = [var_dT + ' [ns]', 'Events'])
    canvas['c_'+tag] = rt.TCanvas('c_'+tag, 'c_'+tag, 800, 600)
    h.Fit('gaus', 'LQR','', np.percentile(dt, 20), np.percentile(dt, 98))
    h = h.DrawCopy('LE')

    l = rt.TLatex()
    l.SetTextSize(0.04);
    l.DrawLatexNDC(0.15, 0.85, 'Unbiased dispersion: {:.1f} #pm {:.1f}'.format(disp_est, disp_unc))
    canvas['c_'+tag].Update()
    canvas['c_'+tag].SaveAs(out_dir + '/'+tagin+'_ch{:02d}/'.format(k)+tag+'.png')

if __name__ == '__main__':
    cebefo_style()

    rt.gErrorIgnoreLevel = 6000

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
    elif len(args.runs_interval)==1 and not args.runs_interval[0].endswith('.txt'):
        N = args.runs_interval[0]
        args.input_file = [args.input_file[0].replace('XXX', str(N))]
        flag = '_'+str(N)
    else:
        if args.runs_interval[0].endswith('.txt'):
            runs = np.sort(np.loadtxt(args.runs_interval[0]).astype(np.int))
            runs = list(runs)
        elif len(args.runs_interval)==2 and args.runs_interval[0] < args.runs_interval[1]:
                runs = range(int(args.runs_interval[0]), int(args.runs_interval[1])+1)
        else:
            runs = args.runs_interval

        file_template = args.input_file[0]
        args.input_file = map(lambda x: file_template.replace('XXX', str(x)), runs)
        if len(runs) > 1:
            flag = '_' + str(runs[0]) + '-' + str(runs[-1])
        else:
            flag = '_' + str(runs[0])

    save_loc = args.save_loc
    if os.path.isdir(save_loc):
        if save_loc[-1] != '/':
            save_loc += '/'
        out_dir = save_loc + 'Run' + str(flag) + '_' + os.path.basename(args.config[:-4])
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.mkdir(out_dir)
        shutil.copy('lib/index.php', out_dir+'/index.php')
    else:
        print 'Save location not existing:', save_loc
        raise Exception(' ')


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
    canvas['pos_sel'] = {}
    canvas['w_pos'] = {}
    canvas['t_res_raw'] = {}
    canvas['space_corr'] = {}
    canvas['t_res_space'] = {}
    canvas['dt_vs_amp'] = {}
    canvas['dtcorr_vs_amp'] = {}
    canvas['dt_corr'] = {}

    rt.gStyle.SetStatY(0.95)
    rt.gStyle.SetStatX(0.85)
    rt.gStyle.SetStatW(0.15)
    rt.gStyle.SetStatH(0.1)
    rt.gStyle.SetHistLineWidth(2)

    '''=========================== AllTracks ==========================='''
    if 'AllTracks' in configurations.plots:
        name = 'h_alltracks'
        title = 'All tracks at z_DUT[0]'

        var = 'y_dut[0]:x_dut[0]'
        dy = configurations.xy_center[1]
        dx = configurations.xy_center[0]
        width = configurations.xy_center[2]

        h = rt.TH2D(name, title, 100, -width+dx, width+dx, 100, -width+dy, width+dy)
        h.SetXTitle('x [mm]')
        h.SetYTitle('y [mm]')

        sel = ''
        if hasattr(configurations, 'TracksCleaning'):
            sel += configurations.TracksCleaning['cuts']
        chain.Project(name, var, sel)

        canvas['AllTracks'] = rt.TCanvas('c_allTracks', 'c_allTracks', 800, 600)
        h.DrawCopy('colz')

        canvas['AllTracks'].Update()
        canvas['AllTracks'].SaveAs(out_dir + '/AllTracks.png')

    print '\n======================= Channels selection loop ==========================='
    for k in configurations.ch_ordered:
        print '---> Channel', k
        conf = configurations.channel[k]

        line = rt.TLine()
        line.SetLineColor(6)
        line.SetLineWidth(2)
        line.SetLineStyle(10)

        selection = '1'

        '''=========================== Amplitude ==========================='''
        if 'Amp' in configurations.plots:
            name = 'h_amp_'+str(k)
            title = 'Amplitude channel '+str(k)

            amp_aux = np.concatenate(list(tree2array( chain, 'amp['+str(k)+']')))

            h = rt.TH1D(name, title, 80, np.percentile(amp_aux, 1), min(990., np.percentile(amp_aux, 99.99)))
            h.SetXTitle('Peak amplitude [mV]')
            h.SetYTitle('Events / {:.1f} mV'.format(h.GetBinWidth(1)))
            chain.Project(name, 'amp['+str(k)+']')

            Range = [0.0, 9999999.0]
            if 'cut' in conf.keys():
                if conf['cut'][0] in ['a', 'A']:
                    if 'min' in conf.keys():
                        Range[0] = float(conf['min'])
                    if 'max' in conf.keys():
                        Range[1] = float(conf['max'])

            x_low, x_up, n_pk = define_range_around_peak(h, [0.2, 0.2], Range)

            canvas['amp'][k] = rt.TCanvas('c_amp_'+str(k), 'c_amp_'+str(k), 800, 600)
            h.DrawCopy('E')

            gr = rt.TGraph(1)
            gr.SetPoint(0, h.GetBinCenter(n_pk), h.GetBinContent(n_pk))
            gr.SetMarkerStyle(23)
            gr.SetMarkerColor(2)
            gr.Draw("P")

            if h.GetMaximum() > 5*h.GetBinContent(n_pk):
                canvas['amp'][k].SetLogy()

            if 'cut' in conf.keys():
                if conf['cut'][0] in ['a', 'A']:
                    selection += ' && (amp[{}]>{} && amp[{}]<{})'.format(k, x_low, k, x_up)
                    line.DrawLine(x_low, 0, x_low, h.GetMaximum())
                    line.DrawLine(x_up, 0, x_up, h.GetMaximum())

            canvas['amp'][k].Update()
            canvas['amp'][k].SaveAs(out_dir + '/Amp_ch{:02d}.png'.format(k))


        '''=========================== Integral ==========================='''
        if 'Int' in configurations.plots:
            canvas['int'][k] = rt.TCanvas('c_int_'+str(k), 'c_int_'+str(k), 800, 600)

            name = 'h_int_'+str(k)
            title = 'Integral channel '+str(k)
            int_aux = -np.concatenate(list(tree2array(chain, 'integral['+str(k)+']', 'integral['+str(k)+'] != 0')))
            h = rt.TH1D(name, title, 100, np.percentile(int_aux, 0.1), np.max(int_aux))
            h.SetXTitle('Integral [pC]')
            h.SetYTitle('Events / {:.1f} pC'.format(h.GetBinWidth(1)))
            chain.Project(name, '-integral['+str(k)+']', '-integral['+str(k)+'] != 0')

            Range = [0.0, 9999999.0]
            if 'cut' in conf.keys():
                if conf['cut'][0] in ['i', 'I']:
                    if 'min' in conf.keys():
                        Range[0] = float(conf['min'])
                    if 'max' in conf.keys():
                        Range[1] = float(conf['max'])

            x_low, x_up, n_pk = define_range_around_peak(h, [0.25, 0.3], Range)

            h.DrawCopy('E')

            gr = rt.TGraph(1)
            gr.SetPoint(0, h.GetBinCenter(n_pk), h.GetBinContent(n_pk))
            gr.SetMarkerStyle(23)
            gr.SetMarkerColor(2)
            gr.Draw("P")

            if(h.GetMaximum() - h.GetMinimum() > 5000):
                canvas['int'][k].SetLogy()

            if 'cut' in conf.keys():
                if conf['cut'][0] in ['i', 'I']:
                    selection += ' && (-integral[{}]>{} && -integral[{}]<{})'.format(k, x_low, k, x_up)
                    line.DrawLine(x_low, 0, x_low, h.GetMaximum())
                    line.DrawLine(x_up, 0, x_up, h.GetMaximum())

            canvas['int'][k].Update()
            canvas['int'][k].SaveAs(out_dir + '/Int_ch{:02d}.png'.format(k))

        '''=========================== Risetime ======================================='''
        if 'Risetime' in configurations.plots:
            print '\tRisetime'
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
                cut = selection
            else:
                cut = selection + '&&' + configurations.channel[conf['idx_ref']]['sel']

            chain.Project(name, 'risetime['+str(k)+']', cut)

            if (h.GetMaximum() - h.GetMinimum() > 50000):
                canvas['risetime'][k].SetLogy()

            ##ploting histogram
            h.DrawCopy('lE')
            canvas['risetime'][k].Update()
            canvas['risetime'][k].SaveAs(out_dir+'/risetime_ch{:02d}.png'.format(k))


        '''=========================== Waveform color chart ==========================='''
        if 'WaveColor' in configurations.plots:
            name = 'h_wave_'+str(k)
            title = 'Waveform color chart channel '+str(k)

            h = rt.TH2D(name, title, 250, 0, 210, 250, -1100, 500)
            h.SetXTitle('Time [ns]')
            h.SetYTitle('Voltage [mV]')

            N_max = chain.GetEntries()
            N_max = max(500, int(N_max*0.1))
            chain.Project(name, 'channel[{}]:time[{}]'.format(k,conf['idx_time'], 'Entry$ < {}'.format(N_max)))

            h.SetStats(0)
            canvas['wave'][k] = rt.TCanvas('c_wave_'+str(k), 'c_wave_'+str(k), 800, 600)
            if(h.GetMaximum() - h.GetMinimum() > 5000):
                canvas['wave'][k].SetLogz()
            h.DrawCopy('colz')
            canvas['wave'][k].Update()
            canvas['wave'][k].SaveAs(out_dir + '/Waveform_ch{:02d}.png'.format(k))

        '''=========================== Track position ==========================='''
        if conf['idx_dut'] >= 0:
            if hasattr(configurations, 'TracksCleaning'):
                selection += ' && ' + configurations.TracksCleaning['cuts']

                if 'ch' in configurations.TracksCleaning.keys():
                    ch = configurations.TracksCleaning['ch']
                    if ch < k:
                        selection += ' && ' + configurations.channel[ch]['sel']

            var = 'y_dut[{}]:x_dut[{}]'.format(conf['idx_dut'], conf['idx_dut'])
            dy = configurations.xy_center[1]
            dx = configurations.xy_center[0]
            width = configurations.xy_center[2]
            N_bins = [100, 100]
            if 'PosRaw' in configurations.plots:
                name = 'h_pos_'+str(k)
                title = 'Track position at z_DUT[' + str(conf['idx_dut']) + '], for channel {} selected event'.format(k)

                h = rt.TH2D(name, title, N_bins[0], -width+dx, width+dx, N_bins[1], -width+dy, width+dy)
                h.SetXTitle('x [mm]')
                h.SetYTitle('y [mm]')

                chain.Project(name, var, selection)

                canvas['pos'][k] = rt.TCanvas('c_pos_'+str(k), 'c_pos_'+str(k), 800, 600)
                h.DrawCopy('colz')
                canvas['pos'][k].Update()
                canvas['pos'][k].SaveAs(out_dir + '/PositionXY_raw_ch{:02d}.png'.format(k))

            if ('PosSel' in configurations.plots) or ('PosSel+' in configurations.plots):
                canvas['pos_sel'][k] = rt.TCanvas('c_pos_sel_'+str(k), 'c_pos_sel_'+str(k), 800, 600)

                name = 'h_pos_sel'
                title = 'Events average selection efficiency'
                h = rt.TH3D(name, title, N_bins[0], -width+dx, width+dx, N_bins[1], -width+dy, width+dy, 2, -0.5, 1.5)

                chain.Project(name, selection + ':' + var, configurations.TracksCleaning['cuts'])
                h = h.Project3DProfile('yx')
                h.SetYTitle('y [mm]')
                h.SetXTitle('x [mm]')
                h.DrawCopy('colz')

                if ('PosSel+' in configurations.plots) and ('shape' in conf.keys()):
                    if not conf['shape'] == 'None':
                        size = conf['shape'].split('x')

                        x = h.GetXaxis().GetBinCenter(1) - h.GetXaxis().GetBinWidth(1)*0.51 + float(size[0])
                        nbx = h.GetXaxis().FindBin(x)

                        y = h.GetYaxis().GetBinCenter(1) - h.GetYaxis().GetBinWidth(1)*0.51 + float(size[1])
                        nby = h.GetYaxis().FindBin(y)

                        arr_h, pos_h = rootTH2_to_np(h, cut=0.2)
                        # print arr_h

                        p_max = 0
                        idx_max = [0, 0]
                        for iy in range(0, arr_h.shape[0]-nby):
                            for ix in range(0, arr_h.shape[1]-nbx):
                                p = np.sum(arr_h[iy:iy+nby, ix:ix+nbx])

                                if p > p_max:
                                    p_max = p
                                    idx_max = [iy, ix]

                        avg_prob = p_max/(nbx*nby)
                        iy, ix = idx_max
                        x_start = h.GetXaxis().GetBinCenter(ix+1) - 0.5*h.GetXaxis().GetBinWidth(ix+1)
                        x_stop = h.GetXaxis().GetBinCenter(ix+nbx) + 0.5*h.GetXaxis().GetBinWidth(ix+nbx)
                        y_start = h.GetYaxis().GetBinCenter(iy+1) - 0.5*h.GetYaxis().GetBinWidth(iy+1)
                        y_stop = h.GetYaxis().GetBinCenter(iy+nby) + 0.5*h.GetYaxis().GetBinWidth(iy+nby)

                        best_theta = 0
                        cx = 0.5*(x_stop + x_start)
                        cy = 0.5*(y_stop + y_start)
                        x = '(x_dut[{}] - {})'.format(conf['idx_dut'], cx)
                        y = '(y_dut[{}] - {})'.format(conf['idx_dut'], cy)
                        var_template = '{c}*{y} - {s}*{x} + {cy}:{c}*{x} + {s}*{y} + {cx}'
                        for theta in np.arange(-0.1, 0.1, 0.005):
                            var_aux = var_template.format(x=x, y=y, c=np.cos(theta), s=np.sin(theta), cy=cy, cx=cx)

                            h = rt.TH3D('h_aux', title, N_bins[0], -width+dx, width+dx, N_bins[1], -width+dy, width+dy, 2, -0.5, 1.5)
                            chain.Project('h_aux', selection + ':' + var_aux, configurations.TracksCleaning['cuts'])
                            h = h.Project3DProfile('yx')
                            arr_h, pos_h = rootTH2_to_np(h, cut=0.2)
                            p = np.sum(arr_h[iy:iy+nby, ix:ix+nbx])
                            if p > p_max:
                                p_max = p
                                best_theta = theta
                        conf['theta_z'] = best_theta
                        conf['x_rot'] = '{c}*{x} + {s}*{y} + {cx}'.format(x=x, y=y, c=np.cos(best_theta), s=np.sin(best_theta), cx=cx)
                        conf['y_rot'] = '{c}*{y} - {s}*{x} + {cy}'.format(x=x, y=y, c=np.cos(best_theta), s=np.sin(best_theta), cy=cy)

                        var = conf['y_rot']+':'+conf['x_rot']
                        h = rt.TH3D('h_aux', title, N_bins[0], -width+dx, width+dx, N_bins[1], -width+dy, width+dy, 2, -0.5, 1.5)
                        chain.Project('h_aux', selection + ':' + var, configurations.TracksCleaning['cuts'])
                        h = h.Project3DProfile('yx')
                        h.SetYTitle('y [mm]')
                        h.SetXTitle('x [mm]')
                        h.DrawCopy('colz')

                        line.SetLineStyle(9)
                        line.SetLineColor(6)
                        line.SetLineWidth(2)
                        line.DrawLine(x_start, y_start, x_stop, y_start)
                        line.DrawLine(x_start, y_stop, x_stop, y_stop)
                        line.DrawLine(x_start, y_start, x_start, y_stop)
                        line.DrawLine(x_stop, y_start, x_stop, y_stop)

                        conf['space_sel'] = {}
                        conf['space_sel']['limits'] = [x_start, x_stop, y_start, y_stop]
                        conf['space_sel']['cut'] = '({y} < {h} && {y} > {l}'.format(y=conf['y_rot'], l=y_start, h=y_stop)
                        conf['space_sel']['cut'] += ' && {x} < {h} && {x} > {l})'.format(x=conf['x_rot'], l=x_start, h=x_stop)


                Set_2D_colz_graphics()
                canvas['pos_sel'][k].Update()
                canvas['pos_sel'][k].SaveAs(out_dir + '/PositionXY_sel_ch{:02d}.png'.format(k))



            if 'PosWeight' in configurations.plots:
                name = 'h_weight_pos_'+str(k)
                title = 'Track position at z_DUT[' + str(conf['idx_dut']) + '] weighted with channel {} selected event'.format(k)
                h_w = rt.TH2D(name, title, N_bins[0], -width+dx, width+dx, N_bins[1], -width+dy, width+dy)
                h_w.SetXTitle('x [mm]')
                h_w.SetYTitle('y [mm]')
                h_w.SetZTitle('Average Amplitude [mV]')

                h = rt.TH2D('h_amp_aux'+str(k), title, N_bins[0], -width+dx, width+dx, N_bins[1], -width+dy, width+dy)
                chain.Project('h_amp_aux'+str(k), var, selection)

                weights = '('+ selection +') * amp[' + str(k) + ']'
                chain.Project(name, var, weights)

                h_w.Divide(h)

                h_w.GetZaxis().SetRangeUser(h_w.GetMinimum(10), h_w.GetMaximum())

                canvas['w_pos'][k] = rt.TCanvas('c_w_pos_'+str(k), 'c_w_pos_'+str(k), 800, 600)
                h_w.DrawCopy('colz')
                Set_2D_colz_graphics()
                canvas['w_pos'][k].Update()
                canvas['w_pos'][k].SaveAs(out_dir + '/PositionXY_amp_weight_ch{:02d}.png'.format(k))

        '''=========================== End Selections ==========================='''
        conf['sel'] =  selection
        # print conf['sel']
        # if 'space_sel' in conf.keys():
            # print conf['space_sel']['cut']

    '''End of channels selection loop'''

    print '\n\n======================= Timing performances loop ==========================='

    headout_dir = out_dir
    best_result = {}
    for k in configurations.ch_ordered:
        conf = configurations.channel[k]
        if (not 'time' in conf['type']) or conf['idx_ref'] < 0:
            continue
        print '---> Channel', k

        results = []

        best_result[k] = Bauble()
        best_result[k].dT = [-1, -1]
        best_result[k].var = [None, None]
        best_result[k].AmpCorr = False

        ref_type = configurations.channel[conf['idx_ref']]['type']

        for v_ref, v_time in itertools.product(configurations.Var[ref_type], configurations.Var['time']):
            print 'Running on:', v_time, v_ref

            time_var_chref = v_ref+'[{}]'.format(conf['idx_ref'])
            time_var = v_time + '[{}]'.format(k)

            out_dir = '{}/{}_{}'.format(headout_dir, v_time, v_ref)
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
                shutil.copy('lib/index.php', out_dir + '/index.php')

            var_dT = time_var + ' - ' + v_ref+'[{}]'.format(conf['idx_ref'])

            selection = conf['sel'] +' && ' + configurations.channel[conf['idx_ref']]['sel'] + ' && {} != 0'.format(time_var) +' && ' + conf['space_sel']['cut']

            if 'TimeRes2DAmp' in configurations.plots:
                for kk in configurations.ch_ordered:
                    if 'amp'+str(k) in configurations.channel[kk]['type']:
                        selection += ' && ' + configurations.channel[kk]['sel']

            if chain.GetEntries(selection) < 5:
                print 'Not enought stat ({})'.format(chain.GetEntries(selection))
                continue

            '''=========================== Raw time resolution ==========================='''
            if 'TimeResRaw' in configurations.plots:
                delta_t =  tree2array(chain, var_dT, selection).flatten()
                if ( len(delta_t) ==0):
                    print 'Empty delta'
                    continue

                if os.uname()[1].startswith('lxplus'):
                    median = np.percentile(delta_t, 50)[0]
                    width = np.abs(np.percentile(delta_t, 10) - np.percentile(delta_t, 90))[0]
                else:
                    median = np.percentile(delta_t, 50)
                    width = np.abs(np.percentile(delta_t, 10) - np.percentile(delta_t, 90))

                name = 'h_delta_t_raw_'+str(k)
                if width == 0:
                    width = np.std(delta_t)
                if width == 0 or np.isnan(width):
                    width = 0.1
                if np.isnan(median):
                    continue
                title = 'Time resolution for channel '+str(k)+', width 10-90 = {:.2f} ns'.format(width)
                h = create_TH1D(delta_t, name, title,
                                    binning = [ None, median-2*width, median+2*width],
                                    axis_title = [var_dT + ' [ns]', 'Events'])

                rt.gStyle.SetStatX(1.)
                canvas['t_res_raw'][k] = rt.TCanvas('c_t_res_raw_'+str(k), 'c_t_res_raw_'+str(k), 800, 600)
                h.Fit('gaus', 'LQR','', median-width, median+width)
                h = h.DrawCopy('LE')
                canvas['t_res_raw'][k].Update()
                canvas['t_res_raw'][k].SaveAs(out_dir + '/TimeResolution_raw_ch{:02d}.png'.format(k))

            '''=========================== Time resolution 2D ==========================='''
            if 'TimeResRaw2D' in configurations.plots:
                if not 'shape' in conf.keys():
                    print 'Please specify a shape for ch', k, 'in order to print TimeResRaw2D'
                    continue

                b = [var_dT, conf['x_rot'], conf['y_rot']]
                if 'TimeRes2DAmp' in configurations.plots:
                    for kk in configurations.ch_ordered:
                        if 'amp'+str(k) in configurations.channel[kk]['type']:
                            b += ['amp[{}]'.format(kk)]
                    if len(b) == 3:
                        b += ['amp[{}]'.format(k)]

                data_raw = tree2array(chain, b, selection).view(np.recarray)
                data = np.zeros((data_raw.shape[0],3))
                if 'TimeRes2DAmp' in configurations.plots:
                    data = np.zeros((data_raw.shape[0],4))
                    for iev, d in enumerate(data_raw):
                        data[iev] = [d[0][0], d[1], d[2], d[3]]
                else:
                    for iev, d in enumerate(data_raw):
                        data[iev] = [d[0][0], d[1], d[2]]

                if ( len(data) ==0):
                    print 'Empty delta'
                    continue

                bin_width = 1.
                lim = np.array(conf['space_sel']['limits']).astype(np.float)
                size = np.array(conf['shape'].split('x')).astype(np.float)

                x_sec, x_step = np.linspace(lim[0], lim[1], 1+int(size[0]/bin_width), retstep=True)
                y_sec, y_step = np.linspace(lim[2], lim[3], 1+int(size[1]/bin_width), retstep=True)

                if lim[1]-lim[0]+2*x_step >= lim[3]-lim[2]+2*y_step:
                    b = [x_sec.shape[0]+1, lim[0]-x_step, lim[1]+x_step]
                    space_diff = lim[1]-lim[0]+2*x_step - (lim[3]-lim[2])
                    adding_steps = int(space_diff/y_step)
                    if adding_steps%2 == 1:
                        adding_steps +=1
                    b += [y_sec.shape[0]-1+adding_steps, lim[2]-y_step*adding_steps/2, lim[3]+y_step*adding_steps/2]
                else:
                    b = [y_sec.shape[0]+1, lim[2]-y_step, lim[3]+y_step]
                    space_diff = lim[3]-lim[2]+2*y_step - (lim[1]-lim[0])
                    adding_steps = int(space_diff/x_step)
                    if adding_steps%2 == 1:
                        adding_steps +=1
                    b = [x_sec.shape[0]-1+adding_steps, lim[0]-x_step*adding_steps/2, lim[1]+x_step*adding_steps/2] + b

                name = 'h_TimeResRaw2D_ch{:02d}'.format(k)
                title = 'Raw time res breakdown ch{:02d}'.format(k)
                h_2D_res_raw = rt.TH2D(name, title, b[0], b[1], b[2], b[3], b[4], b[5])
                ResRaw = []

                if 'TimeRes2DAmp' in configurations.plots:
                    name = 'h_TimeRes2DAmp_ch{:02d}'.format(k)
                    title = 'Time resolution amplitude corrected breakdown ch{:02d}'.format(k)
                    h_2D_res_amp = rt.TH2D(name, title, b[0], b[1], b[2], b[3], b[4], b[5])
                    ResAmp = []


                it = itertools.product(zip(x_sec[:-1], x_sec[1:]), zip(y_sec[:-1], y_sec[1:]))
                os.mkdir(out_dir + '/TimeResRaw2D_ch{:02d}'.format(k))
                shutil.copy('lib/index.php', out_dir + '/TimeResRaw2D_ch{:02d}'.format(k)+'/index.php')
                if 'TimeRes2DAmp' in configurations.plots:
                    os.mkdir(out_dir + '/TimeRes2DAmp_ch{:02d}'.format(k))
                    shutil.copy('lib/index.php', out_dir +  '/TimeRes2DAmp_ch{:02d}'.format(k)+'/index.php')

                for ib, ((xd, xu), (yd, yu)) in enumerate(it):
                    selx = np.logical_and(data[:,1]>xd, data[:,1]<xu)
                    sely = np.logical_and(data[:,2]>yd, data[:,2]<yu)
                    sel = np.logical_and(selx, sely)
                    if np.sum(sel) < 10:
                        print '[WARNING] Low stat in ch', k, 'pos [{:.1f}, {:.1f}, {:.1f}, {:.1f}]'.format(xd, xu, yd, yu)
                    aux_d = data[sel]
                    dt = aux_d[:,0]

                    fill_TimeResHisto(k, dt,
                                      h_2D_res_raw,
                                      ResRaw, 'TimeResRaw2D',
                                      'Raw',
                                      canvas)


                    if 'TimeRes2DAmp' in configurations.plots:
                        dt_sel = np.logical_and(dt>median-2*width, dt<median+2*width)
                        amp = aux_d[:,3]
                        inputs = np.column_stack((np.ones_like(amp), amp, amp**2))
                        coeff, r, rank, s = np.linalg.lstsq(inputs[dt_sel], dt[dt_sel], rcond=None)
                        b = [None, None, None, None, median-2*width, median+2*width]
                        h_TvsA = create_TH2D(np.column_stack((amp, dt)),
                                             'h_TvsA', 'Channel '+str(k),
                                             binning=b,
                                             axis_title=['Amp [mV]', var_dT + ' [ns]']
                                             )
                        canvas['dt_vs_amp'][k] = rt.TCanvas('dt_vs_amp'+str(k), 'dt_vs_amp'+str(k), 800, 600)
                        h_TvsA.DrawCopy('colz')
                        prof = h_TvsA.ProfileX('prof_amp')
                        prof.SetLineColor(6)
                        prof.SetLineWidth(2)
                        prof.DrawCopy('SAMEE1')

                        f = rt.TF1('amp_fit'+str(k),'[0]+[1]*x+[2]*x^2', np.min(amp), np.max(amp))
                        for j,a in enumerate(coeff):
                            f.SetParameter(j, a)
                        f.SetLineColor(6)
                        f.SetLineStyle(9)
                        f.DrawCopy('SAMEL')
                        canvas['dt_vs_amp'][k].Update()
                        canvas['dt_vs_amp'][k].SaveAs(out_dir +  '/TimeRes2DAmp_ch{:02d}'.format(k)+'/TvsAmp_{}.png'.format(ib))

                        dt_corr = dt - np.dot(inputs, coeff)

                        fill_TimeResHisto(k, dt_corr,
                                          h_2D_res_amp,
                                          ResAmp, 'TimeRes2DAmp',
                                          'Amplitude corrected',
                                          canvas)


                ResRaw = np.array(ResRaw)
                tag = 'TimeResRaw2D_ch{:02d}'.format(k)
                canvas['c_'+tag] = rt.TCanvas('c_'+tag, 'c_'+tag, 800, 600)
                Set_2D_colz_graphics()
                h_2D_res_raw.GetZaxis().SetRangeUser(0.9*np.min(ResRaw[:,0]), 1.1*np.max(ResRaw[:,0]))
                h_2D_res_raw.SetStats(0)
                h_2D_res_raw.SetYTitle('y [mm]')
                h_2D_res_raw.SetXTitle('x [mm]')
                h_2D_res_raw.SetZTitle(var_dT+' [ps]')
                h_2D_res_raw.DrawCopy('colz')
                rt.gStyle.SetPaintTextFormat(".1f");
                h_2D_res_raw.DrawCopy('TEXT SAME ERROR')

                avg_time_res = np.average(ResRaw[:,0], weights=1./np.square(ResRaw[:,1]))
                d_avg_time_res = 1./np.sqrt(np.sum(1./np.square(ResRaw[:,1])))
                results.append([avg_time_res, d_avg_time_res, v_time, v_ref, False])
                if best_result[k].dT[0] < 0 or  best_result[k].dT[0] > avg_time_res:
                    best_result[k].dT = [avg_time_res, d_avg_time_res]
                    best_result[k].var = [v_time, v_ref]
                    best_result[k].AmpCorr = False

                msg = 'Avg raw resolution: {:.1f} +/- {:.1f} ps'.format(avg_time_res, d_avg_time_res)
                print msg
                l = rt.TLatex()
                l.SetTextSize(0.04);
                l.DrawLatexNDC(0.15, 0.89, msg.replace('+/-', '#pm'))
                canvas['c_'+tag].Update()
                canvas['c_'+tag].SaveAs(out_dir + '/TimeResRaw2D_ch{:02d}.png'.format(k))

                if 'TimeRes2DAmp' in configurations.plots:
                    ResAmp = np.array(ResAmp)
                    tag = 'TimeRes2DAmp_ch{:02d}'.format(k)
                    canvas['c_'+tag] = rt.TCanvas('c_'+tag, 'c_'+tag, 800, 600)
                    Set_2D_colz_graphics()
                    h_2D_res_amp.GetZaxis().SetRangeUser(0.9*np.min(ResAmp[:,0]), 1.1*np.max(ResAmp[:,0]))
                    h_2D_res_amp.SetStats(0)
                    h_2D_res_amp.SetYTitle('y [mm]')
                    h_2D_res_amp.SetXTitle('x [mm]')
                    h_2D_res_amp.SetZTitle(var_dT+' [ps]')
                    h_2D_res_amp.DrawCopy('colz')
                    rt.gStyle.SetPaintTextFormat(".1f");
                    h_2D_res_amp.DrawCopy('TEXT SAME ERROR')

                    avg_time_res = np.average(ResAmp[:,0], weights=1./np.square(ResAmp[:,1]))
                    d_avg_time_res = 1./np.sqrt(np.sum(1./np.square(ResAmp[:,1])))
                    results.append([avg_time_res, d_avg_time_res, v_time, v_ref, True])
                    if best_result[k].dT[0] < 0 or  best_result[k].dT[0] > avg_time_res:
                        best_result[k].dT = [avg_time_res, d_avg_time_res]
                        best_result[k].var = [v_time, v_ref]
                        best_result[k].AmpCorr = True
                    msg = 'Avg resolution after amp correction: {:.1f} +/- {:.1f} ps'.format(avg_time_res, d_avg_time_res)
                    print msg
                    l = rt.TLatex()
                    l.SetTextSize(0.04);
                    l.DrawLatexNDC(0.15, 0.89, msg.replace('+/-', '#pm'))
                    canvas['c_'+tag].Update()
                    canvas['c_'+tag].SaveAs(out_dir + '/TimeRes2DAmp_ch{:02d}.png'.format(k))

        with open(headout_dir+'/TimeResolution_ch{}.txt'.format(k),'w') as file:
            ln = '#avg_time_res, d_avg_time_res, var_time, var_ref, AmpCorrection\n'
            file.write(ln)
            for r in results:
                ln = '{:.2f}  {:.2f}  {}  {}  {}\n'.format(r[0], r[1], r[2], r[3], r[4])
                file.write(ln)

    print '\n\n======================= Summary =============================='
    table =  PrettyTable(['Ch', 'Best Resolution [ps]', 'Var ref', 'Var timr', 'Amp corrected'])
    for k, res in best_result.iteritems():
        row = [str(k), '{:.2f} +/- {:.2f}'.format(res.dT[0],res.dT[1])]
        row += [res.var[1], res.var[0], 'Yes' if res.AmpCorr else 'No']
        table.add_row(row)

    print table
    table_txt = table.get_string()
    with open(headout_dir+'/SummaryTable.txt','w') as file:
        file.write(table_txt)
