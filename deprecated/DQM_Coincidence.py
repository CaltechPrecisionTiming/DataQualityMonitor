import numpy as np
import os, re, shutil
import argparse
import ROOT as rt
from root_numpy import tree2array
from lib.histo_utilities import create_TH1D, create_TH2D
from lib.cebefo_style import cebefo_style


#command line arguments returns total number of arguments
def parsing():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str, help="input root file", nargs='+') # input a file name: python DQM_SiPM.py filename -C
    parser.add_argument("-C", "--config", type=str, default='config/Na22_config.txt', help="Config file") # optional argument, type -C to add config file
#unsure how files are saved? what format'''
    parser.add_argument("-S", "--save_loc", type=str, default='./out_plots/', help="Saving location")
    # if absent from command line then true else false
    parser.add_argument("-B", "--batch", default=True, action='store_false', help="Root batch mode")

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
            #print('i: ' + str(i))
            #print("l[0]: " + l[0])
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
                #print('l[0]: ' + l[0] + 'n: ' + str(n))
                self.ch_ordered.append(n)
                l = l[1:]
                self.channel[n] = {}
                for k, val in zip(self.labels, l):
                    if 'idx' in k:
                        self.channel[n][k] = int(val)
                    else:
                        self.channel[n][k] = val


def define_range_around_peak(h, perc = [0.4, 0.4], m = 0): # h is a histogram
    SS = rt.TSpectrum()
    n_pks = SS.Search(h, 1, "", 0.2)
    i_max = 0
    max = -1
    y = SS.GetPositionY()
    for i in range(n_pks):
        if y[i] > max and SS.GetPositionX()[i] > m:
            max = y[i]
            i_max = i
    n_pk = h.FindBin(SS.GetPositionX()[i_max])


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
    return x_low, x_up

if __name__ == '__main__':
    cebefo_style()

    rt.gErrorIgnoreLevel = 6000

    args = parsing()

    if args.batch:
        rt.gROOT.SetBatch()

    configurations = Config(args.config) # constructor input is config file
# '''dont need flag variable? can delete? should root file have trees and leaves?'''
    if args.runs_interval==None:
        # ==========search does not match?==========
        aux = re.search(r'Run[0-9]+', args.input_file[0]) # + matches one or more
        flag = aux.group(0)[3:] # error nonetype does not have group
        aux = re.search(r'Run[0-9]+', args.input_file[-1])
        flag += '_'+aux.group(0)[3:]
    elif len(args.runs_interval)==2:
        deduced_file_list = []
        for i in range(args.runs_interval[0], args.runs_interval[1]+1): # if input -N 2 2 then the range is from 2-3, only run once
            aux = args.input_file[0].replace('XXX', str(i)) # file does not have 'XXX', no need to replace
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
        print 'Directory flag: ', flag
        raise Exception(' ')

    save_loc = args.save_loc
    if os.path.isdir(save_loc):
        if save_loc[-1] != '/':
            save_loc += '/'
        out_dir = save_loc + 'Run' + str(flag) + '_plots'
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.mkdir(out_dir)
    else:
        print 'Save location not existing:', save_loc
        raise Exception(' ')


#'''match treename, treename should be pulse'''
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

# 2D array
    canvas = {}
    canvas['amp'] = {}
    canvas['int'] = {}
    canvas['risetime'] = {}
    canvas['baseline_RMS'] = {}
    canvas['IL_2'] = {}
    canvas['IL_5'] = {}
    canvas['IL_8'] = {}
    canvas['IL_10'] = {}
    canvas['IL_15'] = {}
    canvas['IL_20'] = {}
    canvas['IL_25'] = {}
    canvas['IL_30'] = {}
    canvas['IL_35'] = {}
    canvas['IL_40'] = {}
    canvas['IL_45'] = {}
    canvas['IL_50'] = {}

    rt.gStyle.SetStatY(0.98);
    rt.gStyle.SetStatX(1.);
    rt.gStyle.SetStatW(0.2);
    rt.gStyle.SetStatH(0.1);
    rt.gStyle.SetHistLineWidth(2);


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
            print '\tAmplitude'
            name = 'h_amp_'+str(k)
            title = 'Amplitude channel '+str(k) # str(k) is the channel number

            amp_aux = np.concatenate(list(tree2array( chain, 'amp['+str(k)+']')))

            h = rt.TH1D(name, title, 80, np.percentile(amp_aux, 1), min(900., np.percentile(amp_aux, 99)))
            h.SetXTitle('Peak amplitude [mV]')
            h.SetYTitle('Events / {:.1f} mV'.format(h.GetBinWidth(1)))
            chain.Project(name, 'amp['+str(k)+']')

            x_low, x_up = define_range_around_peak(h, [0.35, 0.3], 20)

            canvas['amp'][k] = rt.TCanvas('c_amp_'+str(k), 'c_amp_'+str(k), 800, 600)
            h.DrawCopy('E')
            line.DrawLine(x_low, 0, x_low, h.GetMaximum())
            line.DrawLine(x_up, 0, x_up, h.GetMaximum())

            canvas['amp'][k].Update()
            canvas['amp'][k].SaveAs(out_dir + '/Amp_ch'+str(k)+'.png')

            selection += ' && (amp[{}]>{} && amp[{}]<{})'.format(k, x_low, k, x_up)

        '''=========================== Integral ==========================='''
        if 'Int' in configurations.plots:
            print '\tIntegral'

            canvas['int'][k] = rt.TCanvas('c_int_'+str(k), 'c_int_'+str(k), 800, 600)

            name = 'h_int_'+str(k)
            title = 'Integral channel '+str(k)
            int_aux = -np.concatenate(list(tree2array(chain, 'integral['+str(k)+']', 'integral['+str(k)+'] != 0')))
            h = rt.TH1D(name, title, 100, np.percentile(int_aux, 0.1), np.percentile(int_aux, 99.9))
            h.SetXTitle('Integral [pC]')
            h.SetYTitle('Events / {:.1f} pC'.format(h.GetBinWidth(1)))
            chain.Project(name, '-integral['+str(k)+']', '-integral['+str(k)+'] != 0')

            x_low, x_up = define_range_around_peak(h, [0.25, 0.3], 1)

            if conf['idx_ref'] >= 0:
                res = h.Fit('landau','LQSR', '', x_low, x_up)

            if(h.GetMaximum() - h.GetMinimum() > 5000):
                canvas['int'][k].SetLogy()
            h.DrawCopy('E')
            line.DrawLine(x_low, 0, x_low, h.GetMaximum())
            line.DrawLine(x_up, 0, x_up, h.GetMaximum())

            canvas['int'][k].Update()
            canvas['int'][k].SaveAs(out_dir + '/Int_ch'+str(k)+'.png')

            selection += ' && (-integral[{}]>{} && -integral[{}]<{})'.format(k, x_low, k, x_up)


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
                cut = selection + '&&' + configurations.channel[conf['idx_ref']]['sel'] #sel???

            chain.Project(name, 'risetime['+str(k)+']', cut)

            if (h.GetMaximum() - h.GetMinimum() > 50000):
                canvas['risetime'][k].SetLogy()

            ##ploting histogram
            h.DrawCopy('lE')
            canvas['risetime'][k].Update()
            canvas['risetime'][k].SaveAs(out_dir+'/risetime_ch'+str(k)+'.png')

        #'''=========================== End Selections ==========================='''
            conf['sel'] =  selection

        '''=======================================Baseline RMS============================='''
        if 'Baseline_RMS' in configurations.plots:
            print '\tBaseline_RMS'
            canvas['baseline_RMS'][k] = rt.TCanvas('c_baseline_RMS_'+str(k), 'c_baseline_RMS_'+str(k), 800, 600)
            name   = 'h_baseline_RMS_' + str(k)
            title  = 'Baseline_RMS ch' + str(k)
            baseline_RMS_aux = np.concatenate(list(tree2array( chain, 'baseline_RMS['+str(k)+']')))
            h = rt.TH1D(name, title, 80, np.percentile(baseline_RMS_aux, 1), min(900., np.percentile(baseline_RMS_aux, 99)))
            h.SetXTitle('Baseline RMS [mV]')
            h.SetYTitle('Events')
            chain.Project(name, 'baseline_RMS['+str(k)+']')
            h.DrawCopy('E')
            canvas['baseline_RMS'][k].Update()
            canvas['baseline_RMS'][k].SaveAs(out_dir + '/Baseline_RMS_ch'+str(k)+'.png')

        '''===========================================IL_2======================================'''
        if 'IL_2' in configurations.plots and k > 0:
            print '\tIL_2[1] - IL_2[0]'
            canvas['IL_2'][k] = rt.TCanvas('c_IL_2_'+str(k), 'c_IL_2_'+str(k), 800, 600)
            name   = 'h_IL_2_' + str(k)
            title  = 'IL_2[1] - IL_2[0]'
            IL_2_chk = np.concatenate(list(tree2array( chain, 'IL_2['+str(k)+']')))
            IL_2_ch0 = np.concatenate(list(tree2array( chain, 'IL_2['+str(0)+']')))
            IL_2_chk_M_ch0 = []
            for i in range(len(IL_2_chk)):
                IL_2_chk_M_ch0.append(IL_2_chk[i] - IL_2_ch0[i])
            median = np.percentile(IL_2_chk_M_ch0, 50) # 50th percentile, median
            #print('median ' + str(median))
            width = np.percentile(IL_2_chk_M_ch0, 90) - np.percentile(IL_2_chk_M_ch0, 10)
            #print('width ' + str(width))
            h = rt.TH1D(name, title, 80, median - 2.* width, median + 2.* width)
            chain.Project(name, 'IL_2['+str(k)+'] - IL_2[0]')
            #for i in range(len(IL_2_chk_M_ch0)):
            #    h.Fill(IL_2_chk_M_ch0[i])
            f1 = rt.TF1("multi_gaus","gaus(0)", median - 0.5 * width, median + 0.5 * width)
            f1.SetParameter(0, 3000)
            f1.SetParameter(1, h.GetMean())
            f1.SetParameter(2, width * 0.25)
            h.Fit("multi_gaus", "LR", "")
            h.SetXTitle('IL_2 [ns]')
            h.SetYTitle('Events')
            h.DrawCopy('E')
            #print h.GetMean(), h.GetRMS()
            canvas['IL_2'][k].Update()
            canvas['IL_2'][k].SaveAs(out_dir + '/IL_2_ch'+str(k)+'.png')

        '''===========================================IL_5======================================'''
        if 'IL_5' in configurations.plots and k > 0:
            print '\tIL_5[1] - IL_5[0]'
            canvas['IL_5'][k] = rt.TCanvas('c_IL_5_'+str(k), 'c_IL_5_'+str(k), 800, 600)
            name   = 'h_IL_5_' + str(k)
            title  = 'IL_5[1] - IL_5[0]'
            IL_5_chk = np.concatenate(list(tree2array( chain, 'IL_5['+str(k)+']')))
            IL_5_ch0 = np.concatenate(list(tree2array( chain, 'IL_5['+str(0)+']')))
            IL_5_chk_M_ch0 = []
            for i in range(len(IL_5_chk)):
                IL_5_chk_M_ch0.append(IL_5_chk[i] - IL_5_ch0[i])
            median = np.percentile(IL_5_chk_M_ch0, 50) # 50th percentile, median
            width = np.abs(np.percentile(IL_5_chk_M_ch0, 10) - np.percentile(IL_5_chk_M_ch0, 90))
            h = rt.TH1D(name, title, 80, median - 2.*width, median + 2.*width)
            chain.Project(name, 'IL_5['+str(k)+'] - IL_5[0]')
            #for i in range(len(IL_5_chk_M_ch0)):
            #    h.Fill(IL_5_chk_M_ch0[i])
            f1 = rt.TF1("multi_gaus","gaus(0)", median - 0.5 * width, median + 0.5 * width)
            f1.SetParameter(0, 3000)
            f1.SetParameter(1, h.GetMean())
            f1.SetParameter(2, width * 0.25)
            h.Fit("multi_gaus", "LR", "")
            h.SetXTitle('IL_2 [ns]')
            h.SetYTitle('Events')
            h.DrawCopy('E')
            #print h.GetMean(), h.GetRMS()
            canvas['IL_5'][k].Update()
            canvas['IL_5'][k].SaveAs(out_dir + '/IL_5_ch'+str(k)+'.png')

        '''===========================================IL_8======================================'''
        if 'IL_8' in configurations.plots and k > 0:
                print '\tIL_8[1] - IL_8[0]'
                canvas['IL_8'][k] = rt.TCanvas('c_IL_8_'+str(k), 'c_IL_8_'+str(k), 800, 600)
                name   = 'h_IL_8_' + str(k)
                title  = 'IL_8[1] - IL_8[0]'
                IL_8_chk = np.concatenate(list(tree2array( chain, 'IL_8['+str(k)+']')))
                IL_8_ch0 = np.concatenate(list(tree2array( chain, 'IL_8['+str(0)+']')))
                IL_8_chk_M_ch0 = []
                for i in range(len(IL_8_chk)):
                    IL_8_chk_M_ch0.append(IL_8_chk[i] - IL_8_ch0[i])
                median = np.percentile(IL_8_chk_M_ch0, 50) # 50th percentile, median
                width = np.abs(np.percentile(IL_8_chk_M_ch0, 10) - np.percentile(IL_8_chk_M_ch0, 90))
                h = rt.TH1D(name, title, 80, median - 2. * width, median + 2. * width)
                chain.Project(name, 'IL_8['+str(k)+'] - IL_8[0]')
                #for i in range(len(IL_8_chk_M_ch0)):
                #    h.Fill(IL_8_chk_M_ch0[i])
                f1 = rt.TF1("multi_gaus","gaus(0)", median - 0.5 * width, median + 0.5 * width)
                f1.SetParameter(0, 3000)
                f1.SetParameter(1, h.GetMean())
                f1.SetParameter(2, width * 0.25)
                h.Fit("multi_gaus", "LR", "")
                h.SetXTitle('IL_2 [ns]')
                h.SetYTitle('Events')
                h.DrawCopy('E')
                #print h.GetMean(), h.GetRMS()
                canvas['IL_8'][k].Update()
                canvas['IL_8'][k].SaveAs(out_dir + '/IL_8_ch'+str(k)+'.png')

        '''===========================================IL_10======================================'''
        if 'IL_10' in configurations.plots and k > 0:
                print '\tIL_10[1] - IL_10[0]'
                canvas['IL_10'][k] = rt.TCanvas('c_IL_10_'+str(k), 'c_IL_10_'+str(k), 800, 600)
                name   = 'h_IL_10_' + str(k)
                title  = 'IL_10[1] - IL_10[0]'
                IL_10_chk = np.concatenate(list(tree2array( chain, 'IL_10['+str(k)+']')))
                IL_10_ch0 = np.concatenate(list(tree2array( chain, 'IL_10['+str(0)+']')))
                IL_10_chk_M_ch0 = []
                for i in range(len(IL_10_chk)):
                    IL_10_chk_M_ch0.append(IL_10_chk[i] - IL_10_ch0[i])
                median = np.percentile(IL_10_chk_M_ch0, 50) # 50th percentile, median
                width = np.abs(np.percentile(IL_10_chk_M_ch0, 10) - np.percentile(IL_10_chk_M_ch0, 90))
                h = rt.TH1D(name, title, 80, median - 2. * width, median + 2. * width)
                chain.Project(name, 'IL_10['+str(k)+'] - IL_10[0]')
                #for i in range(len(IL_10_chk_M_ch0)):
                #    h.Fill(IL_10_chk_M_ch0[i])
                f1 = rt.TF1("multi_gaus","gaus(0)", median - 0.5 * width, median + 0.5 * width)
                f1.SetParameter(0, 3000)
                f1.SetParameter(1, h.GetMean())
                f1.SetParameter(2, width * 0.25)
                h.Fit("multi_gaus", "LR", "")
                h.SetXTitle('IL_2 [ns]')
                h.SetYTitle('Events')
                h.DrawCopy('E')
                #print h.GetMean(), h.GetRMS()
                canvas['IL_10'][k].Update()
                canvas['IL_10'][k].SaveAs(out_dir + '/IL_10_ch'+str(k)+'.png')

        '''===========================================IL_15======================================'''
        if 'IL_15' in configurations.plots and k > 0:
                print '\tIL_15[1] - IL_15[0]'
                canvas['IL_15'][k] = rt.TCanvas('c_IL_15_'+str(k), 'c_IL_15_'+str(k), 800, 600)
                name   = 'h_IL_15_' + str(k)
                title  = 'IL_15[1] - IL_15[0]'
                IL_15_chk = np.concatenate(list(tree2array( chain, 'IL_15['+str(k)+']')))
                IL_15_ch0 = np.concatenate(list(tree2array( chain, 'IL_15['+str(0)+']')))
                IL_15_chk_M_ch0 = []
                for i in range(len(IL_15_chk)):
                    IL_15_chk_M_ch0.append(IL_15_chk[i] - IL_15_ch0[i])
                median = np.percentile(IL_15_chk_M_ch0, 50) # 50th percentile, median
                width = np.abs(np.percentile(IL_15_chk_M_ch0, 10) - np.percentile(IL_15_chk_M_ch0, 90))
                h = rt.TH1D(name, title, 80, median - 2.* width, median + 2.* width)
                chain.Project(name, 'IL_15['+str(k)+'] - IL_15[0]')
                #for i in range(len(IL_15_chk_M_ch0)):
                #    h.Fill(IL_15_chk_M_ch0[i])
                f1 = rt.TF1("multi_gaus","gaus(0)", median - 0.5 * width, median + 0.5 * width)
                f1.SetParameter(0, 3000)
                f1.SetParameter(1, h.GetMean())
                f1.SetParameter(2, width * 0.25)
                h.Fit("multi_gaus", "LR", "")
                h.SetXTitle('IL_2 [ns]')
                h.SetYTitle('Events')
                h.DrawCopy('E')
                #print h.GetMean(), h.GetRMS()
                canvas['IL_15'][k].Update()
                canvas['IL_15'][k].SaveAs(out_dir + '/IL_15_ch'+str(k)+'.png')
        '''===========================================IL_20======================================'''
        if 'IL_20' in configurations.plots and k > 0:
                print '\tIL_20[1] - IL_20[0]'
                canvas['IL_20'][k] = rt.TCanvas('c_IL_20_'+str(k), 'c_IL_20_'+str(k), 800, 600)
                name   = 'h_IL_20_' + str(k)
                title  = 'IL_20[1] - IL_20[0]'
                IL_20_chk = np.concatenate(list(tree2array( chain, 'IL_20['+str(k)+']')))
                IL_20_ch0 = np.concatenate(list(tree2array( chain, 'IL_20['+str(0)+']')))
                IL_20_chk_M_ch0 = []
                for i in range(len(IL_20_chk)):
                    IL_20_chk_M_ch0.append(IL_20_chk[i] - IL_20_ch0[i])
                median = np.percentile(IL_20_chk_M_ch0, 50) # 50th percentile, median
                width = np.abs(np.percentile(IL_20_chk_M_ch0, 10) - np.percentile(IL_20_chk_M_ch0, 90))
                h = rt.TH1D(name, title, 80, median - 2.*width, median + 2.*width)
                chain.Project(name, 'IL_20['+str(k)+'] - IL_20[0]')
                #for i in range(len(IL_20_chk_M_ch0)):
                #    h.Fill(IL_20_chk_M_ch0[i])
                f1 = rt.TF1("multi_gaus","gaus(0)", median - 0.5 * width, median + 0.5 * width)
                f1.SetParameter(0, 3000)
                f1.SetParameter(1, h.GetMean())
                f1.SetParameter(2, width * 0.25)
                h.Fit("multi_gaus", "LR", "")
                h.SetXTitle('IL_2 [ns]')
                h.SetYTitle('Events')
                h.DrawCopy('E')
                #print h.GetMean(), h.GetRMS()
                canvas['IL_20'][k].Update()
                canvas['IL_20'][k].SaveAs(out_dir + '/IL_20_ch'+str(k)+'.png')

        '''===========================================IL_25======================================'''
        if 'IL_25' in configurations.plots and k > 0:
                print '\tIL_25[1] - IL_25[0]'
                canvas['IL_25'][k] = rt.TCanvas('c_IL_25_'+str(k), 'c_IL_25_'+str(k), 800, 600)
                name   = 'h_IL_25_' + str(k)
                title  = 'IL_25[1] - IL_25[0]'
                IL_25_chk = np.concatenate(list(tree2array( chain, 'IL_25['+str(k)+']')))
                IL_25_ch0 = np.concatenate(list(tree2array( chain, 'IL_25['+str(0)+']')))
                IL_25_chk_M_ch0 = []
                for i in range(len(IL_25_chk)):
                    IL_25_chk_M_ch0.append(IL_25_chk[i] - IL_25_ch0[i])
                median = np.percentile(IL_25_chk_M_ch0, 50) # 50th percentile, median
                width = np.abs(np.percentile(IL_25_chk_M_ch0, 10) - np.percentile(IL_25_chk_M_ch0, 90))
                h = rt.TH1D(name, title, 80, median - 2.*width, median + 2.*width)
                chain.Project(name, 'IL_25['+str(k)+'] - IL_25[0]')
                #for i in range(len(IL_25_chk_M_ch0)):
                #    h.Fill(IL_25_chk_M_ch0[i])
                f1 = rt.TF1("multi_gaus","gaus(0)", median - 0.5 * width, median + 0.5 * width)
                f1.SetParameter(0, 3000)
                f1.SetParameter(1, h.GetMean())
                f1.SetParameter(2, width * 0.25)
                h.Fit("multi_gaus", "LR", "")
                h.SetXTitle('IL_2 [ns]')
                h.SetYTitle('Events')
                h.DrawCopy('E')
                #print h.GetMean(), h.GetRMS()
                canvas['IL_25'][k].Update()
                canvas['IL_25'][k].SaveAs(out_dir + '/IL_25_ch'+str(k)+'.png')

        '''===========================================IL_30======================================'''
        if 'IL_30' in configurations.plots and k > 0:
                print '\tIL_30[1] - IL_30[0]'
                canvas['IL_30'][k] = rt.TCanvas('c_IL_30_'+str(k), 'c_IL_30_'+str(k), 800, 600)
                name   = 'h_IL_30_' + str(k)
                title  = 'IL_30[1] - IL_30[0]'
                IL_30_chk = np.concatenate(list(tree2array( chain, 'IL_30['+str(k)+']')))
                IL_30_ch0 = np.concatenate(list(tree2array( chain, 'IL_30['+str(0)+']')))
                IL_30_chk_M_ch0 = []
                for i in range(len(IL_30_chk)):
                    IL_30_chk_M_ch0.append(IL_30_chk[i] - IL_30_ch0[i])
                median = np.percentile(IL_30_chk_M_ch0, 50) # 50th percentile, median
                width = np.abs(np.percentile(IL_30_chk_M_ch0, 10) - np.percentile(IL_30_chk_M_ch0, 90))
                h = rt.TH1D(name, title, 80, median - 2.*width, median + 2.*width)
                chain.Project(name, 'IL_30['+str(k)+'] - IL_30[0]')
                #for i in range(len(IL_30_chk_M_ch0)):
                #    h.Fill(IL_30_chk_M_ch0[i])
                f1 = rt.TF1("multi_gaus","gaus(0)", median - 0.5 * width, median + 0.5 * width)
                f1.SetParameter(0, 3000)
                f1.SetParameter(1, h.GetMean())
                f1.SetParameter(2, width * 0.25)
                h.Fit("multi_gaus", "LR", "")
                h.SetXTitle('IL_2 [ns]')
                h.SetYTitle('Events')
                h.DrawCopy('E')
                #print h.GetMean(), h.GetRMS()
                canvas['IL_30'][k].Update()
                canvas['IL_30'][k].SaveAs(out_dir + '/IL_30_ch'+str(k)+'.png')

        '''===========================================IL_35======================================'''
        if 'IL_35' in configurations.plots and k > 0:
                print '\tIL_35[1] - IL_35[0]'
                canvas['IL_35'][k] = rt.TCanvas('c_IL_35_'+str(k), 'c_IL_35_'+str(k), 800, 600)
                name   = 'h_IL_35_' + str(k)
                title  = 'IL_35[1] - IL_35[0]'
                IL_35_chk = np.concatenate(list(tree2array( chain, 'IL_35['+str(k)+']')))
                IL_35_ch0 = np.concatenate(list(tree2array( chain, 'IL_35['+str(0)+']')))
                IL_35_chk_M_ch0 = []
                for i in range(len(IL_35_chk)):
                    IL_35_chk_M_ch0.append(IL_35_chk[i] - IL_35_ch0[i])
                median = np.percentile(IL_35_chk_M_ch0, 50) # 50th percentile, median
                width = np.abs(np.percentile(IL_35_chk_M_ch0, 10) - np.percentile(IL_35_chk_M_ch0, 90))
                h = rt.TH1D(name, title, 80, median - 2.*width, median + 2.*width)
                chain.Project(name, 'IL_35['+str(k)+'] - IL_35[0]')
                #for i in range(len(IL_35_chk_M_ch0)):
                #    h.Fill(IL_35_chk_M_ch0[i])
                f1 = rt.TF1("multi_gaus","gaus(0)", median - 0.5 * width, median + 0.5 * width)
                f1.SetParameter(0, 3000)
                f1.SetParameter(1, h.GetMean())
                f1.SetParameter(2, width * 0.25)
                h.Fit("multi_gaus", "LR", "")
                h.SetXTitle('IL_2 [ns]')
                h.SetYTitle('Events')
                h.DrawCopy('E')
                #print h.GetMean(), h.GetRMS()
                canvas['IL_35'][k].Update()
                canvas['IL_35'][k].SaveAs(out_dir + '/IL_35_ch'+str(k)+'.png')

        '''===========================================IL_40======================================'''
        if 'IL_40' in configurations.plots and k > 0:
                print '\tIL_40[1] - IL_40[0]'
                canvas['IL_40'][k] = rt.TCanvas('c_IL_40_'+str(k), 'c_IL_40_'+str(k), 800, 600)
                name   = 'h_IL_40_' + str(k)
                title  = 'IL_40[1] - IL_40[0]'
                IL_40_chk = np.concatenate(list(tree2array( chain, 'IL_40['+str(k)+']')))
                IL_40_ch0 = np.concatenate(list(tree2array( chain, 'IL_40['+str(0)+']')))
                IL_40_chk_M_ch0 = []
                for i in range(len(IL_40_chk)):
                    IL_40_chk_M_ch0.append(IL_40_chk[i] - IL_40_ch0[i])
                median = np.percentile(IL_40_chk_M_ch0, 50) # 50th percentile, median
                width = np.abs(np.percentile(IL_40_chk_M_ch0, 10) - np.percentile(IL_40_chk_M_ch0, 90))
                h = rt.TH1D(name, title, 80, median - 2.*width, median + 2.*width)
                chain.Project(name, 'IL_40['+str(k)+'] - IL_40[0]')
                #for i in range(len(IL_40_chk_M_ch0)):
                #    h.Fill(IL_40_chk_M_ch0[i])
                f1 = rt.TF1("multi_gaus","gaus(0)", median - 0.5 * width, median + 0.5 * width)
                f1.SetParameter(0, 3000)
                f1.SetParameter(1, h.GetMean())
                f1.SetParameter(2, width * 0.25)
                h.Fit("multi_gaus", "LR", "")
                h.SetXTitle('IL_2 [ns]')
                h.SetYTitle('Events')
                h.DrawCopy('E')
                #print h.GetMean(), h.GetRMS()
                canvas['IL_40'][k].Update()
                canvas['IL_40'][k].SaveAs(out_dir + '/IL_40_ch'+str(k)+'.png')
        '''===========================================IL_45======================================'''
        if 'IL_45' in configurations.plots and k > 0:
                print '\tIL_45[1] - IL_45[0]'
                canvas['IL_45'][k] = rt.TCanvas('c_IL_45_'+str(k), 'c_IL_45_'+str(k), 800, 600)
                name   = 'h_IL_45_' + str(k)
                title  = 'IL_45[1] - IL_45[0]'
                IL_45_chk = np.concatenate(list(tree2array( chain, 'IL_45['+str(k)+']')))
                IL_45_ch0 = np.concatenate(list(tree2array( chain, 'IL_45['+str(0)+']')))
                IL_45_chk_M_ch0 = []
                for i in range(len(IL_45_chk)):
                    IL_45_chk_M_ch0.append(IL_45_chk[i] - IL_45_ch0[i])
                median = np.percentile(IL_45_chk_M_ch0, 50) # 50th percentile, median
                width = np.abs(np.percentile(IL_45_chk_M_ch0, 10) - np.percentile(IL_45_chk_M_ch0, 90))
                h = rt.TH1D(name, title, 80, median - 2.*width, median + 2.*width)
                chain.Project(name, 'IL_45['+str(k)+'] - IL_45[0]')
                #for i in range(len(IL_45_chk_M_ch0)):
                #    h.Fill(IL_45_chk_M_ch0[i])
                f1 = rt.TF1("multi_gaus","gaus(0)", median - 0.5 * width, median + 0.5 * width)
                f1.SetParameter(0, 3000)
                f1.SetParameter(1, h.GetMean())
                f1.SetParameter(2, width * 0.25)
                h.Fit("multi_gaus", "LR", "")
                h.SetXTitle('IL_2 [ns]')
                h.SetYTitle('Events')
                h.DrawCopy('E')
                #print h.GetMean(), h.GetRMS()
                canvas['IL_45'][k].Update()
                canvas['IL_45'][k].SaveAs(out_dir + '/IL_45_ch'+str(k)+'.png')
        '''===========================================IL_50======================================'''
        if 'IL_50' in configurations.plots and k > 0:
                print '\tIL_50[1] - IL_50[0]'
                canvas['IL_50'][k] = rt.TCanvas('c_IL_50_'+str(k), 'c_IL_50_'+str(k), 800, 600)
                name   = 'h_IL_50_' + str(k)
                title  = 'IL_50[1] - IL_50[0]'
                IL_50_chk = np.concatenate(list(tree2array( chain, 'IL_50['+str(k)+']')))
                IL_50_ch0 = np.concatenate(list(tree2array( chain, 'IL_50['+str(0)+']')))
                IL_50_chk_M_ch0 = []
                for i in range(len(IL_50_chk)):
                    IL_50_chk_M_ch0.append(IL_50_chk[i] - IL_50_ch0[i])
                median = np.percentile(IL_50_chk_M_ch0, 50) # 50th percentile, median
                width = np.abs(np.percentile(IL_50_chk_M_ch0, 10) - np.percentile(IL_50_chk_M_ch0, 90))
                h = rt.TH1D(name, title, 80, median - 2.*width, median + 2.*width)
                chain.Project(name, 'IL_50['+str(k)+'] - IL_50[0]')
                #for i in range(len(IL_50_chk_M_ch0)):
                #    h.Fill(IL_50_chk_M_ch0[i])
                f1 = rt.TF1("multi_gaus","gaus(0)", median - 0.5 * width, median + 0.5 * width)
                f1.SetParameter(0, 3000)
                f1.SetParameter(1, h.GetMean())
                f1.SetParameter(2, width * 0.25)
                h.Fit("multi_gaus", "LR", "")
                h.SetXTitle('IL_2 [ns]')
                h.SetYTitle('Events')
                h.DrawCopy('E')
                #print h.GetMean(), h.GetRMS()
                canvas['IL_50'][k].Update()
                canvas['IL_50'][k].SaveAs(out_dir + '/IL_50_ch'+str(k)+'.png')
