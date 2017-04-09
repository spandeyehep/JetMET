#!/usr/bin/env python
''' make result script from directory of L3res output 
'''
#
# Standard imports and batch mode
#
import ROOT
ROOT.gROOT.SetBatch(True)
import itertools
import os

from math                                import sqrt, cos, sin, pi, atan2
from RootTools.core.standard             import *
from JetMET.tools.user                   import plot_directory as user_plot_directory
from JetMET.tools.helpers                import getObjFromFile

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging" )
argParser.add_argument('--version',            action='store',      default='V6',            help='JEC version as postfix to 23Sep2016' )
argParser.add_argument('--inputDir',           action='store',      default='/afs/hephy.at/user/r/rschoefbeck/www/JEC/L3res/V6/DYnJets/inclusive/ee_log/ptll30-njet1p/',      help="input directory")
argParser.add_argument('--outputFile',         action='store',      help="output file. If not present, use <plotDirectory>/L3res.root")
args = argParser.parse_args()

#
# Logger
#
import JetMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

# pt and thresholds
from JetMET.JEC.L3res.thresholds import ptll_bins, abs_eta_bins

if not args.outputFile:
    args.outputFile = os.path.join( args.inputDir, "L3res.root" )


def getProfiles( fname ):
    f = ROOT.TFile(fname)
    assert not f.IsZombie()
    f.cd()
    canvas = f.Get(f.GetListOfKeys()[0].GetName())
    top, bottom = canvas.GetListOfPrimitives()
    mc, data = None, None
    for o in top.GetListOfPrimitives():
#        print o.GetName(), fname.replace('.root','')+'_all_mc_combined'
        if o.GetName().startswith(os.path.basename(fname).replace('.root','')+'_all_mc_combined'): 
            mc   = o
        if o.GetName().startswith(os.path.basename(fname).replace('.root','')+'_data'):
            data = o

    ROOT.gDirectory.cd('PyROOT:/')
    mc_c, data_c    = None, None
    if mc:   mc_c   = mc.Clone()
    if data: data_c = data.Clone()

    f.Close()

    return mc_c, data_c

#for era in ['Run2016BCD', 'Run2016EF', 'Run2016GH']:

method_name = {'mpf':'MPF_CHS', 'ptbal':'PtBal_CHS', 'mpfNoType1':'MPF-notypeI_CHS'}
stuff = []
for abs_eta_bin in abs_eta_bins:
    ptll_vs_ptll_file = os.path.join( args.inputDir, 'ptll_profile_ptll_for_eta_%3.2f_%3.2f.root'%( abs_eta_bin ) )
   
    ptll_mc, ptll_data = getProfiles(ptll_vs_ptll_file)
    for method in ['mpf', 'ptbal', 'mpfNoType1']:
        response_file = os.path.join( args.inputDir, 'r_%s_%s_profile_ptll_for_eta_%3.2f_%3.2f.root'%( ( method, args.version ) + abs_eta_bin ) ) 
        r_mc, r_data = getProfiles(response_file)

        n_bins = r_mc.GetNbinsX()

        mc    = ROOT.TGraphErrors(n_bins)
        mc.   SetName("MC_%s_a30_eta%s_%s_L1L2L3"%   ( method_name[method], str(int(10*abs_eta_bin[0])).zfill(2), str(int(10*abs_eta_bin[1])).zfill(2) ))
        mc.SetTitle( "%s;%s;%s"%( mc.GetName(), r_mc.GetXaxis().GetTitle(), r_mc.GetYaxis().GetTitle()) )
        data  = ROOT.TGraphErrors(n_bins)
        data. SetName("Data_%s_a30_eta%s_%s_L1L2L3"% ( method_name[method], str(int(10*abs_eta_bin[0])).zfill(2), str(int(10*abs_eta_bin[1])).zfill(2) ))
        data.SetTitle( "%s;%s;%s"%( data.GetName(), r_data.GetXaxis().GetTitle(), r_data.GetYaxis().GetTitle()) )
        ratio = ROOT.TGraphErrors(n_bins)
        ratio.SetName("Ratio_%s_a30_eta%s_%s_L1L2L3"%( method_name[method], str(int(10*abs_eta_bin[0])).zfill(2), str(int(10*abs_eta_bin[1])).zfill(2) ))
        ratio.SetTitle( "Ratio %s;Data/MC %s;%s"%( data.GetName(), r_data.GetXaxis().GetTitle(), r_data.GetYaxis().GetTitle()) )

        for i_bin in range(1,n_bins+1):
            pt_mean_mc      = ptll_mc.GetBinContent(i_bin)
            pt_error_mc     = ptll_mc.GetBinError(i_bin)
            r_mean_mc       = r_mc.GetBinContent(i_bin)
            r_error_mc      = r_mc.GetBinError(i_bin)
            pt_mean_data    = ptll_data.GetBinContent(i_bin)
            pt_error_data   = ptll_data.GetBinError(i_bin)
            r_mean_data     = r_data.GetBinContent(i_bin)
            r_error_data    = r_data.GetBinError(i_bin)

            mc.   SetPoint(i_bin,       pt_mean_mc,     r_mean_mc)
            mc.   SetPointError(i_bin,  pt_error_mc,    r_error_mc)

            data. SetPoint(i_bin,       pt_mean_data,   r_mean_data)
            data. SetPointError(i_bin,  pt_error_data,  r_error_data)

            if r_mean_mc>0 and r_mean_data>0:
                ratio.SetPoint(i_bin, pt_mean_mc, r_mean_data/r_mean_mc)
                ratio.SetPointError(i_bin, pt_error_mc, r_mean_data/r_mean_mc*sqrt(r_error_data**2/r_mean_data**2 + r_error_mc**2/r_mean_mc**2) )

        for tg in [mc, data, ratio]:
            for p in reversed(range(tg.GetN()+1)):
                if tg.GetX()[p]==0 or tg.GetY()[p]==0:
                    tg.RemovePoint(p)

        stuff.append( mc )
        stuff.append( data )
        stuff.append( ratio )

f_out = ROOT.TFile( args.outputFile, "recreate")
f_out.cd()
for o in stuff:
    o.Write()

f_out.Close() 
