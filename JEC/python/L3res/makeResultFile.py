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
argParser.add_argument('--inputDir',           action='store',      default='/afs/hephy.at/user/r/rschoefbeck/www/JEC/L3res_small/V6/DYnJets/Run2016G/mumu_log/ptll30-njet1p/',      help="input directory")
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

def fillTGraphs( r_data, r_mc, ptll_data, ptll_mc, ratio, data, mc) : 
    n_bins = r_mc.GetNbinsX()

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


obj_name = {'mpf':'MPF_CHS', 'ptbal':'PtBal_CHS', 'mpfNoType1':'MPF-notypeI_CHS'}

logger.info( "Reading %s", args.inputDir )

stuff = []

# Eta binned TGraphs vs. pt
for abs_eta_bin in abs_eta_bins:
    ptll_vs_ptll_file = os.path.join( args.inputDir, ( 'ptll_profile_ptll_for_eta_%4.3f_%4.3f'%( abs_eta_bin ) ).replace( '.','')+'.root' )
   
    ptll_mc, ptll_data = getProfiles(ptll_vs_ptll_file)
    for method in ['mpf', 'ptbal', 'mpfNoType1']:
        for alpha in ['10', '15', '20', '30']:
            response_file = os.path.join( args.inputDir, ('r_%s_%s_a%s_profile_ptll_for_eta_%4.3f_%4.3f'%( ( method, args.version, alpha ) + abs_eta_bin )).replace( '.','')+'.root' ) 

            r_mc, r_data    = getProfiles(response_file)
            n_bins = r_mc.GetNbinsX()

            mc    = ROOT.TGraphErrors(n_bins)
            mc.   SetName("MC_%s_a%s_eta_%s_%s_L1L2L3"%   ( obj_name[method], alpha, str(int(1000*abs_eta_bin[0])).zfill(4), str(int(1000*abs_eta_bin[1])).zfill(4) ))
            mc.SetTitle( "%s;%s;%s"%( mc.GetName(), r_mc.GetXaxis().GetTitle(), r_mc.GetYaxis().GetTitle()) )
            data  = ROOT.TGraphErrors(n_bins)
            data. SetName("Data_%s_a%s_eta_%s_%s_L1L2L3"% ( obj_name[method], alpha, str(int(1000*abs_eta_bin[0])).zfill(4), str(int(1000*abs_eta_bin[1])).zfill(4) ))
            data.SetTitle( "%s;%s;%s"%( data.GetName(), r_data.GetXaxis().GetTitle(), r_data.GetYaxis().GetTitle()) )
            ratio = ROOT.TGraphErrors(n_bins)
            ratio.SetName("Ratio_%s_a%s_eta_%s_%s_L1L2L3"%( obj_name[method], alpha, str(int(1000*abs_eta_bin[0])).zfill(4), str(int(1000*abs_eta_bin[1])).zfill(4) ))
            ratio.SetTitle( "Ratio %s;Data/MC %s;%s"%( data.GetName(), r_data.GetXaxis().GetTitle(), r_data.GetYaxis().GetTitle()) )

            fillTGraphs( r_data, r_mc, ptll_data, ptll_mc, ratio, data, mc) 

            stuff.append( mc )
            stuff.append( data )
            stuff.append( ratio )

#Mass plots
ptll_vs_ptll_file  = os.path.join( args.inputDir, 'ptll_profile_ptll_for_eta_0000_1300.root')
ptll_mc, ptll_data = getProfiles(ptll_vs_ptll_file)

dl_mass_file       = os.path.join( args.inputDir, 'dl_mass_profile_pt.root' )
mass_mc, mass_data = getProfiles(dl_mass_file)

mc    = ROOT.TGraphErrors(n_bins)
mc.   SetName("MC_ZMass_L1L2L3")
mc.SetTitle( "%s;%s;%s"%( mc.GetName(), mass_mc.GetXaxis().GetTitle(), mass_mc.GetYaxis().GetTitle()) )
data  = ROOT.TGraphErrors(n_bins)
data. SetName("Data_ZMass_L1L2L3")
data.SetTitle( "%s;%s;%s"%( data.GetName(), mass_data.GetXaxis().GetTitle(), mass_data.GetYaxis().GetTitle()) )
ratio = ROOT.TGraphErrors(n_bins)
ratio.SetName("Ratio_ZMass_L1L2L3")
ratio.SetTitle( "Ratio %s;Data/MC %s;%s"%( data.GetName(), mass_data.GetXaxis().GetTitle(), mass_data.GetYaxis().GetTitle()) )

fillTGraphs( mass_data, mass_mc, ptll_data, ptll_mc, ratio, data, mc) 

stuff.append( mc )
stuff.append( data )
stuff.append( ratio )

f_out = ROOT.TFile( args.outputFile, "recreate")
f_out.cd()
for o in stuff:
    o.Write()

f_out.Close()
logger.info( "Written %s", args.outputFile ) 
