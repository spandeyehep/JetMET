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
argParser.add_argument('--inputDir',           action='store',      default='/afs/hephy.at/user/r/rschoefbeck/www/JEC/L3res_small/V6/DYnJets/inclusive/mumu_log/ptll30-njet1p/',      help="input directory")
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


def getPTllProfiles( fname ):
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
for abs_eta_bin in abs_eta_bins:
    ptll_vs_ptll_file = os.path.join( args.inputDir, 'ptll_profile_ptll_for_eta_%3.2f_%3.2f.root'%( abs_eta_bin ) )
   
    mc_ptll, data_ptll = getPTllProfiles(ptll_vs_ptll_file)
    for method in ['mpf', 'ptbal', 'mpfNoType1']:
        response_file = os.path.join( args.inputDir, 'r_%s_%s_profile_ptll_for_eta_%3.2f_%3.2f.root'%( ( method, args.version ) + abs_eta_bin ) ) 
        
        
