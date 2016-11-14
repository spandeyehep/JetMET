''' Make simple jet trees for rereco jet response comparison 
'''

# Standard imports
import sys
import os
import logging
import ROOT

#RootTools
from RootTools.core.standard import *

# argParser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel', 
      action='store',
      nargs='?',
      choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],
      default='INFO',
      help="Log level for logging"
)

argParser.add_argument('--era', 
      action='store',
      #nargs=1,
      type=str,
      choices=['Run2016B', 'Run2016C', 'Run2016D', 'Run2016E', 'Run2016F', 'Run2016F_early', 'Run2016F_late', 'Run2016G'],
      default='Run2016B',
      help="era"
)

argParser.add_argument('--input_path', 
      action='store',
      #nargs=1,
      type=str,
      default='/scratch/rschoefbeck/cmgTuples/jetTuples/',
      help="era"
)

argParser.add_argument('--run', 
      action='store',
      type = int,
      default=-1,
#      default=280187,
      help="run"
)

argParser.add_argument('--small', 
      action='store_true',
      # default=True,
      help="small?"
)

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

prefix = "clean"
if args.small:
    prefix+='_small'
all_files = [f for f in os.listdir(args.input_path) if os.path.isfile(os.path.join(args.input_path, f))]

files = filter( lambda f: args.era.replace('_early','').replace('_late','') in f, all_files )
if args.run>0:
    files = filter( lambda f: str(args.run) in f, files)

if len(files)==0:
    print "No files found for era %s and run %i." %( args.era, args.run )
    sys.exit(0)
else:
    print "Found %i files for era %s and run %i." %( len(files), args.era, args.run )

abs_files = [os.path.join(args.input_path, f) for f in files]

jets = Sample.fromFiles('jets', files = abs_files, treeName = 'jets' )

if args.era == 'Run2016F_early':
    jets.setSelectionString('run<278802')
elif args.era == 'Run2016F_late':
    jets.setSelectionString('run>=278802')

if args.small:
    jets.reduceFiles( to = 1 )

import numpy as np
from math import exp, log
pt_max = 10000
pt_min = 20
n_bin = 200
pt_thresholds = [exp(x) for x in np.arange(log(pt_min),log(pt_max), (log(pt_max)-log(pt_min))/n_bin) ]

pt_2D = jets.get2DHistoFromDraw("rawPt_rereco:rawPt_prompt", selectionString = "id_prompt>=1&&abs(muEF_prompt-muEF_rereco)<0.1", binning = ( pt_thresholds, pt_thresholds ), binningIsExplicit = True)

name = '_'.join([prefix, args.era, 'pt_2D' ]) 
#p_pt_2D = Plot2D.fromHisto(name = name, histos = [[pt_2D]], texX = "prompt p_{T}", texY = "rereco p_{T}" )
c1 = ROOT.TCanvas()
pt_2D.SetTitle("")
pt_2D.Draw('COLZ')
pt_2D.GetYaxis().SetTitle("rereco p_{T}")
c1.SetLogy()
pt_2D.GetXaxis().SetTitle("prompt p_{T}")
c1.SetLogx()
c1.SetLogz()
ROOT.gStyle.SetOptStat(0)
c1.Print(os.path.join( "/afs/hephy.at/user/r/rschoefbeck/www/80X_JetHT_rereco/", name+'.png' ))


