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
      choices=['Run2016B', 'Run2016C', 'Run2016D', 'Run2016E', 'Run2016F', 'Run2016G'],
      default='Run2016G',
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

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

all_files = [f for f in os.listdir(args.input_path) if os.path.isfile(os.path.join(args.input_path, f))]

files = filter( lambda f: args.era in f, all_files)
if args.run>0:
    files = filter( lambda f: str(args.run) in f, files)

if len(files)==0:
    print "No files found for era %s and run %i." %( args.era, args.run )
    sys.exit(0)
else:
    print "Found %i files for era %s and run %i." %( len(files), args.era, args.run )

abs_files = [os.path.join(args.input_path, f) for f in files]

jets = Sample.fromFiles('jets', files = abs_files, treeName = 'jets' )


thresholds = [10**(x/10.) for x in range(11,36)] 
profile = jets.get1DHistoFromDraw("rawPt_rereco/rawPt_prompt:0.5*(rawPt_rereco+rawPt_prompt)", binning = thresholds, binningIsExplicit = True, isProfile=True)

profiles = [ profile ]
#prefix=preprefix+"maxN_%s_"%maxN if maxN is not None and maxN>0 else ""
histos = [ [h.ProjectionX()] for h in profiles ]
prefix = "test"
name = '_'.join([prefix, args.era, 'run_%i'%args.run ]) 
jetResponsePlot = Plot.fromHisto(name = name, histos = histos, texX = "avg. p_{T}", texY = "response rereco/prompt" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = True, yRange=(0.95,1.07))

