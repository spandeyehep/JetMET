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

argParser.add_argument('--small', 
      action='store_true',
      default=False,
      help="small?"
)

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

prefix = "prompt"
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

pt_thresholds = [10**(x/10.) for x in range(11,36)] 
profile_pt_eta_inc = jets.get1DHistoFromDraw("rawPt_rereco/rawPt_prompt:rawPt_prompt", selectionString = "id_prompt>=1&&abs(muEF_rereco-muEF_prompt)<0.1", binning = pt_thresholds, binningIsExplicit = True, isProfile=True)
profile_pt_eta_barrel = jets.get1DHistoFromDraw("rawPt_rereco/rawPt_prompt:rawPt_prompt", selectionString = "id_prompt>=1&&abs(muEF_rereco-muEF_prompt)<0.1&&abs(eta_prompt)<1.3", binning = pt_thresholds, binningIsExplicit = True, isProfile=True)
profile_pt_eta_endcap = jets.get1DHistoFromDraw("rawPt_rereco/rawPt_prompt:rawPt_prompt", selectionString = "id_prompt>=1&&abs(muEF_rereco-muEF_prompt)<0.1&&abs(eta_prompt)>1.3 && abs(eta_prompt)<3", binning = pt_thresholds, binningIsExplicit = True, isProfile=True)
profile_pt_eta_hf = jets.get1DHistoFromDraw("rawPt_rereco/rawPt_prompt:rawPt_prompt", selectionString = "id_prompt>=1&&abs(muEF_rereco-muEF_prompt)<0.1&&abs(eta_prompt)>3", binning = pt_thresholds, binningIsExplicit = True, isProfile=True)

profiles = [ profile_pt_eta_inc, profile_pt_eta_barrel, profile_pt_eta_endcap, profile_pt_eta_hf]
histos = [ [h.ProjectionX()] for h in profiles ]
histos[0][0].style = styles.errorStyle( ROOT.kBlack, markerSize = 0.5)
histos[0][0].legendText = "all"
histos[1][0].style = styles.errorStyle( ROOT.kRed, markerSize = 0.5)
histos[1][0].legendText = "barrel |#eta| < 1.3 "
histos[2][0].style = styles.errorStyle( ROOT.kBlue, markerSize = 0.5)
histos[2][0].legendText = "endcap 1.3 < |#eta| < 3"
histos[3][0].style = styles.errorStyle( ROOT.kGreen, markerSize = 0.5)
histos[3][0].legendText = "forward |#eta| > 3"

name = '_'.join([prefix, args.era, 'pt' ]) 
jetResponsePlot = Plot.fromHisto(name = name, histos = histos, texX = "prompt p_{T}", texY = "response rereco/prompt" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/80X_JetHT_rereco/", ratio = None, logY = False, logX = True, yRange=(0.95,1.07))

#eta_thresholds = [0,1,2,3,4,5] 
eta_thresholds = [ i/10. for i in range(52) ]
profile_eta_pt_inc   = jets.get1DHistoFromDraw("rawPt_rereco/rawPt_prompt:abs(eta_prompt)", selectionString = "id_prompt>=1&&abs(muEF_rereco-muEF_prompt)<0.1&&rawPt_prompt>20", binning = eta_thresholds, binningIsExplicit = True, isProfile=True)
profile_eta_pt_50    = jets.get1DHistoFromDraw("rawPt_rereco/rawPt_prompt:abs(eta_prompt)", selectionString = "id_prompt>=1&&abs(muEF_rereco-muEF_prompt)<0.1&&rawPt_prompt<50&&rawPt_prompt>20", binning = eta_thresholds, binningIsExplicit = True, isProfile=True)
profile_eta_pt_100   = jets.get1DHistoFromDraw("rawPt_rereco/rawPt_prompt:abs(eta_prompt)", selectionString = "id_prompt>=1&&abs(muEF_rereco-muEF_prompt)<0.1&&rawPt_prompt<100&&rawPt_prompt>50", binning = eta_thresholds, binningIsExplicit = True, isProfile=True)
profile_eta_pt_300   = jets.get1DHistoFromDraw("rawPt_rereco/rawPt_prompt:abs(eta_prompt)", selectionString = "id_prompt>=1&&abs(muEF_rereco-muEF_prompt)<0.1&&rawPt_prompt<300&&rawPt_prompt>100", binning = eta_thresholds, binningIsExplicit = True, isProfile=True)
profile_eta_pt_gr300 = jets.get1DHistoFromDraw("rawPt_rereco/rawPt_prompt:abs(eta_prompt)", selectionString = "id_prompt>=1&&abs(muEF_rereco-muEF_prompt)<0.1&&rawPt_prompt>300", binning = eta_thresholds, binningIsExplicit = True, isProfile=True)

profiles = [ profile_eta_pt_inc, profile_eta_pt_50, profile_eta_pt_100, profile_eta_pt_300, profile_eta_pt_gr300]
histos = [ [h.ProjectionX()] for h in profiles ]
histos[0][0].style = styles.errorStyle( ROOT.kBlack, markerSize = 0.5)
histos[0][0].legendText = "all"
histos[1][0].style = styles.errorStyle( ROOT.kRed, markerSize = 0.5)
histos[1][0].legendText = "20 < p_{T} < 50"
histos[2][0].style = styles.errorStyle( ROOT.kBlue, markerSize = 0.5)
histos[2][0].legendText = "50 < p_{T} < 100"
histos[3][0].style = styles.errorStyle( ROOT.kGreen, markerSize = 0.5)
histos[3][0].legendText = "100 < p_{T} < 300"
histos[4][0].style = styles.errorStyle( ROOT.kMagenta, markerSize = 0.5)
histos[4][0].legendText = "p_{T} > 300"

name = '_'.join([prefix, args.era, 'eta' ]) 
jetResponsePlot = Plot.fromHisto(name = name, histos = histos, texX = "|prompt #eta|", texY = "response rereco/prompt" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/80X_JetHT_rereco/", ratio = None, logY = False, logX = False, yRange=(0.95,1.12), legend = (0.15,0.9-0.05*len(histos),0.47,0.93))

