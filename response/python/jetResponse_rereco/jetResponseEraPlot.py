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

argParser.add_argument('--input_path', 
      action='store',
      #nargs=1,
      type=str,
      default='/scratch/rschoefbeck/cmgTuples/jetTuples/',
      help="era"
)

argParser.add_argument('--small', 
      action='store_true',
      default=False,
      help="small?"
)

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

prefix = "prompt_era"
if args.small:
    prefix+='_small'

eras = ["Run2016B", "Run2016C", "Run2016D", "Run2016E", "Run2016F_early", "Run2016F_late", "Run2016G"]
files =  {era: [f for f in os.listdir(args.input_path) if os.path.isfile(os.path.join(args.input_path, f)) if era.replace('_late','').replace('_early','') in f] for era in eras }

abs_files = {era: [os.path.join(args.input_path, f) for f in files[era]] for era in eras }

jets = {era:Sample.fromFiles('jets', files = abs_files[era], treeName = 'jets' ) for era in eras }

jets['Run2016F_early'].setSelectionString('run<278802')
jets['Run2016F_late'].setSelectionString('run>=278802')

if args.small:
    for era in eras:
        jets[era].reduceFiles( to = 1 )

pt_thresholds = [10**(x/10.) for x in range(11,36)] 
profile_pt = {era:jets[era].get1DHistoFromDraw("rawPt_rereco/rawPt_prompt:rawPt_prompt", selectionString = "id_prompt>=1&&abs(muEF_rereco-muEF_prompt)<0.1", binning = pt_thresholds, binningIsExplicit = True, isProfile=True) for era in eras}

colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange]

profiles = [ profile_pt[era] for era in eras]
histos = [ [h.ProjectionX()] for h in profiles ]
for i_era, era in enumerate(eras):
    histos[i_era][0].style = styles.errorStyle( colors[i_era], markerSize = 0.5)
    histos[i_era][0].legendText = era

name = '_'.join([prefix, 'pt' ]) 
jetResponsePlot = Plot.fromHisto(name = name, histos = histos, texX = "prompt p_{T}", texY = "response rereco/prompt" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/80X_JetHT_rereco/", ratio = None, logY = False, logX = True, yRange=(0.95,1.07))

#eta_thresholds = [0,1,2,3,4,5] 
eta_thresholds = [ i/10. for i in range(52) ]
profile_eta   = {era:jets[era].get1DHistoFromDraw("rawPt_rereco/rawPt_prompt:abs(eta_prompt)", selectionString = "id_prompt>=1&&abs(muEF_rereco-muEF_prompt)<0.1&&rawPt_prompt>20", binning = eta_thresholds, binningIsExplicit = True, isProfile=True) for era in eras}

profiles = [ profile_eta[era] for era in eras]
histos = [ [h.ProjectionX()] for h in profiles ]
for i_era, era in enumerate(eras):
    histos[i_era][0].style = styles.errorStyle( colors[i_era], markerSize = 0.5)
    histos[i_era][0].legendText = era

name = '_'.join([prefix, 'eta' ]) 
jetResponsePlot = Plot.fromHisto(name = name, histos = histos, texX = "|prompt #eta|", texY = "response rereco/prompt" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/80X_JetHT_rereco/", ratio = None, logY = False, logX = False, yRange=(0.95,1.12), legend = (0.15,0.9-0.05*len(histos),0.47,0.93))

nvtx_thresholds = range(40) 
profile_nvtx = {era:jets[era].get1DHistoFromDraw("rawPt_rereco/rawPt_prompt:nVert", selectionString = "id_prompt>=1&&abs(muEF_rereco-muEF_prompt)<0.1&&rawPt_prompt>20", binning = nvtx_thresholds, binningIsExplicit = True, isProfile=True) for era in eras}

profiles = [ profile_nvtx[era] for era in eras]
histos = [ [h.ProjectionX()] for h in profiles ]
for i_era, era in enumerate(eras):
    histos[i_era][0].style = styles.errorStyle( colors[i_era], markerSize = 0.5)
    histos[i_era][0].legendText = era

name = '_'.join([prefix, 'nvtx' ]) 
jetResponsePlot = Plot.fromHisto(name = name, histos = histos, texX = "vertex multiplicity", texY = "response rereco/prompt" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/80X_JetHT_rereco/", ratio = None, logY = False, logX = False, yRange=(0.95,1.07))
