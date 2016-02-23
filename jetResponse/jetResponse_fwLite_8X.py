''' FWLiteReader example: Loop over a sample and write some data to a histogram.
'''
# Standard imports
import os
import logging
import ROOT
import array

#RootTools
from RootTools.core.Variable import Variable, ScalarType, VectorType
from RootTools.core.logger import get_logger
from RootTools.fwlite.FWLiteSample import FWLiteSample
from RootTools.fwlite.FWLiteReader import FWLiteReader
from RootTools.plot.Plot import Plot
from RootTools.plot.styles import lineStyle
import RootTools.plot.plotting as plotting

#Helper
from helpers import deltaR2

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

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

#Histograms
thresholds = [10**(x/10.) for x in range(11,36)]
jetResponse = ROOT.TProfile("response", "response", len(thresholds)-1, array.array('d', thresholds) )
jetResponse.texName = "JetHT 260627" 
jetResponse.style = lineStyle(ROOT.kBlue)

import files 
# 8X mAOD, assumes eos mount in home directory 
maxN = -1
#dirname = "~/eos/cms/store/relval/CMSSW_8_0_0_pre6/JetHT/MINIAOD/80X_dataRun2_v4_RelVal_jetHT2015HLHT-v1/10000/"
JetHT_8X_260627 = FWLiteSample.fromFiles("JetHT_8X_260627", files = ["root://eoscms.cern.ch/"+s for s in files.JetHT_8X_260627], maxN = maxN)
maxN = -1
#JetHT_16Dec_260627 = FWLiteSample.fromFiles("JetHT_16Dec_260627", files = [os.path.expanduser("~/eos/cms"+s) for s in JetHT_16DecRereco_260627], maxN = maxN)
JetHT_16Dec_260627 = FWLiteSample.fromFiles("JetHT_16Dec_260627", files = ["root://eoscms.cern.ch/"+s for s in files.JetHT_16DecRereco_260627], maxN = maxN)

products = {
    'jets':{'type':'vector<pat::Jet>', 'label':("slimmedJets")} 
    }

r1 = JetHT_8X_260627.fwliteReader( products = products )
r2 = JetHT_16Dec_260627.fwliteReader( products = products  )

r1.start()
runs_1 = set()
position_r1 = {}
while r1.run( readProducts = False ):
    if r1.evt[0] == 260627:
        position_r1[r1.evt] = r1.position-1

r2.start()
runs_2 = set()
position_r2 = {}
while r2.run( readProducts = False ):
    if r2.evt[0] == 260627:
        position_r2[r2.evt] = r2.position-1

logger.info( "Have %i events in first samle and %i in second", len(position_r1), len(position_r2) )

# Fast intersect
intersec = set(position_r1.keys()).intersection(set(position_r2.keys()))
positions = [(position_r1[i], position_r2[i]) for i in intersec]

# Without sorting, there is a jump between files with almost every event -> extremly slow
positions.sort()
logger.info("Have %i events in common.", len(intersec))

#Looping over common events
for i, p in enumerate(positions):
    p1,p2 = p
    r1.goToPosition(p1)
    r2.goToPosition(p2)
    if i%10000==0: logger.info("At %i/%i of common events.", i, len(positions))
        
    jets1_ = [ j.correctedJet("Uncorrected") for j in r1.products['jets'] ]
    jets2_ = [ j.correctedJet("Uncorrected") for j in r2.products['jets'] ]

    jets1 = [{'pt':j.pt(), 'eta':j.eta(), 'phi':j.phi()} for j in jets1_]
    jets2 = [{'pt':j.pt(), 'eta':j.eta(), 'phi':j.phi()} for j in jets2_]

    for c in zip(jets1, jets2):
        if deltaR2(*c)<0.2**2 and abs(c[0]['eta'])<1.3 and abs(c[1]['eta'])<1.3:
            jetResponse.Fill(c[1]['pt'], c[0]['pt']/c[1]['pt'] )
#    print "jets1",jets1
#    print "jets2",jets2

## Make plot
prefix=""
jetResponsePlot = Plot.fromHisto(name = prefix+"pT_ratio_800pre6_76X", histos = [[jetResponse]], texX = "raw Jet p_{T} 76X", texY = "response ratio 800pre6/76X" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", logX = True, logY = False, sorting = False, yRange = (0.965, 1.045) )


