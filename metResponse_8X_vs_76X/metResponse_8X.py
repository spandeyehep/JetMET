''' FWLiteReader example: Loop over a sample and write some data to a histogram.
'''
# Standard imports
import os
import logging
import ROOT
import array

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

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

#Histograms
pt_thresholds = [10**(x/10.) for x in range(11,36)]

metScatter = ROOT.TH2F("METScatter", "METScatter", 100,0,100,100,0,100 )
MExy_8X = ROOT.TH1F("MExy_8X", "MExy_8X", 100,-100,100 )
MExy_8X.texName = "MExy 800pre6"
MExy_8X.style = styles.lineStyle(ROOT.kBlack)
MExy_16Dec = ROOT.TH1F("MExy_16Dec", "MExy_16Dec", 100,-100,100 )
MExy_16Dec.texName = "MExy 76X"
MExy_16Dec.style = styles.lineStyle(ROOT.kBlue)

import files 
maxN = 25
# 8X mAOD, assumes eos mount in home directory 
#dirname = "~/eos/cms/store/relval/CMSSW_8_0_0_pre6/SingleMuon/MINIAOD/80X_dataRun2_v4_RelVal_jetHT2015HLHT-v1/10000/"
SingleMuon_8X_260627    = FWLiteSample.fromDAS("SingleMuon_8X_260627", "/SingleMuon/CMSSW_8_0_0_pre6-80X_dataRun2_v4_RelVal_sigMu2015HLHT-v1/MINIAOD", instance = "global", maxN = maxN)
SingleMuon_16Dec_260627 = FWLiteSample.fromFiles("SingleMuon_16Dec_260627", files = ["root://eoscms.cern.ch/"+s for s in files.SingleMuon_16DecRereco_260627], maxN = maxN)
logger.info("Loading of samples is done.")

products = {
    'METs':{'type':'vector<pat::MET>', 'label':("slimmedMETs")},
    'Muons':{'type':'vector<pat::Muon>', 'label':("slimmedMuons")}
    }

r1 = SingleMuon_8X_260627.fwliteReader( products = products )
r2 = SingleMuon_16Dec_260627.fwliteReader( products = products  )

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

# Withou sorting, there is a jump between files with almost every event -> extremly slow
positions.sort()
logger.info("Have %i events in common.", len(intersec))

#Looping over common events
for i, p in enumerate(positions):
    p1,p2 = p
    r1.goToPosition(p1)
    r2.goToPosition(p2)
    if i%10000==0: logger.info("At %i/%i of common events.", i, len(positions))

    muons = r1.products['Muons']
    if muons.size()>=2: 
        pZ =  (muons[0].p4() + muons[1].p4())
        if not abs(91.2 - pZ.mass())<15.: continue

        MExy_8X.Fill(r1.products['METs'][0].px())
        MExy_8X.Fill(r1.products['METs'][0].py())
        MExy_16Dec.Fill(r2.products['METs'][0].px())
        MExy_16Dec.Fill(r2.products['METs'][0].py())
        metScatter.Fill(r1.products['METs'][0].pt(), r2.products['METs'][0].pt())

prefix=""        
MExy_plot = Plot.fromHisto(name = prefix+"MExy", histos =  [[MExy_16Dec],[MExy_8X]] , texX = "ME_{x,y}" )  
plotting.draw(MExy_plot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", logX = False, logY = True)

c1 = ROOT.TCanvas()
metScatter.Draw("COLZ")
l=ROOT.TLine(0,0,100,100)
l.Draw()
c1.Print("/afs/hephy.at/user/r/rschoefbeck/www/etc/metScatter.png")
