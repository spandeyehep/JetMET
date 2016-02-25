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
from RootTools.plot.styles import lineStyle, fillStyle
import RootTools.plot.plotting as plotting
import RootTools.core.helpers as helpers

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
pt_thresholds = [10**(x/10.) for x in range(11,36)]

ratio_jetResponse = ROOT.TProfile("response", "response", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
ratio_jetResponse.style = lineStyle(ROOT.kBlue)
ratio_jetResponse.legendText = "JetHT 260627" 

ratio_chargedEmEnergy = ROOT.TProfile("chargedEmEnergy", "chargedEmEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
ratio_chargedEmEnergy.style = lineStyle(ROOT.kBlue)

ratio_chargedHadronEnergy = ROOT.TProfile("chargedHadronEnergy", "chargedHadronEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
ratio_chargedHadronEnergy.style = lineStyle(ROOT.kRed)

ratio_chargedHadronMultiplicity = ROOT.TProfile("chargedHadronMultiplicity", "chargedHadronMultiplicity", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
ratio_chargedHadronMultiplicity.style = lineStyle(ROOT.kGreen)

ratio_neutralEmEnergy = ROOT.TProfile("neutralEmEnergy", "neutralEmEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
ratio_neutralEmEnergy.style = lineStyle(ROOT.kBlack)

ratio_neutralHadronEnergy = ROOT.TProfile("neutralHadronEnergy", "neutralHadronEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
ratio_neutralHadronEnergy.style = lineStyle(ROOT.kMagenta)

ratio_neutralHadronMultiplicity = ROOT.TProfile("neutralHadronMultiplicity", "neutralHadronMultiplicity", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
ratio_neutralHadronMultiplicity.style = lineStyle(ROOT.kMagenta)

pt_8X_chargedHadronEnergy = ROOT.TProfile("pt_8X_chargedHadronEnergy", "pt_8X_chargedHadronEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_8X_chargedHadronEnergy.style = lineStyle( ROOT.kBlack )

pt_8X_neutralHadronEnergy = ROOT.TProfile("pt_8X_neutralHadronEnergy", "pt_8X_neutralHadronEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_8X_neutralHadronEnergy.style = lineStyle( ROOT.kBlack )

pt_8X_photonEnergy = ROOT.TProfile("pt_8X_photonEnergy", "pt_8X_photonEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_8X_photonEnergy.style = lineStyle( ROOT.kBlack )

pt_8X_electronEnergy = ROOT.TProfile("pt_8X_electronEnergy", "pt_8X_electronEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_8X_electronEnergy.style = lineStyle( ROOT.kBlack )

pt_8X_muonEnergy = ROOT.TProfile("pt_8X_muonEnergy", "pt_8X_muonEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_8X_muonEnergy.style = lineStyle( ROOT.kBlack )

pt_8X_HFHadronEnergy = ROOT.TProfile("pt_8X_HFHadronEnergy", "pt_8X_HFHadronEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_8X_HFHadronEnergy.style = lineStyle( ROOT.kBlack )

pt_8X_HFEMEnergy = ROOT.TProfile("pt_8X_HFEMEnergy", "pt_8X_HFEMEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_8X_HFEMEnergy.style = lineStyle( ROOT.kBlack )

pt_76X_chargedHadronEnergy = ROOT.TProfile("pt_76X_chargedHadronEnergy", "pt_76X_chargedHadronEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_76X_chargedHadronEnergy.legendText = "charged hadrons"
pt_76X_chargedHadronEnergy.style = fillStyle(ROOT.kRed - 4, lineColor = None)

pt_76X_neutralHadronEnergy = ROOT.TProfile("pt_76X_neutralHadronEnergy", "pt_76X_neutralHadronEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_76X_neutralHadronEnergy.legendText = "neutral hadrons"
pt_76X_neutralHadronEnergy.style = fillStyle(ROOT.kGreen - 3, lineColor = None)

pt_76X_photonEnergy = ROOT.TProfile("pt_76X_photonEnergy", "pt_76X_photonEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_76X_photonEnergy.legendText = "photons"
pt_76X_photonEnergy.style = fillStyle(ROOT.kBlue - 7, lineColor = None)

pt_76X_electronEnergy = ROOT.TProfile("pt_76X_electronEnergy", "pt_76X_electronEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_76X_electronEnergy.style = fillStyle(ROOT.kCyan - 2, lineColor = None)

pt_76X_muonEnergy = ROOT.TProfile("pt_76X_muonEnergy", "pt_76X_muonEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_76X_muonEnergy.legendText = "leptons"
pt_76X_muonEnergy.style = fillStyle(ROOT.kCyan - 2, lineColor = None)

pt_76X_HFHadronEnergy = ROOT.TProfile("pt_76X_HFHadronEnergy", "pt_76X_HFHadronEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_76X_HFHadronEnergy.legendText = "HF hadrons"
pt_76X_HFHadronEnergy.style = fillStyle(ROOT.kMagenta, lineColor = None)

pt_76X_HFEMEnergy = ROOT.TProfile("pt_76X_HFEMEnergy", "pt_76X_HFEMEnergy", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
pt_76X_HFEMEnergy.legendText = "HF EM"
pt_76X_HFEMEnergy.style = fillStyle(ROOT.kGray, lineColor = None)


eta_thresholds = [ x*5.2/10. for x in range(-10,11) ]

eta_8X_chargedHadronEnergy = ROOT.TProfile("eta_8X_chargedHadronEnergy", "eta_8X_chargedHadronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_8X_chargedHadronEnergy.style = lineStyle( ROOT.kBlack )

eta_8X_neutralHadronEnergy = ROOT.TProfile("eta_8X_neutralHadronEnergy", "eta_8X_neutralHadronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_8X_neutralHadronEnergy.style = lineStyle( ROOT.kBlack )

eta_8X_photonEnergy = ROOT.TProfile("eta_8X_photonEnergy", "eta_8X_photonEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_8X_photonEnergy.style = lineStyle( ROOT.kBlack )

eta_8X_electronEnergy = ROOT.TProfile("eta_8X_electronEnergy", "eta_8X_electronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_8X_electronEnergy.style = lineStyle( ROOT.kBlack )

eta_8X_muonEnergy = ROOT.TProfile("eta_8X_muonEnergy", "eta_8X_muonEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_8X_muonEnergy.style = lineStyle( ROOT.kBlack )

eta_8X_HFHadronEnergy = ROOT.TProfile("eta_8X_HFHadronEnergy", "eta_8X_HFHadronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_8X_HFHadronEnergy.style = lineStyle( ROOT.kBlack )

eta_8X_HFEMEnergy = ROOT.TProfile("eta_8X_HFEMEnergy", "eta_8X_HFEMEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_8X_HFEMEnergy.style = lineStyle( ROOT.kBlack )

eta_76X_chargedHadronEnergy = ROOT.TProfile("eta_76X_chargedHadronEnergy", "eta_76X_chargedHadronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_76X_chargedHadronEnergy.legendText = "charged hadrons"
eta_76X_chargedHadronEnergy.style = fillStyle(ROOT.kRed - 4, lineColor = None)

eta_76X_neutralHadronEnergy = ROOT.TProfile("eta_76X_neutralHadronEnergy", "eta_76X_neutralHadronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_76X_neutralHadronEnergy.legendText = "neutral hadrons"
eta_76X_neutralHadronEnergy.style = fillStyle(ROOT.kGreen - 3, lineColor = None)

eta_76X_photonEnergy = ROOT.TProfile("eta_76X_photonEnergy", "eta_76X_photonEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_76X_photonEnergy.legendText = "photons"
eta_76X_photonEnergy.style = fillStyle(ROOT.kBlue - 7, lineColor = None)

eta_76X_electronEnergy = ROOT.TProfile("eta_76X_electronEnergy", "eta_76X_electronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_76X_electronEnergy.style = fillStyle(ROOT.kCyan - 2, lineColor = None)

eta_76X_muonEnergy = ROOT.TProfile("eta_76X_muonEnergy", "eta_76X_muonEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_76X_muonEnergy.legendText = "leptons"
eta_76X_muonEnergy.style = fillStyle(ROOT.kCyan - 2, lineColor = None)

eta_76X_HFHadronEnergy = ROOT.TProfile("eta_76X_HFHadronEnergy", "eta_76X_HFHadronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_76X_HFHadronEnergy.legendText = "HF hadrons"
eta_76X_HFHadronEnergy.style = fillStyle(ROOT.kMagenta, lineColor = None)

eta_76X_HFEMEnergy = ROOT.TProfile("eta_76X_HFEMEnergy", "eta_76X_HFEMEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_76X_HFEMEnergy.legendText = "HF EM"
eta_76X_HFEMEnergy.style = fillStyle(ROOT.kGray, lineColor = None)

eta_lowPt_8X_chargedHadronEnergy = ROOT.TProfile("eta_lowPt_8X_chargedHadronEnergy", "eta_lowPt_8X_chargedHadronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_8X_chargedHadronEnergy.style = lineStyle( ROOT.kBlack )

eta_lowPt_8X_neutralHadronEnergy = ROOT.TProfile("eta_lowPt_8X_neutralHadronEnergy", "eta_lowPt_8X_neutralHadronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_8X_neutralHadronEnergy.style = lineStyle( ROOT.kBlack )

eta_lowPt_8X_photonEnergy = ROOT.TProfile("eta_lowPt_8X_photonEnergy", "eta_lowPt_8X_photonEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_8X_photonEnergy.style = lineStyle( ROOT.kBlack )

eta_lowPt_8X_electronEnergy = ROOT.TProfile("eta_lowPt_8X_electronEnergy", "eta_lowPt_8X_electronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_8X_electronEnergy.style = lineStyle( ROOT.kBlack )

eta_lowPt_8X_muonEnergy = ROOT.TProfile("eta_lowPt_8X_muonEnergy", "eta_lowPt_8X_muonEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_8X_muonEnergy.style = lineStyle( ROOT.kBlack )

eta_lowPt_8X_HFHadronEnergy = ROOT.TProfile("eta_lowPt_8X_HFHadronEnergy", "eta_lowPt_8X_HFHadronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_8X_HFHadronEnergy.style = lineStyle( ROOT.kBlack )

eta_lowPt_8X_HFEMEnergy = ROOT.TProfile("eta_lowPt_8X_HFEMEnergy", "eta_lowPt_8X_HFEMEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_8X_HFEMEnergy.style = lineStyle( ROOT.kBlack )

eta_lowPt_76X_chargedHadronEnergy = ROOT.TProfile("eta_lowPt_76X_chargedHadronEnergy", "eta_lowPt_76X_chargedHadronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_76X_chargedHadronEnergy.legendText = "charged hadrons"
eta_lowPt_76X_chargedHadronEnergy.style = fillStyle(ROOT.kRed - 4, lineColor = None)

eta_lowPt_76X_neutralHadronEnergy = ROOT.TProfile("eta_lowPt_76X_neutralHadronEnergy", "eta_lowPt_76X_neutralHadronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_76X_neutralHadronEnergy.legendText = "neutral hadrons"
eta_lowPt_76X_neutralHadronEnergy.style = fillStyle(ROOT.kGreen - 3, lineColor = None)

eta_lowPt_76X_photonEnergy = ROOT.TProfile("eta_lowPt_76X_photonEnergy", "eta_lowPt_76X_photonEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_76X_photonEnergy.legendText = "photons"
eta_lowPt_76X_photonEnergy.style = fillStyle(ROOT.kBlue - 7, lineColor = None)

eta_lowPt_76X_electronEnergy = ROOT.TProfile("eta_lowPt_76X_electronEnergy", "eta_lowPt_76X_electronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_76X_electronEnergy.style = fillStyle(ROOT.kCyan - 2, lineColor = None)

eta_lowPt_76X_muonEnergy = ROOT.TProfile("eta_lowPt_76X_muonEnergy", "eta_lowPt_76X_muonEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_76X_muonEnergy.legendText = "leptons"
eta_lowPt_76X_muonEnergy.style = fillStyle(ROOT.kCyan - 2, lineColor = None)

eta_lowPt_76X_HFHadronEnergy = ROOT.TProfile("eta_lowPt_76X_HFHadronEnergy", "eta_lowPt_76X_HFHadronEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_76X_HFHadronEnergy.legendText = "HF hadrons"
eta_lowPt_76X_HFHadronEnergy.style = fillStyle(ROOT.kMagenta, lineColor = None)

eta_lowPt_76X_HFEMEnergy = ROOT.TProfile("eta_lowPt_76X_HFEMEnergy", "eta_lowPt_76X_HFEMEnergy", len(eta_thresholds)-1, array.array('d', eta_thresholds) )
eta_lowPt_76X_HFEMEnergy.legendText = "HF EM"
eta_lowPt_76X_HFEMEnergy.style = fillStyle(ROOT.kGray, lineColor = None)


import files 
maxN = 15
# 8X mAOD, assumes eos mount in home directory 
#dirname = "~/eos/cms/store/relval/CMSSW_8_0_0_pre6/JetHT/MINIAOD/80X_dataRun2_v4_RelVal_jetHT2015HLHT-v1/10000/"
JetHT_8X_260627 = FWLiteSample.fromFiles("JetHT_8X_260627", files = ["root://eoscms.cern.ch/"+s for s in files.JetHT_8X_260627], maxN = maxN)
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

    jets1 = [{'pt':j.pt(), 'eta':j.eta(), 'phi':j.phi(), 'j':j} for j in jets1_]
    jets2 = [{'pt':j.pt(), 'eta':j.eta(), 'phi':j.phi(), 'j':j} for j in jets2_]

    for c in zip(jets1, jets2):
        if deltaR2(*c)<0.2**2:

            if abs(c[0]['eta'])<1.3 and abs(c[1]['eta'])<1.3:
                ratio_jetResponse.Fill(c[1]['pt'], c[0]['pt']/c[1]['pt'] )

                if c[1]['j'].chargedEmEnergy()>0:       
                    ratio_chargedEmEnergy.Fill(c[1]['pt'], c[0]['j'].chargedEmEnergy()/c[1]['j'].chargedEmEnergy() )
                if c[1]['j'].chargedHadronEnergy() >0: 
                    ratio_chargedHadronEnergy.Fill(c[1]['pt'], c[0]['j'].chargedHadronEnergy()/c[1]['j'].chargedHadronEnergy() )
                if c[1]['j'].chargedHadronMultiplicity() >0: 
                    ratio_chargedHadronMultiplicity.Fill(c[1]['pt'], c[0]['j'].chargedHadronMultiplicity()/c[1]['j'].chargedHadronMultiplicity() )
                if c[1]['j'].neutralEmEnergy() >0: 
                    ratio_neutralEmEnergy.Fill(c[1]['pt'], c[0]['j'].neutralEmEnergy()/c[1]['j'].neutralEmEnergy() )
                if c[1]['j'].neutralHadronEnergy() >0: 
                    ratio_neutralHadronEnergy.Fill(c[1]['pt'], c[0]['j'].neutralHadronEnergy()/c[1]['j'].neutralHadronEnergy() )
                if c[1]['j'].neutralHadronMultiplicity() >0: 
                    ratio_neutralHadronMultiplicity.Fill(c[1]['pt'], c[0]['j'].neutralHadronMultiplicity()/c[1]['j'].neutralHadronMultiplicity() )

            if abs(c[0]['eta'])<1.3 and abs(c[1]['eta'])<1.3:
                pt_8X_chargedHadronEnergy.Fill(c[1]['pt'], c[1]['j'].chargedHadronEnergy())
                pt_8X_neutralHadronEnergy.Fill(c[1]['pt'], c[1]['j'].neutralHadronEnergy() - c[1]['j'].HFHadronEnergy())
                pt_8X_photonEnergy.Fill(c[1]['pt'], c[1]['j'].photonEnergy())
                pt_8X_electronEnergy.Fill(c[1]['pt'], c[1]['j'].electronEnergy())
                pt_8X_muonEnergy.Fill(c[1]['pt'], c[1]['j'].muonEnergy())
                pt_8X_HFHadronEnergy.Fill(c[1]['pt'], c[1]['j'].HFHadronEnergy())
                pt_8X_HFEMEnergy.Fill(c[1]['pt'], c[1]['j'].HFEMEnergy())
                pt_76X_chargedHadronEnergy.Fill(c[1]['pt'], c[0]['j'].chargedHadronEnergy())
                pt_76X_neutralHadronEnergy.Fill(c[1]['pt'], c[0]['j'].neutralHadronEnergy() - c[0]['j'].HFHadronEnergy())
                pt_76X_photonEnergy.Fill(c[1]['pt'], c[0]['j'].photonEnergy())
                pt_76X_electronEnergy.Fill(c[1]['pt'], c[0]['j'].electronEnergy())
                pt_76X_muonEnergy.Fill(c[1]['pt'], c[0]['j'].muonEnergy())
                pt_76X_HFHadronEnergy.Fill(c[1]['pt'], c[0]['j'].HFHadronEnergy())
                pt_76X_HFEMEnergy.Fill(c[1]['pt'], c[0]['j'].HFEMEnergy())
            if abs(c[0]['pt'])>30 and abs(c[1]['pt'])>30:
                eta_8X_chargedHadronEnergy.Fill(c[1]['eta'], c[1]['j'].chargedHadronEnergy())
                eta_8X_neutralHadronEnergy.Fill(c[1]['eta'], c[1]['j'].neutralHadronEnergy() - c[1]['j'].HFHadronEnergy())
                eta_8X_photonEnergy.Fill(c[1]['eta'], c[1]['j'].photonEnergy())
                eta_8X_electronEnergy.Fill(c[1]['eta'], c[1]['j'].electronEnergy())
                eta_8X_muonEnergy.Fill(c[1]['eta'], c[1]['j'].muonEnergy())
                eta_8X_HFHadronEnergy.Fill(c[1]['eta'], c[1]['j'].HFHadronEnergy())
                eta_8X_HFEMEnergy.Fill(c[1]['eta'], c[1]['j'].HFEMEnergy())
                eta_76X_chargedHadronEnergy.Fill(c[1]['eta'], c[0]['j'].chargedHadronEnergy())
                eta_76X_neutralHadronEnergy.Fill(c[1]['eta'], c[0]['j'].neutralHadronEnergy() - c[0]['j'].HFHadronEnergy())
                eta_76X_photonEnergy.Fill(c[1]['eta'], c[0]['j'].photonEnergy())
                eta_76X_electronEnergy.Fill(c[1]['eta'], c[0]['j'].electronEnergy())
                eta_76X_muonEnergy.Fill(c[1]['eta'], c[0]['j'].muonEnergy())
                eta_76X_HFHadronEnergy.Fill(c[1]['eta'], c[0]['j'].HFHadronEnergy())
                eta_76X_HFEMEnergy.Fill(c[1]['eta'], c[0]['j'].HFEMEnergy())
            if not (abs(c[0]['pt'])>30 and abs(c[1]['pt'])>30):
                eta_lowPt_8X_chargedHadronEnergy.Fill(c[1]['eta'], c[1]['j'].chargedHadronEnergy())
                eta_lowPt_8X_neutralHadronEnergy.Fill(c[1]['eta'], c[1]['j'].neutralHadronEnergy() - c[1]['j'].HFHadronEnergy())
                eta_lowPt_8X_photonEnergy.Fill(c[1]['eta'], c[1]['j'].photonEnergy())
                eta_lowPt_8X_electronEnergy.Fill(c[1]['eta'], c[1]['j'].electronEnergy())
                eta_lowPt_8X_muonEnergy.Fill(c[1]['eta'], c[1]['j'].muonEnergy())
                eta_lowPt_8X_HFHadronEnergy.Fill(c[1]['eta'], c[1]['j'].HFHadronEnergy())
                eta_lowPt_8X_HFEMEnergy.Fill(c[1]['eta'], c[1]['j'].HFEMEnergy())
                eta_lowPt_76X_chargedHadronEnergy.Fill(c[1]['eta'], c[0]['j'].chargedHadronEnergy())
                eta_lowPt_76X_neutralHadronEnergy.Fill(c[1]['eta'], c[0]['j'].neutralHadronEnergy() - c[0]['j'].HFHadronEnergy())
                eta_lowPt_76X_photonEnergy.Fill(c[1]['eta'], c[0]['j'].photonEnergy())
                eta_lowPt_76X_electronEnergy.Fill(c[1]['eta'], c[0]['j'].electronEnergy())
                eta_lowPt_76X_muonEnergy.Fill(c[1]['eta'], c[0]['j'].muonEnergy())
                eta_lowPt_76X_HFHadronEnergy.Fill(c[1]['eta'], c[0]['j'].HFHadronEnergy())
                eta_lowPt_76X_HFEMEnergy.Fill(c[1]['eta'], c[0]['j'].HFEMEnergy())

## Make plot
prefix="maxN_"+str(maxN)+'_'
ratio_jetResponse
jetResponsePlot = Plot.fromHisto(name = prefix+"pT_ratio_800pre6_76X", histos = [[ratio_jetResponse]], texX = "raw Jet p_{T} 76X", texY = "response ratio 800pre6/76X" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", logX = True, logY = False, sorting = False, yRange = (0.965, 1.045) )

energyFractionsPlot = Plot.fromHisto(name = prefix+"Efrac_ratio_800pre6_76X", histos = [[ratio_chargedHadronEnergy],[ratio_neutralEmEnergy],[ratio_neutralHadronEnergy] ], texX = "raw Jet p_{T} 76X", texY = "response ratio 800pre6/76X" )
plotting.draw(energyFractionsPlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", logX = True, logY = False, sorting = False, yRange = (0.965, 1.08) )

pt_8X = [pt_8X_HFEMEnergy, pt_8X_HFHadronEnergy, pt_8X_electronEnergy, pt_8X_muonEnergy, pt_8X_neutralHadronEnergy, pt_8X_photonEnergy, pt_8X_chargedHadronEnergy]
pt_8X_frac = []
for i, h in enumerate(pt_8X):
    h_new = h.ProjectionX().Clone()
    h_new.__dict__.update(h.__dict__) #Copy attributes
    pt_8X_frac.append( h_new )
    for j, h2 in enumerate(pt_8X[i+1:]):
        pt_8X_frac[-1].Add(h2.ProjectionX()) 

for h in reversed(pt_8X_frac):
    h.Divide(pt_8X_frac[0])

pt_76X = [pt_76X_HFEMEnergy, pt_76X_HFHadronEnergy, pt_76X_electronEnergy, pt_76X_muonEnergy, pt_76X_neutralHadronEnergy, pt_76X_photonEnergy, pt_76X_chargedHadronEnergy]
pt_76X_frac = []
for i, h in enumerate(pt_76X):
    h_new = h.ProjectionX().Clone()
    h_new.__dict__.update(h.__dict__) #Copy attributes
    pt_76X_frac.append( h_new )
    for j, h2 in enumerate(pt_76X[i+1:]):
        pt_76X_frac[-1].Add(h2.ProjectionX()) 

for h in reversed(pt_76X_frac):
    h.Divide(pt_76X_frac[0])

pt_plot = Plot.fromHisto(name = prefix+"pt", histos =  [[h] for h in pt_76X_frac+pt_8X_frac] , texX = "Jet p_{T} 76X", texY = "energy fractions" )  
plotting.draw(pt_plot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", logX = True, logY = False, sorting = False, yRange = (0,1), legend=(0.15,0.13,0.5,0.5))

eta_8X = [eta_8X_HFEMEnergy, eta_8X_HFHadronEnergy, eta_8X_electronEnergy, eta_8X_muonEnergy, eta_8X_neutralHadronEnergy, eta_8X_photonEnergy, eta_8X_chargedHadronEnergy]
eta_8X_frac = []
for i, h in enumerate(eta_8X):
    h_new = h.ProjectionX().Clone()
    h_new.__dict__.update(h.__dict__) #Copy attributes
    eta_8X_frac.append( h_new )
    for j, h2 in enumerate(eta_8X[i+1:]):
        eta_8X_frac[-1].Add(h2.ProjectionX()) 

for h in reversed(eta_8X_frac):
    h.Divide(eta_8X_frac[0])

eta_76X = [eta_76X_HFEMEnergy, eta_76X_HFHadronEnergy, eta_76X_electronEnergy, eta_76X_muonEnergy, eta_76X_neutralHadronEnergy, eta_76X_photonEnergy, eta_76X_chargedHadronEnergy]
eta_76X_frac = []
for i, h in enumerate(eta_76X):
    h_new = h.ProjectionX().Clone()
    h_new.__dict__.update(h.__dict__) #Copy attributes
    eta_76X_frac.append( h_new )
    for j, h2 in enumerate(eta_76X[i+1:]):
        eta_76X_frac[-1].Add(h2.ProjectionX()) 

for h in reversed(eta_76X_frac):
    h.Divide(eta_76X_frac[0])

eta_plot = Plot.fromHisto(name = prefix+"eta", histos =  [[h] for h in eta_76X_frac+eta_8X_frac] , texX = "raw Jet #eta 76X", texY = "energy fractions" )  
plotting.draw(eta_plot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", logX = False, logY = False, sorting = False, yRange = (0,1), legend=(0.40,0.13,0.7,0.4))

eta_lowPt_8X = [eta_lowPt_8X_HFEMEnergy, eta_lowPt_8X_HFHadronEnergy, eta_lowPt_8X_electronEnergy, eta_lowPt_8X_muonEnergy, eta_lowPt_8X_neutralHadronEnergy, eta_lowPt_8X_photonEnergy, eta_lowPt_8X_chargedHadronEnergy]
eta_lowPt_8X_frac = []
for i, h in enumerate(eta_lowPt_8X):
    h_new = h.ProjectionX().Clone()
    h_new.__dict__.update(h.__dict__) #Copy attributes
    eta_lowPt_8X_frac.append( h_new )
    for j, h2 in enumerate(eta_lowPt_8X[i+1:]):
        eta_lowPt_8X_frac[-1].Add(h2.ProjectionX()) 

for h in reversed(eta_lowPt_8X_frac):
    h.Divide(eta_lowPt_8X_frac[0])

eta_lowPt_76X = [eta_lowPt_76X_HFEMEnergy, eta_lowPt_76X_HFHadronEnergy, eta_lowPt_76X_electronEnergy, eta_lowPt_76X_muonEnergy, eta_lowPt_76X_neutralHadronEnergy, eta_lowPt_76X_photonEnergy, eta_lowPt_76X_chargedHadronEnergy]
eta_lowPt_76X_frac = []
for i, h in enumerate(eta_lowPt_76X):
    h_new = h.ProjectionX().Clone()
    h_new.__dict__.update(h.__dict__) #Copy attributes
    eta_lowPt_76X_frac.append( h_new )
    for j, h2 in enumerate(eta_lowPt_76X[i+1:]):
        eta_lowPt_76X_frac[-1].Add(h2.ProjectionX()) 

for h in reversed(eta_lowPt_76X_frac):
    h.Divide(eta_lowPt_76X_frac[0])

eta_lowPt_plot = Plot.fromHisto(name = prefix+"eta_lowPt_8X", histos =  [[h] for h in eta_lowPt_76X_frac+eta_lowPt_8X_frac] , texX = "raw Jet #eta 76X", texY = "energy fractions" )  
plotting.draw(eta_lowPt_plot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", logX = False, logY = False, sorting = False, yRange = (0,1), legend=(0.40,0.13,0.7,0.4))
