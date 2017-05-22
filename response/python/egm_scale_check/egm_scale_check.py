''' FWLiteReader example: Loop over a sample and write some data to a histogram.
'''
# Standard imports
import os
import logging
import ROOT
import array

#RootTools
from RootTools.core.standard import *

#Helper
import JetMET.tools.helpers as helpers
from math import sqrt

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

from FWCore.PythonUtilities.LumiList import LumiList
# Apply golden JSON
json = '$CMSSW_BASE/src/CMGTools/TTHAnalysis/data/json/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' 
lumiList = LumiList(os.path.expandvars(json))
logger.info( "Loaded json %s", json )

max_events = 1000 
max_files = 1
##data

Run2016B = FWLiteSample.fromDAS("Run2016B", "/SingleElectron/Run2016B-03Feb2017_ver1-v1/MINIAOD", maxN = 2*max_files)
#Run2016C = FWLiteSample.fromDAS("Run2016C", "/SingleElectron/Run2016C-03Feb2017-v1/MINIAOD", maxN = 2*max_files)
Run2016D = FWLiteSample.fromDAS("Run2016D", "/SingleElectron/Run2016D-03Feb2017-v1/MINIAOD", maxN = 2*max_files)
Run2016E = FWLiteSample.fromDAS("Run2016E", "/SingleElectron/Run2016E-03Feb2017-v1/MINIAOD", maxN = 2*max_files)
#Run2016F = FWLiteSample.fromDAS("Run2016F", "/SingleElectron/Run2016F-03Feb2017-v1/MINIAOD", maxN = 2*max_files)
Run2016G = FWLiteSample.fromDAS("Run2016G", "/SingleElectron/Run2016G-03Feb2017-v1/MINIAOD", maxN = 2*max_files)
Run2016H = FWLiteSample.fromDAS("Run2016H", "/SingleElectron/Run2016H-03Feb2017_ver3-v1/MINIAOD", maxN = 2*max_files)

preprefix="gOverE_"

data_samples = [Run2016B, Run2016D, Run2016E, Run2016G, Run2016H]
Run2016B.color=ROOT.kBlack
Run2016D.color=ROOT.kRed
Run2016E.color=ROOT.kGreen
Run2016G.color=ROOT.kMagenta
Run2016H.color=ROOT.kOrange
for s in data_samples:
    s.isData = True

DY_Pt250to400 = FWLiteSample.fromDAS("DY_Pt250to400", "/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM", maxN = max_files)
DY_Pt250to400.isData = False
DY_Pt250to400.color = ROOT.kBlue
DY_Pt400to650 = FWLiteSample.fromDAS("DY_Pt400to650", "/DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM", maxN = max_files)
DY_Pt400to650.isData = False
DY_Pt400to650.color = ROOT.kGreen
DY_Pt650toInf = FWLiteSample.fromDAS("DY_Pt650toInf", "/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM", maxN = max_files)
DY_Pt650toInf.isData = False
DY_Pt650toInf.color = ROOT.kGreen

DY = FWLiteSample.combine("DY_PtAbove250", [DY_Pt250to400, DY_Pt400to650, DY_Pt650toInf])
DY.isData = False
DY.color  = ROOT.kBlue

#/DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM
#/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM

#mc_/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia_20M/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM

samples = [ DY ] + data_samples


# binning 
thresholds = [10**(x/10.) for x in range(11,36)] 
thresholds_r  = [ 0.95 + i*(1.07-0.95)/30. for i in range(30) ] 

gOverE = {}
gOverE_2D = {}
for sample in samples:
    gOverE[sample.name]          = ROOT.TProfile("gOverE", "gOverE", len(thresholds)-1, array.array('d', thresholds) )
    gOverE[sample.name].texName  = "#gamma/e ( "+sample.name+" )"
    gOverE[sample.name].style    = styles.lineStyle(sample.color)
    #gOverE['BeforeGSFix']           = ROOT.TProfile("gOverEBeforeGSFix", "gOverEBeforeGSFix", len(thresholds)-1, array.array('d', thresholds) )
    #gOverE['BeforeGSFix'].texName   = "#gamma/e (before GS Fix )"
    #gOverE['BeforeGSFix'].style     = styles.lineStyle(ROOT.kBlue)

    gOverE_2D[sample.name]        = ROOT.TH2D("gOverE_2D", "gOverE_2D", len(thresholds)-1, array.array('d', thresholds), len(thresholds_r)-1, array.array('d', thresholds_r) )
    #gOverE_2D['BeforeGSFix']= ROOT.TH2D("goverEBeforeGSFix_2D", "goverEBeforeGSFix_2D", len(thresholds)-1, array.array('d', thresholds), len(thresholds_r)-1, array.array('d', thresholds_r) )

    products = {
        'electrons':                 {'type':'vector<pat::Electron>', 'label':("slimmedElectrons")},
        'photons':                   {'type':'vector<pat::Photon>', 'label':("slimmedPhotons")},
    #    'electronsBeforeGSFix':      {'type':'vector<pat::Electron>', 'label':("slimmedElectronsBeforeGSFix")},
    #    'photonsBeforeGSFix':        {'type':'vector<pat::Photon>', 'label':("slimmedPhotonsBeforeGSFix")},
        }

    r = sample.fwliteReader( products = products )

    r.start()
    runs = set()
    position_r = {}
    count=0

    while r.run( ):
        if max_events is not None and max_events>0 and count>=max_events:break
        if sample.isData and ( not lumiList.contains(r.evt[0], r.evt[1])) : continue
            
        electrons = [ e for e in r.products['electrons'] if e.energy()>50 ]
        photons   = [ p for p in r.products['photons'] if p.energy()>100 ]
        for e in electrons:

            if abs(e.eta())<1.3 : continue
            if not e.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight"): continue 

            e_sc = {'eta':e.superCluster().eta(), 'phi':e.superCluster().phi()}

            for p in photons:
                p_sc = {'eta':p.superCluster().eta(), 'phi':p.superCluster().phi()}
                if helpers.deltaR2(e_sc, p_sc) == 0.:
                    gOverE[sample.name].Fill( p.pt(), p.pt()/e.pt() )
                    gOverE_2D[sample.name].Fill( p.pt(), p.pt()/e.pt() )

        count += 1

# Make plot
#profiles = [jetResponse_M2_0_100, jetResponse_M2_5_500, jetResponse_M2_0_500, jetResponse_M2_5_100, jetResponse_M0, jetResponse_M21p, jetResponse_M23p]
profiles = [ gOverE[sample.name] for sample in samples ]
prefix   = preprefix+"max_events_%s_"%max_events if max_events is not None and max_events>0 else preprefix
histos = [ [h.ProjectionX()] for h in profiles ]
for i, h in enumerate(histos):
    h[0].__dict__.update(profiles[i].__dict__)
#    h[0].Divide(jetResponse_M2_0_100.ProjectionX())

jetResponsePlot = Plot.fromHisto(name = prefix+"gOverE", histos = histos, texX = "photon p_{T}", texY = "p_{T}(#gamma)/p_{T}(e)" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = {'num':0,'den':1, 'yRange':(0.98,1.02)}, logY = False, logX = True, yRange=(0.95,1.07))

for sample in samples:
    jetResponsePlot2D = Plot2D.fromHisto(name = prefix+"gOverE_2D_"+sample.name, histos = [[gOverE_2D[sample.name]]], texX = "photon p_{T}", texY = "p_{T}(#gamma)/p_{T}(e)" )
    plotting.draw2D(jetResponsePlot2D, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", logY = False, logX = True)
