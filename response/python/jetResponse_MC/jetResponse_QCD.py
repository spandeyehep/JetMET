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

max_events = 4000000
max_files = 30 

#new       = FWLiteSample.fromDAS("new"      , "/RelValTTbar_13/CMSSW_8_0_21-PU25ns_80X_mcRun2_asymptotic_2016_TrancheIV_v6_Tr4GT_v6-v1/MINIAODSIM")
#old       = FWLiteSample.fromDAS("old"      , "/RelValTTbar_13/CMSSW_8_0_20-PU25ns_80X_mcRun2_asymptotic_2016_TrancheIV_v2_Tr4GT_v2-v1/MINIAODSIM")
new       = FWLiteSample.fromDAS("new"      , "/QCD_Pt-15to7000_TuneCUETP8M1_FlatP6_13TeV_pythia8/RunIISummer16MiniAODv2-NoPU_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", maxN = max_files)
old       = FWLiteSample.fromDAS("old"      , "/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISummer16MiniAODv2-NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v4-v2/MINIAODSIM", maxN = max_files)

preprefix = "QCD_v2_vs_v6"

useGenPt = True

# define TProfiles
thresholds = [10**(x/10.) for x in range(11,36)] 

eta_thresholds = [0, 1.3, 2.5, 3]
color={0:ROOT.kBlack, 1.3:ROOT.kBlue, 2.5:ROOT.kRed, 3:ROOT.kMagenta}
resp = {}
for s in [old.name, new.name]:
    resp[s] = {} 
    for i_eta_th, eta_th in enumerate(eta_thresholds):
        resp[s][eta_th] = ROOT.TProfile("response", "response", len(thresholds)-1, array.array('d', thresholds) )
        resp[s][eta_th].style = styles.lineStyle(color[eta_th], dashed = (s==old.name) )
        if s==new.name:
            resp[s][eta_th].legendText = "%2.1f<=#eta"%eta_th
            if eta_th!=eta_thresholds[-1]: resp[s][eta_th].legendText += "<%2.1f"%eta_thresholds[i_eta_th+1]
            resp[s][eta_th].legendText += " (old)" if s==old.name else "" 

products = {
    'jets':      {'type':'vector<pat::Jet>', 'label':("slimmedJets")},
    'genInfo':   {'type':' GenEventInfoProduct', 'label': "generator"},
    }

for sample in [new, old]:
    r1 = sample.fwliteReader( products = products )
    r1.start()
    i=0
    while r1.run():
        #pt_hat = r1.products['genInfo'].binningValues()[0]
        #if i%10000==0: logger.info("At %i", i)
        #if not (pt_hat > pt_hat_min and pt_hat <pt_hat_max): continue
            
        id_jets = [ j.correctedJet("Uncorrected") for j in r1.products['jets'] if helpers.jetID( j )]
        for j in id_jets:
            gj = j.genJet()
            if gj:
                for eta_th in reversed(eta_thresholds):
                    if abs(gj.eta())>eta_th:
                        resp[sample.name][eta_th].Fill( gj.pt(), j.pt() / gj.pt() )
                        #print eta_th, gj.eta(), j.pt(), gj.pt()
                        break

        i+=1
        if i>max_events: break

        
# Make plot
profiles = [resp["old"][t] for t in eta_thresholds] + [resp["new"][t] for t in eta_thresholds]
#profiles = [ jetResponse_NJC]
prefix=preprefix+("_max_events_%s_"%max_events if max_events is not None and max_events>0 else "_")
histos = [ [h.ProjectionX()] for h in profiles ]
for i, h in enumerate(histos):
    h[0].__dict__.update(profiles[i].__dict__)

jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponse_QCD", histos = histos, texX = "gen Jet p_{T}" , texY = "relative response" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = True, yRange=(0.7,1.2),  legend = [0.50,0.92-0.05*4,0.92,0.88])

# Ratios
histos = [ [h.ProjectionX()] for h in profiles ]
for i, h in enumerate(histos):
    h[0].__dict__.update(profiles[i].__dict__)

for i in range(len(eta_thresholds)):
    histos[i+4][0].Divide( histos[i][0] )

histos = histos[4:]

jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponse_QCD_ratio", histos = histos, texX = "gen Jet p_{T}" , texY = "relative response: v6/v2" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = True, yRange=(0.9,1.1),  legend = [0.50,0.92-0.05*4,0.92,0.88])
