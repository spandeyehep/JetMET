''' FWLiteReader based closure plots for L1L2L3.
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

max_events = 1
max_files = 1

spring16   = FWLiteSample.fromDAS("spring16"      , "/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM", maxN = max_files)
moriond17  = FWLiteSample.fromDAS("moriond17"     , "/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", maxN = max_files)

#pt_hat_min = 300
#pt_hat_max = 350
#preprefix = "relval_ttbar_pt_hat_%i_%i" % ( pt_hat_min, pt_hat_max)
preprefix = "relval_QCD_HT_100to200"

# define TProfiles
thresholds = [10**(x/10.) for x in range(11,36)] 

eta_thresholds = [0, 1.3, 2.5, 3]
color={0:ROOT.kBlack, 1.3:ROOT.kBlue, 2.5:ROOT.kRed, 3:ROOT.kMagenta}
resp = {}
for s in [spring16.name, moriond17.name]:
    resp[s] = {} 
    for i_eta_th, eta_th in enumerate(eta_thresholds):
        resp[s][eta_th] = ROOT.TProfile("response", "response", len(thresholds)-1, array.array('d', thresholds) )
        resp[s][eta_th].style = styles.lineStyle(color[eta_th], dashed = (s==spring16.name) )
        if s==moriond17.name:
            resp[s][eta_th].legendText = "%2.1f<=#eta"%eta_th
            if eta_th!=eta_thresholds[-1]: resp[s][eta_th].legendText += "<%2.1f"%eta_thresholds[i_eta_th+1]
            resp[s][eta_th].legendText += " (Spring16)" if s==spring16.name else " (Moriond17)" 

products = {
    'jets':      {'type':'vector<pat::Jet>', 'label':("slimmedJets")},
#    'genInfo':   {'type':' GenEventInfoProduct', 'label': "generator"},

#    'pfClusters':{'type':"vector<reco::PFCluster>", 'label':("particleFlowClusterHBHE")}, 
#    'pfRecHits': {'type':"vector<reco::PFRecHit>", 'label':("particleFlowRecHitHBHE")},
    }

for sample in [spring16, moriond17]:
    r1 = sample.fwliteReader( products = products )
    r1.start()
    i=0
    while r1.run():
#        pt_hat = r1.products['genInfo'].binningValues()[0]
#        if i%10000==0: logger.info("At %i", i)
#        if not (pt_hat > pt_hat_min and pt_hat <pt_hat_max): continue
            
        #id_jets = [ j.correctedJet("Uncorrected") for j in r1.products['jets'] if helpers.jetID( j )]
        id_jets = [ j for j in r1.products['jets'] if helpers.jetID( j )]
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

        

## Make plot
#profiles = [resp["spring16"][t] for t in eta_thresholds] + [resp["moriond17"][t] for t in eta_thresholds]
##profiles = [ jetResponse_NJC]
#prefix=preprefix+"_max_events_%s_"%max_events if max_events is not None and max_events>0 else ""
#histos = [ [h.ProjectionX()] for h in profiles ]
#for i, h in enumerate(histos):
#    h[0].__dict__.update(profiles[i].__dict__)
#
#jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponse_QCD", histos = histos, texX = "gen Jet p_{T}" , texY = "corrected response" )
#plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = True, yRange=(0.9,1.2))
