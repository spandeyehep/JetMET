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

max_events = 500000
max_files = 30

#'/scratch/rschoefbeck/flat_HCal_M2/flat_pi500_2017_realistic_PF.root'
#'/scratch/rschoefbeck/flat_HCal_M2/flat_pi500_2017_realistic_PF_PhotoStatUnc.root'


#preprefix = "cluster_pt100_pi-20-30_10k_2017"
preprefix = "cluster_pt100_pi-pi500_10k_2017"
#preprefix = "cluster_pt100_pi-90-100_10k_2017"

files_def = [ 
#'/scratch/rschoefbeck/flat_HCal_M2/flat_pi-20-30_10k_2017_realistic_PF.root'
'/scratch/rschoefbeck/flat_HCal_M2/flat_pi500_10k_2017_realistic_PF.root'
#'/scratch/rschoefbeck/flat_HCal_M2/flat_pi-90-100_10k_2017_realistic_PF.root'
]
files_fix = [ 
#'/scratch/rschoefbeck/flat_HCal_M2/flat_pi-20-30_10k_2017_realistic_PF_PhotoStatUnc.root'
'/scratch/rschoefbeck/flat_HCal_M2/flat_pi500_10k_2017_realistic_PF_PhotoStatUnc.root'
#'/scratch/rschoefbeck/flat_HCal_M2/flat_pi-90-100_10k_2017_realistic_PF_PhotoStatUnc.root'
]

s_def   = Sample.fromFiles("def"      , files = files_def, maxN = max_files, treeName = "s")
s_fix   = Sample.fromFiles("fix"      , files = files_fix, maxN = max_files, treeName = "s")

for name, cut in [ ["HPD", "abs(eta)<1.3"], ["SiPM", "abs(eta)>1.3"]]:
    s_fix.setSelectionString( "ecal==0&&true>100&&"+cut )
    s_def.setSelectionString( "ecal==0&&true>100&&"+cut )

    prefix=preprefix+"_"+name+"_"+s_def.name

    h_def = s_def.get1DHistoFromDraw("hcal/true", [50,0,2])
    h_fix = s_fix.get1DHistoFromDraw("hcal/true", [50,0,2])

    ## define TProfiles
    #thresholds = [10**(x/10.) for x in range(11,36)] 
    #
    #eta_thresholds = [0, 1.3, 2.5, 3]
    #color={0:ROOT.kBlack, 1.3:ROOT.kBlue, 2.5:ROOT.kRed, 3:ROOT.kMagenta}
    #resp = {}
    #for s in [old.name, new.name]:
    #    resp[s] = {} 
    #    for i_eta_th, eta_th in enumerate(eta_thresholds):
    #        resp[s][eta_th] = ROOT.TProfile("response", "response", len(thresholds)-1, array.array('d', thresholds) )
    #        resp[s][eta_th].style = styles.lineStyle(color[eta_th], dashed = (s==old.name) )
    #        if s==new.name:
    #            resp[s][eta_th].legendText = "%2.1f<=#eta"%eta_th
    #            if eta_th!=eta_thresholds[-1]: resp[s][eta_th].legendText += "<%2.1f"%eta_thresholds[i_eta_th+1]
    #            resp[s][eta_th].legendText += " (old)" if s==old.name else " (new)" 
    #
    #products = {
    #    'jets':      {'type':'vector<pat::Jet>', 'label':("slimmedJets")},
    #    'genInfo':   {'type':' GenEventInfoProduct', 'label': "generator"},
    #
    ##    'pfClusters':{'type':"vector<reco::PFCluster>", 'label':("particleFlowClusterHBHE")}, 
    ##    'pfRecHits': {'type':"vector<reco::PFRecHit>", 'label':("particleFlowRecHitHBHE")},
    #    }

    #for sample in [new, old]:
    #    r1 = sample.fwliteReader( products = products )
    #    r1.start()
    #    i=0
    #    while r1.run():
    #        pt_hat = r1.products['genInfo'].binningValues()[0]
    #        if i%10000==0: logger.info("At %i", i)
    #        if not (pt_hat > pt_hat_min and pt_hat <pt_hat_max): continue
    #            
    #        id_jets = [ j.correctedJet("Uncorrected") for j in r1.products['jets'] if helpers.jetID( j )]
    #        for j in id_jets:
    #            gj = j.genJet()
    #            if gj:
    #                for eta_th in reversed(eta_thresholds):
    #                    if abs(gj.eta())>eta_th:
    #                        resp[sample.name][eta_th].Fill( gj.pt(), j.pt() / gj.pt() )
    #                        #print eta_th, gj.eta(), j.pt(), gj.pt()
    #                        break
    #
    #        i+=1
    #        if i>max_events: break
    #
    #        
    #
    ## Make plot
    #profiles = [resp["old"][t] for t in eta_thresholds] + [resp["new"][t] for t in eta_thresholds]
    ##profiles = [ jetResponse_NJC]
    #histos = [ [h.ProjectionX()] for h in profiles ]
    #for i, h in enumerate(histos):
    #    h[0].__dict__.update(profiles[i].__dict__)

    h_def.legendText = "def (RMS %3.2f #pm %3.2f)"% (h_def.GetRMS()/h_def.GetMean(), h_def.GetRMSError()/h_def.GetMean())
    h_def.style = styles.lineStyle( ROOT.kBlue )
    h_fix.legendText = "fix (RMS %3.2f #pm %3.2f)"% (h_fix.GetRMS()/h_fix.GetMean(), h_fix.GetRMSError()/h_fix.GetMean())
    h_fix.style = styles.lineStyle( ROOT.kRed )
    clusterResponsPlot = Plot.fromHisto(name = prefix+"cluster_test", histos = [ [ h_def], [ h_fix] ], texX = "H/T" , texY = "relative response" )
    plotting.draw(clusterResponsPlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = False ) #, yRange=(0.7,1.2))
