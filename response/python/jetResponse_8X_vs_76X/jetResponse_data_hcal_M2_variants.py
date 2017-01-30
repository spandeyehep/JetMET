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
import helpers

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

##data
JetHT           = FWLiteSample.fromDAS("JetHT", "/JetHT/schoef-crab_JetHT_Run2015D_lumiBased_reduced-4fff70efe810c67b5c65aa7d4a7cd41d/USER", instance="phys03", prefix = "root://hephyse.oeaw.ac.at/")
#JetHT_M2_5_500  = FWLiteSample.fromDAS("JetHT_M2_5_500", "/JetHT/schoef-crab_JetHT_Run2015D_M2_5_500_lumiBased_reduced-8e13882dc7c4566a38618e8b59bae173/USER", instance="phys03", prefix = "root://hephyse.oeaw.ac.at/")
#JetHT_M2_0_500  = FWLiteSample.fromDAS("JetHT_M2_0_500", "/JetHT/schoef-crab_JetHT_Run2015D_lumiBased_reduced_M2_0_500-1698d5c0a2c4014a47bae89464ec1363/USER", instance="phys03", prefix = "root://hephyse.oeaw.ac.at/")
#JetHT_M2_5_100  = FWLiteSample.fromDAS("JetHT_M2_5_100", "/JetHT/schoef-crab_JetHT_Run2015D_lumiBased_reduced_M2_5_100-1844535d82563cf58dccb520d79351d2/USER", instance="phys03", prefix = "root://hephyse.oeaw.ac.at/")
#JetHT_M0        = FWLiteSample.fromDAS("JetHT_M0", "/JetHT/schoef-crab_JetHT_Run2015D_lumiBased_reduced_M0-9d4a632a722feffaafd5b729e0393ea4/USER", instance="phys03", prefix = "root://hephyse.oeaw.ac.at/")
#JetHT_M21p      = FWLiteSample.fromDAS("JetHT_M21p", "/JetHT/schoef-crab_JetHT_Run2015D_lumiBased_reduced_M21p-4396cd947410146fe025d6c0e01e6549/USER", instance="phys03", prefix = "root://hephyse.oeaw.ac.at/")
#JetHT_M23p      = FWLiteSample.fromDAS("JetHT_M23p", "/JetHT/schoef-crab_JetHT_Run2015D_lumiBased_reduced_M23p-dee0be6699db8222b92d05d38956587e/USER", instance="phys03", prefix = "root://hephyse.oeaw.ac.at/")
JetHT_NJC      = FWLiteSample.fromDAS("JetHT_NJC", "/JetHT/schoef-crab_JetHT_Run2015D_lumiBased_reduced_noJetCoreTracking-0db42af97b9a6a6b284e4d8f7e0fb0f2/USER", instance="phys03", prefix = "root://hephyse.oeaw.ac.at/")

preprefix="NJC_"

# define TProfiles
thresholds = [10**(x/10.) for x in range(11,36)] 

#jetResponse_M2_0_100 = ROOT.TProfile("response", "response", len(thresholds)-1, array.array('d', thresholds) )
#jetResponse_M2_0_100.texName = "760 M2(0,100)"
#jetResponse_M2_0_100.style = styles.lineStyle(ROOT.kBlack)
#jetResponse_M2_5_500 = ROOT.TProfile("response_M2_5_500", "response_M2_5_500", len(thresholds)-1, array.array('d', thresholds) )
#jetResponse_M2_5_500.texName = "M2(5,500)"
#jetResponse_M2_5_500.style = styles.lineStyle(ROOT.kOrange)
#jetResponse_M2_0_500 = ROOT.TProfile("response_M2_0_500", "response_M2_0_500", len(thresholds)-1, array.array('d', thresholds) )
#jetResponse_M2_0_500.texName = "M2(0,500)"
#jetResponse_M2_0_500.style = styles.lineStyle(ROOT.kRed - 7)
#jetResponse_M2_5_100 = ROOT.TProfile("response_M2_5_100", "response_M2_5_100", len(thresholds)-1, array.array('d', thresholds) )
#jetResponse_M2_5_100.texName = "M2(5,100)"
#jetResponse_M2_5_100.style = styles.lineStyle(ROOT.kRed)
#jetResponse_M0 = ROOT.TProfile("response_M0", "response_M0", len(thresholds)-1, array.array('d', thresholds) )
#jetResponse_M0.texName = "M0"
#jetResponse_M0.style = styles.lineStyle(ROOT.kBlue)
#jetResponse_M21p = ROOT.TProfile("response_M21p", "response_M21p", len(thresholds)-1, array.array('d', thresholds) )
#jetResponse_M21p.texName = "1-pulse"
#jetResponse_M21p.style = styles.lineStyle(ROOT.kGreen) 
#jetResponse_M23p = ROOT.TProfile("response_M23p", "response_M23p", len(thresholds)-1, array.array('d', thresholds) )
#jetResponse_M23p.texName = "3-pulse"
#jetResponse_M23p.style = styles.lineStyle(ROOT.kMagenta)
jetResponse_NJC = ROOT.TProfile("response_NJC", "response_NJC", len(thresholds)-1, array.array('d', thresholds) )
jetResponse_NJC.texName = "no Jet core tracking"
jetResponse_NJC.style = styles.lineStyle(ROOT.kOrange)

products = {
    'jets':      {'type':'vector<pat::Jet>', 'label':("slimmedJets")},
#    'pfClusters':{'type':"vector<reco::PFCluster>", 'label':("particleFlowClusterHBHE")}, 
#    'pfRecHits': {'type':"vector<reco::PFRecHit>", 'label':("particleFlowRecHitHBHE")},
    }

maxN = 100000

for response, sample1, sample2 in [\
#    [jetResponse_M2_0_100, JetHT, JetHT],
#    [jetResponse_M2_5_500, JetHT_M2_5_500, JetHT],
#    [jetResponse_M2_0_500, JetHT_M2_0_500, JetHT],
#    [jetResponse_M2_5_100, JetHT_M2_5_100, JetHT],
#    [jetResponse_M0, JetHT_M0, JetHT],
#    [jetResponse_M21p, JetHT_M21p, JetHT],
#    [jetResponse_M23p, JetHT_M23p, JetHT],
    [jetResponse_NJC, JetHT_NJC, JetHT],
    ]:
    r1 = sample1.fwliteReader( products = products )
    r2 = sample2.fwliteReader( products = products )

    r1.start()
    runs_1 = set()
    position_r1 = {}
    count=0
    while r1.run( readProducts = False ):
            if r1.evt[0]==260431: 
                position_r1[r1.evt] = r1.position-1
                count+=1
            if maxN is not None and maxN>0 and count>=maxN:break

    r2.start()
    runs_2 = set()
    position_r2 = {}
    count=0
    while r2.run( readProducts = False ):
            if r2.evt[0]==260431: 
                position_r2[r2.evt] = r2.position-1
                count+=1
            if maxN is not None and maxN>0 and count>=maxN:break

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
            if helpers.deltaR2(*c)<0.2**2:
                if not ( helpers.jetID(c[0]['j']) and helpers.jetID(c[1]['j']) ): continue

                if abs(c[0]['eta'])<1.3 and abs(c[1]['eta'])<1.3:
                    response.Fill(c[1]['pt'], c[0]['pt']/c[1]['pt'] )


# Make plot
#profiles = [jetResponse_M2_0_100, jetResponse_M2_5_500, jetResponse_M2_0_500, jetResponse_M2_5_100, jetResponse_M0, jetResponse_M21p, jetResponse_M23p]
profiles = [ jetResponse_NJC ]
prefix=preprefix+"maxN_%s_"%maxN if maxN is not None and maxN>0 else ""
histos = [ [h.ProjectionX()] for h in profiles ]
for i, h in enumerate(histos):
    h[0].__dict__.update(profiles[i].__dict__)
#    h[0].Divide(jetResponse_M2_0_100.ProjectionX())

jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponse_data_eta13", histos = histos, texX = "raw Jet p_{T}", texY = "response ratio wrt. 760" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = True, yRange=(0.95,1.07))
