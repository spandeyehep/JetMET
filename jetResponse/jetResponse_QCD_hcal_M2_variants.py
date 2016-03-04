''' Reader example: Loop over a sample and write some data to a histogram.
'''
# Standard imports
import os
import logging
import ROOT
import array

#RootTools
from RootTools.core.standard import *

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

##data
#JetHT           = FWLiteSample.fromDAS("JetHT", "/JetHT/schoef-crab_JetHT_Run2015D_M2_5_500_lumiBased_reduced-8e13882dc7c4566a38618e8b59bae173/USER", instance="phys03", prefix = "root://hephyse.oeaw.ac.at/")
#JetHT_M2_5_500  = FWLiteSample.fromDAS("JetHT_M2_5_500", "/JetHT/schoef-crab_JetHT_Run2015D_lumiBased_reduced-4fff70efe810c67b5c65aa7d4a7cd41d/USER", instance="phys03", prefix = "root://hephyse.oeaw.ac.at/")
#JetHT_M0        = FWLiteSample.fromDAS("JetHT_M0", "/JetHT/schoef-crab_JetHT_Run2015D_lumiBased_reduced_M0-9d4a632a722feffaafd5b729e0393ea4/USER", instance="phys03", prefix = "root://hephyse.oeaw.ac.at/")
#JetHT_M21p      = FWLiteSample.fromDAS("JetHT_M21p", "/JetHT/schoef-crab_JetHT_Run2015D_lumiBased_reduced_M21p-4396cd947410146fe025d6c0e01e6549/USER", instance="phys03", prefix = "root://hephyse.oeaw.ac.at/")
#JetHT_M23p      = FWLiteSample.fromDAS("JetHT_M23p", "/JetHT/schoef-crab_JetHT_Run2015D_lumiBased_reduced_M23p-dee0be6699db8222b92d05d38956587e/USER", instance="phys03", prefix = "root://hephyse.oeaw.ac.at/")

##sim
QCD_Pt_15to3000 = Sample.fromCMGOutput("QCD_Pt_15to3000", baseDirectory = "/data/rschoefbeck/cmgTuples/QCD_HCal/QCD_Pt_15to3000", treeFilename='tree.root')
QCD_Pt_15to3000_M2_5_500 = Sample.fromCMGOutput("QCD_Pt_15to3000_M2_5_500", baseDirectory = "/data/rschoefbeck/cmgTuples/QCD_HCal/QCD_Pt_15to3000_M2_5_500", treeFilename='tree.root')
QCD_Pt_15to3000_M0 = Sample.fromCMGOutput("QCD_Pt_15to3000_M0", baseDirectory = "/data/rschoefbeck/cmgTuples/QCD_HCal/", chunkString="QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8_M0", treeFilename='tree.root')
QCD_Pt_15to3000_M21p = Sample.fromCMGOutput("QCD_Pt_15to3000_M21p", baseDirectory = "/data/rschoefbeck/cmgTuples/QCD_HCal/", chunkString="QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8_M21p", treeFilename='tree.root')
QCD_Pt_15to3000_M23p = Sample.fromCMGOutput("QCD_Pt_15to3000_M23p", baseDirectory = "/data/rschoefbeck/cmgTuples/QCD_HCal/", chunkString="QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8_M23p", treeFilename='tree.root')

# define TProfiles
thresholds = [10**(x/10.) for x in range(11,36)] 
jetResponse_M2_0_100 = ROOT.TProfile("response", "response", len(thresholds)-1, array.array('d', thresholds) )
jetResponse_M2_0_100.texName = "760 M2(0,100)"
jetResponse_M2_0_100.style = styles.lineStyle(ROOT.kBlack)
jetResponse_M2_5_500 = ROOT.TProfile("response_M2_5_500", "response_M2_5_500", len(thresholds)-1, array.array('d', thresholds) )
jetResponse_M2_5_500.texName = "M2(5,500)"
jetResponse_M2_5_500.style = styles.lineStyle(ROOT.kRed)
jetResponse_M0 = ROOT.TProfile("response_M0", "response_M0", len(thresholds)-1, array.array('d', thresholds) )
jetResponse_M0.texName = "M0"
jetResponse_M0.style = styles.lineStyle(ROOT.kBlue)
jetResponse_M21p = ROOT.TProfile("response_M21p", "response_M21p", len(thresholds)-1, array.array('d', thresholds) )
jetResponse_M21p.texName = "1-pulse"
jetResponse_M21p.style = styles.lineStyle(ROOT.kGreen) 
jetResponse_M23p = ROOT.TProfile("response_M23p", "response_M23p", len(thresholds)-1, array.array('d', thresholds) )
jetResponse_M23p.texName = "3-pulse"
jetResponse_M23p.style = styles.lineStyle(ROOT.kMagenta)

maxN = 100000

variables =  map( Variable.fromString, ['evt/l','run/I', 'lumi/I', 'Jet[pt/F,rawPt/F,phi/F,eta/F]', 'nJet/I'] )

for response, sample1, sample2 in [\
    [jetResponse_M2_0_100, QCD_Pt_15to3000, QCD_Pt_15to3000],
    [jetResponse_M2_5_500, QCD_Pt_15to3000_M2_5_500, QCD_Pt_15to3000],
    [jetResponse_M0, QCD_Pt_15to3000_M0, QCD_Pt_15to3000],
    [jetResponse_M21p, QCD_Pt_15to3000_M21p, QCD_Pt_15to3000],
    [jetResponse_M23p, QCD_Pt_15to3000_M23p, QCD_Pt_15to3000],
    ]:

    # Read all run/lumi/event for first tree
    r1 = sample1.treeReader( variables = variables , selectionString = "1")
    r1.start()
    position_r1 = {}
    while r1.run():
        position_r1[(r1.data.run, r1.data.lumi, r1.data.evt)] = r1.eList.GetEntry(r1.position-1)
        if maxN>0 and r1.position>=maxN:break

    # Read all run/lumi/event for second tree
    r2 = sample2.treeReader( variables = variables , selectionString = "1" )
    r2.start()
    position_r2={}
    while r2.run():
        position_r2[(r2.data.run, r2.data.lumi, r2.data.evt)] = r2.eList.GetEntry(r2.position-1)
        if maxN>0 and r2.position>=maxN:break

    logger.info( "Have %i events in first samle and %i in second", len(position_r1), len(position_r2) )

    # Fast intersect
    intersec = set(position_r1.keys()).intersection(set(position_r2.keys()))
    positions = [(position_r1[i], position_r2[i]) for i in intersec]
    positions.sort()
    logger.info("Have %i events in common.", len(intersec))

    #Loop over pairs
    for nEvent, pos in enumerate(positions):
        if nEvent%10000==0: 
            logger.info("At event %i/%i", nEvent, len(intersec))

        p1, p2 = pos
        r1.goToPosition(p1)
        r2.goToPosition(p2)

        n = min([r1.data.nJet, r2.data.nJet])
        jets1 = [{'rawPt':r1.data.Jet_rawPt[i], 'eta':r1.data.Jet_eta[i], 'phi':r1.data.Jet_phi[i]} for i in xrange(n) if abs(r1.data.Jet_eta[i])<1.3]
        jets2 = [{'rawPt':r2.data.Jet_rawPt[i], 'eta':r2.data.Jet_eta[i], 'phi':r2.data.Jet_phi[i]} for i in xrange(n) if abs(r2.data.Jet_eta[i])<1.3]
        for c in zip(jets1, jets2):
            if deltaR2(*c)<0.2**2:
                response.Fill(c[1]['rawPt'],c[0]['rawPt']/c[1]['rawPt'])


# Make plot
profiles = [jetResponse_M2_0_100, jetResponse_M2_5_500, jetResponse_M0, jetResponse_M21p, jetResponse_M23p]
prefix="maxN_%s_"%maxN if maxN is not None and maxN>0 else ""
histos = [ [h.ProjectionX()] for h in profiles ]
for i, h in enumerate(histos):
    h[0].__dict__.update(profiles[i].__dict__)
#    h[0].Divide(jetResponse_M2_0_100.ProjectionX())

jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponse_QC_Deta13", histos = histos, texX = "raw Jet p_{T}", texY = "response ratio wrt. 760" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = True, yRange=(0.95,1.07))
