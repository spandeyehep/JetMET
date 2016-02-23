#Standard imports
import sys
import logging
import ROOT
import array
from math import log

#RootTools
from RootTools.core.Sample import Sample
from RootTools.core.Variable import Variable, ScalarType, VectorType
from RootTools.core.TreeReader import TreeReader
from RootTools.core.logger import get_logger
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

#make samples. Samples are statisticall compatible.
maxN = -1
jetHT_M2_5_500  = Sample.fromCMGOutput("JetHT_M2_5_500", baseDirectory = "/data/rschoefbeck/cmgTuples/JetHT_HCal/JetHT_schoef-crab_JetHT_Run2015D_M2_5_500/", treeFilename='tree.root', maxN=maxN)
jetHT           = Sample.fromCMGOutput("JetHT", baseDirectory = "/data/rschoefbeck/cmgTuples/JetHT_HCal/JetHT_schoef-crab_JetHT_Run2015D/", treeFilename='tree.root', maxN=maxN)
maxN = -1
QCD_Pt_15to3000_M2_5_500 = Sample.fromCMGOutput("QCD_Pt_15to3000_M2_5_500", baseDirectory = "/data/rschoefbeck/cmgTuples/QCD_HCal/QCD_Pt_15to3000_M2_5_500", treeFilename='tree.root', maxN=maxN)
QCD_Pt_15to3000 = Sample.fromCMGOutput("QCD_Pt_15to3000", baseDirectory = "/data/rschoefbeck/cmgTuples/QCD_HCal/QCD_Pt_15to3000", treeFilename='tree.root', maxN=maxN)

variables =  map( Variable.fromString, ['evt/l','run/I', 'lumi/I', 'Jet[pt/F,rawPt/F,phi/F,eta/F]', 'nJet/I'] )

# define TProfiles
thresholds = [10**(x/10.) for x in range(11,36)]
jetResponse_data = ROOT.TProfile("response_data", "response_data", len(thresholds)-1, array.array('d', thresholds) )
jetResponse_data.texName = "JetHT 260431" 
jetResponse_data.style = lineStyle(ROOT.kBlue)
jetResponse_mc = ROOT.TProfile("response_mc", "response_mc", len(thresholds)-1, array.array('d', thresholds) )
jetResponse_mc.texName = "QCD flat"
jetResponse_mc.style = lineStyle(ROOT.kRed)

# Loop over data and MC
for jetResponse, sample1, sample2 in [\
    [jetResponse_data, jetHT_M2_5_500, jetHT],
    [jetResponse_mc, QCD_Pt_15to3000_M2_5_500, QCD_Pt_15to3000],
      ]:
    # Read all run/lumi/event for first tree
    r1 = sample1.treeReader( variables = variables , selectionString = "1")
    r1.start()
    position_r1 = {}
    while r1.run():
        position_r1[(r1.data.run, r1.data.lumi, r1.data.evt)] = r1.eList.GetEntry(r1.position-1)

    # Read all run/lumi/event for second tree
    r2 = sample2.treeReader( variables = variables , selectionString = "1" )
    r2.start()
    position_r2={}
    while r2.run():
        position_r2[(r2.data.run, r2.data.lumi, r2.data.evt)] = r2.eList.GetEntry(r2.position-1)

    logger.info( "Have %i events in first samle and %i in second", len(position_r1), len(position_r2) )

    # Fast intersect
    intersec = set(position_r1.keys()).intersection(set(position_r2.keys()))
    positions = [(position_r1[i], position_r2[i]) for i in intersec]
    logger.info("Have %i events in common.", len(intersec))

    # Make event list from common events
    event_list1 = ROOT.TEventList("l1","l1")
    event_list2 = ROOT.TEventList("l2","l2")
    #Loop over pairs
    for nEvent, pos in enumerate(positions):
        if nEvent%10000==0: 
            logger.info("At event %i/%i", nEvent, len(intersec))

        p1, p2 = pos
        event_list1.Enter(p1)
        event_list2.Enter(p2)

    # Setthose 
    r1.setEventList(event_list1)
    r2.setEventList(event_list2)

    r1.start()
    r2.start()
    while r1.run() and r2.run():
        n = min([r1.data.nJet, r2.data.nJet])
        jets1 = [{'rawPt':r1.data.Jet_rawPt[i], 'eta':r1.data.Jet_eta[i], 'phi':r1.data.Jet_eta[i]} for i in xrange(n) if abs(r1.data.Jet_eta[i])<1.3]
        jets2 = [{'rawPt':r2.data.Jet_rawPt[i], 'eta':r2.data.Jet_eta[i], 'phi':r2.data.Jet_eta[i]} for i in xrange(n) if abs(r2.data.Jet_eta[i])<1.3]
        for c in zip(jets1, jets2):
            if deltaR2(*c)<0.2**2:
                jetResponse.Fill(c[1]['rawPt'],c[0]['rawPt']/c[1]['rawPt'])

# Make plot
prefix="test"
jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponse_eta13", histos = [[jetResponse_mc], [jetResponse_data]], texX = "raw Jet p_{T}", texY = "response ratio M2(5,500)/M2(0,100)" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = {'yRange':(0.985,1.025)}, logX = True, logY = False, sorting = False, yRange = (0.965, 1.045) )
