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

max_events = -1
max_files = -1

eta_bins = [0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191]

pt_thresholds = [ 10**(x/10.) for x in range(11,36) ] 
#eta_thresholds = [0, 1.3, 2.5, 3.2]
eta_thresholds = [2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839]
color={0:ROOT.kBlack, 1.3:ROOT.kBlue, 2.5:ROOT.kRed, 3.2:ROOT.kMagenta }
color={2.500:ROOT.kBlack, 2.650:ROOT.kBlue, 2.853:ROOT.kRed, 2.964:ROOT.kMagenta, 3.139:ROOT.kGreen, 3.489:ROOT.kOrange, 3.839:ROOT.kAzure }
resp={}
for i_eta_th, eta_th in enumerate( eta_thresholds ):
    resp[eta_th] = ROOT.TProfile("responses", "response", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
    resp[eta_th].style = styles.lineStyle(color[eta_th] )
    resp[eta_th].legendText = "%4.3f<=#eta"%eta_th
    if eta_th!=eta_thresholds[-1]: resp[eta_th].legendText += "<%4.3f"%eta_thresholds[i_eta_th+1]

#resp_eta = ROOT.TProfile("response_eta", "response_eta", 52,-5.2,5.2 )
resp_eta = ROOT.TProfile("response_eta", "response_eta", len(eta_bins)-1, array.array('d', eta_bins) )

JetHT_22Feb2017 = Sample.fromDirectory("JetHT_22Feb2017", "/scratch/rschoefbeck/cmgTuples/legacy_JetHT/JetHT_Run2016H_22Feb2017/", maxN = max_files, treeName = "tree")
JetHT_03Feb2016 = Sample.fromDirectory("JetHT_03Feb2016", "/scratch/rschoefbeck/cmgTuples/legacy_JetHT/JetHT_Run2016H-03Feb2017_ver2-v1/JetHT_Run2016H-03Feb2017_ver2-v1", maxN = max_files, treeName = "tree")

jet_str = "pt/F,eta/F,rawPt/F,phi/F,id/I"
read_variables = ["evt/l", "run/I", "lumi/I", "nJet/I", "nVert/I", "Jet[%s]"%jet_str]

r1 = JetHT_22Feb2017.treeReader( variables = map( TreeVariable.fromString, read_variables) , selectionString = "run==281693")
r2 = JetHT_03Feb2016.treeReader( variables = map( TreeVariable.fromString, read_variables) , selectionString = "run==281693")

from math import pi

jetVars = ['pt', 'eta', 'rawPt', 'phi', 'id']

def getJets(c, jetVars=jetVars ):
    return [ {var:getattr(c, "Jet_"+var)[i] for var in jetVars} for i in range( getattr( c, "nJet" ) )  ]

r1.start()
position_r1 = {}
count=0
while r1.run():
    position_r1[(r1.event.run, r1.event.lumi, r1.event.evt) ] = r1.position-1
    count+=1
    if max_events is not None and max_events>0 and count>=max_events:break

r2.start()
position_r2 = {}
count=0
while r2.run():
    position_r2[(r2.event.run, r2.event.lumi, r2.event.evt) ] = r2.position-1
    count+=1
    if max_events is not None and max_events>0 and count>=max_events:break

logger.info( "Have %i events in first samle and %i in second", len(position_r1), len(position_r2) )

# Fast intersect
intersec = set(position_r1.keys()).intersection(set(position_r2.keys()))
positions = [(position_r1[i], position_r2[i]) for i in intersec]

if len(positions)==0:
    logger.error( "Found no common events." )
    sys.exit(0)

# Without sorting, there is a jump between files with almost every event -> extremly slow
positions.sort()
logger.info("Have %i events in common.", len(intersec))

#Looping over common events
for i, p in enumerate(positions):
    p1,p2 = p
    r1.goToPosition(p1)
    r2.goToPosition(p2)
    if i%10000==0: logger.info("At %i/%i of common events.", i, len(positions))
    jets1 = filter( lambda j:j['id'], getJets( r1.event ) ) 
    jets2 = filter( lambda j:j['id'], getJets( r2.event ) )
    for c in zip(jets1, jets2):
        if c[1]['rawPt']<=20: continue
        if helpers.deltaR2(*c)<0.2**2:
            resp_eta.Fill( c[0]['eta'], c[0]['rawPt']/c[1]['rawPt'] )
            for eta_th in reversed(eta_thresholds):
                if abs(c[0]['eta'])>eta_th:
                    resp[eta_th].Fill( c[0]['rawPt'], c[0]['rawPt']/c[1]['rawPt'] )
                    break
# Make eta plot
profiles = [resp[eta_th] for eta_th in eta_thresholds ]
prefix= "response" + ("_max_events_%s_"%max_events if max_events is not None and max_events>0 else "" )
histos = [ [h.ProjectionX()] for h in profiles ]
for i, h in enumerate(histos):
    h[0].__dict__.update(profiles[i].__dict__)

jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponseRatio_JetHT", histos = histos, texX = "raw jet p_{T} (legacy)" , texY = "response ratio legacy/reminiaod" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = True, yRange=(0.95,1.25))

jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponseRatio_eta_JetHT", histos = [[resp_eta.ProjectionX()]], texX = " Jet #eta (legacy)" , texY = "response ratio legacy/reminiaod" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = False, yRange=(0.85,1.2))
