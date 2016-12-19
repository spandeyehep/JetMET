''' FWLiteReader example: Loop over a sample and write some data to a histogram.
'''
# Standard imports
import os
import logging
import ROOT
from math import cos, sin, atan2, sqrt

#pdgToName
from pdgToName import pdgToName

#RootTools
from RootTools.core.standard import *

#StopsDilepton
from StopsDilepton.tools.mcTools import pdgToName
from StopsDilepton.tools.helpers import deltaR

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

# Helper
def toDict(p):
    return {'pt':p.pt(), 'eta':p.eta(), 'phi':p.phi(), 'pdgId':p.pdgId()}

def select_in_cone( p, collection, dR = 0.3 ):
    return [q for q in collection if deltaR(toDict(p), toDict(q)) < dR  ]

def vecSumPt(particles):
    px = sum(p.px() for p in particles)
    py = sum(p.py() for p in particles)
    return sqrt(px**2+py**2)

def bold(s):
    return '\033[1m'+s+'\033[0m'

## 8X mAOD, assumes eos mount in home directory 
## from Directory
dirname = "/data/rschoefbeck/pickEvents/StopsDilepton/" 
prompt = FWLiteSample.fromFiles("prompt", files = [ \
    "root://eoscms.cern.ch//store/user/zdemirag/DoubleMuon/crab_pickEventsPrompt/161219_193844/0000/pickevents_prompt_%i.root"%i for i in range(1,6) 
    ])

rereco = FWLiteSample.fromFiles("rereco", files = [ \
    "root://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/schoef/DoubleMuon/crab_pickEvents/161219_192421/0000/pickevents_%i.root"%i for i in range(1,24) + range(25, 116) 
    ])

products = {
    'pfCands':{'type':'vector<pat::PackedCandidate>', 'label':"packedPFCandidates"},
    'pfJets':{'type':'vector<pat::Jet>', 'label': ("slimmedJets")},
    'pfMet':{'type':'vector<pat::MET>','label':( "slimmedMETs" )},
    'electrons':{'type':'vector<pat::Electron>','label':( "slimmedElectrons" )},
    'muons':{'type':'vector<pat::Muon>', 'label':("slimmedMuons") },
}

# Make reader for FWLite products
r1 = prompt.fwliteReader( products = products )
r2 = rereco.fwliteReader( products = products )

# Align the datasets
# Run over 1st dataset
position_r1 = {}
r1.start()
while r1.run( readProducts = False ):
    position_r1[r1.evt] = r1.position-1

# Run over 2nd dataset
position_r2 = {}
r2.start()
while r2.run( readProducts = False ):
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
    if 2259466578 in r1.evt: 
        logger.info( "Found evt %i:%i:%i. MET: prompt %3.2f rereco: %3.2f"% (r1.evt[0], r1.evt[1], r1.evt[2], r1.products['pfMet'][0].pt(), r2.products['pfMet'][0].pt() ) )
        break
        
cands_prompt = list(r1.products['pfCands'])
cands_rereco = list(r2.products['pfCands'])
cands_prompt.sort( key = lambda p: -p.pt() )
cands_rereco.sort( key = lambda p: -p.pt() )

for p in cands_rereco[:10]:
    if p.pt()<10: continue
    c_prompt = select_in_cone( p, cands_prompt )
    c_rereco = select_in_cone( p, cands_rereco )

    logger.info( "Rereco particle %s pdgId %i pt %3.2f eta %3.2f phi %3.2f", pdgToName(p.pdgId()), p.pdgId(), p.pt(), p.eta(), p.phi() )
    for q in c_prompt:
        logger.info( "  Cone 0.3 prompt %s pdgId %i pt %3.2f eta %3.2f phi %3.2f", pdgToName(q.pdgId()), q.pdgId(), q.pt(), q.eta(), q.phi() )
    logger.info( "  Total sumPt %3.2f", vecSumPt(c_prompt) ) 

    for q in c_rereco:
        logger.info( "  Cone 0.3 rereco %s pdgId %i pt %3.2f eta %3.2f phi %3.2f", pdgToName(q.pdgId()), q.pdgId(), q.pt(), q.eta(), q.phi() )
    logger.info( "  Total sumPt %3.2f", vecSumPt(c_rereco) ) 
         
