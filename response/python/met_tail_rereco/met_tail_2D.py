''' FWLiteReader example: Loop over a sample and write some data to a histogram.
'''
# Standard imports
import os
import logging
import ROOT
import itertools
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

h_dr_prompt = ROOT.TH1D( "dr_prompt", "dr", 200,0,2)
h_dr_rereco = ROOT.TH1D( "dr_rereco", "dr", 200,0,2)
h_dr_2D     = ROOT.TH2D( "dr_2D", "dr", 100,0,2,200,0,2)

from zeynep import passing, prompt_files, rereco_files
prompt = FWLiteSample.fromFiles("prompt", files = prompt_files)
rereco = FWLiteSample.fromFiles("rereco", files = rereco_files)

#evt = 2259466578 # the first duplicated, Robert, Dec. 19th 
#evt = 1821418817 # Zeynep Dec. 20th
#

products = {
    'pfCands':{'type':'vector<pat::PackedCandidate>', 'label':"packedPFCandidates"},
    'pfJets':{'type':'vector<pat::Jet>', 'label': ("slimmedJets")},
    'pfMet':{'type':'vector<pat::MET>','label':( "slimmedMETs" )},
    'electrons':{'type':'vector<pat::Electron>','label':( "slimmedElectrons" )},
    'muons':{'type':'vector<pat::Muon>', 'label':("slimmedMuons") },
}

# Make reader for FWLite products
r_prompt = prompt.fwliteReader( products = products )
r_rereco = rereco.fwliteReader( products = products )

# Align the datasets
# Run over 1st dataset
position_r_prompt = {}
r_prompt.start()
while r_prompt.run( readProducts = False ):
    position_r_prompt[r_prompt.evt] = r_prompt.position-1

# Run over 2nd dataset
position_r_rereco = {}
r_rereco.start()
while r_rereco.run( readProducts = False ):
    position_r_rereco[r_rereco.evt] = r_rereco.position-1

logger.info( "Have %i events in first samle and %i in second", len(position_r_prompt), len(position_r_rereco) )

# Fast intersect
intersec = set(position_r_prompt.keys()).intersection(set(position_r_rereco.keys()))
positions = [(position_r_prompt[i], position_r_rereco[i]) for i in intersec]

# Withou sorting, there is a jump between files with almost every event -> extremly slow
positions.sort()
logger.info("Have %i events in common.", len(intersec))

def minDR( particles ):
    if len( particles )<2: return 999.
    return min( deltaR( toDict(c[0]), toDict(c[1]) ) for c in itertools.combinations( particles, 2 ) )

#Looping over common events
for i, p in enumerate(positions):
    p1,p2 = p
    r_prompt.goToPosition(p1)
    r_rereco.goToPosition(p2)
    if r_rereco.evt[-1] in passing: 
        logger.info( "Taking %i:%i:%i, because it's not filtered"%r_prompt.evt )
        muons_prompt = filter( lambda p: abs(p.pdgId())==13, r_prompt.products['pfCands'])
        min_dr_prompt = minDR( muons_prompt )
        h_dr_prompt.Fill(min_dr_prompt)
        muons_rereco = filter( lambda p: abs(p.pdgId())==13, r_rereco.products['pfCands'])
        min_dr_rereco = minDR( muons_rereco )
        h_dr_rereco.Fill(min_dr_rereco)

        h_dr_2D.Fill(min_dr_prompt, min_dr_rereco)
    else:
        logger.info( "Skipping %i:%i:%i, because it's filtered"%r_prompt.evt )

h_dr_prompt.style = styles.lineStyle( ROOT.kBlue )
h_dr_prompt.legendText = "prompt"
h_dr_rereco.style = styles.lineStyle( ROOT.kRed )
h_dr_rereco.legendText = "rereco"
        
p_dr = Plot.fromHisto(name = "min_dr", histos =  [[h_dr_prompt],[h_dr_rereco]] , texX = "minDR(PF mu)" )  
plotting.draw(p_dr, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", logX = False, logY = False)

c1 = ROOT.TCanvas()
h_dr_2D.Draw("COLZ")
h_dr_2D.GetXaxis().SetTitle("prompt")
h_dr_2D.GetYaxis().SetTitle("rereco")
c1.SetPadRightMargin(0.15)
c1.Print( "/afs/hephy.at/user/r/rschoefbeck/www/etc/min_dr_2D.png" )
