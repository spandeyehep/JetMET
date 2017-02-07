''' FWLite example
'''
# Standard imports
import ROOT
from DataFormats.FWLite import Events, Handle
from PhysicsTools.PythonAnalysis import *

small = True

# Use 'edmDumpEventContent <file>' to list all products. Then, make a simple dictinary as below for the products you want to read.

# example file
events = Events(['root://eoscms.cern.ch//store/relval/CMSSW_8_0_25/JetHT/RECO/2016_12_21_06_38_PRnewco_80X_dataRun2_2016LegacyRepro_Candidate_v2-v2/10000/007FD8F2-E8CA-E611-8B61-0025905B85D2.root'])

# miniAOD
#products = {
#    'pfCands':{'type':'vector<pat::PackedCandidate>', 'label':"packedPFCandidates"},
#    'pfJets':{'type':'vector<pat::Jet>', 'label': ("slimmedJets")},
#    'pfMet':{'type':'vector<pat::MET>','label':( "slimmedMETs" )},
#    'electrons':{'type':'vector<pat::Electron>','label':( "slimmedElectrons" )},
#    'muons':{'type':'vector<pat::Muon>', 'label':("slimmedMuons") },
#}

# RECO
edmCollections = { 
    'pfRecHitsHBHE':{ 'label':("particleFlowRecHitHBHE"), 'type':"vector<reco::PFRecHit>"},
    'pfMet':        { 'label':('pfMet'), 'type':'vector<reco::PFMET>'},
    'caloRecHits':  { 'label':("reducedHcalRecHits"), 'type':'edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >'},
    'pf':           { 'label':('particleFlow'), 'type':'vector<reco::PFCandidate>'},
   }

# add handles
for k, v in edmCollections.iteritems():
    v['handle'] = Handle(v['type'])

nevents = 1 if small else events.size()

for i in range(nevents):
  events.to(i)

  eaux  = events.eventAuxiliary()

  # run/lumi/event
  run   = eaux.run()
  event = eaux.event()
  lumi  = eaux.luminosityBlock()

  #read all products as specifed in edmCollections
  products = {}
  for k, v in edmCollections.iteritems():
    events.getByLabel(v['label'], v['handle'])
    products[v['name']] = v['handle'].product()

  print run,lumi,event
  
  #printRecHits
  for i, rh in enumerate(products["caloRecHits"]):
    print "caloRechit n %i E %3.2f"%(i, rh.energy())
  for i, rh in enumerate(products["pfRecHitsHBHE"]):
    print "pfRecHit   n %i E %3.2f"%(i, rh.energy())
