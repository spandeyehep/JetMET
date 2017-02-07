''' FWLite example
'''
# Standard imports
import ROOT
from DataFormats.FWLite import Events, Handle
from PhysicsTools.PythonAnalysis import *

small = True

# example file
events = Events(['root://eoscms.cern.ch//store/relval/CMSSW_8_0_25/JetHT/RECO/2016_12_21_06_38_PRnewco_80X_dataRun2_2016LegacyRepro_Candidate_v2-v2/10000/007FD8F2-E8CA-E611-8B61-0025905B85D2.root'])

# Use 'edmDumpEventContent <file>' to list all products. Then, make a simple dictinary as below for the products you want to read.
# These are the PF rechit collections:
#vector<reco::PFRecHit>                "particleFlowRecHitECAL"    "Cleaned"         "reRECO"   
#vector<reco::PFRecHit>                "particleFlowRecHitHBHE"    "Cleaned"         "reRECO"   
#vector<reco::PFRecHit>                "particleFlowRecHitHF"      "Cleaned"         "reRECO"   
#vector<reco::PFRecHit>                "particleFlowRecHitHO"      "Cleaned"         "reRECO"   
#vector<reco::PFRecHit>                "particleFlowRecHitPS"      "Cleaned"         "reRECO" 

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
    'pfMet':        { 'label':('pfMet'), 'type':'vector<reco::PFMET>'},
    #'pfRecHitsHBHE':{ 'label':("particleFlowRecHitHBHE"), 'type':"vector<reco::PFRecHit>"},
    #'caloRecHits':  { 'label':("reducedHcalRecHits"), 'type':'edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >'},
    'clusterHCAL':  {  'label': "particleFlowClusterHCAL", "type":"vector<reco::PFCluster>"},
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
    products[k] = v['handle'].product()

  print run,lumi,event
  
  #print RecHits
  for i, cl in enumerate(products["clusterHCAL"]):
    print "cluster   n %i E %3.2f"%(i, cl.energy())
  #for i, rh in enumerate(products["caloRecHits"]):
  #  print "caloRechit n %i E %3.2f"%(i, rh.energy())
