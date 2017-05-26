#!/usr/bin/env python

# standard imports
import ROOT
import sys
import os
import copy
import subprocess
import shutil

from operator import mul
from math import sqrt, atan2, sin, cos

# RootTools
from RootTools.core.standard import *

# User specific
import JetMET.tools.user as user

#from StopsDilepton.tools.objectSelection import getMuons, getElectrons, muonSelector, eleSelector, getGoodLeptons, getGoodAndOtherLeptons,  getGoodBJets, getGoodJets, isBJet, jetId, isBJet, getGoodPhotons, getGenPartsAll, multiIsoWPInt

# central configuration
targetLumi = 1000 #pb-1 Which lumi to normalize to

def get_parser():
    ''' Argument parser for post-processing module.
    '''
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")

    argParser.add_argument('--logLevel',
        action='store',
        nargs='?',
        choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],
        default='INFO',
        help="Log level for logging"
        )

    argParser.add_argument('--overwrite',
        action='store_true',
#        default = True,
        help="Overwrite existing output files, bool flag set to True  if used")

    argParser.add_argument('--sample',
        action='store',
        nargs=1,
        type=str,
        default='JetHT_Run2016B',
        help="List of samples to be processed as defined in L2res_master"
        )

    argParser.add_argument('--eventsPerJob',
        action='store',
        nargs='?',
        type=int,
        default=300000,
        help="Maximum number of events per job (Approximate!)."
        )

    argParser.add_argument('--nJobs',
        action='store',
        nargs='?',
        type=int,
        default=1,
        help="Maximum number of simultaneous jobs."
        )
    argParser.add_argument('--job',
        action='store',
        nargs='*',
        type=int,
        default=[],
        help="Run only jobs i"
        )

    argParser.add_argument('--minNJobs',
        action='store',
        nargs='?',
        type=int,
        default=1,
        help="Minimum number of simultaneous jobs."
        )

    argParser.add_argument('--targetDir',
        action='store',
        nargs='?',
        type=str,
        default=user.data_output_directory,
        help="Name of the directory the post-processed files will be saved"
        )

    argParser.add_argument('--processingEra',
        action='store',
        nargs='?',
        type=str,
        default='V1',
        help="Name of the processing era"
        )

    argParser.add_argument('--skim',
        action='store',
        nargs='?',
        type=str,
        default='default',
        help="Skim conditions to be applied for post-processing"
        )

    argParser.add_argument('--small',
        action='store_true',
#        default=True,
        help="Run the file on a small sample (for test purpose), bool flag set to True if used",
        #default = True
        )

    return argParser

options = get_parser().parse_args()

# Logging
import JetMET.tools.logger as logger
logFile = '/tmp/%s_%s_%s_njob%s.txt'%(options.skim, '_'.join(options.samples), os.environ['USER'], str(0 if options.nJobs==1 else options.job[0]))
logger  = logger.get_logger(options.logLevel, logFile = logFile)

import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(options.logLevel, logFile = None )

# Flags 
defSkim = options.skim.lower().startswith('default')

# Skim condition
skimConds = []
if defSkim:
    skimConds.append( "0.5*(Jet_pt[0]+Jet_pt[1])>50" )

#Samples: Load samples
maxN = 1 if options.small else None
from JetMET.JEC.samples.L2res_master import L2res_master
sample = Sample.combine( options.sample, samples = [ Sample.fromCMGOutput("%i"%i_sample, data_path = path, maxN = maxN) for i_sample, path in enumerate( L2res_master[options.sample] ) ] )
logger.debug("Reading from CMG tuple %s which are %i files.", options.sample, len(sample.files) )
    
isData = 'Run2016' in sample.name 
isMC   =  not isData 

#if isMC:
#    from StopsDilepton.tools.puReweighting import getReweightingFunction
#    mcProfile = "Summer16"
#    # nTrueIntReweighting
#    nTrueInt36fb_puRW        = getReweightingFunction(data="PU_2016_36000_XSecCentral", mc=mcProfile)
#    nTrueInt36fb_puRWDown    = getReweightingFunction(data="PU_2016_36000_XSecDown",    mc=mcProfile)
#    nTrueInt36fb_puRWUp      = getReweightingFunction(data="PU_2016_36000_XSecUp",      mc=mcProfile)
        


## systematic variations
#addSystematicVariations = (not isData) and (not options.skipSystematicVariations)
#if addSystematicVariations:
#    # B tagging SF
#    from StopsDilepton.tools.btagEfficiency import btagEfficiency
#    btagEff = btagEfficiency( fastSim = False )

## LHE cut (DY samples)
#if options.LHEHTCut>0:
#    sample.name+="_lheHT"+str(options.LHEHTCut)
#    logger.info( "Adding upper LHE cut at %f", options.LHEHTCut )
#    skimConds.append( "lheHTIncoming<%f"%options.LHEHTCut )

directory  = os.path.join(options.targetDir, options.processingEra) 

output_directory = os.path.join( directory, options.skim, sample.name )

if os.path.exists(output_directory) and options.overwrite:
    if options.nJobs > 1:
        logger.warning( "NOT removing directory %s because nJobs = %i", output_directory, options.nJobs )
    else:
        logger.info( "Output directory %s exists. Deleting.", output_directory )
        shutil.rmtree(output_directory)

try:    #Avoid trouble with race conditions in multithreading
    os.makedirs(output_directory)
    logger.info( "Created output directory %s.", output_directory )
except:
    pass

#branches to be kept for data and MC
branchKeepStrings_DATAMC = [\
    "run", "lumi", "evt", "isData", "rho", "nVert",
    "met_pt", "met_phi","met_Jet*", "met_Unclustered*", "met_sumEt", "met_rawPt","met_rawPhi", "met_rawSumEt", "met_caloPt", "met_caloPhi",
#        "metNoHF_pt", "metNoHF_phi",
#        "puppiMet_pt","puppiMet_phi","puppiMet_sumEt","puppiMet_rawPt","puppiMet_rawPhi","puppiMet_rawSumEt",
    "Flag_*","HLT_*",
#        "nDiscJet", "DiscJet_*",
#        "nJetFailId", "JetFailId_*",
    "nLepGood", "LepGood_*",
    "nLepOther", "LepOther_*",
#        "nTauGood", "TauGood_*",
]
#branches to be kept for MC samples only
branchKeepStrings_MC = [\
    "nTrueInt", "genWeight", "xsec", "met_gen*", "lheHTIncoming",
#        "ngenPartAll","genPartAll_*","ngenLep","genLep_*"
]


if isMC:
    if isTiny or isSmall:
        jetMCInfo = ['mcPt/F', 'hadronFlavour/I','mcMatchId/I']
    else:
        jetMCInfo = ['mcMatchFlav/I', 'partonId/I', 'partonMotherId/I', 'mcPt/F', 'mcFlavour/I', 'hadronFlavour/I', 'mcMatchId/I']
        if not (options.susySignal):
            jetMCInfo.append('partonFlavour/I')
else:
    jetMCInfo = []

if isData:
    lumiScaleFactor=None
    branchKeepStrings = branchKeepStrings_DATAMC + branchKeepStrings_DATA
    from FWCore.PythonUtilities.LumiList import LumiList
    # Apply golden JSON
    json = '$CMSSW_BASE/src/CMGTools/TTHAnalysis/data/json/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
    lumiList = LumiList(os.path.expandvars(json))
    logger.info( "Loaded json %s", json )
else:
    lumiScaleFactor = targetLumi/float(sample.normalization) 
    branchKeepStrings = branchKeepStrings_DATAMC + branchKeepStrings_MC

jetVars = ['pt/F', 'rawPt/F', 'eta/F', 'phi/F', 'id/I', 'btagCSV/F', 'area/F'] + jetCorrInfo + jetMCInfo
jetVarNames = [x.split('/')[0] for x in jetVars]

read_variables = map(TreeVariable.fromString, ['met_pt/F', 'met_phi/F', 'run/I', 'lumi/I', 'evt/l', 'nVert/I'] )

new_variables = [ 'weight/F']
if isMC:
    read_variables+= [TreeVariable.fromString('nTrueInt/F')]

read_variables += [\
#    TreeVariable.fromString('nLepGood/I'),
#    VectorTreeVariable.fromString('LepGood[pt/F,eta/F,etaSc/F,phi/F,pdgId/I,tightId/I,miniRelIso/F,relIso03/F,sip3d/F,mediumMuonId/I,mvaIdSpring15/F,lostHits/I,convVeto/I,dxy/F,dz/F,jetPtRelv2/F,jetPtRatiov2/F,eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz/I]'),
#    TreeVariable.fromString('nJet/I'),
#    VectorTreeVariable.fromString('Jet[%s]'% ( ','.join(jetVars) ) )
]

new_variables += [\
    'JetGood[%s]'% ( ','.join(jetVars) )
]

if isData: new_variables.extend( ['jsonPassed/I'] )

# Define a reader
reader = sample.treeReader( \
    variables = read_variables ,
    selectionString = "&&".join(skimConds)
    )

def filler( event ):
    # shortcut
    r = reader.event
    elif isMC:
        event.weight = lumiScaleFactor*r.genWeight if lumiScaleFactor is not None else 1

    # lumi lists and vetos
    if isData:
        #event.vetoPassed  = vetoList.passesVeto(r.run, r.lumi, r.evt)
        event.jsonPassed  = lumiList.contains(r.run, r.lumi)
        # store decision to use after filler has been executed
        event.jsonPassed_ = event.jsonPassed

    if isMC:
        event.reweightPU     = nTrueInt_puRW       ( r.nTrueInt ) 
        event.reweightPUDown = nTrueInt_puRWDown   ( r.nTrueInt ) 
        event.reweightPUUp   = nTrueInt_puRWUp     ( r.nTrueInt ) 

# Create a maker. Maker class will be compiled. This instance will be used as a parent in the loop
treeMaker_parent = TreeMaker(
    sequence  = [ filler ],
    variables = [ TreeVariable.fromString(x) for x in new_variables ],
    treeName = "Events"
    )

# Split input in ranges
if options.nJobs>1:
    eventRanges = reader.getEventRanges( nJobs = options.nJobs )
else:
    eventRanges = reader.getEventRanges( maxNEvents = options.eventsPerJob, minJobs = options.minNJobs )

logger.info( "Splitting into %i ranges of %i events on average.",  len(eventRanges), (eventRanges[-1][1] - eventRanges[0][0])/len(eventRanges) )

#Define all jobs
jobs = [(i, eventRanges[i]) for i in range(len(eventRanges))]

filename, ext = os.path.splitext( os.path.join(output_directory, sample.name + '.root') )

clonedEvents = 0
convertedEvents = 0
outputLumiList = {}
for ievtRange, eventRange in enumerate( eventRanges ):

    if len(options.job)>0 and not ievtRange in options.job: continue

    logger.info( "Processing range %i/%i from %i to %i which are %i events.",  ievtRange, len(eventRanges), eventRange[0], eventRange[1], eventRange[1]-eventRange[0] )

    # Check whether file exists
    outfilename = filename+'_'+str(ievtRange)+ext
    if os.path.isfile(outfilename):
        logger.info( "Output file %s found.", outfilename)
        if not checkRootFile(outfilename, checkForObjects=["Events"]):
            logger.info( "File %s is broken. Overwriting.", outfilename)
        elif not options.overwrite:
            logger.info( "Skipping.")
            continue
        else:
            logger.info( "Overwriting.")

    tmp_directory = ROOT.gDirectory
    outputfile = ROOT.TFile.Open(outfilename, 'recreate')
    tmp_directory.cd()

    if options.small: 
        logger.info("Running 'small'. Not more than 10000 events") 
        nMaxEvents = eventRange[1]-eventRange[0]
        eventRange = ( eventRange[0], eventRange[0] +  min( [nMaxEvents, 10000] ) )

    # Set the reader to the event range
    reader.setEventRange( eventRange )

    clonedTree = reader.cloneTree( branchKeepStrings, newTreename = "Events", rootfile = outputfile )
    clonedEvents += clonedTree.GetEntries()
    # Clone the empty maker in order to avoid recompilation at every loop iteration
    maker = treeMaker_parent.cloneWithoutCompile( externalTree = clonedTree )

    maker.start()
    # Do the thing
    reader.start()

    while reader.run():
        maker.run()
        if isData:
            if maker.event.jsonPassed_:
                if reader.event.run not in outputLumiList.keys():
                    outputLumiList[reader.event.run] = set([reader.event.lumi])
                else:
                    if reader.event.lumi not in outputLumiList[reader.event.run]:
                        outputLumiList[reader.event.run].add(reader.event.lumi)

    convertedEvents += maker.tree.GetEntries()
    maker.tree.Write()
    outputfile.Close()
    logger.info( "Written %s", outfilename)

  # Destroy the TTree
    maker.clear()

logger.info( "Converted %i events of %i, cloned %i",  convertedEvents, reader.nEvents , clonedEvents )

# Storing JSON file of processed events
if isData:
    jsonFile = filename+'.json'
    LumiList( runsAndLumis = outputLumiList ).writeJSON(jsonFile)
    logger.info( "Written JSON file %s",  jsonFile )

logger.info("Copying log file to %s", output_directory )
copyLog = subprocess.call(['cp', logFile, output_directory] )
if copyLog:
    logger.info( "Copying log from %s to %s failed", logFile, output_directory)
else:
    logger.info( "Successfully copied log file" )
    os.remove(logFile)
    logger.info( "Removed temporary log file" )

