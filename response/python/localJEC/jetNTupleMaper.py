#Standard imports
import sys
import logging
import ROOT
import array
from math import log
import os
import pickle

#RootTools
from RootTools.core.Sample import Sample
from RootTools.core.Variable import Variable, ScalarType, VectorType
from RootTools.core.TreeReader import TreeReader
from RootTools.core.TreeMaker import TreeMaker
from RootTools.core.logger import get_logger
from RootTools.plot.Plot import Plot
from RootTools.plot.styles import lineStyle
import RootTools.plot.plotting as plotting

#Helper
from JetMET.tools.helpers import deltaR2, getObjFromFiles

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

argParser.add_argument('--outputDir',
    dest='outputDir',
    action='store',
    nargs='?',
    type=str,
    default='/data/rschoefbeck/JetMET/localJEC/',
    help="output directory"
    )
args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

#make samples. 
maxN = -1
#qcd_AllChGood_Flat0to50 = Sample.fromCMGOutput("qcd_AllChGood_Flat0to50", "/data/rschoefbeck/cmgTuples/MC25ns_v2_0l/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_schoef-crab_QCD_Pt-15to7000_AllChlgoodAsymptFlat0to50bx25_AllChannelsGood_v0-v2_RunIISpring15MiniAODv2-74X", treeFilename='tree.root', maxN=maxN)
#qcd_Flat0to50 = Sample.fromCMGOutput("qcd_Flat0to50", "/data/rschoefbeck/cmgTuples/MC25ns_v2_0l/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_schoef-crab_QCD_Pt-15to7000_AsymptFlat0to50bx25Reco_MCRUN2_74_V9-v3_RunIISpring15MiniAODv2-74X_3", treeFilename='tree.root', maxN=maxN)
qcd_AllChGood_noPU = Sample.fromCMGOutput("qcd_AllChGood_noPU", "/data/rschoefbeck/cmgTuples/MC25ns_v2_0l/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_schoef-crab_QCD_Pt-15to7000_AllChlgoodAsymptNoPUbx25_AllChannelsGood_v0-v2_RunIISpring15MiniAODv2-74X", treeFilename='tree.root', maxN=maxN)
qcd_noPU = Sample.fromCMGOutput("qcd_noPU", "/data/rschoefbeck/cmgTuples/MC25ns_v2_0l/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_schoef-crab_QCD_Pt-15to7000_AsymptNoPUbx25Reco_MCRUN2_74_V9-v3_RunIISpring15MiniAODv2-74X_3", treeFilename='tree.root', maxN=maxN)

# Select samples
sample_ideal, sample_real = qcd_AllChGood_noPU, qcd_noPU

# Pickle file which will hold the commonf events
positions_filename = os.path.join(args.outputDir, "common_positions_%s_%s.pkl"%(sample_ideal.name, sample_real.name) )

# Make common event lists
if not os.path.isfile(positions_filename):
    logger.info( "File %s not found, therefore finding positions of common events.", positions_filename)
    variables =  map( Variable.fromString, ['evt/l','run/I', 'lumi/I'] )

    reader_ideal = sample_ideal.treeReader( variables = variables , selectionString = "1")
    reader_ideal.start()
    position_reader_ideal = {}
    while reader_ideal.run():
        position_reader_ideal[(reader_ideal.data.run, reader_ideal.data.lumi, reader_ideal.data.evt)] = reader_ideal.eList.GetEntry(reader_ideal.position-1)

    # Read all run/lumi/event for second tree
    reader_real = sample_real.treeReader( variables = variables , selectionString = "1" )
    reader_real.start()
    position_reader_real={}
    while reader_real.run():
        position_reader_real[(reader_real.data.run, reader_real.data.lumi, reader_real.data.evt)] = reader_real.eList.GetEntry(reader_real.position-1)

    logger.info( "Have %i events in first samle and %i in second", len(position_reader_ideal), len(position_reader_real) )

    # Fast intersect
    intersec = list(set(position_reader_ideal.keys()).intersection(set(position_reader_real.keys())))
    intersec.sort()
    positions = [(position_reader_ideal[i], position_reader_real[i]) for i in iter(intersec)]
    logger.info("Have %i events in common.", len(intersec))
    
    # Save positions
    pickle.dump(positions, file(positions_filename, "w") )
else:
    logger.info( "File %s found, loading positions.", positions_filename)
    positions = pickle.load( file(positions_filename) )

# Make n-tuples per jet
jetComponents = ['pt/F', 'eta/F', 'phi/F', 'area/F', 'ptd/F', 'axis2/F', 'mult/I', 'partonId/I', 'partonMotherId/I', \
                 'nLeptons/F', 'id/I', 'btagCSV/F', 'btagCMVA/F', 'rawPt/F', 'mcPt/F', 'mcFlavour/I', 'hadronFlavour/I', 'mcMatchId/I', 
                 'corr_JECUp/F', 'corr_JECDown/F', 'corr/F', 'mass/F', 'mcMatchFlav/I', 
                 'chHEF/F', 'neHEF/F', 'phEF/F', 'eEF/F', 'muEF/F', 'HFHEF/F', 'HFEMEF/F', 
                 'chHMult/I', 'neHMult/I', 'phMult/I', 'eMult/I', 'muMult/I', 'HFHMult/I', 'HFEMMult/I', 'charge/F']

jetComponentNames = [c.split('/')[0] for c in jetComponents]

jet_string = 'Jet[%s]'%(','.join(jetComponents))
variables =  map( Variable.fromString, ['evt/l','run/I', 'lumi/I', jet_string, 'nJet/I'] )

#Maker (only scalar variables since we store per jet)
jet_variables =  map( Variable.fromString, [x.replace('/','_ideal/') for x in jetComponents] + [x.replace('/','_real/') for x in jetComponents] ) 

# Wrapper for threading
def wrapper(job):

    nJob, positions = job

    # Checking whether the jobs was already done
    outputFilename = os.path.join(args.outputDir, "jetTrees", "tree_%s_%s_%i.root"%(sample_ideal.name, sample_real.name, nJob) )
    if  os.path.exists(outputFilename):
        logger.info( "Found file %s -> Skipping job.", outputFilename)
        return

    # Filler for data struct of maker
    def jet_filler(struct, jet_ideal, jet_real):
        for var in jetComponentNames:
            setattr(struct, var+'_ideal', jet_ideal[var])
            setattr(struct, var+'_real', jet_real[var])

    # Matching criterion for jets
    def match(j1, j2):
        return deltaR2(j1,j2)<0.2**2 

    # readers for the two samples
    reader_ideal = sample_ideal.treeReader( variables = variables , selectionString = "1")
    reader_real  = sample_real.treeReader( variables = variables , selectionString = "1")

    # Maker
    jetTreeMaker  =    TreeMaker( filler = None, variables = jet_variables, treeName = "jets") 

    # Getting started
    reader_ideal.start()
    reader_real.start()
    jetTreeMaker.start()

    # Looping over common positions
    for i_pos, pos in enumerate(positions):
        if i_pos%10000 == 0: logger.info( "At position %6i/%6i.", i_pos, len(positions) )
        reader_ideal.goToPosition(pos[0])
        reader_real.goToPosition(pos[1])

        # Get jets of event
        jets_ideal = [{c:getattr(reader_ideal.data, "Jet_"+c)[i] for c in jetComponentNames} for i in xrange(reader_ideal.data.nJet)]
        jets_real =  [{c:getattr(reader_real.data, "Jet_"+c)[i] for c in jetComponentNames} for i in xrange(reader_real.data.nJet)]

        # find matched pairs    
        pairs = [(j1, j2) for j1 in jets_ideal for j2 in jets_real if match(j1, j2)]
        
        # Fill for each jet
        for jet_ideal, jet_real in pairs:
            jet_filler(jetTreeMaker.data, jet_ideal, jet_real)
            jetTreeMaker.run()

    # Write tree from memory to file
    if not os.path.exists(os.path.dirname(outputFilename)):
        os.makedirs(os.path.dirname(outputFilename))

    outputFile = ROOT.TFile(outputFilename, "recreate")
    jetTreeMaker.tree.Write()
    outputFile.Close()
    logger.info( "Written file %s.", outputFilename)
    # Clean up (delete TTree)
    jetTreeMaker.clear()

    return

# Threading
nPos = 100000
jobs = [ (i, positions[i:i+nPos]) for i in xrange(0, len(positions), nPos)]

from multiprocessing import Pool
pool = Pool(processes=10)

results = pool.map(wrapper, jobs)
pool.close()
pool.join()
