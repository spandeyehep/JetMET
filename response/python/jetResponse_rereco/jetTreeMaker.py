''' Make simple jet trees for rereco jet response comparison 
'''
# Standard imports
import sys
import os
import logging
import ROOT

#RootTools
from RootTools.core.standard import *

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

argParser.add_argument('--era', 
      action='store',
      #nargs=1,
      type=str,
      choices=['Run2016B', 'Run2016C', 'Run2016D', 'Run2016E', 'Run2016F', 'Run2016G'],
      default='Run2016G',
      help="era"
)

argParser.add_argument('--output_path', 
      action='store',
      #nargs=1,
      type=str,
      default='/scratch/rschoefbeck/cmgTuples/jetTuples/',
      help="era"
)

argParser.add_argument('--overwrite', 
      action='store_true',
      default=False,
)

argParser.add_argument('--run', 
      action='store',
      type = int,
      default=-1,
#      default=280187,
      help="run"
)

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

if args.era == 'Run2016B':
    s_prompt = 'JetHT_Run2016B_PromptReco_v2'
    s_rereco = 'JetHT_Run2016B_23Sep2016_v3'
elif args.era == 'Run2016C':
    s_prompt = 'JetHT_Run2016C_PromptReco_v2'
    s_rereco = 'JetHT_Run2016C_23Sep2016_v1'
elif args.era == 'Run2016D':
    s_prompt = 'JetHT_Run2016D_PromptReco_v2'
    s_rereco = 'JetHT_Run2016D_23Sep2016_v1'
elif args.era == 'Run2016E':
    s_prompt = 'JetHT_Run2016E_PromptReco_v2'
    s_rereco = 'JetHT_Run2016E_23Sep2016_v1'
elif args.era == 'Run2016F':
    s_prompt = 'JetHT_Run2016F_PromptReco_v1'
    s_rereco = 'JetHT_Run2016F_23Sep2016_v1'
elif args.era == 'Run2016G':
    s_prompt = 'JetHT_Run2016G_PromptReco_v1'
    s_rereco = 'JetHT_Run2016G_23Sep2016_v1'

from helpers import fromHeppySample
prompt = fromHeppySample( s_prompt, data_path = "/scratch/rschoefbeck/cmgTuples/80X_JetHT/" )
rereco = fromHeppySample( s_rereco, data_path = "/scratch/rschoefbeck/cmgTuples/80X_JetHT/" )

if args.run<0:
    from FWCore.PythonUtilities.LumiList import LumiList
    json = '$CMSSW_BASE/src/CMGTools/TTHAnalysis/data/json/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
    print "Find runs that are in datasets %s and %s and json %s"%(prompt.heppy.dataset, rereco.heppy.dataset, json)
    lumiList = LumiList(os.path.expandvars(json))

    def _dasPopen(dbs):
        if 'LSB_JOBID' in os.environ:
            raise RuntimeError, "Trying to do a DAS query while in a LXBatch job (env variable LSB_JOBID defined)\nquery was: %s" % dbs
        if 'X509_USER_PROXY' in os.environ:
            dbs += " --key {0} --cert {0}".format(os.environ['X509_USER_PROXY'])
        logger.info('DAS query\t: %s',  dbs)
        return os.popen(dbs)

    dbs='das_client --query="run dataset=%s instance=prod/%s" --limit %i'%(prompt.heppy.dataset, 'global', 0)
    prompt_runs = [int(r) for r in _dasPopen(dbs).readlines()]
    dbs='das_client --query="run dataset=%s instance=prod/%s" --limit %i'%(rereco.heppy.dataset, 'global', 0)
    rereco_runs =[int(r) for r in  _dasPopen(dbs).readlines()]

    runs = []
    for str_run in lumiList.getRuns():
        run = int(str_run)
        if run in prompt_runs and run in rereco_runs:
            runs.append(run)

    print "Now running %i jobs: %r"%( len(runs), runs )

    import subprocess
    def wrapper( run_ ):
        subprocess.call(["python", "jetTreeMaker.py", ("--era=%s"%args.era), ("--run=%i"%run_) ])

    from multiprocessing import Pool
    pool = Pool( 10 )
    results = pool.map(wrapper, runs)
    pool.close()

    sys.exit(0)

outputFilename = os.path.join( args.output_path, 'jet_%s_%i.root'%(args.era, args.run) )
if os.path.exists(outputFilename) and not args.overwrite:
    print "Found file %s. Skipping." % outputFilename
    sys.exit(0)

jet_str = "pt/F,eta/F,rawPt/F,phi/F,btagCSV/F,chHEF/F,neHEF/F,phEF/F,eEF/F,muEF/F,HFHEF/F,HFEMEF/F,chHMult/I,neHMult/I,phMult/I,eMult/I,muMult/I,HFHMult/I,HFEMMult/I,id/I"
read_variables = ["evt/l", "run/I", "lumi/I", "nJet/I", "nVert/I", "Jet[%s]"%jet_str]

r_prompt = prompt.treeReader( variables = map( TreeVariable.fromString, read_variables) , selectionString = "run==%i"%args.run)
r_rereco = rereco.treeReader( variables = map( TreeVariable.fromString, read_variables) , selectionString = "run==%i"%args.run)

from math import pi

maxN = -1 #20000

from helpers import getVarValue, getObjDict, deltaR2
jetVars = ['pt', 'eta', 'btagCSV', 'rawPt', 'phi', 'chHEF', 'neHEF', 'phEF', 'eEF', 'muEF', 'HFHEF', 'HFEMEF', 'chHMult', 'neHMult', 'phMult', 'eMult', 'muMult', 'HFHMult', 'HFEMMult', 'id']

def getJets(c, jetVars=jetVars, jetColl="Jet"):
    return [getObjDict(c, jetColl+'_', jetVars, i) for i in range(int(getVarValue(c, 'n'+jetColl)))]

r_prompt.start()
position_r_prompt = {}
count=0
while r_prompt.run():
    position_r_prompt[(r_prompt.event.run, r_prompt.event.lumi, r_prompt.event.evt) ] = r_prompt.position-1
    count+=1
    if maxN is not None and maxN>0 and count>=maxN:break

r_rereco.start()
position_r_rereco = {}
count=0
while r_rereco.run():
    position_r_rereco[(r_rereco.event.run, r_rereco.event.lumi, r_rereco.event.evt) ] = r_rereco.position-1
    count+=1
    if maxN is not None and maxN>0 and count>=maxN:break

logger.info( "Have %i events in first samle and %i in second", len(position_r_prompt), len(position_r_rereco) )

# Fast intersect
intersec = set(position_r_prompt.keys()).intersection(set(position_r_rereco.keys()))
positions = [(position_r_prompt[i], position_r_rereco[i]) for i in intersec]

if len(positions)==0:
    print "Found no common events in era %s and run %i. Quit." %( args.era, args.run )
    sys.exit(0)


# Without sorting, there is a jump between files with almost every event -> extremly slow
positions.sort()
logger.info("Have %i events in common.", len(intersec))

new_variables = [ "evt/l", "run/I", "lumi/I", "nVert/I" ] + jet_str.replace('/','_prompt/').split(',') +  jet_str.replace('/','_rereco/').split(',')

# Maker
jetTreeMaker  =    TreeMaker( sequence = [], variables = map( TreeVariable.fromString, new_variables ), treeName = "jets")

# Filler for data struct of maker
def jet_filler(struct, jet_prompt, jet_rereco):
    for var in jetVars:
        setattr(struct, var+'_prompt', jet_prompt[var])
        setattr(struct, var+'_rereco', jet_rereco[var])

#Looping over common events
jetTreeMaker.start()
for i, p in enumerate(positions):
    p1,p2 = p
    r_prompt.goToPosition(p1)
    r_rereco.goToPosition(p2)
    if i%10000==0: logger.info("At %i/%i of common events.", i, len(positions))

    jets_rereco = filter( lambda j:j['id']>0, getJets( r_rereco.event ) )
    jets_prompt = filter( lambda j:j['id']>0, getJets( r_prompt.event ) ) 

    counter = 0 
    for jet_prompt, jet_rereco in zip(jets_prompt, jets_rereco):
        if deltaR2( jet_prompt, jet_rereco )<0.2**2:
            jet_filler(jetTreeMaker.event, jet_prompt, jet_rereco)
            for attr in ['evt', 'run', 'lumi', 'nVert']:
                setattr(jetTreeMaker.event, attr,getattr(r_prompt.event, attr))
            jetTreeMaker.run()

if not os.path.exists( args.output_path ):
    os.makedirs( args.output_path )

outputFile = ROOT.TFile(outputFilename, "recreate")
jetTreeMaker.tree.Write()
outputFile.Close()
logger.info( "Written file %s", outputFilename)
