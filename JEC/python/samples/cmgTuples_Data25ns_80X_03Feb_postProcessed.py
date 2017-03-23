import copy, os, sys
from RootTools.core.Sample import Sample 
import ROOT

# Logging
import logging
logger = logging.getLogger(__name__)

# Data directory
try:    data_directory = sys.modules['__main__'].data_directory
except: from JetMET.tools.user import data_directory

# Take post processing directory if defined in main module
try:    postProcessing_directory = sys.modules['__main__'].postProcessing_directory
except: postProcessing_directory = 'postProcessed_80X_v36/dilepTiny'

logger.info("Loading data samples from directory %s", os.path.join(data_directory, postProcessing_directory))

dirs = {}
for (run, version) in [('B','_v2'),('C',''),('D',''),('E',''),('F',''),('G',''),('H','_v2'),('H','_v3')]:
  runTag = 'Run2016' + run + '_03Feb2017' + version
  dirs["DoubleEG_Run2016"   + run + version + "_backup"] = ["DoubleEG_"   + runTag + "_Trig_ee",   "SingleElectron_" + runTag + "_Trig_e_for_ee"]
  dirs["DoubleMuon_Run2016" + run + version + "_backup"] = ["DoubleMuon_" + runTag + "_Trig_mumu", "SingleMuon_"     + runTag + "_Trig_mu_for_mumu"]
  dirs["MuonEG_Run2016"     + run + version + "_backup"] = ["MuonEG_"     + runTag + "_Trig_mue",  "SingleElectron_" + runTag + "_Trig_e_for_mue", "SingleMuon_" + runTag + "_Trig_mu_for_mue"]

def merge(pd, totalRunName, listOfRuns):
  dirs[pd + '_' + totalRunName + '_backup'] = []
  for run in listOfRuns: dirs[pd + '_' + totalRunName + '_backup'].extend(dirs[pd + '_' + run + '_backup'])

for pd in ['DoubleEG','DoubleMuon','MuonEG']:
  merge(pd, 'Run2016BCD',    ['Run2016B_v2', 'Run2016C', 'Run2016D'])
  merge(pd, 'Run2016BCDEFG', ['Run2016BCD', 'Run2016E', 'Run2016F', 'Run2016G'])
  merge(pd, 'Run2016EF', ['Run2016E', 'Run2016F'])  #E Fearly
  merge(pd, 'Run2016GH', ['Run2016F', 'Run2016G', 'Run2016H_v2', 'Run2016H_v3']) #Flate+GH
  merge(pd, 'Run2016',       ['Run2016BCDEFG', 'Run2016H_v2', 'Run2016H_v3'])

for key in dirs:
  dirs[key] = [ os.path.join( data_directory, postProcessing_directory, dir) for dir in dirs[key]]


def getSample(pd, runName, lumi):
  sample      = Sample.fromDirectory(name=(pd + '_' + runName + '_backup'), treeName="Events", texName=(pd + ' (' + runName + ')'), directory=dirs[pd + '_' + runName + '_backup'])
  sample.lumi = lumi
  return sample

DoubleEG_Run2016BCD_backup      = getSample('DoubleEG',   'Run2016BCD',    (5.744+2.573+4.248)*1000)
DoubleMuon_Run2016BCD_backup    = getSample('DoubleMuon', 'Run2016BCD',    (5.744+2.573+4.248)*1000)  
MuonEG_Run2016BCD_backup        = getSample('MuonEG',     'Run2016BCD',    (5.743+2.573+4.248)*1000) 

DoubleEG_Run2016BCDEFG_backup   = getSample('DoubleEG',   'Run2016BCDEFG', (5.744+2.573+4.248+4.009+3.101+7.540)*1000)
DoubleMuon_Run2016BCDEFG_backup = getSample('DoubleMuon', 'Run2016BCDEFG', (5.744+2.573+4.248+4.009+3.101+7.540)*1000)
MuonEG_Run2016BCDEFG_backup     = getSample('MuonEG',     'Run2016BCDEFG', (5.743+2.573+4.248+4.009+3.101+7.540)*1000)

#Fearly 2.666, Flate 0.397 according to brilcalc
DoubleEG_Run2016EF_backup   = getSample('DoubleEG',   'Run2016EF', (4.009+3.101-0.397)*1000)
DoubleMuon_Run2016EF_backup = getSample('DoubleMuon', 'Run2016EF', (4.009+3.101-0.397)*1000)
MuonEG_Run2016EF_backup     = getSample('MuonEG',     'Run2016EF', (4.009+3.101-0.397)*1000)

DoubleEG_Run2016GH_backup   = getSample('DoubleEG',   'Run2016GH', (0.397+7.540+8.329+0.210)*1000)
DoubleMuon_Run2016GH_backup = getSample('DoubleMuon', 'Run2016GH', (0.397+7.540+8.329+0.210)*1000)
MuonEG_Run2016GH_backup     = getSample('MuonEG',     'Run2016GH', (0.397+7.540+8.329+0.210)*1000)

DoubleEG_Run2016_backup         = getSample('DoubleEG',   'Run2016',       (5.744+2.573+4.248+4.009+3.101+7.540+8.329+0.210)*1000)
DoubleMuon_Run2016_backup       = getSample('DoubleMuon', 'Run2016',       (5.744+2.573+4.248+4.009+3.101+7.540+8.329+0.210)*1000)
MuonEG_Run2016_backup           = getSample('MuonEG',     'Run2016',       (5.743+2.573+4.248+4.009+3.101+7.540+8.327+0.210)*1000)

allSamples_Data25ns = []
allSamples_Data25ns += [DoubleMuon_Run2016BCD_backup,    DoubleEG_Run2016BCD_backup,    MuonEG_Run2016BCD_backup]
allSamples_Data25ns += [DoubleMuon_Run2016BCDEFG_backup, DoubleEG_Run2016BCDEFG_backup, MuonEG_Run2016BCDEFG_backup]
allSamples_Data25ns += [DoubleMuon_Run2016_backup,       DoubleEG_Run2016_backup,       MuonEG_Run2016_backup]

for s in allSamples_Data25ns:
  s.color   = ROOT.kBlack
  s.isData  = True
