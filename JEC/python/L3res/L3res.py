#!/usr/bin/env python
''' Analysis script to remake Claudia's PU plots
'''
#
# Standard imports and batch mode
#
import ROOT
ROOT.gROOT.SetBatch(True)
import itertools
import os

from math                                import sqrt, cos, sin, pi
from RootTools.core.standard             import *
from JetMET.tools.user                   import plot_directory
from JetMET.tools.helpers                import deltaPhi, deltaR

# Object selection
from JetMET.tools.objectSelection        import getFilterCut, getJets, jetVars
#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?', )
argParser.add_argument('--plot_directory',     action='store',      default='JEC/L3res')
args = argParser.parse_args()

#
# Logger
#
import JetMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)


## JEC on the fly
from JetMET.JetCorrector.jetCorrectors_Summer16 import jetCorrector_data, jetCorrector_mc
#
## pT_corr = pT_raw*L1(pT_raw)*L2L3(pT_raw*L1(pT_raw))*L2L3Res(pT_raw*L1(pT_raw)*L2L3(pT_raw*L1(pT_raw)))
jetCorrector_L1L2L3MC      = jetCorrector_mc.reduceLevels(correctionLevels   = ['L1FastJet', 'L2Relative', 'L3Absolute'] ) 
jetCorrector_L1L2L3Data    = jetCorrector_data.reduceLevels(correctionLevels = ['L1FastJet', 'L2Relative', 'L3Absolute'] ) 

if args.small: args.plot_directory += "_small"
#
# Make samples, will be searched for in the postProcessing directory
#
data_directory = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
postProcessing_directory = "postProcessed_80X_v36/dilepTiny/"
from JetMET.JEC.samples.cmgTuples_Summer16_mAODv2_postProcessed import *
data_directory = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
postProcessing_directory = "postProcessed_80X_v36/dilepTiny/"
from JetMET.JEC.samples.cmgTuples_Data25ns_80X_03Feb_postProcessed import *

selection       = 'ptll10-btb-njet1p'
selectionString = 'dl_pt>10&&cos(dl_phi-JetGood_phi[0])<-0.5&&nJetGood>=1'

#
# Text on the plots
#
def drawObjects( dataMCScale, lumi_scale ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Preliminary'), 
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) )
    ]
    return [tex.DrawLatex(*l) for l in lines] 


# Formatting for 1D plots
def draw1DPlots(plots, mode, dataMCScale):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, args.plot_directory, mode + ("_log" if log else ""), selection)
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
      
      plotting.draw(plot,
        plot_directory = plot_directory_,
        ratio = {'yRange':(0.6,1.4)},
        logX = False, logY = log, sorting = True,
        yRange = (0.03, "auto") if log else (0.001, "auto"),
        scaling = {},
        legend = (0.50,0.88-0.04*sum(map(len, plot.histos)),0.9,0.88),
        drawObjects = drawObjects( dataMCScale , lumi_scale )
      )

#Formatting for 1D profiles
def draw1DProfiles(plots, mode, dataMCScale):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, args.plot_directory, mode + ("_log" if log else ""), selection)
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot

      plotting.draw(plot,
        plot_directory = plot_directory_,
        ratio = {'yRange':(0.8,1.2)},
        logX = True, logY = log, sorting = True,
        yRange = (0.03, "auto") if log else (0.61, 1.41),
        scaling = {},
        legend = (0.50,0.88-0.04*sum(map(len, plot.histos)),0.9,0.88),
        drawObjects = drawObjects( dataMCScale , lumi_scale )
      )
#
# Read variables and sequences
#
read_variables = ["run/I", "weight/F", "l1_eta/F" , "l1_phi/F", "l2_eta/F", "l2_phi/F", 
                  "JetGood[pt/F,eta/F,phi/F,area/F,btagCSV/F,rawPt/F]", 
                  "dl_mass/F", "dl_eta/F", "dl_pt/F","dl_phi/F",
                  "met_chsPt/F", "met_chsPhi/F", "metSig/F", "ht/F", "nBTag/I", "nJetGood/I", 'nVert/I']
sequence = []

# extra lepton stuff
read_variables += [
#     "nLepGood/I", "LepGood[pt/F,eta/F,etaSc/F,phi/F,pdgId/I,tightId/I,miniRelIso/F,relIso03/F,sip3d/F,mediumMuonId/I,mvaIdSpring15/F,lostHits/I,convVeto/I,dxy/F,dz/F,jetPtRelv2/F,jetPtRatiov2/F,eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz/I]",
#     "nLepOther/I", "LepOther[pt/F,eta/F,etaSc/F,phi/F,pdgId/I,tightId/I,miniRelIso/F,relIso03/F,sip3d/F,mediumMuonId/I,mvaIdSpring15/F,lostHits/I,convVeto/I,dxy/F,dz/F,jetPtRelv2/F,jetPtRatiov2/F,eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz/I]",
    ]

jetMCBranches = [ "mcMatchFlav/I", "partonId/I", "mcFlavour/I", "partonFlavour/I", "hadronFlavour/I", "mcMatchId/I", "mcPt/F" ]
jetMCVars = [ s.split('/')[0] for s in jetMCBranches]

# mcMatchFlav   Flavour of associated parton from hard scatter (if any) for Jets failing id after jet-lepton cleaning 
# partonId  parton flavour (manually matching to status 23 particles) for Jets failing id after jet-lepton cleaning 
# mcFlavour parton flavour (physics definition, i.e. including b's from shower) for Jets failing id after jet-lepton cleaning 
# partonFlavour purely parton-based flavour for Jets failing id after jet-lepton cleaning 
# hadronFlavour hadron flavour (ghost matching to B/C hadrons) for Jets failing id after jet-lepton cleaning 
# mcMatchId Match to source from hard scatter (pdgId of heaviest particle in chain, 25 for H, 6 for t, 23/24 for W/Z), zero if non-prompt or fake for Jets failing id after jet-lepton cleaning 

jetVars += ['rawPt']

nan_jet = {key:float('nan') for key in jetVars}

def makeL3ResObservables( event, sample ):
    good_jets = getJets( event, jetColl="JetGood", jetVars = jetVars)

    # leading jet
    event.leading_jet =  good_jets[0]
    # subleading jet
    event.subleading_jet = good_jets[1] if len(good_jets)>=2 else nan_jet
    event.alpha = event.subleading_jet['rawPt'] / event.dl_pt

    event.r_ptbal = event.leading_jet['rawPt'] / event.dl_pt
    event.r_mpf   = 1. + event.met_chsPt * cos(event.met_chsPhi - event.dl_phi) / event.dl_pt
    

sequence.append( makeL3ResObservables )

log_pt_thresholds = [10**(x/10.) for x in range(11,36)]

z_window = 10
def getLeptonSelection( mode ):
  if   mode=="mumu": return "nGoodMuons==2&&nGoodElectrons==0&&isOS&&isMuMu&&abs(dl_mass-91.2)<%f" % z_window
  elif mode=="ee":   return "nGoodMuons==0&&nGoodElectrons==2&&isOS&&isEE&&abs(dl_mass-91.2)<%f" % z_window

weight_mc   = lambda event, sample: event.weight*event.reweightLeptonSF*event.reweightDilepTriggerBackup*event.reweightPU36fb
weight_data = lambda event, sample: event.weight

#
# Loop over channels
#
yields     = {}
allPlots   = {}
allModes   = ['mumu']#,'ee']
for index, mode in enumerate(allModes):
  yields[mode] = {}
  if   mode=="mumu": data = DoubleMuon_Run2016_backup
  elif mode=="ee":   data = DoubleEG_Run2016_backup
  if   mode=="mumu": data.texName = "data (2 #mu)"
  elif mode=="ee":   data.texName = "data (2 e)"

  data.setSelectionString([getFilterCut(isData=True), getLeptonSelection(mode)])
  data.name           = "data"
  data.style          = styles.errorStyle(ROOT.kBlack)
  data.weight         = weight_data 

  lumi_scale          = data.lumi/1000

  #mc             = [DY_HT_LO] + [ Top_pow, TTZ_LO, TTXNoZ, multiBoson]
  #DY_sample = DY
  DY_sample = DY_HT_LO
  TTJets_sample = Top
  other_mc_samples = [TTZ_LO, TTXNoZ, multiBoson]
  other_mc  = Sample.combine( name = "other_mc", texName = "VV/VVV/TTX/tZq/tWZ", samples = other_mc_samples, color = ROOT.kMagenta )
  mc        = [DY_sample, TTJets_sample, other_mc]
  all_mc   = Sample.combine( name = "all_mc",   texName = "simulation", samples = [DY_sample, Top_pow, TTZ_LO, TTXNoZ, multiBoson], color = ROOT.kBlue )
  all_mc.style = styles.lineStyle( all_mc.color, errors = True )    

  for sample in mc: sample.style = styles.fillStyle(sample.color)

  for sample in mc:
    sample.scale          = lumi_scale
    sample.read_variables = ['reweightDilepTriggerBackup/F','reweightLeptonSF/F','reweightPU36fb/F', 'nTrueInt/F']
    #sample.read_variables +=["JetGood[%s]"%( ",".join( jetMCBranches ) )]
    #sample.weight         = lambda event, sample: event.reweightLeptonSF*event.reweightLeptonHIPSF*event.reweightDilepTriggerBackup*nTrueInt27fb_puRW(event.nTrueInt)*event.reweightBTag_SF
    sample.weight         = weight_mc 
    sample.setSelectionString([getFilterCut(isData=False),  getLeptonSelection(mode)])

  stack          = Stack(mc,  data )
  stack_profile  = Stack([all_mc], data )

  if args.small:
        for sample in stack.samples + stack_profile.samples:
            sample.reduceFiles( to = 1 )

  # Use some defaults
  Plot.setDefaults(stack = stack, weight = None,   selectionString = selectionString, addOverFlowBin='upper')
  
  plots                 = []
  profiles1D = []

  plots.append(Plot(
    name = 'yield', texX = 'yield', texY = 'Number of Events',
    attribute = lambda event, sample: 0.5 + index,
    binning=[3, 0, 3],
  ))

  plots.append(Plot(
    texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events / 30 GeV',
    name = 'jet1_pt', attribute = lambda event, sample: event.JetGood_pt[0],
    binning=[600/30,0,600],
  ))

  plots.append( 
    Plot(
    texX = '#eta(leading jet) (GeV)', texY = 'Number of Events',
    name = 'jet1_eta', attribute = lambda event, sample: event.JetGood_eta[0],
    binning=[104,-5.2,5.2],
  ))

  plots.append(Plot(
    texX = 'p_{T}(subleading jet) (GeV)', texY = 'Number of Events / 30 GeV',
    name = 'jet1_pt', attribute = lambda event, sample: event.JetGood_pt[1],
    binning=[600/30,0,600],
  ))

  plots.append( 
    Plot(
    texX = '#eta(subleading jet) (GeV)', texY = 'Number of Events',
    name = 'jet1_eta', attribute = lambda event, sample: event.JetGood_eta[1],
    binning=[104,-5.2,5.2],
  ))

  plots.append( 
    Plot(
    texX = '#eta(subleading jet) (GeV)', texY = 'Number of Events',
    name = 'jet1_eta', attribute = lambda event, sample: event.JetGood_eta[1],
    binning=[104,-5.2,5.2],
  ))

  plots.append( 
    Plot(
    texX = '#alpha', texY = 'Number of Events',
    name = 'alpha', attribute = lambda event, sample: event.alpha,
    binning=[50,0,1],
  ))

  plots.append( 
    Plot(
    texX = '#Delta R(j_{1}, j_{2})', texY = 'Number of Events',
    name = 'dRj1j2', attribute = lambda event, sample: deltaR( event.leading_jet, event.subleading_jet),
    binning=[60,0,6],
  ))

  plots.append(Plot(
    texX = 'number of jets', texY = 'Number of Events',
    attribute = TreeVariable.fromString('nJetGood/I'),
    binning=[14,0,14],
  ))

  plots.append(Plot(
    texX = 'number of medium b-tags (CSVM)', texY = 'Number of Events',
    attribute = TreeVariable.fromString('nBTag/I'),
    binning=[8,0,8],
  ))

  plots.append(Plot(
    texX = 'H_{T} (GeV)', texY = 'Number of Events / 25 GeV',
    attribute = TreeVariable.fromString( "ht/F" ),
    binning=[500/25,0,600],
  ))

  plots.append(Plot(
    texX = 'm(ll) of leading dilepton (GeV)', texY = 'Number of Events / 4 GeV',
    attribute = TreeVariable.fromString( "dl_mass/F" ),
    binning=[200/4,0,200],
  ))

  plots.append(Plot(
    texX = 'p_{T}(ll) (GeV)', texY = 'Number of Events / 10 GeV',
    attribute = TreeVariable.fromString( "dl_pt/F" ),
    binning=[20,0,400],
  ))

  plots.append(Plot(
      texX = '#eta(ll) ', texY = 'Number of Events',
      name = 'dl_eta', attribute = lambda event, sample: event.dl_eta, 
      read_variables = ['dl_eta/F'],
      binning=[52,-5.2,5.2],
  ))

  plots.append(Plot(
    texX = '#phi(ll)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "dl_phi/F" ),
    binning=[10,-pi,pi],
  ))

  plots.append(Plot(
      texX = 'chs-E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
      attribute = TreeVariable.fromString( "met_chsPt/F" ),
      binning=[400/20,0,400],
  ))

  plots.append(Plot(
      name = "R_ptbal",
      texX = 'R_{pt-bal}', texY = 'Number of Events',
      attribute = lambda event, sample: event.r_ptbal,
      binning=[100,0,2],
  ))
  plots.append(Plot(
      name = "R_mpf",
      texX = 'R_{MPF}', texY = 'Number of Events',
      attribute = lambda event, sample: event.r_mpf,
      binning=[100,0,2],
  ))

  profiles1D.append(Plot(
    name = 'r_ptbal_profile_eta', texX = '#eta (jet)', texY = 'R_{pt-bal}',
    histo_class = ROOT.TProfile,
    stack = stack_profile,
    attribute = (
        lambda event, sample: event.dl_pt,
        lambda event, sample: event.r_ptbal,
    ),
    binning = Binning.fromThresholds(log_pt_thresholds),
  ))

  profiles1D.append(Plot(
    name = 'r_mpf_profile_eta', texX = '#eta (jet)', texY = 'R_{MPF}',
    histo_class = ROOT.TProfile,
    stack = stack_profile,
    attribute = (
        lambda event, sample: event.dl_pt,
        lambda event, sample: event.r_mpf,
    ),
    binning = Binning.fromThresholds(log_pt_thresholds),
  ))

  plotting.fill( plots + profiles1D , read_variables = read_variables, sequence = sequence )

  # Get normalization yields from yield histogram
  for plot in plots:
    if plot.name == "yield":
      for i, l in enumerate(plot.histos):
        for j, h in enumerate(l):
          yields[mode][plot.stack[i][j].name] = h.GetBinContent( h.FindBin( 0.5 + index ) )
          h.GetXaxis().SetBinLabel(1, "#mu#mu")
          h.GetXaxis().SetBinLabel(2, "e#mu")
          h.GetXaxis().SetBinLabel(3, "ee")

  yields[mode]["MC"] = sum(yields[mode][s.name] for s in mc)
  dataMCScale        = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')

  draw1DPlots( plots, mode, dataMCScale )
  draw1DProfiles( profiles1D, mode, dataMCScale )
