#!/usr/bin/env python
''' Analysis script for L1 closure plots
'''
#
# Standard imports and batch mode
#
import ROOT, os
ROOT.gROOT.SetBatch(True)
import itertools

from math                                import sqrt, cos, sin, pi
from RootTools.core.standard             import *
from JetMET.tools.user                   import plot_directory
from JetMET.tools.helpers                import deltaPhi
from JetMET.tools.helpers                import deltaR

# Object selection
from JetMET.tools.objectSelection        import getFilterCut, getJets
#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?', )
argParser.add_argument('--plot_directory',     action='store',      default='L1closure')
args = argParser.parse_args()

#
# Logger
#
import JetMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)


# JEC on the fly
from JetMET.JetCorrector.jetCorrectors_Spring16 import jetCorrector_data, jetCorrector_mc

# pT_corr = pT_raw*L1(pT_raw)*L2L3(pT_raw*L1(pT_raw))*L2L3Res(pT_raw*L1(pT_raw)*L2L3(pT_raw*L1(pT_raw)))
jetCorrector_L1MC          = jetCorrector_mc.fromLevels(correctionLevels   = ['L1FastJet'] )
jetCorrector_L1Data        = jetCorrector_data.fromLevels(correctionLevels = ['L1FastJet'] )
jetCorrector_L1L2L3MC      = jetCorrector_mc.fromLevels(correctionLevels   = ['L1FastJet', 'L2Relative', 'L3Absolute'] ) 
jetCorrector_L1L2L3Data    = jetCorrector_data.fromLevels(correctionLevels = ['L1FastJet', 'L2Relative', 'L3Absolute'] ) 
jetCorrector_L1L2L3ResData = jetCorrector_data.fromLevels(correctionLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'] )

if args.small:                        args.plot_directory += "_small"
#
# Make samples, will be searched for in the postProcessing directory
#
data_directory = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
postProcessing_directory = "postProcessed_80X_v26/dilep/"
from JetMET.JEC.samples.cmgTuples_Spring16_mAODv2_postProcessed import *
data_directory = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
postProcessing_directory = "postProcessed_80X_v26/dilep"
from JetMET.JEC.samples.cmgTuples_Data25ns_80X_23Sep_postProcessed import *

selection       = 'ptll30-btb-njet1p-nbtag0'
selectionString = 'dl_pt>30&&cos(dl_phi-JetGood_phi[0])<-0.5&&nJetGood>=1&&nBTag==0'

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

def drawPlots(plots, mode, dataMCScale):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, args.plot_directory, mode + ("_log" if log else ""), selection)
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot

      plotting.draw(plot,
	    plot_directory = plot_directory_,
	    ratio = {'yRange':(0.1,1.9)},
	    logX = False, logY = log, sorting = True,
	    yRange = (0.03, "auto") if log else (0.001, "auto"),
	    scaling = {},
	    legend = (0.50,0.88-0.04*sum(map(len, plot.histos)),0.9,0.88),
	    drawObjects = drawObjects( dataMCScale , lumi_scale )
      )

#
# Read variables and sequences
#
read_variables = ["weight/F", "l1_eta/F" , "l1_phi/F", "l2_eta/F", "l2_phi/F", "JetGood[pt/F,eta/F,phi/F,area/F,btagCSV/F,rawPt/F]", "dl_mass/F", "dl_eta/F", "dl_mt2ll/F", "dl_mt2bb/F", "dl_mt2blbl/F",
                  "met_pt/F", "met_phi/F", "metSig/F", "ht/F", "nBTag/I", "nJetGood/I", "rho/F"]
sequence = []

# extra lepton stuff
read_variables += [
     "nLepGood/I", "LepGood[pt/F,eta/F,etaSc/F,phi/F,pdgId/I,tightId/I,miniRelIso/F,relIso03/F,sip3d/F,mediumMuonId/I,mvaIdSpring15/F,lostHits/I,convVeto/I,dxy/F,dz/F,jetPtRelv2/F,jetPtRatiov2/F,eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz/I]",
     "nLepOther/I", "LepOther[pt/F,eta/F,etaSc/F,phi/F,pdgId/I,tightId/I,miniRelIso/F,relIso03/F,sip3d/F,mediumMuonId/I,mvaIdSpring15/F,lostHits/I,convVeto/I,dxy/F,dz/F,jetPtRelv2/F,jetPtRatiov2/F,eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz/I]",
    ]

def makeJetBalancing( event, sample ):
    good_jets = getJets( event, jetColl="JetGood")

    # leading jet
    rawPt = good_jets[0]['rawPt']
    eta   = good_jets[0]['eta']
    event.area  = good_jets[0]['area']

    # compute correction factors
    if sample.isData:
        corr    = jetCorrector_L1L2L3ResData.correction(rawPt, eta, event.area, event.rho) 
        corr_L1 = jetCorrector_L1Data.correction(rawPt, eta, event.area, event.rho)
    else:  
        corr    = jetCorrector_L1L2L3MC.correction(rawPt, eta, event.area, event.rho) 
        corr_L1 = jetCorrector_L1MC.correction(rawPt, eta, event.area, event.rho)

    event.deltaPU = rawPt*corr*(1-1./corr_L1)
    event.deltaPUPerArea    =   event.deltaPU/event.area

    if event.rho>0:
        event.deltaPUPerRho     =   event.deltaPU/event.rho 
        event.deltaPUPerAreaRho =   event.deltaPU/(event.rho*event.area)
    else:
        event.deltaPUPerRho     =  float('nan') 
        event.deltaPUPerAreaRho =  float('nan') 

    # balance raw jet
    event.r_ptbal_raw  = rawPt / event.dl_pt
    # balance jet corrected for everything EXCEPT L1, however, L2L3+L2L3res are evaluated at the L1 corrected pT
    event.r_ptbal_noL1  = ( corr / corr_L1 ) / event.dl_pt
    # balance L1 corrected jet 
    event.r_ptbal_L1   = rawPt*corr_L1 / event.dl_pt
    # balance fully corrected jet
    event.r_ptbal_corr = corr / event.dl_pt

    return

sequence.append( makeJetBalancing )

z_window = 7
def getLeptonSelection( mode ):
  if   mode=="mumu": return "nGoodMuons==2&&nGoodElectrons==0&&isOS&&isMuMu&&abs(dl_mass-91.2)<%f" % z_window
  elif mode=="ee":   return "nGoodMuons==0&&nGoodElectrons==2&&isOS&&isEE&&abs(dl_mass-91.2)<%f" % z_window

#
# Loop over channels
#
yields     = {}
allPlots   = {}
allModes   = ['mumu']#,'ee']
for index, mode in enumerate(allModes):
  yields[mode] = {}
  if   mode=="mumu": data_sample = DoubleMuon_Run2016_backup
  elif mode=="ee":   data_sample = DoubleEG_Run2016_backup
  if   mode=="mumu": data_sample.texName = "data (2 #mu)"
  elif mode=="ee":   data_sample.texName = "data (2 e)"

  data_sample.setSelectionString([getFilterCut(isData=True), getLeptonSelection(mode)])
  data_sample.name           = "data"
  data_sample.read_variables = ["evt/I","run/I"]
  data_sample.style          = styles.errorStyle(ROOT.kBlack)
  lumi_scale                 = data_sample.lumi/1000

  weight_ = lambda event, sample: event.weight

  multiBosonList = [multiBoson]
  #mc             = [DY_HT_LO] + [ Top_pow, TTZ_LO, TTXNoZ] + multiBosonList
  mc             = [DY] + [ Top_pow, TTZ_LO, TTXNoZ] + multiBosonList 

  for sample in mc: sample.style = styles.fillStyle(sample.color)

  for sample in mc:
    sample.scale          = lumi_scale
    sample.read_variables = ['reweightLeptonHIPSF/F','reweightDilepTriggerBackup/F','reweightLeptonSF/F','reweightPU36fb/F', 'nTrueInt/F']
   #sample.weight         = lambda event, sample: event.reweightLeptonSF*event.reweightLeptonHIPSF*event.reweightDilepTriggerBackup*nTrueInt27fb_puRW(event.nTrueInt)*event.reweightBTag_SF
    sample.weight         = lambda event, sample: event.reweightLeptonSF*event.reweightDilepTriggerBackup*event.reweightPU36fb
    sample.setSelectionString([getFilterCut(isData=False),  getLeptonSelection(mode)])

  stack = Stack(mc, data_sample)

  if args.small:
        for sample in stack.samples:
            sample.reduceFiles( to = 1 )

  # Use some defaults
  Plot.setDefaults(stack = stack, weight = weight_, selectionString = selectionString, addOverFlowBin='upper')
  
  plots = []

  plots.append(Plot(
    name = 'yield', texX = 'yield', texY = 'Number of Events',
    attribute = lambda event, sample: 0.5 + index,
    binning=[3, 0, 3],
  ))

  plots.append(Plot(
    name = 'nVtxs', texX = 'vertex multiplicity', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "nVert/I" ),
    binning=[50,0,50],
  ))

  plots.append(Plot(
    name = 'rho', texX = 'rho', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "rho/F" ),
    binning=[50,0,50],
  ))

  plots.append(Plot(
    name = 'leading_jet_area', texX = 'area', texY = 'Number of Events',
    attribute = lambda event, sample: event.area,
    binning=[50,0,1],
  ))

  plots.append(Plot(
    name = 'leading_jet_radius', texX = 'area', texY = 'Number of Events',
    attribute = lambda event, sample: sqrt(event.area/pi),
    binning=[50,0,1],
  ))

  plots.append(Plot(
      texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
      attribute = TreeVariable.fromString( "met_pt/F" ),
      binning=[400/20,0,400],
  ))

  plots.append(Plot(
      texX = '#phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
      attribute = TreeVariable.fromString( "met_phi/F" ),
      binning=[10,-pi,pi],
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
    binning=[40, 91.2 - z_window, 91.2 + z_window],
  ))

  plots.append(Plot(
    texX = 'p_{T}(ll) (GeV)', texY = 'Number of Events / 10 GeV',
    attribute = TreeVariable.fromString( "dl_pt/F" ),
    binning=[20,0,400],
  ))

  plots.append(Plot(
      texX = '#eta(ll) ', texY = 'Number of Events',
      name = 'dl_eta', attribute = lambda event, sample: abs(event.dl_eta), read_variables = ['dl_eta/F'],
      binning=[10,0,3],
  ))

  plots.append(Plot(
    texX = '#phi(ll)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "dl_phi/F" ),
    binning=[10,-pi,pi],
  ))

  plots.append(Plot(
    texX = 'Cos(#Delta#phi(ll, E_{T}^{miss}))', texY = 'Number of Events',
    name = 'cosZMetphi',
    attribute = lambda event, sample: cos( event.dl_phi - event.met_phi ), 
    read_variables = ["met_phi/F", "dl_phi/F"],
    binning = [40,-1,1],
  ))

  plots.append(Plot(
    texX = 'p_{T}(l_{1}) (GeV)', texY = 'Number of Events / 5 GeV',
    attribute = TreeVariable.fromString( "l1_pt/F" ),
    binning=[20,0,300],
  ))

  plots.append(Plot(
    texX = '#eta(l_{1})', texY = 'Number of Events',
    name = 'l1_eta', attribute = lambda event, sample: abs(event.l1_eta), read_variables = ['l1_eta/F'],
    binning=[15,0,3],
  ))

  plots.append(Plot(
    texX = '#phi(l_{1})', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "l1_phi/F" ),
    binning=[10,-pi,pi],
  ))

  plots.append(Plot(
    texX = 'p_{T}(l_{2}) (GeV)', texY = 'Number of Events / 5 GeV',
    attribute = TreeVariable.fromString( "l2_pt/F" ),
    binning=[20,0,300],
  ))

  plots.append(Plot(
    texX = '#eta(l_{2})', texY = 'Number of Events',
    name = 'l2_eta', attribute = lambda event, sample: abs(event.l2_eta), read_variables = ['l2_eta/F'],
    binning=[15,0,3],
  ))

  plots.append(Plot(
    texX = '#phi(l_{2})', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "l2_phi/F" ),
    binning=[10,-pi,pi],
  ))

  plots.append(Plot(
    texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events / 30 GeV',
    name = 'jet1_pt', attribute = lambda event, sample: event.JetGood_pt[0],
    binning=[600/30,0,600],
  ))

  plots.append(Plot(
    texX = '#eta(leading jet) (GeV)', texY = 'Number of Events',
    name = 'jet1_eta', attribute = lambda event, sample: abs(event.JetGood_eta[0]),
    binning=[10,0,3],
  ))

  plots.append(Plot(
    texX = '#phi(leading jet) (GeV)', texY = 'Number of Events',
    name = 'jet1_phi', attribute = lambda event, sample: event.JetGood_phi[0],
    binning=[10,-pi,pi],
  ))


  plots.append(Plot(
    name = 'cosZJet1phi',
    texX = 'Cos(#Delta#phi(Z, leading jet))', texY = 'Number of Events',
    attribute = lambda event, sample: cos( event.dl_phi - event.JetGood_phi[0] ) ,
    read_variables =  ["dl_phi/F", "JetGood[phi/F]"],
    binning = [10,-1,1],
  ))

  plots.append(Plot(
    texX = 'p_{T}(2nd leading jet) (GeV)', texY = 'Number of Events / 30 GeV',
    name = 'jet2_pt', attribute = lambda event, sample: event.JetGood_pt[1],
    binning=[600/30,0,600],
  ))

  plots.append(Plot(
    texX = '#eta(2nd leading jet) (GeV)', texY = 'Number of Events',
    name = 'jet2_eta', attribute = lambda event, sample: abs(event.JetGood_eta[1]),
    binning=[10,0,3],
  ))

  plots.append(Plot(
    texX = '#phi(2nd leading jet) (GeV)', texY = 'Number of Events',
    name = 'jet2_phi', attribute = lambda event, sample: event.JetGood_phi[1],
    binning=[10,-pi,pi],
  ))

  plots.append(Plot(
    name = 'cosMetJet2phi',
    texX = 'Cos(#Delta#phi(E_{T}^{miss}, second jet))', texY = 'Number of Events',
    attribute = lambda event, sample: cos( event.met_phi - event.JetGood_phi[1] ) , 
    read_variables = ["met_phi/F", "JetGood[phi/F]"],
    binning = [10,-1,1],
  ))
    
  plots.append(Plot(
    name = 'cosMetJet2phi_smallBinning',
    texX = 'Cos(#Delta#phi(E_{T}^{miss}, second jet))', texY = 'Number of Events',
    attribute = lambda event, sample: cos( event.met_phi - event.JetGood_phi[1] ) , 
    read_variables = ["met_phi/F", "JetGood[phi/F]"],
    binning = [20,-1,1],
  ))

  plots.append(Plot(
    name = 'cosZJet2phi',
    texX = 'Cos(#Delta#phi(Z, 2nd leading jet))', texY = 'Number of Events',
    attribute = lambda event, sample: cos( event.dl_phi - event.JetGood_phi[0] ),
    read_variables = ["dl_phi/F", "JetGood[phi/F]"],
    binning = [10,-1,1],
  ))

  plots.append(Plot(
    name = 'cosJet1Jet2phi',
    texX = 'Cos(#Delta#phi(leading jet, 2nd leading jet))', texY = 'Number of Events',
    attribute = lambda event, sample: cos( event.JetGood_phi[1] - event.JetGood_phi[0] ) ,
    read_variables =  ["JetGood[phi/F]"],
    binning = [10,-1,1],
  ))

  plotting.fill(plots, read_variables = read_variables, sequence = sequence)

  # Get normalization yields from yield histogram
  for plot in plots:
    if plot.name == "yield":
      for i, l in enumerate(plot.histos):
        for j, h in enumerate(l):
          yields[mode][plot.stack[i][j].name] = h.GetBinContent(h.FindBin(0.5+index))
          h.GetXaxis().SetBinLabel(1, "#mu#mu")
          h.GetXaxis().SetBinLabel(2, "e#mu")
          h.GetXaxis().SetBinLabel(3, "ee")

  yields[mode]["MC"] = sum(yields[mode][s.name] for s in mc)
  dataMCScale        = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')

  drawPlots(plots, mode, dataMCScale)
