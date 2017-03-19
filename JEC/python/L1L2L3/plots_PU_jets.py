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
from JetMET.tools.helpers                import deltaPhi
from JetMET.tools.helpers                import deltaR

# Object selection
from JetMET.tools.objectSelection        import getFilterCut, getJets, jetVars
#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?', )
argParser.add_argument('--plot_directory',     action='store',      default='JEC/PUjets')
args = argParser.parse_args()

#
# Logger
#
import JetMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)


## JEC on the fly
#from JetMET.JetCorrector.jetCorrectors_Spring16 import jetCorrector_data, jetCorrector_mc
#
## pT_corr = pT_raw*L1(pT_raw)*L2L3(pT_raw*L1(pT_raw))*L2L3Res(pT_raw*L1(pT_raw)*L2L3(pT_raw*L1(pT_raw)))
#jetCorrector_L1MC          = jetCorrector_mc.reduceLevels(correctionLevels   = ['L1FastJet'] )
#jetCorrector_L1Data        = jetCorrector_data.reduceLevels(correctionLevels = ['L1FastJet'] )
#jetCorrector_L1L2L3MC      = jetCorrector_mc.reduceLevels(correctionLevels   = ['L1FastJet', 'L2Relative', 'L3Absolute'] ) 
#jetCorrector_L1L2L3Data    = jetCorrector_data.reduceLevels(correctionLevels = ['L1FastJet', 'L2Relative', 'L3Absolute'] ) 
#jetCorrector_L1L2L3ResData = jetCorrector_data.reduceLevels(correctionLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'] )

if args.small: args.plot_directory += "_small"
#
# Make samples, will be searched for in the postProcessing directory
#
data_directory = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
postProcessing_directory = "postProcessed_80X_v26/dilep/"
from JetMET.JEC.samples.cmgTuples_Spring16_mAODv2_postProcessed import *
data_directory = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
postProcessing_directory = "postProcessed_80X_v26/dilep"
from JetMET.JEC.samples.cmgTuples_Data25ns_80X_23Sep_postProcessed import *

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
# Formatting for 1D plots
def draw1DPlotsFlavourPairs(plot_pairs, mode, dataMCScale):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, args.plot_directory, mode + ("_log" if log else ""), selection)
    for plot_inc, plot_flavours in plot_pairs:
      if not max(l[0].GetMaximum() for l in plot_inc.histos): continue # Empty plot
      
      # Legend inclusive plot
      h_plot_inc = plot_inc.histos  
      for si, ss in enumerate( h_plot_inc ):
        for  sj in range(len(ss)):
            h_plot_inc[si][sj].legendText = plot_inc.stack[si][sj].texName

      # Legend and style per contribution
      h_plot_flavours = plot_flavours.histos
      for ih, h in enumerate(h_plot_flavours):
        h[0].style = jet_flavours[ih]['style']
        h[0].legendText = jet_flavours[ih]['name']
      
      plot_ = Plot.fromHisto( plot_inc.name, h_plot_inc + h_plot_flavours, texX = plot_inc.texX, texY = plot_inc.texY )
      plot_.stack = None
      plotting.draw( plot_,
	    plot_directory = plot_directory_,
	    ratio = {'yRange':(0.1,1.9)},
	    logX = False, logY = log, sorting = True,
	    yRange = (0.03, "auto") if log else (0.001, "auto"),
	    scaling = {},
	    legend = (0.5,0.97-0.035*sum(map(len, plot_.histos)),0.9,0.88),
	    drawObjects = drawObjects( dataMCScale , lumi_scale )
      )

#
# Read variables and sequences
#
read_variables = ["run/I", "weight/F", "l1_eta/F" , "l1_phi/F", "l2_eta/F", "l2_phi/F", 
                  "JetGood[pt/F,eta/F,phi/F,area/F,btagCSV/F,rawPt/F]", 
                  "dl_mass/F", "dl_eta/F", "dl_pt/F",
                  "met_pt/F", "met_phi/F", "metSig/F", "ht/F", "nBTag/I", "nJetGood/I", "rho/F", 'nVert/I']
sequence = []

# extra lepton stuff
read_variables += [
     "nLepGood/I", "LepGood[pt/F,eta/F,etaSc/F,phi/F,pdgId/I,tightId/I,miniRelIso/F,relIso03/F,sip3d/F,mediumMuonId/I,mvaIdSpring15/F,lostHits/I,convVeto/I,dxy/F,dz/F,jetPtRelv2/F,jetPtRatiov2/F,eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz/I]",
     "nLepOther/I", "LepOther[pt/F,eta/F,etaSc/F,phi/F,pdgId/I,tightId/I,miniRelIso/F,relIso03/F,sip3d/F,mediumMuonId/I,mvaIdSpring15/F,lostHits/I,convVeto/I,dxy/F,dz/F,jetPtRelv2/F,jetPtRatiov2/F,eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz/I]",
    ]

jetMCBranches = [ "mcMatchFlav/I", "partonId/I", "mcFlavour/I", "partonFlavour/I", "hadronFlavour/I", "mcMatchId/I", "mcPt/F" ]
jetMCVars = [ s.split('/')[0] for s in jetMCBranches]

# mcMatchFlav   Flavour of associated parton from hard scatter (if any) for Jets failing id after jet-lepton cleaning 
# partonId  parton flavour (manually matching to status 23 particles) for Jets failing id after jet-lepton cleaning 
# mcFlavour parton flavour (physics definition, i.e. including b's from shower) for Jets failing id after jet-lepton cleaning 
# partonFlavour purely parton-based flavour for Jets failing id after jet-lepton cleaning 
# hadronFlavour hadron flavour (ghost matching to B/C hadrons) for Jets failing id after jet-lepton cleaning 
# mcMatchId Match to source from hard scatter (pdgId of heaviest particle in chain, 25 for H, 6 for t, 23/24 for W/Z), zero if non-prompt or fake for Jets failing id after jet-lepton cleaning 

def jetHasMCMatch( jet ):
    return jet['mcPt'] > 10 #and abs(1-jet['pt']/jet['mcPt'])

def makeJetMatchingInfo( event, sample ):
    if sample.isData:
        good_jets = getJets( event, jetColl="JetGood", jetVars = jetVars)
    else:
        good_jets = getJets( event, jetColl="JetGood", jetVars = jetVars + jetMCVars)

    # leading jet
    event.rawPt = good_jets[0]['rawPt']
    event.eta   = good_jets[0]['eta']
    event.area  = good_jets[0]['area']

    if not sample.isData:
        event.isQuark = jetHasMCMatch( good_jets[0] ) and abs(good_jets[0]['partonFlavour'])>0 and abs(good_jets[0]['partonFlavour'])<6
        event.isGluon = jetHasMCMatch( good_jets[0] ) and good_jets[0]['partonFlavour'] == 21
        event.isPU    = good_jets[0]['mcPt']==0 and good_jets[0]['partonFlavour'] == 0
        event.isOther = not (event.isQuark or event.isGluon or event.isPU)
        #print sample.name, good_jets[0]['pt'], good_jets[0]['partonFlavour'], good_jets[0]['mcPt'], event.isQuark, event.isGluon, event.isPU, event.isOther
    else:
        event.isQuark = None 
        event.isGluon = None 
        event.isPU    = None 
        event.isOther = None 
    return

sequence.append( makeJetMatchingInfo )

z_window = 7
def getLeptonSelection( mode ):
  if   mode=="mumu": return "nGoodMuons==2&&nGoodElectrons==0&&isOS&&isMuMu&&abs(dl_mass-91.2)<%f" % z_window
  elif mode=="ee":   return "nGoodMuons==0&&nGoodElectrons==2&&isOS&&isEE&&abs(dl_mass-91.2)<%f" % z_window

weight_mc   = lambda event, sample: event.weight*event.reweightLeptonSF*event.reweightDilepTriggerBackup*event.reweightPU36fb
weight_data = lambda event, sample: event.weight

# weights for leading jet flavours
weight_isQuark  = lambda event, sample: event.isQuark
weight_isGluon  = lambda event, sample: event.isGluon
weight_isPU     = lambda event, sample: event.isPU
weight_isOther  = lambda event, sample: event.isOther

jet_flavours = [ \
    {'name':'quark'    , 'style':styles.lineStyle( ROOT.kBlack, width = 2)      , 'weight': weight_isQuark},
    {'name':'gluon'    , 'style':styles.lineStyle( ROOT.kRed, width = 2 )    , 'weight': weight_isGluon},
    {'name':'PU'       , 'style':styles.lineStyle( ROOT.kBlue, width = 2 )   , 'weight': weight_isPU},
    {'name':'other'    , 'style':styles.lineStyle( ROOT.kMagenta, width = 2 ) , 'weight': weight_isOther},
]

def firstJetEta( ptlow, pthigh ):
    return lambda event, sample: event.JetGood_eta[0] if event.JetGood_pt[0]>=ptlow and event.JetGood_pt[0]<pthigh else None
def firstJetPt ( absetalow, absetahigh ):
    return lambda event, sample: event.JetGood_pt[0] if abs(event.JetGood_eta[0])>=absetalow and abs(event.JetGood_eta[0])<absetahigh else None
def firstJetRawPt ( absetalow, absetahigh ):
    return lambda event, sample: event.JetGood_rawPt[0] if abs(event.JetGood_eta[0])>=absetalow and abs(event.JetGood_eta[0])<absetahigh else None

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
  TTJets_sample = TTJets
  other_mc_samples = [TTZ_LO, TTXNoZ, multiBoson]
  other_mc  = Sample.combine( name = "other_mc", texName = "VV/VVV/TTX/tZq/tWZ", samples = other_mc_samples, color = ROOT.kMagenta )
  all_mc    = Sample.combine( name = "all_mc",   texName = "simulation", samples = [DY_sample, TTJets] + other_mc_samples, color = ROOT.kBlue )
  mc        = [DY_sample, Top_pow, other_mc]

  for sample in mc: sample.style = styles.fillStyle(sample.color)
  all_mc.style = styles.lineStyle( all_mc.color, errors = True )    

  for sample in mc + [other_mc, all_mc]:
    sample.scale          = lumi_scale
    sample.read_variables = ['reweightLeptonHIPSF/F','reweightDilepTriggerBackup/F','reweightLeptonSF/F','reweightPU36fb/F', 'nTrueInt/F']
    sample.read_variables +=["JetGood[%s]"%( ",".join( jetMCBranches ) )]
   #sample.weight         = lambda event, sample: event.reweightLeptonSF*event.reweightLeptonHIPSF*event.reweightDilepTriggerBackup*nTrueInt27fb_puRW(event.nTrueInt)*event.reweightBTag_SF
    sample.weight         = weight_mc 
    sample.setSelectionString([getFilterCut(isData=False),  getLeptonSelection(mode)])

  stack          = Stack(mc,  data )

  jet_flavour_stack = Stack( *[[all_mc] for flavour in jet_flavours ] ) 

  if args.small:
        for sample in stack.samples + jet_flavour_stack.samples:
            sample.reduceFiles( to = 1 )

  # Use some defaults
  Plot.setDefaults(stack = stack, weight = None,   selectionString = selectionString, addOverFlowBin='upper')
  
  plots                 = []
  plot_flavour_pairs    = []

  plots.append(Plot(
    name = 'yield', texX = 'yield', texY = 'Number of Events',
    attribute = lambda event, sample: 0.5 + index,
    binning=[3, 0, 3],
  ))

  plot_flavour_pairs.append( [ 
    Plot(
    texX = '#eta(leading jet) (GeV)', texY = 'Number of Events',
    name = 'jet1_eta', attribute = lambda event, sample: event.JetGood_eta[0],
    binning=[104,-5.2,5.2],
  ), Plot(
    texX = '#eta(leading jet) (GeV)', texY = 'Number of Events',
    name = 'jet1_eta_flavours', attribute = lambda event, sample: event.JetGood_eta[0],
    stack  = jet_flavour_stack,
    weight = [ [flavour['weight'] ] for flavour in jet_flavours ],
    binning=[104,-5.2,5.2],
  )])

  for ptlow, pthigh in [ (20,30), (30,40), (40,50), (50,100), (100,1000)]:
      plot_flavour_pairs.append( [ 
        Plot(
        texX = '#eta(leading jet) (GeV)', texY = 'Number of Events',
        name = 'jet1_eta_pt_%i_%i'%(ptlow, pthigh), attribute = firstJetEta( ptlow, pthigh ),
        binning=[104,-5.2,5.2],
      ), Plot(
        texX = '#eta(leading jet) (GeV)', texY = 'Number of Events',
        name = 'jet1_eta_flavours_pt_%i_%i'%(ptlow, pthigh), attribute = firstJetEta( ptlow, pthigh ),
        stack  = jet_flavour_stack,
        weight = [ [flavour['weight'] ] for flavour in jet_flavours ],
        binning=[104,-5.2,5.2],
      )])
  for etalow, etahigh in [ (0,1), (1,2), (2,2.8), (2.8,3), (3,3.2), (3.2,4), (4,5.2)]:
      plot_flavour_pairs.append( [ 
        Plot(
        texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events',
        name = 'jet1_pt_eta_%i_%i'%(10*etalow, 10*etahigh), attribute = firstJetPt( etalow, etahigh ),
        binning=[100,30,230],
      ), Plot(
        texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events',
        name = 'jet1_pt_flavours_eta_%i_%i'%(10*etalow, 10*etahigh), attribute = firstJetPt( etalow, etahigh ),
        stack  = jet_flavour_stack,
        weight = [ [flavour['weight'] ] for flavour in jet_flavours ],
        binning=[100,30,230],
      )])
  for etalow, etahigh in [ (0,1), (1,2), (2,2.8), (2.8,3), (3,3.2), (3.2,4), (4,5.2)]:
      plot_flavour_pairs.append( [ 
        Plot(
        texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events',
        name = 'jet1_rawPt_eta_%i_%i'%(10*etalow, 10*etahigh), attribute = firstJetRawPt( etalow, etahigh ),
        binning=[100,30,230],
      ), Plot(
        texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events',
        name = 'jet1_rawPt_flavours_eta_%i_%i'%(10*etalow, 10*etahigh), attribute = firstJetRawPt( etalow, etahigh ),
        stack  = jet_flavour_stack,
        weight = [ [flavour['weight'] ] for flavour in jet_flavours ],
        binning=[100,30,230],
      )])


  plotting.fill( plots + sum(plot_flavour_pairs,[]), read_variables = read_variables, sequence = sequence )

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
  draw1DPlotsFlavourPairs( plot_flavour_pairs, mode, dataMCScale )
