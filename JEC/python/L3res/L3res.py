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

from math                                import sqrt, cos, sin, pi, atan2
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
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?' )#, default = True)
argParser.add_argument('--plot_directory',     action='store',      default='JEC/L3res')
args = argParser.parse_args()

#
# Logger
#
import JetMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

# JEC on the fly, tarball configuration
from JetMET.JetCorrector.JetCorrector import JetCorrector, correction_levels_data, correction_levels_mc

# config
Summer16_23Sep2016_DATA = \
[(1, 'Summer16_23Sep2016BCDV5_DATA'),
 (276831, 'Summer16_23Sep2016EFV5_DATA'),
 (278802, 'Summer16_23Sep2016GV5_DATA'),
 (280919, 'Summer16_23Sep2016HV5_DATA')]

Summer16_23Sep2016_MC = [(1, 'Summer16_23Sep2016V5_MC') ]

correction_levels_data  = [ 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' ]
correction_levels_mc    = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]

## pT_corr = pT_raw*L1(pT_raw)*L2L3(pT_raw*L1(pT_raw))*L2L3Res(pT_raw*L1(pT_raw)*L2L3(pT_raw*L1(pT_raw)))

jetCorrector_data = JetCorrector.fromTarBalls( Summer16_23Sep2016_DATA, correctionLevels = correction_levels_data )
jetCorrector_mc   = JetCorrector.fromTarBalls( Summer16_23Sep2016_MC, correctionLevels = correction_levels_mc )

if args.small: args.plot_directory += "_small"
#
# Make samples, will be searched for in the postProcessing directory
#
data_directory = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
postProcessing_directory = "postProcessed_80X_v37/dilepTiny/"
from JetMET.JEC.samples.cmgTuples_Summer16_mAODv2_postProcessed import *
data_directory = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
postProcessing_directory = "postProcessed_80X_v37/dilepTiny/"
from JetMET.JEC.samples.cmgTuples_Data25ns_80X_03Feb_postProcessed import *

selection       = 'ptll10-btb-njet1p'
selectionString = 'dl_pt>10&&cos(dl_phi-JetGood_phi[0])<-0.5&&Sum$(JetGood_pt>20 && JetGood_id)>=1' #&&(nJetGood==1||JetGood_pt[1]/dl_pt<0.3)'

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
def draw1DProfiles(plots, mode, dataMCScale, logX = False):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, args.plot_directory, mode + ("_log" if log else ""), selection)
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot

      plotting.draw(plot,
        plot_directory = plot_directory_,
        ratio = {'yRange': (0.8,1.2) },
        logX = 'profile_pt' in plot.name, logY = log, sorting = True,
        yRange =  (85,95) if plot.name.startswith('dl_mass') else (0.3, 1.5),# if log else (0.61, 1.41),
        scaling = {},
        legend = (0.50,0.88-0.04*sum(map(len, plot.histos)),0.9,0.88),
        drawObjects = drawObjects( dataMCScale , lumi_scale )
      )
#
# Read variables and sequences
#
read_variables = ["run/I", "weight/F", "l1_eta/F" , "l1_phi/F", "l1_dxy/F", "l2_eta/F", "l2_phi/F", "l2_dxy/F", "rho/F", 
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

# List all correction levels
corr_levels = ['raw', 'corr'] 

nan_jet = {key:float('nan') for key in jetVars + ['pt_%s'%corr_level for corr_level in corr_levels]}

def makeL3ResObservables( event, sample ):
    good_jets = filter( lambda j:j['pt']>20, getJets( event, jetColl="JetGood", jetVars = jetVars) )

    for j in good_jets:
        # Raw correction level
        j['pt_raw']  = j['rawPt']
        # 'Corr' correction level: L1L2L3 L2res
        if sample.isData:
            j['pt_corr'] =  jetCorrector_data.correction( j['rawPt'], j['eta'], j['area'], event.rho, event.run ) * j['rawPt'] 
        else:
            j['pt_corr'] =  jetCorrector_mc.  correction( j['rawPt'], j['eta'], j['area'], event.rho, event.run ) * j['rawPt'] 

    # compute type-1 MET shifts for chs met
    type1_met_shifts = \
        { corr_level: 
                {'px' :sum( ( j['rawPt'] - j['pt_%s'%corr_level] )*cos(j['phi']) for j in good_jets), 
                 'py' :sum( ( j['rawPt'] - j['pt_%s'%corr_level] )*sin(j['phi']) for j in good_jets) } 
          for corr_level in corr_levels }
    
    # leading jet
    event.leading_jet    =  good_jets[0]
    # subleading jet
    event.subleading_jet = good_jets[1] if len(good_jets)>=2 else nan_jet

    # alpha 
    event.alpha = event.subleading_jet['pt_corr'] / event.dl_pt #FIXME

    # alpha cut flag
    event.alpha_passed = ( event.nJetGood == 1 ) or (event.alpha < 0.3)

    for corr_level in corr_levels:

        # chs MET 
        chs_MEx_corr = event.met_chsPt*cos(event.met_chsPhi) + type1_met_shifts[corr_level]['px']
        chs_MEy_corr = event.met_chsPt*sin(event.met_chsPhi) + type1_met_shifts[corr_level]['py']

        chs_MEt_corr    = sqrt(  chs_MEx_corr**2 + chs_MEy_corr**2 )
        chs_MEphi_corr  = atan2( chs_MEy_corr, chs_MEx_corr )

        setattr( event, "met_chsPt_type1_%s"%corr_level,  chs_MEt_corr )
        setattr( event, "met_chsPhi_type1_%s"%corr_level, chs_MEphi_corr )

        # PT-bal
        setattr( event, "r_ptbal_%s"%corr_level,  event.leading_jet['pt_%s'%corr_level] / event.dl_pt )
        # MPF 
        setattr( event, "r_mpf_%s"%corr_level,  1. + chs_MEt_corr * cos(chs_MEphi_corr - event.dl_phi) / event.dl_pt )


sequence.append( makeL3ResObservables )

log_pt_thresholds = [10**(x/10.) for x in range(11,36)]

# Functor to retrieve attributes from the event object (Inline defined lambdas are wrongly bound)
def att_getter( arg, key = None):
    if key is None:
        def _f( event, sample ):
            return getattr( event, arg )
        return _f
    else:
        def _f( event, sample ):
            return getattr( event, arg )[key]
        return _f

def make_weight( dl_pt_bin = ( -1, -1), abs_eta_bin = None, require_alpha_passed = True):

    def _w( event, sample ):
        return \
            ( (not require_alpha_passed ) or event.alpha_passed ) \
            and  ( ( dl_pt_bin[0]<0 ) or (event.dl_pt > dl_pt_bin[0]) ) \
            and  ( ( dl_pt_bin[1]<0 ) or (event.dl_pt < dl_pt_bin[1]) ) \
            and  ( ( abs_eta_bin is None ) or (abs(event.leading_jet['eta']) >= abs_eta_bin[0] and abs(event.leading_jet['eta']) < abs_eta_bin[1]) ) 
    return _w

dl_pt_bins = [ (-1,-1), (20,30), (30,40), (40,50), (50, 100), (100, 200), (200, 500), (500, -1 )]
abs_eta_bins = [ (0, 1.3), (1.3, 2.5), (2.5, 3), (3, 5 )]

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
allModes   = ['mumu','ee']
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

  DY_sample     = DY_HT_LO
  TTJets_sample = Top

  other_mc_samples  = [TTZ_LO, TTXNoZ, multiBoson]
  all_mc_samples    = [DY_sample, TTJets_sample] + other_mc_samples

  other_mc          = Sample.combine( name = "other_mc", texName = "VV/VVV/TTX/tZq/tWZ", samples = other_mc_samples, color = ROOT.kMagenta )
  all_mc_combined   = Sample.combine( name = "all_mc_combined",   texName = "simulation", samples = all_mc_samples , color = ROOT.kBlue )
  all_mc_combined.style = styles.lineStyle( all_mc_combined.color, errors = True )    


  for sample in all_mc_samples: sample.style = styles.fillStyle(sample.color)

  for sample in all_mc_samples + [other_mc, all_mc_combined]:
    sample.scale          = lumi_scale
    sample.read_variables = ['reweightDilepTriggerBackup/F','reweightLeptonSF/F','reweightPU36fb/F', 'nTrueInt/F']
    sample.weight         = weight_mc 
    sample.setSelectionString([getFilterCut(isData=False),  getLeptonSelection(mode)])

  mc = [DY_sample, TTJets_sample, other_mc]
  stack          = Stack( mc,  data )
  stack_profile  = Stack( [all_mc_combined], data )

  if args.small:
        for sample in stack.samples + stack_profile.samples:
            sample.reduceFiles( to = 1 )

  # Use some defaults
  Plot.setDefaults( stack = stack, \
                    weight = lambda event, sample: event.alpha_passed,   
                    selectionString = selectionString, 
                    addOverFlowBin = None
    )
  
  plots      = []
  profiles1D = []

  plots.append(Plot(
    name = 'yield', texX = 'yield', texY = 'Number of Events',
    attribute = lambda event, sample: 0.5 + index,
    binning=[3, 0, 3],
  ))

  plots.append( 
    Plot(
    name = 'jet1_eta', 
    texX = '#eta(leading jet) (GeV)', texY = 'Number of Events',
    attribute = lambda event, sample: event.leading_jet['eta'],
    binning=[104,-5.2,5.2],
  ))

  plots.append( 
    Plot(
    name = 'jet2_eta', 
    texX = '#eta(subleading jet) (GeV)', texY = 'Number of Events',
    attribute = lambda event, sample: event.subleading_jet['eta'],
    binning=[104,-5.2,5.2],
  ))

  plots.append( 
    Plot(
    name = 'alpha', 
    texX = '#alpha', texY = 'Number of Events',
    attribute = lambda event, sample: event.alpha,
    binning=[50,0,1],
  ))

  plots.append( 
    Plot(
    name = 'dRj1j2', 
    texX = '#Delta R(j_{1}, j_{2})', texY = 'Number of Events',
    attribute = lambda event, sample: deltaR( event.leading_jet, event.subleading_jet),
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
      name = 'dl_eta', 
      texX = '#eta(ll) ', texY = 'Number of Events',
      attribute = lambda event, sample: event.dl_eta, 
      read_variables = ['dl_eta/F'],
      binning=[52,-5.2,5.2],
  ))

  plots.append(Plot(
    texX = '#phi(ll)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "dl_phi/F" ),
    binning=[10,-pi,pi],
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
    texX = 'p_{T}(l_{1}) (GeV)', texY = 'Number of Events / 15 GeV',
    attribute = TreeVariable.fromString( "l1_pt/F" ),
    binning=[20,0,300],
  ))

  plots.append(Plot(
    texX = 'd_{xy}(l_{1})', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "l1_dxy/F" ),
    binning=[40,-0.2,0.2],
  ))

  plots.append(Plot(
    texX = '#eta(l_{1})', texY = 'Number of Events',
    name = 'l1_eta', attribute = lambda event, sample: event.l1_eta, read_variables = ['l1_eta/F'],
    binning=[60,-3,3],
  ))

  plots.append(Plot(
    texX = '#phi(l_{1})', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "l1_phi/F" ),
    binning=[10,-pi,pi],
  ))

  plots.append(Plot(
    texX = 'p_{T}(l_{2}) (GeV)', texY = 'Number of Events / 15 GeV',
    attribute = TreeVariable.fromString( "l2_pt/F" ),
    binning=[20,0,300],
  ))

  plots.append(Plot(
    texX = 'd_{xy}(l_{2})', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "l2_dxy/F" ),
    binning=[40,-0.2,0.2],
  ))

  plots.append(Plot(
    texX = '#eta(l_{2})', texY = 'Number of Events',
    name = 'l2_eta', attribute = lambda event, sample: event.l2_eta, read_variables = ['l2_eta/F'],
    binning=[60,-3,3],
  ))

  plots.append(Plot(
    texX = '#phi(l_{2})', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "l2_phi/F" ),
    binning=[10,-pi,pi],
  ))

  profiles1D.append(Plot(
    name = 'dl_mass_profile_pt', texX = 'p_{T}(ll) (GeV)', texY = 'm(ll) (GeV)',
    histo_class = ROOT.TProfile,
    stack = stack_profile,
    attribute = (
        lambda event, sample: event.dl_pt,
        lambda event, sample: event.dl_mass,
    ),
    binning = Binning.fromThresholds(log_pt_thresholds),
  ))

  for corr_level in corr_levels:

      plots.append(Plot(
        name = 'jet1_pt_%s'%corr_level, 
        texX = 'p_{T}(leading %s jet) (GeV)'%corr_level, texY = 'Number of Events / 30 GeV',
        attribute = att_getter( "leading_jet", 'pt_%s'%corr_level ),
        binning=[600/30,0,600],
      ))

      plots.append(Plot(
        name = 'jet2_pt_%s'%corr_level, 
        texX = 'p_{T}(subleading %s jet) (GeV)'%corr_level, texY = 'Number of Events / 30 GeV',
        attribute = att_getter( "subleading_jet", 'pt_%s'%corr_level ),
        binning=[600/30,0,600],
      ))

      plots.append(Plot(
          name = "chsMET_%s" %corr_level,
          texX = 'chs-E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
          attribute = att_getter( "met_chsPt_type1_%s"%corr_level ),
          binning=[400/20,0,400],
      ))

      for method in ['ptbal', 'mpf']: 

          # inclusive response
          plots.append(Plot(
              name = "R_%s_%s"%(method, corr_level),
              texX = '%s R_{%s}'%(corr_level, method), texY = 'Number of Events',
              attribute = att_getter( "r_%s_%s"%(method, corr_level) ),
              binning=[100,0,3],
          ))

          # response profile wrt nvert
          profiles1D.append(Plot(
            name = 'r_%s_%s_profile_nvtx'%(method, corr_level), 
            texX = 'vertex multiplicity', 
            texY = '%s R_{%s}'%(corr_level, method),
            histo_class = ROOT.TProfile,
            stack = stack_profile,
            attribute = (
                att_getter( "nVert" ),
                att_getter( "r_%s_%s"%(method, corr_level) ),
            ),
            binning = [50,0,50],
            weight = make_weight( dl_pt_bin = (30, -1), abs_eta_bin = (0, 1.3) ),
          ))

          # response profile vs dl_pt
          for abs_eta_bin in abs_eta_bins:
              profiles1D.append(Plot(
                name = 'r_%s_%s_profile_ptll_for_eta_%3.2f_%3.2f'%( ( method, corr_level ) + abs_eta_bin ), 
                texX = 'p_{T}(ll) (GeV)', 
                texY = '%s R_{%s}'%(corr_level, method),
                histo_class = ROOT.TProfile,
                stack = stack_profile,
                attribute = (
                    "dl_pt",
                    att_getter( "r_%s_%s"%(method, corr_level) ),
                ),
                binning = Binning.fromThresholds(log_pt_thresholds),
                weight = make_weight( abs_eta_bin = abs_eta_bin ),
              ))

          # response profile vs eta
          for dl_pt_bin in dl_pt_bins:
              profiles1D.append(Plot(
                name = 'r_%s_%s_profile_jet_eta_for_dlpt_%i_%i'%( ( method, corr_level ) +  dl_pt_bin ), 
                texX = '#eta (jet)', 
                texY = '%s R_{%s}'%( corr_level, method),
                histo_class = ROOT.TProfile,
                stack = stack_profile,
                attribute = (
                    att_getter( "leading_jet", 'eta' ),
                    att_getter( "r_%s_%s"%(method, corr_level) ),
                ),
                binning = [26,-5.2,5.2],
                weight = make_weight( dl_pt_bin = dl_pt_bin ),
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
