#!/usr/bin/env python
''' Analysis script for L3 residuals (Z balancing) 
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
from JetMET.tools.user                   import plot_directory as user_plot_directory
from JetMET.tools.helpers                import deltaPhi, deltaR

# Object selection
from JetMET.tools.objectSelection        import getFilterCut, getJets, jetVars

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging" )
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')#, default = True)
argParser.add_argument('--plot_directory',     action='store',      default='JEC/L3res_new', help="subdirectory for plots")
args = argParser.parse_args()

#
# Logger
#
import JetMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

# JEC on the fly, tarball configuration
from JetMET.JetCorrector.JetCorrector import JetCorrector

## JetCorrector config
#Summer16_03Feb2017_DATA = \
#[(1,      'Summer16_03Feb2017BCD_%s_DATA'%args.version),
# (276831, 'Summer16_03Feb2017EF_%s_DATA'%args.version),
# (278802, 'Summer16_03Feb2017G_%s_DATA'%args.version),
# (280919, 'Summer16_03Feb2017H_%s_DATA'%args.version)]
#
#Summer16_03Feb2017_MC = [(1, 'Summer16_03Feb2017_V1_MC') ]
#
#correction_levels_data  = [ 'L1FastJet', 'L2Relative', 'L3Absolute' , 'L2L3Residual' ] if not args.noRes else [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]
#correction_levels_mc    = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]

#jetCorrector_data    = JetCorrector.fromTarBalls( Summer16_03Feb2017_DATA, correctionLevels = correction_levels_data )
#jetCorrector_mc      = JetCorrector.fromTarBalls( Summer16_03Feb2017_MC,   correctionLevels = correction_levels_mc )
#
#jetCorrector_RC_data = JetCorrector.fromTarBalls( Summer16_03Feb2017_DATA, correctionLevels = [ 'L1RC'] )
#jetCorrector_RC_mc   = JetCorrector.fromTarBalls( Summer16_03Feb2017_MC,   correctionLevels = [ 'L1RC'] )

#

from JetMET.JEC.samples.L2res_skim import *

#selection       = 'ptll%s-njet1p' % args.minptll
#selectionString = 'dl_pt>%s&&Sum$(JetGood_pt>10 && JetGood_id)>=1' % args.minptll 
#if args.btb:
#    selection+='-btb'
#    selectionString += '&&cos(dl_phi - JetGood_phi[0])<-0.5'
#
#selectionString = 'Sum$(JetGood_pt>10 && JetGood_id)>=1' #FIXME
## pt and thresholds
#from JetMET.JEC.L3res.thresholds import ptll_thresholds, ptll_bins, abs_eta_bins, coarse_ptll_bins, coarse_abs_eta_bins, L2res_abs_eta_bins
#
#
#all_abs_eta_bins = abs_eta_bins
#if args.addL2resEtaBins:
#    all_abs_eta_bins += L2res_abs_eta_bins
#
##if args.small:
##    ptll_bins = ptll_bins[:1]
##    all_abs_eta_bins = all_abs_eta_bins[:1]
#
## Text on the plots
##
#tex = ROOT.TLatex()
#tex.SetNDC()
#tex.SetTextSize(0.04)
#tex.SetTextAlign(11) # align right
#
#def drawObjects( dataMCScale, lumi_scale ):
#    lines = [
#      (0.15, 0.95, args.era), 
#      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) )
#    ]
#    return [tex.DrawLatex(*l) for l in lines] 
#
#plot_directory = os.path.join( user_plot_directory, args.plot_directory, args.version, args.dy, args.era )
#
## Formatting for 1D plots
#def draw1DPlots(plots, mode, dataMCScale):
#  for log in [False, True]:
#    plot_directory_ = os.path.join(plot_directory, mode + ("_log" if log else ""), selection)
#    for plot in plots:
#      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
#      p_drawObjects = map( lambda l:tex.DrawLatex(*l), getattr(plot, "drawObjects", [] ) )
#
#      if hasattr( plot, "subdir"):
#        plot_directory__ = os.path.join( plot_directory_, plot.subdir)
#      else:
#        plot_directory__ = plot_directory_  
#
#      plotting.draw(plot,
#        plot_directory = plot_directory__,
#        ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
#        logX = False, logY = log, sorting = True,
#        yRange         = (0.03, "auto") if log else (0.001, "auto"),
#        scaling        = {0:1} if len(plot.stack)==2 else {},
#        legend         = (0.50,0.88-0.04*sum(map(len, plot.histos)),0.9,0.88),
#        drawObjects    = drawObjects( dataMCScale , lumi_scale ) + p_drawObjects
#      )
#
##Formatting for 1D profiles
#def draw1DProfiles(plots, mode, dataMCScale):
#  for log in [False, True]:
#    plot_directory_ = os.path.join(plot_directory, mode + ("_log" if log else ""), selection)
#    for plot in plots:
#
#      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
#
#      p_drawObjects = map( lambda l:tex.DrawLatex(*l), getattr(plot, "drawObjects", [] ) )
# 
#      if plot.name.startswith('dl_mass'):
#          y_range     = (89, 93)
#          ratio_range = (0.99, 1.01)
#          logY = False
#      elif plot.name.startswith( 'ptll_profile_ptll'):
#          y_range     = (ptll_thresholds[0], ptll_thresholds[-1])
#          ratio_range = (0.97, 1.03)
#          logY = True
#      else:
#          y_range     = (0.3, 1.5)
#          ratio_range = (0.8, 1.2)
#          logY = False
#
#      if hasattr( plot, "subdir"):
#        plot_directory__ = os.path.join( plot_directory_, plot.subdir)
#      else:
#        plot_directory__ = plot_directory_  
#      
#      plotting.draw(plot,
#        plot_directory = plot_directory__,
#        ratio          = {'yRange': ratio_range } if not plot.name.startswith('r_gen_') else None,
#        logX           = 'profile_pt' in plot.name, 
#        logY           = logY, 
#        sorting        = False,
#        yRange         = y_range,
#        scaling        = {},
#        legend         = (0.50,0.88-0.04*sum(map(len, plot.histos)),0.9,0.88),
#        drawObjects    = drawObjects( dataMCScale , lumi_scale ) + p_drawObjects, 
#      )
#
##Formatting for 2D plots
#def draw2DPlots(plots, mode, dataMCScale):
#  for log in [False, True]:
#    plot_directory_ = os.path.join(plot_directory, mode + ("_log" if log else ""), selection)
#    for plot in plots:
#
#      p_drawObjects = map( lambda l:tex.DrawLatex(*l), getattr(plot, "drawObjects", [] ) )
#
#      if hasattr( plot, "subdir"):
#        plot_directory__ = os.path.join( plot_directory_, plot.subdir)
#      else:
#        plot_directory__ = plot_directory_  
#
#      plotting.draw2D(plot,
#        plot_directory = plot_directory__,
#        logX = False, 
#        #logX = '_vs_dl_pt' in plot.name, 
#        logY = False, 
#        #yRange =  (85,95) if plot.name.startswith('dl_mass') else (0.3, 1.5),# if log else (0.61, 1.41),
#        drawObjects = drawObjects( dataMCScale , lumi_scale ) + p_drawObjects, 
#      )
#
## decoration
#def dl_pt_string( b ):
#    if b[0]<0 and b[1]<0: return ""
#    res = "p_{T}(ll)"
#    if b[0]>0:
#        res = "%i"%(b[0]) + " #leq " + res
#    if b[1]>0:
#        res = res + "< %i"%(b[1]) 
#    return res
#
#def abs_eta_string( b ):
#    if b[0]<0 and b[1]<0: return ""
#    res = "|#eta(jet)|"
#    if b[0]>0:
#        res = "%4.3f"%(b[0]) + " #leq " + res
#    if b[1]>0:
#        res = res + "< %4.3f"%(b[1])
#    return res 
#
##
## Read variables and sequences
##
#read_variables = ["evt/l", "run/I", "weight/F", "l1_eta/F" , "l1_phi/F", "l1_dxy/F", "l2_eta/F", "l2_phi/F", "l2_dxy/F", "rho/F", 
#                  "dl_mass/F", "dl_eta/F", "dl_pt/F","dl_phi/F",
#                  "met_chsPt/F", "met_chsPhi/F", "metSig/F", "ht/F", "nBTag/I", "nJetGood/I", 'nVert/I']
#sequence = []
#
## extra lepton stuff
#read_variables += [
##     "nLepGood/I", "LepGood[pt/F,eta/F,etaSc/F,phi/F,pdgId/I,tightId/I,miniRelIso/F,relIso03/F,sip3d/F,mediumMuonId/I,mvaIdSpring15/F,lostHits/I,convVeto/I,dxy/F,dz/F,jetPtRelv2/F,jetPtRatiov2/F,eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz/I]",
##     "nLepOther/I", "LepOther[pt/F,eta/F,etaSc/F,phi/F,pdgId/I,tightId/I,miniRelIso/F,relIso03/F,sip3d/F,mediumMuonId/I,mvaIdSpring15/F,lostHits/I,convVeto/I,dxy/F,dz/F,jetPtRelv2/F,jetPtRatiov2/F,eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz/I]",
#    ]
#
#jetMCBranches = [ "mcMatchFlav/I", "partonId/I", "mcFlavour/I", "partonFlavour/I", "hadronFlavour/I", "mcMatchId/I", "mcPt/F" ]
#jetMCVars = [ s.split('/')[0] for s in jetMCBranches ]
#
## mcMatchFlav   Flavour of associated parton from hard scatter (if any) for Jets failing id after jet-lepton cleaning 
## partonId  parton flavour (manually matching to status 23 particles) for Jets failing id after jet-lepton cleaning 
## mcFlavour parton flavour (physics definition, i.e. including b's from shower) for Jets failing id after jet-lepton cleaning 
## partonFlavour purely parton-based flavour for Jets failing id after jet-lepton cleaning 
## hadronFlavour hadron flavour (ghost matching to B/C hadrons) for Jets failing id after jet-lepton cleaning 
## mcMatchId Match to source from hard scatter (pdgId of heaviest particle in chain, 25 for H, 6 for t, 23/24 for W/Z), zero if non-prompt or fake for Jets failing id after jet-lepton cleaning 
#
#jetVars += [ 'rawPt','mcPt' ]
#
#
#null_jet = {key:float('nan') for key in jetVars}
#null_jet['pt']         = 0
#null_jet['pt_corr']    = 0
#null_jet['pt_corr_RC'] = 0
#
#alphas = ["30", "20", "15", "10"]
#
#def makeL3ResObservables( event, sample ):
#    good_jets = getJets( event, jetColl="JetGood", jetVars = jetVars)
#
#    for j in good_jets:
#        # 'Corr' correction level: L1L2L3 L2res
#        if sample.isData:
#            jet_corr_factor    =  jetCorrector_data.   correction( j['rawPt'], j['eta'], j['area'], event.rho, event.run )
#            jet_corr_factor_RC =  jetCorrector_RC_data.correction( j['rawPt'], j['eta'], j['area'], event.rho, event.run )
#        else:
#            jet_corr_factor    =  jetCorrector_mc.     correction( j['rawPt'], j['eta'], j['area'], event.rho, event.run )  
#            jet_corr_factor_RC =  jetCorrector_RC_mc.  correction( j['rawPt'], j['eta'], j['area'], event.rho, event.run )  
#        
#        # corrected jet
#        j['pt_corr']    =  jet_corr_factor * j['rawPt'] 
#
#        # noL1 -> divide out L1FastJet, remove 
#        if args.noL1: 
#            if sample.isData:
#                jet_corr_factor_L1 =  jetCorrector_L1_data.correction( j['rawPt'], j['eta'], j['area'], event.rho, event.run ) 
#            else: 
#                jet_corr_factor_L1 =  jetCorrector_L1_mc.  correction( j['rawPt'], j['eta'], j['area'], event.rho, event.run ) 
#            # noL1 -> divide out L1FastJet, remove 
#            j['pt_corr']    =  j['pt_corr']/jet_corr_factor_L1 
#            # no L1RC if 'noL1'
#            j['pt_corr_RC'] =  j['rawPt'] 
#        else:
#            # L1RC 
#            j['pt_corr_RC'] =  jet_corr_factor_RC * j['rawPt'] 
#
#
#    # compute type-1 MET shifts for chs met L1L2L3 - L1RC (if 'noL1', then L1FastJets is divided out and L1RC is not applied )
#    type1_met_shifts = \
#                {'px' :sum( ( j['pt_corr_RC'] - j['pt_corr'] )*cos(j['phi']) for j in good_jets), 
#                 'py' :sum( ( j['pt_corr_RC'] - j['pt_corr'] )*sin(j['phi']) for j in good_jets) } 
#
#    # leading jet
#    event.leading_jet    = good_jets[0]
#    # subleading jet
#    event.subleading_jet = good_jets[1] if len(good_jets)>=2 else null_jet
#
#    # alpha 
#    event.alpha = event.subleading_jet['pt_corr'] / event.dl_pt
#    # alpha cut flag
#    event.alpha_30_passed = ( event.alpha < 0.3)
#    event.alpha_20_passed = ( event.alpha < 0.2)
#    event.alpha_15_passed = ( event.alpha < 0.15)
#    event.alpha_10_passed = ( event.alpha < 0.1)
#
#    # chs MET 
#    chs_MEx_corr = event.met_chsPt*cos(event.met_chsPhi) + type1_met_shifts['px']
#    chs_MEy_corr = event.met_chsPt*sin(event.met_chsPhi) + type1_met_shifts['py']
#
#    chs_MEt_corr    = sqrt(  chs_MEx_corr**2 + chs_MEy_corr**2 )
#    chs_MEphi_corr  = atan2( chs_MEy_corr, chs_MEx_corr )
##    if sample.isData: #FIXME
##        print "evt", event.evt
##        for i, j in enumerate(good_jets):
##            print "jet", i, "pt(raw)", j['rawPt'], "pt(L1L2L3)", j['pt_corr'], "pt(L1RC)", j['pt_corr_RC'], "phi", cos(j['phi'])
##            print "cont. to type1 from jet ", i, ( j['pt_corr_RC'] - j['pt_corr'] )*cos(j['phi']), "py", ( j['pt_corr_RC'] - j['pt_corr'] )*sin(j['phi'])
##        print "type1 shifts", type1_met_shifts['px'], type1_met_shifts['py'] 
##        print "raw chs met: pt", event.met_chsPt, 'phi', event.met_chsPhi
##        print "type1 chs met: px", chs_MEx_corr, 'py', chs_MEy_corr 
##        print "             : pt",chs_MEt_corr,"phi",chs_MEphi_corr 
##        print 
#    setattr( event, "met_chsPt_type1",  chs_MEt_corr )
#    setattr( event, "met_chsPhi_type1", chs_MEphi_corr )
#
#    # PT-bal
#    event.r_ptbal      = event.leading_jet['pt_corr'] / event.dl_pt
#    # PT-bal raw
#    event.r_ptbalRaw   = event.leading_jet['rawPt'] / event.dl_pt 
#    # MPF 
#    event.r_mpf        = 1. + chs_MEt_corr * cos(chs_MEphi_corr - event.dl_phi) / event.dl_pt
#    # MPF no type-1 
#    event.r_mpfNoType1 = 1. + event.met_chsPt * cos(event.met_chsPhi - event.dl_phi) / event.dl_pt
#
#    # gen 
#    if not sample.isData and event.leading_jet['mcPt']>0:
#        event.r_gen    = event.leading_jet['pt_corr']/event.leading_jet['mcPt']
#        #event.r_gen    = event.leading_jet['pt']/event.leading_jet['mcPt']
#    else:
#        event.r_gen    = None 
#
#sequence.append( makeL3ResObservables )
#
## Functor to retrieve attributes from the event object (Inline defined lambdas are wrongly bound)
#def att_getter( arg, key = None):
#    if key is None:
#        def _f( event, sample ):
#            return getattr( event, arg )
#        return _f
#    else:
#        def _f( event, sample ):
#            return getattr( event, arg )[key]
#        return _f
#
#def make_weight( ptll_bin = ( -1, -1), abs_eta_bin = None, alpha_passed = "alpha_30_passed", is_finite = []):
#
#    def _w( event, sample ):
#        return \
#            ( ( alpha_passed is None) or getattr(event, alpha_passed) ) \
#            and  ( ( ptll_bin[0]<0 ) or (event.dl_pt > ptll_bin[0]) ) \
#            and  ( ( ptll_bin[1]<0 ) or (event.dl_pt < ptll_bin[1]) ) \
#            and  ( ( abs_eta_bin is None ) or (abs(event.leading_jet['eta']) >= abs_eta_bin[0] and abs(event.leading_jet['eta']) < abs_eta_bin[1]) ) \
#            and  ( False not in [getattr(event, att)<float('inf') for att in is_finite] )
#    return _w
#
#
#z_window = 10
#def getLeptonSelection( mode ):
#  if   mode=="mumu": return "nGoodMuons==2&&nGoodElectrons==0&&isOS&&isMuMu&&abs(dl_mass-91.2)<%f" % z_window
#  elif mode=="ee":   return "nGoodMuons==0&&nGoodElectrons==2&&isOS&&isEE&&abs(dl_mass-91.2)<%f" % z_window
#
#weight_mc   = lambda event, sample: event.weight*event.reweightLeptonSF*event.reweightDilepTriggerBackup*event.reweightPU36fb
#weight_data = lambda event, sample: event.weight
#
#if   args.mode=="mumu": 
#    if args.era == 'inclusive':
#        data = DoubleMuon_Run2016_backup
#        runrange = None
#    elif args.era == "Run2016BCD":
#        data = DoubleMuon_Run2016BCD_backup
#        runrange = None
#    elif args.era == "Run2016EFearly":
#        data = DoubleMuon_Run2016EFearly_backup #F early
#        runrange = "run<=278801"
#    elif args.era == "Run2016FlateGH":
#        data = DoubleMuon_Run2016FlateGH_backup # contains also F for Flate
#        runrange = "run>=278802" 
#    elif args.era == "Run2016FlateG":
#        data = DoubleMuon_Run2016FlateG_backup # contains also F for Flate
#        runrange = "run>=278802" 
#    elif args.era == "Run2016H":
#        data = DoubleMuon_Run2016H_backup
#        runrange = None 
#
#    data.texName = "data (2 #mu)"
#    index = 0
#
#elif args.mode=="ee":
#    if args.era == 'inclusive':
#        data = DoubleEG_Run2016_backup
#        runrange = None
#    elif args.era == "Run2016BCD":
#        data = DoubleEG_Run2016BCD_backup
#        runrange = None
#    elif args.era == "Run2016EFearly":
#        data = DoubleEG_Run2016EFearly_backup # F early
#        runrange = "run<=278801"
#    elif args.era == "Run2016FlateGH":
#        data = DoubleEG_Run2016FlateGH_backup # contains also F for Flate
#        runrange = "run>=278802" 
#    elif args.era == "Run2016FlateG":
#        data = DoubleEG_Run2016FlateG_backup # contains also F for Flate
#        runrange = "run>=278802" 
#    elif args.era == "Run2016H":
#        data = DoubleEG_Run2016H_backup
#        runrange = None
#
#    data.texName = "data (2 e)"
#    index = 1
#
##data = test_data #FIXME 
##data.lumi=1
#
#yields     = {}
#allPlots   = {}
#yields[args.mode] = {}
#
#data_selection = [getFilterCut(isData=True), getLeptonSelection(args.mode)]
#if runrange is not None:
#    data_selection.append( runrange )
#
#logger.info( "Set data selectionstring: %s", "&&".join(data_selection) )
#
##data.setSelectionString( data_selection ) #FIXME
#data.name           = "data"
#data.style          = styles.errorStyle(ROOT.kBlack)
#data.weight         = weight_data 
#data.read_variables = [ "JetGood[pt/F,eta/F,phi/F,area/F,btagCSV/F,rawPt/F]" ]
##data.read_variables += [ "photon[pt/F,eta/F,phi/F]", "nPhotonGood/I" ]
#
#lumi_scale          = data.lumi/1000
#
#if args.dy == 'DY_HT_LO':
#    DY_sample      = DY_HT_LO
#    dy_legend_text = "DY(HT)"
#elif args.dy == 'DYnJets':
#    DY_sample      = DYnJets
#    dy_legend_text = "DY(Njet)"
#
#DY_sample.legendText = dy_legend_text
#
#TTJets_sample = TTLep_pow
#
##other_mc_samples  = [TTZ_LO, TTXNoZ, multiBoson]
#all_mc_samples    = [DY_sample, TTJets_sample]# + other_mc_samples
#
##other_mc          = Sample.combine( name = "other_mc", texName = "VV/VVV/TTX/tZq/tWZ", samples = other_mc_samples, color = ROOT.kMagenta )
#all_mc_combined   = Sample.combine( name = "all_mc_combined",   texName = "%s + rest"%dy_legend_text, samples = all_mc_samples , color = ROOT.kBlue )
#all_mc_combined.style = styles.lineStyle( all_mc_combined.color, errors = True )    
#
#
#for sample in all_mc_samples: sample.style = styles.fillStyle(sample.color)
##other_mc.style = styles.fillStyle(ROOT.kMagenta)
#
#for sample in all_mc_samples + [all_mc_combined]: #+[other_mc]
#    sample.scale          = lumi_scale
#    sample.read_variables = ['reweightDilepTriggerBackup/F','reweightLeptonSF/F','reweightPU36fb/F', 'nTrueInt/F', "JetGood[pt/F,eta/F,phi/F,area/F,btagCSV/F,rawPt/F,mcPt/F]"]
#    sample.weight         = weight_mc 
#    sample.setSelectionString([getFilterCut(isData=False),  getLeptonSelection(args.mode)])
#
#mc = [DY_sample, TTJets_sample] #, other_mc]
#stack             = Stack( mc,  data )
#stack_gen         = Stack( mc )
#stack_profile     = Stack( [all_mc_combined], data )
#stack_profile_gen = Stack( [all_mc_combined] )
#
#if args.small:
#    for sample in stack.samples + stack_profile.samples:
#        sample.reduceFiles( to = 1 )
#
## Use some defaults
#Plot.setDefaults( stack = stack, \
#                weight = lambda event, sample: event.alpha_30_passed,   
#                selectionString = selectionString, 
#                addOverFlowBin = None
#)
#
#plots      = []
#profiles1D = []
#plots2D    = []
#
## 1D plots
#plots.append(Plot(
#name = 'yield', texX = 'yield', texY = 'Number of Events',
#attribute = lambda event, sample: 0.5 + index,
#binning=[3, 0, 3],
#))
#plots[-1].subdir = "plots"
#
#if not args.skipOtherPlots:
#
#    plots.append( 
#    Plot(
#    name = 'jet1_eta', 
#    texX = '#eta(leading jet) (GeV)', texY = 'Number of Events',
#    attribute = lambda event, sample: event.leading_jet['eta'],
#    binning=[104,-5.2,5.2],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append( 
#    Plot(
#    name = 'jet2_eta', 
#    texX = '#eta(subleading jet) (GeV)', texY = 'Number of Events',
#    attribute = lambda event, sample: event.subleading_jet['eta'],
#    binning=[104,-5.2,5.2],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append( 
#    Plot(
#    name = 'alpha', 
#    texX = '#alpha', texY = 'Number of Events',
#    attribute = lambda event, sample: event.alpha,
#    weight = lambda event, sample: 1, # Don't use alpha weight for alpha plot
#    binning=[50,0,1],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append( 
#    Plot(
#    name = 'dRj1j2', 
#    texX = '#Delta R(j_{1}, j_{2})', texY = 'Number of Events',
#    attribute = lambda event, sample: deltaR( event.leading_jet, event.subleading_jet),
#    binning=[60,0,6],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = 'number of jets', texY = 'Number of Events',
#    attribute = TreeVariable.fromString('nJetGood/I'),
#    binning=[14,0,14],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = 'number of medium b-tags (CSVM)', texY = 'Number of Events',
#    attribute = TreeVariable.fromString('nBTag/I'),
#    binning=[8,0,8],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = 'H_{T} (GeV)', texY = 'Number of Events / 25 GeV',
#    attribute = TreeVariable.fromString( "ht/F" ),
#    binning=[500/25,0,600],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = 'gen-H_{T} (GeV)', texY = 'Number of Events / 25 GeV',
#    attribute = TreeVariable.fromString( "lheHTIncoming/F" ),
#    stack = Stack( mc ),
#    binning=[500/25,0,600],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = 'm(ll) of leading dilepton (GeV)', texY = 'Number of Events / 4 GeV',
#    attribute = TreeVariable.fromString( "dl_mass/F" ),
#    binning=[200/4,0,200],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = 'p_{T}(ll) (GeV)', texY = 'Number of Events / 70 GeV',
#    attribute = TreeVariable.fromString( "dl_pt/F" ),
#    binning = Binning.fromThresholds(ptll_thresholds),
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#      name = 'dl_eta', 
#      texX = '#eta(ll) ', texY = 'Number of Events',
#      attribute = lambda event, sample: event.dl_eta, 
#      read_variables = ['dl_eta/F'],
#      binning=[52,-5.2,5.2],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = '#phi(ll)', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "dl_phi/F" ),
#    binning=[10,-pi,pi],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    name = 'nVtxs', texX = 'vertex multiplicity', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "nVert/I" ),
#    binning=[50,0,50],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    name = 'rho', texX = 'rho', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "rho/F" ),
#    binning=[50,0,50],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = 'p_{T}(l_{1}) (GeV)', texY = 'Number of Events / 15 GeV',
#    attribute = TreeVariable.fromString( "l1_pt/F" ),
#    binning=[20,0,300],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = 'd_{xy}(l_{1})', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "l1_dxy/F" ),
#    binning=[40,-0.2,0.2],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = '#eta(l_{1})', texY = 'Number of Events',
#    name = 'l1_eta', attribute = lambda event, sample: event.l1_eta, read_variables = ['l1_eta/F'],
#    binning=[60,-3,3],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = '#phi(l_{1})', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "l1_phi/F" ),
#    binning=[10,-pi,pi],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = 'p_{T}(l_{2}) (GeV)', texY = 'Number of Events / 15 GeV',
#    attribute = TreeVariable.fromString( "l2_pt/F" ),
#    binning=[20,0,300],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = 'd_{xy}(l_{2})', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "l2_dxy/F" ),
#    binning=[40,-0.2,0.2],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = '#eta(l_{2})', texY = 'Number of Events',
#    name = 'l2_eta', attribute = lambda event, sample: event.l2_eta, read_variables = ['l2_eta/F'],
#    binning=[60,-3,3],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    texX = '#phi(l_{2})', texY = 'Number of Events',
#    attribute = TreeVariable.fromString( "l2_phi/F" ),
#    binning=[10,-pi,pi],
#    ))
#    plots[-1].subdir = "plots"
#
#    # 2D plots
#    plots2D.append(Plot2D(
#    name = 'alpha_vs_dl_pt_data', 
#    texX = 'p_{T}(ll) (GeV)', 
#    texY = '#alpha',
#    stack = Stack(data),
#    selectionString = selectionString, 
#    attribute = (
#        lambda event, sample: event.dl_pt,
#        lambda event, sample: event.alpha,
#    ),
#    binning=[50,0,200,50,0,1],
#    ))
#    plots2D[-1].subdir = "plots"
#
#    plots2D.append(Plot2D(
#    name = 'alpha_vs_dl_pt_mc', 
#    texX = 'p_{T}(ll) (GeV)', 
#    texY = '#alpha',
#    stack = Stack(mc),
#    selectionString = selectionString, 
#    attribute = (
#        lambda event, sample: event.dl_pt,
#        lambda event, sample: event.alpha,
#    ),
#    binning=[50,0,200,50,0,1],
#    weight = None,
#    ))
#    plots2D[-1].subdir = "plots"
#
#    plots2D.append(Plot2D(
#    name = 'DY_Rmpf_vs_Rptbal', 
#    texX = '%s R_{mpf}'%args.version,
#    texY = '%s R_{ptbal}'%args.version,
#    stack = Stack(DY_sample),
#    selectionString = selectionString, 
#    attribute = (
#        att_getter( "r_mpf" ),
#        att_getter( "r_ptbal" )
#    ),
#    binning=[50,-.5,2,50,0,2],
#    weight = lambda event, sample: event.alpha_30_passed,
#    ))
#    plots2D[-1].subdir = "plots"
#
#    plots.append(Plot(
#    name = 'jet1_pt_%s'%args.version, 
#    texX = 'p_{T}(leading %s jet) (GeV)'%args.version, texY = 'Number of Events / 30 GeV',
#    attribute = att_getter( "leading_jet", 'pt_corr' ),
#    binning=[600/30,0,600],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#    name = 'jet2_pt_%s'%args.version, 
#    texX = 'p_{T}(subleading %s jet) (GeV)'%args.version, texY = 'Number of Events / 30 GeV',
#    attribute = att_getter( "subleading_jet", 'pt_corr' ),
#    binning=[600/30,0,600],
#    ))
#    plots[-1].subdir = "plots"
#
#    plots.append(Plot(
#      name = "chsMET_%s" %args.version,
#      texX = 'chs-E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
#      attribute = att_getter( "met_chsPt_type1" ),
#      binning=[400/20,0,400],
#    ))
#    plots[-1].subdir = "plots"
#
## pt profile vs dl_pt
#for abs_eta_bin in all_abs_eta_bins:
#  profiles1D.append(Plot(
#    name = ( 'ptll_profile_ptll_for_eta_%4.3f_%4.3f'%( abs_eta_bin )).replace('.',''), 
#    texX = 'p_{T}(ll) (GeV)', 
#    texY = 'p_{T}(ll) (GeV)', 
#    histo_class = ROOT.TProfile,
#    stack = stack_profile,
#    attribute = (
#        "dl_pt",
#        "dl_pt",
#    ),
#    binning = Binning.fromThresholds(ptll_thresholds),
#    weight = make_weight( abs_eta_bin = abs_eta_bin ),
#  ))
#  profiles1D[-1].drawObjects = [(0.5, 0.76, abs_eta_string(abs_eta_bin))]
#
#profiles1D.append(Plot(
#name = 'dl_mass_profile_pt', texX = 'p_{T}(ll) (GeV)', texY = 'm(ll) (GeV)',
#histo_class = ROOT.TProfile,
#stack = stack_profile,
#attribute = (
#    lambda event, sample: event.dl_pt,
#    lambda event, sample: event.dl_mass,
#),
#binning = Binning.fromThresholds(ptll_thresholds),
#weight = make_weight( abs_eta_bin = (0, 1.3), alpha_passed = "alpha_30_passed" ),
#))
#profiles1D[-1].drawObjects = [(0.2, 0.8, abs_eta_string((0, 1.3))), (0.2,0.75, "#alpha<0.3") ]
#
#
#for_comparison_eta = {abs_eta_bin:{alpha:{} for alpha in alphas } for abs_eta_bin in all_abs_eta_bins}
#for_comparison_pt  = {ptll_bin:   {alpha:{} for alpha in alphas } for ptll_bin    in ptll_bins}
#
#methods = ['ptbal', 'mpf',  'mpfNoType1']
#if not args.skipOtherPlots: methods += ['ptbalRaw', 'gen']
#
#for method in methods: 
#
#  stack_profile_ = stack_profile_gen if method=='gen' else stack_profile
#  stack_         = stack_gen if method=='gen' else stack
#
#  # response profile vs dl_pt
#  for alpha in [ "30", "20", "15", "10"]: 
#    for abs_eta_bin in all_abs_eta_bins:
#      profiles1D.append(Plot(
#        name = ('r_%s_%s_a%s_profile_ptll_for_eta_%4.3f_%4.3f'%( ( method, args.version, alpha) + abs_eta_bin )).replace('.',''), 
#        texX = 'p_{T}(ll) (GeV)', 
#        texY = '%s R_{%s}'%(args.version, method),
#        histo_class = ROOT.TProfile,
#        stack = stack_profile_,
#        attribute = (
#            "dl_pt",
#            att_getter( "r_%s"%method ),
#        ),
#        binning = Binning.fromThresholds(ptll_thresholds),
#        weight = make_weight( abs_eta_bin = abs_eta_bin, alpha_passed = "alpha_%s_passed"%alpha,  is_finite = [ "r_%s"%method ] ),
#      ))
#      profiles1D[-1].drawObjects = [(0.5, 0.76, abs_eta_string(abs_eta_bin)), (0.5, 0.71, "#alpha<%2.1f"% ( float(alpha)/100.)) ]
#      for_comparison_eta[abs_eta_bin][alpha][method] = profiles1D[-1]
#
#    if not args.skipOtherPlots:
#
#      # inclusive response
#      plots.append(Plot(
#          name = "R_%s_%s"%(method, args.version),
#          texX = '%s R_{%s}'%(args.version, method), texY = 'Number of Events',
#          attribute = att_getter( "r_%s"%method ),
#          stack = stack_,
#          binning=[100,0,3],
#          weight = make_weight( ptll_bin = (30, -1), abs_eta_bin = (0, 1.3), is_finite = [ "r_%s"%method ] ),
#      ))
#      plots[-1].subdir = "response_plots"
#
#      # response profile vs eta
#      for ptll_bin in ptll_bins:
#        profiles1D.append(Plot(
#          name = 'r_%s_%s_profile_jet_eta_for_dlpt_%i_%i'%( ( method, args.version ) +  ptll_bin ), 
#          texX = '#eta (jet)', 
#          texY = '%s R_{%s}'%( args.version, method),
#          histo_class = ROOT.TProfile,
#          stack = stack_profile_,
#          attribute = (
#              att_getter( "leading_jet", 'eta' ),
#              att_getter( "r_%s"%method ),
#          ),
#          binning = [26,-5.2,5.2],
#          weight = make_weight( ptll_bin = ptll_bin, alpha_passed = "alpha_%s_passed"%alpha, is_finite = [ "r_%s"%method ] ),
#        ))
#        profiles1D[-1].drawObjects = [(0.5, 0.76, dl_pt_string(ptll_bin)), (0.5, 0.71, "#alpha<%2.1f"% ( float(alpha)/100.)) ]
#        for_comparison_pt[ptll_bin][alpha][method] = profiles1D[-1]
#
#      # response profile wrt nvert, binned in pT and eta
#      for abs_eta_bin in coarse_abs_eta_bins:
#          for ptll_bin in coarse_ptll_bins:
#              profiles1D.append(Plot(
#                name = ('r_%s_%s_profile_nvtx_dlpt_%i_%i_eta_%4.3f_%4.3f'% ( (method, args.version) + ptll_bin + abs_eta_bin )).replace('.',''), 
#                texX = 'vertex multiplicity', 
#                texY = '%s R_{%s}'%(args.version, method),
#                histo_class = ROOT.TProfile,
#                stack = stack_profile_,
#                attribute = (
#                    att_getter( "nVert" ),
#                    att_getter( "r_%s"%method ),
#                ),
#                binning = [50,0,50],
#                weight = make_weight( ptll_bin = ptll_bin, abs_eta_bin = abs_eta_bin, is_finite = [ "r_%s"%method ] ),
#              ))
#              profiles1D[-1].drawObjects = [(0.5, 0.76, dl_pt_string(ptll_bin)), (0.5, 0.71, abs_eta_string(abs_eta_bin)), (0.5, 0.66, "#alpha<0.3") ]
#              profiles1D[-1].subdir = "response_nvtx"
#
#      # response, binned in pT and eta 
#      for abs_eta_bin in abs_eta_bins:
#          for ptll_bin in ptll_bins:
#              plots.append(Plot(
#                  name = ( "R_%s_%s_dlpt_%i_%i_eta_%4.3f_%4.3f"%( (method, args.version) + ptll_bin + abs_eta_bin )).replace('.',''),
#                  texX = '%s R_{%s}'%(args.version, method), texY = 'Number of Events',
#                  attribute = att_getter( "r_%s"%method ),
#                  stack = stack_,
#                  binning=[100,0,3],
#                  weight = make_weight( ptll_bin = ptll_bin, abs_eta_bin = abs_eta_bin, is_finite = [ "r_%s"%method ] ),
#              ))
#              plots[-1].drawObjects = [(0.5, 0.76, dl_pt_string(ptll_bin)), (0.5, 0.71, abs_eta_string(abs_eta_bin)), (0.5, 0.66, "#alpha<0.3") ]
#              plots[-1].subdir = "response_plots"
#
#
#plotting.fill( plots + profiles1D + plots2D , read_variables = read_variables, sequence = sequence, max_events = 50000 if args.small else -1) #FIXME
#
## Get normalization yields from yield histogram
#for plot in plots:
#    if plot.name == "yield":
#      for i, l in enumerate(plot.histos):
#        for j, h in enumerate(l):
#          yields[args.mode][plot.stack[i][j].name] = h.GetBinContent( h.FindBin( 0.5 + index ) )
#          h.GetXaxis().SetBinLabel(1, "#mu#mu")
#          h.GetXaxis().SetBinLabel(2, "e#mu")
#          h.GetXaxis().SetBinLabel(3, "ee")
#
#yields[args.mode]["MC"] = sum(yields[args.mode][s.name] for s in mc)
#dataMCScale        = yields[args.mode]["data"]/yields[args.mode]["MC"] if yields[args.mode]["MC"] != 0 else float('nan')
#
#draw1DPlots(    plots,      args.mode, dataMCScale )
#draw1DProfiles( profiles1D, args.mode, dataMCScale )
#draw2DPlots(    plots2D,    args.mode, dataMCScale )
#
#if not args.skipOtherPlots:
#    for alpha in alphas:
#        for abs_eta_bin in abs_eta_bins:
#          cp = for_comparison_eta[abs_eta_bin][alpha]['mpf']
#
#          cp.histos += for_comparison_eta[abs_eta_bin][alpha]['mpfNoType1'].histos[:1]
#          cp.stack  += for_comparison_eta[abs_eta_bin][alpha]['mpfNoType1'].stack[:1]
#          cp.histos[-1][0].style      = styles.lineStyle( ROOT.kBlue, dashed = True ) 
#          cp.histos[-1][0].legendText = "R_{mpf} raw" 
#
#          cp.histos += for_comparison_eta[abs_eta_bin][alpha]['gen'].histos 
#          cp.stack  += for_comparison_eta[abs_eta_bin][alpha]['gen'].stack 
#          cp.histos[-1][0].style      = styles.lineStyle( ROOT.kRed ) 
#          cp.histos[-1][0].legendText = "generated" 
#
#          cp.histos += for_comparison_eta[abs_eta_bin][alpha]['ptbal'].histos[:1]
#          cp.stack  += for_comparison_eta[abs_eta_bin][alpha]['ptbal'].stack[:1]
#          cp.histos[-1][0].style      = styles.lineStyle( ROOT.kGreen ) 
#          cp.histos[-1][0].legendText = "R_{ptbal}" 
#
#          cp.histos += for_comparison_eta[abs_eta_bin][alpha]['ptbalRaw'].histos[:1]
#          cp.stack  += for_comparison_eta[abs_eta_bin][alpha]['ptbalRaw'].stack[:1]
#          cp.histos[-1][0].style      = styles.lineStyle( ROOT.kGreen, dashed = True ) 
#          cp.histos[-1][0].legendText = "raw R_{ptbal}" 
#
#          cp.name = ('comparison_%s_a%s_profile_pt_for_eta_%4.3f_%4.3f'% ( (args.version, alpha) +  abs_eta_bin ) ).replace('.','')
#          cp.drawObjects = [(0.2, 0.8, abs_eta_string(abs_eta_bin)), (0.2,0.75, "#alpha<%2.1f"% ( float(alpha)/100.) )]
#          draw1DProfiles(  [cp] , args.mode, dataMCScale )
#
#        for ptll_bin in ptll_bins:
#          cp = for_comparison_pt[ptll_bin][alpha]['mpf']
#
#          cp.histos += for_comparison_pt[ptll_bin][alpha]['mpfNoType1'].histos[:1]
#          cp.stack  += for_comparison_pt[ptll_bin][alpha]['mpfNoType1'].stack[:1]
#          cp.histos[-1][0].style      = styles.lineStyle( ROOT.kBlue, dashed=True ) 
#          cp.histos[-1][0].legendText = "R_{mpf} raw"
#
#          cp.histos += for_comparison_pt[ptll_bin][alpha]['gen'].histos 
#          cp.stack += for_comparison_pt[ptll_bin][alpha]['gen'].stack 
#          cp.histos[-1][0].style      = styles.lineStyle( ROOT.kBlack ) 
#          cp.histos[-1][0].legendText = "generated"
#
#          cp.histos += for_comparison_pt[ptll_bin][alpha]['ptbal'].histos[:1]
#          cp.stack  += for_comparison_pt[ptll_bin][alpha]['ptbal'].stack[:1]
#          cp.histos[-1][0].style      = styles.lineStyle( ROOT.kGreen ) 
#          cp.histos[-1][0].legendText = "R_{ptbal}" 
#
#          cp.histos += for_comparison_pt[ptll_bin][alpha]['ptbalRaw'].histos[:1]
#          cp.stack  += for_comparison_pt[ptll_bin][alpha]['ptbalRaw'].stack[:1]
#          cp.histos[-1][0].style      = styles.lineStyle( ROOT.kGreen, dashed = True ) 
#          cp.histos[-1][0].legendText = "raw R_{ptbal}" 
#
#          cp.name = 'comparison_%s_a%s_profile_jet_eta_for_dlpt_%i_%i'% ( (args.version, alpha) +  ptll_bin ) 
#          cp.drawObjects = [(0.2, 0.8, dl_pt_string(ptll_bin)), (0.2,0.75, "#alpha<%2.1f"% ( float(alpha)/100.) )]
#
#          draw1DProfiles(  [cp] , args.mode, dataMCScale )
