#!/usr/bin/env python
''' Analysis script for L1 closure plots
'''
#
# Standard imports and batch mode
#
import ROOT
ROOT.gROOT.SetBatch(True)
import os
import array

from RootTools.core.standard             import *
from JetMET.tools.user                   import plot_directory

# Object selection
from JetMET.tools.objectSelection        import getJets, jetVars

# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action  = 'store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',              default = False,         action ='store_true',    help='Run only on a small subset of the data?', )
argParser.add_argument('--plot_directory',     action  = 'store',      default='JEC/mcEta')
args = argParser.parse_args()

#
# Logger
#
import JetMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small: args.plot_directory += "_small"

QCD_flat = Sample.fromDPMDirectory( "QCD_flat", "/dpm/oeaw.ac.at/home/cms/store/user/schoef/cmgTuples/80X_1l_31/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISummer16MiniAODv2-NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v4-v2_80X_1l_31/170218_094140/0000", treeName = "tree", maxN = 1 if args.small else None )

official_eta_binning = [ -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191 ]
fine_eta_binning = [208,-5.2,5.2]
hyperfine_eta_binning = [ -5.2+10.4*i/208. for i in range(209) ] 
prefix = "hyperfine_eta_binning"
eta_binning = hyperfine_eta_binning

##plots wrt to jet reco eta
#for name, var, texY, yRange in [
#    ["response_recoEta", "Jet_rawPt/Jet_mcPt", "p_{T}(raw)/p_{T}(gen)", (0.7,1.2)],
#    ["corrResponse_recoEta", "Jet_pt/Jet_mcPt", "p_{T}(corr)/p_{T}(gen)", (0.9,1.2)],
##    ["deltaEta", "Jet_eta-Jet_mcEta", "#eta(reco)-#eta(gen)", (-0.2,0.2)],
#    ]:
#    deltaEta_profile        = QCD_flat.get1DHistoFromDraw( var+":Jet_eta", fine_eta_binning, "Jet_mcPt>20",                isProfile = True )
#    deltaEta_profile_20_30  = QCD_flat.get1DHistoFromDraw( var+":Jet_eta", fine_eta_binning, "Jet_mcPt>20&&Jet_mcPt<30",   isProfile = True )
#    deltaEta_profile_30_50  = QCD_flat.get1DHistoFromDraw( var+":Jet_eta", fine_eta_binning, "Jet_mcPt>30&&Jet_mcPt<50",   isProfile = True )
#    deltaEta_profile_50_100 = QCD_flat.get1DHistoFromDraw( var+":Jet_eta", fine_eta_binning, "Jet_mcPt>50&&Jet_mcPt<100",  isProfile = True )
#    deltaEta_profile_100    = QCD_flat.get1DHistoFromDraw( var+":Jet_eta", fine_eta_binning, "Jet_mcPt>100",               isProfile = True )
#
#    histos = [ [p.ProjectionX()] for p in [deltaEta_profile, deltaEta_profile_20_30, deltaEta_profile_30_50, deltaEta_profile_50_100, deltaEta_profile_100] ] 
#    histos[0][0].style = styles.lineStyle(ROOT.kBlack)
#    histos[1][0].style = styles.lineStyle(ROOT.kBlue)
#    histos[2][0].style = styles.lineStyle(ROOT.kRed)
#    histos[3][0].style = styles.lineStyle(ROOT.kGreen)
#    histos[4][0].style = styles.lineStyle(ROOT.kMagenta)
#    histos[0][0].legendText = "p_{T} > 20 GeV" 
#    histos[1][0].legendText = "20 < p_{T} < 30 GeV" 
#    histos[2][0].legendText = "30 < p_{T} < 50 GeV" 
#    histos[3][0].legendText = "50 < p_{T} < 100 GeV" 
#    histos[4][0].legendText = "100 < p_{T}" 
#
#    deltaEta_plot = Plot.fromHisto(name = name+'_'+prefix, histos = histos, texX = "#eta(reco)", texY = "texY" )
#    plotting.draw( deltaEta_plot,  plot_directory = os.path.join( plot_directory, args.plot_directory ), ratio = None, logY = False, logX = False, yRange = yRange)

# plots wrt to gen eta
for name, var, texY, yRange in [
    ["response", "Jet_rawPt/Jet_mcPt", "p_{T}(raw)/p_{T}(gen)", (0.7,1.2)],
#    ["corrResponse", "Jet_pt/Jet_mcPt", "p_{T}(corr)/p_{T}(gen)", (0.9,1.2)],
#    ["deltaEta", "Jet_eta-Jet_mcEta", "#eta(reco)-#eta(gen)", (-0.2,0.2)],
    ]:
    deltaEta_profile        = QCD_flat.get1DHistoFromDraw( var+":Jet_mcEta", eta_binning, binningIsExplicit = True, selectionString = "Jet_mcPt>20",                isProfile = True )
    deltaEta_profile_20_30  = QCD_flat.get1DHistoFromDraw( var+":Jet_mcEta", eta_binning, binningIsExplicit = True, selectionString = "Jet_mcPt>20&&Jet_mcPt<30",   isProfile = True )
    deltaEta_profile_30_50  = QCD_flat.get1DHistoFromDraw( var+":Jet_mcEta", eta_binning, binningIsExplicit = True, selectionString = "Jet_mcPt>30&&Jet_mcPt<50",   isProfile = True )
    deltaEta_profile_50_100 = QCD_flat.get1DHistoFromDraw( var+":Jet_mcEta", eta_binning, binningIsExplicit = True, selectionString = "Jet_mcPt>50&&Jet_mcPt<100",  isProfile = True )
    deltaEta_profile_100    = QCD_flat.get1DHistoFromDraw( var+":Jet_mcEta", eta_binning, binningIsExplicit = True, selectionString = "Jet_mcPt>100",               isProfile = True )

    histos = [ [p.ProjectionX()] for p in [deltaEta_profile, deltaEta_profile_20_30, deltaEta_profile_30_50, deltaEta_profile_50_100, deltaEta_profile_100] ]

    if name == "response":
        jet_respons_genEta = {(20,30):deltaEta_profile_20_30.ProjectionX(), (30,50):deltaEta_profile_30_50.ProjectionX(), (50,100):deltaEta_profile_50_100.ProjectionX(), (100,10**6):deltaEta_profile_100.ProjectionX()}

    histos[0][0].style = styles.lineStyle(ROOT.kBlack)
    histos[1][0].style = styles.lineStyle(ROOT.kBlue)
    histos[2][0].style = styles.lineStyle(ROOT.kRed)
    histos[3][0].style = styles.lineStyle(ROOT.kGreen)
    histos[4][0].style = styles.lineStyle(ROOT.kMagenta)
    histos[0][0].legendText = "p_{T} > 20 GeV" 
    histos[1][0].legendText = "20 < p_{T} < 30 GeV" 
    histos[2][0].legendText = "30 < p_{T} < 50 GeV" 
    histos[3][0].legendText = "50 < p_{T} < 100 GeV" 
    histos[4][0].legendText = "100 < p_{T}" 

    deltaEta_plot = Plot.fromHisto(name = name+"_"+prefix, histos = histos, texX = "#eta(gen)", texY = "texY" )
    plotting.draw( deltaEta_plot,  plot_directory = os.path.join( plot_directory, args.plot_directory ), ratio = None, logY = False, logX = False, yRange = yRange)

# Profiles for response ratio of gen-eta and reco-eta
migration_response_ratio         = ROOT.TProfile("ratio", "ratio", len(eta_binning)-1, array.array('d', eta_binning)) 
migration_response_ratio_20_30   = ROOT.TProfile("ratio_20_30", "ratio_20_30", len(eta_binning)-1, array.array('d', eta_binning)) 
migration_response_ratio_30_50   = ROOT.TProfile("ratio_30_50", "ratio_30_50", len(eta_binning)-1, array.array('d', eta_binning)) 
migration_response_ratio_50_100  = ROOT.TProfile("ratio_50_100", "ratio_50_100", len(eta_binning)-1, array.array('d', eta_binning)) 
migration_response_ratio_100     = ROOT.TProfile("ratio_100", "ratio_100", len(eta_binning)-1, array.array('d', eta_binning)) 
# Profiles for assessing bin mugration
migration_fraction         = ROOT.TProfile("fraction", "fraction", len(eta_binning)-1, array.array('d', eta_binning)) 
migration_fraction_20_30   = ROOT.TProfile("fraction_20_30", "fraction_20_30", len(eta_binning)-1, array.array('d', eta_binning)) 
migration_fraction_30_50   = ROOT.TProfile("fraction_30_50", "fraction_30_50", len(eta_binning)-1, array.array('d', eta_binning)) 
migration_fraction_50_100  = ROOT.TProfile("fraction_50_100", "fraction_50_100", len(eta_binning)-1, array.array('d', eta_binning)) 
migration_fraction_100     = ROOT.TProfile("fraction_100", "fraction_100", len(eta_binning)-1, array.array('d', eta_binning)) 

pt_bins = jet_respons_genEta.keys()
def getPtBin( pt ):
    for ptb in pt_bins:
        if pt>=ptb[0] and (ptb[1]<0 or pt<ptb[1]):
            return ptb
    
r = QCD_flat.treeReader( variables = map( TreeVariable.fromString, ["Jet[mcPt/F,rawPt/F,mcEta/F,eta/F]", 'nJet/I'] ) )
r.start()
while r.run():
    jets = getJets( r.event, jetColl="Jet", jetVars = ['mcPt', 'rawPt', 'mcEta', 'eta'] )
    for i in range( r.event.nJet ):

        # find bins
        reco_bin = migration_response_ratio.FindBin( r.event.Jet_eta[i] )
        gen_bin  = migration_response_ratio.FindBin( r.event.Jet_mcEta[i] )
        pt_bin = getPtBin( r.event.Jet_mcPt[i] )

        if not pt_bin: continue

        # compute response ratios
        reco_response = jet_respons_genEta[pt_bin].GetBinContent( reco_bin )
        gen_response  = jet_respons_genEta[pt_bin].GetBinContent( gen_bin )
        try:
            ratio = reco_response/gen_response 
        except ZeroDivisionError:
            continue

        # bin migration?
        mig = (reco_bin == gen_bin)

        # fill
        migration_response_ratio.Fill( r.event.Jet_mcEta[i], ratio )
        migration_fraction.Fill( r.event.Jet_mcEta[i], mig  )
        if r.event.Jet_mcPt[i]>20 and r.event.Jet_mcPt[i]<30:
            migration_response_ratio_20_30.Fill( r.event.Jet_mcEta[i], ratio )
            migration_fraction_20_30.Fill( r.event.Jet_mcEta[i], mig )
        elif r.event.Jet_mcPt[i]>30 and r.event.Jet_mcPt[i]<50:
            migration_response_ratio_30_50.Fill( r.event.Jet_mcEta[i], ratio )
            migration_fraction_30_50.Fill( r.event.Jet_mcEta[i], mig )
        elif r.event.Jet_mcPt[i]>50 and r.event.Jet_mcPt[i]<100:
            migration_response_ratio_50_100.Fill( r.event.Jet_mcEta[i], ratio )
            migration_fraction_50_100.Fill( r.event.Jet_mcEta[i], mig )
        elif r.event.Jet_mcPt[i]>100:
            migration_response_ratio_100.Fill( r.event.Jet_mcEta[i], ratio )
            migration_fraction_100.Fill( r.event.Jet_mcEta[i], mig )

#        print i, r.event.Jet_mcEta[i], reco_response/gen_response

histos = [[h.ProjectionX()] for h in [migration_response_ratio, migration_response_ratio_20_30, migration_response_ratio_30_50, migration_response_ratio_50_100, migration_response_ratio_100]]

histos[0][0].style = styles.lineStyle(ROOT.kBlack)
histos[1][0].style = styles.lineStyle(ROOT.kBlue)
histos[2][0].style = styles.lineStyle(ROOT.kRed)
histos[3][0].style = styles.lineStyle(ROOT.kGreen)
histos[4][0].style = styles.lineStyle(ROOT.kMagenta)
histos[0][0].legendText = "p_{T} > 20 GeV" 
histos[1][0].legendText = "20 < p_{T} < 30 GeV" 
histos[2][0].legendText = "30 < p_{T} < 50 GeV" 
histos[3][0].legendText = "50 < p_{T} < 100 GeV" 
histos[4][0].legendText = "100 < p_{T}" 

migration_response_ratio_plot = Plot.fromHisto(name = "migration_response_ratio_"+prefix, histos = histos, texX = "#eta(gen)", texY = "response ratio" )
plotting.draw( migration_response_ratio_plot,  plot_directory = os.path.join( plot_directory, args.plot_directory ), ratio = None, logY = False, logX = False, yRange = (0.97,1.07))

histos = [[h.ProjectionX()] for h in [migration_fraction, migration_fraction_20_30, migration_fraction_30_50, migration_fraction_50_100, migration_fraction_100]]

histos[0][0].style = styles.lineStyle(ROOT.kBlack)
histos[1][0].style = styles.lineStyle(ROOT.kBlue)
histos[2][0].style = styles.lineStyle(ROOT.kRed)
histos[3][0].style = styles.lineStyle(ROOT.kGreen)
histos[4][0].style = styles.lineStyle(ROOT.kMagenta)
histos[0][0].legendText = "p_{T} > 20 GeV" 
histos[1][0].legendText = "20 < p_{T} < 30 GeV" 
histos[2][0].legendText = "30 < p_{T} < 50 GeV" 
histos[3][0].legendText = "50 < p_{T} < 100 GeV" 
histos[4][0].legendText = "100 < p_{T}" 

migration_response_ratio_plot = Plot.fromHisto(name = "migration_fraction_"+prefix, histos = histos, texX = "#eta(gen)", texY = "fraction of non-migrating jets" )
plotting.draw( migration_response_ratio_plot,  plot_directory = os.path.join( plot_directory, args.plot_directory ), ratio = None, logY = False, logX = False, yRange = (0, 1.4))

y_binning_explicit = [-0.4+(0.8)/100*i for i in range(0,101)]
deltaEta_histo2D        = QCD_flat.get2DHistoFromDraw( "Jet_eta-Jet_mcEta:Jet_mcEta", (eta_binning, y_binning_explicit), binningIsExplicit = True, selectionString = "Jet_mcPt>20")

quantiles = range(10,100,10)
n_quantiles = len(quantiles)

h_quantiles = { q : ROOT.TH1D("q_%i"%q, "q_%i"%q, len(eta_binning)-1, array.array('d', eta_binning)) for q in range(10,100,10) }
s=0
h_binsize_plus  = ROOT.TH1D("plus", "plus", len(eta_binning)-1, array.array('d', eta_binning)) 
h_binsize_minus = ROOT.TH1D("minus", "minus", len(eta_binning)-1, array.array('d', eta_binning)) 

for i in range(len(eta_binning)):
    h_binsize_plus.SetBinContent(i+1, h_binsize_plus.GetBinWidth(i+1))
    h_binsize_minus.SetBinContent(i+1, -h_binsize_minus.GetBinWidth(i+1))
    
    proj_y = deltaEta_histo2D.ProjectionY("proj_%i"%i, i+1,i+1)
    s+=proj_y.Integral()
    result = array.array('d', [ROOT.Double()] * n_quantiles )  
    proj_y.GetQuantiles( n_quantiles, result, array.array('d', [q/100. for q in quantiles]) )
    #print i, deltaEta_histo2D.GetXaxis().GetBinLowEdge(i+1), proj_y.Integral(), proj_y.GetMean(), result
    for iq, q in enumerate(quantiles):
        h_quantiles[q].SetBinContent( i+1, result[iq] )

h_quantiles[10].style = styles.lineStyle(ROOT.kMagenta)
h_quantiles[20].style = styles.lineStyle(ROOT.kGreen)
h_quantiles[30].style = styles.lineStyle(ROOT.kRed)
h_quantiles[40].style = styles.lineStyle(ROOT.kBlue)
h_quantiles[50].style = styles.lineStyle(ROOT.kMagenta)
h_quantiles[60].style = styles.lineStyle(ROOT.kBlue)
h_quantiles[70].style = styles.lineStyle(ROOT.kRed)
h_quantiles[80].style = styles.lineStyle(ROOT.kGreen)
h_quantiles[90].style = styles.lineStyle(ROOT.kMagenta)

histos = [ [h_quantiles[q]] for q in quantiles ] + [[h_binsize_plus] ,[h_binsize_minus] ] 

deltaEta_quantile_plot = Plot.fromHisto(name = "deltaEta_quantile_"+prefix, histos = histos, texX = "#eta(gen)", texY = "#eta(reco)-#eta(gen)" )
plotting.draw( deltaEta_quantile_plot,  plot_directory = os.path.join( plot_directory, args.plot_directory ), ratio = None, logY = False, logX = False, yRange = (-0.2,0.2))
