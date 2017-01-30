import ROOT
from RootTools.core.standard import *

maxN = -1
noPU = Sample.fromDPMDirectory( name = "QCD_flat_noPU", directory = "/dpm/oeaw.ac.at/home/cms/store/user/schoef/cmgTuples/80X_JME/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISummer16MiniAODv2-PUFlat0to70_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v4-v1_80X_JME/170121_182640/0000", maxN = maxN, treeName = "tree")

central_jet = "Sum$(Jet_pt>100&&Jet_pt<120&&Jet_id>=3&&abs(Jet_eta)<1.3)>0"
endcap_jet  = "Sum$(Jet_pt>100&&Jet_pt<120&&Jet_id>=3&&abs(Jet_eta)>2)>0"

binning_eta = [20, -5.2, 5.2]
binning_pt = [100,0,50]
h_pt_vs_eta = noPU.get2DHistoFromDraw( "pf_pt:pf_eta", binning_eta + binning_pt, "&&".join( [ central_jet, endcap_jet ] ) ) 
c1 = ROOT.TCanvas()
h_pt_vs_eta.Draw("COLZ")
c1.SetLogz()
c1.Print("/afs/hephy.at/user/r/rschoefbeck/www/etc/h_pt_vs_eta.png")

