# RootTools
from RootTools.core.Sample import Sample

jets = Sample.fromDirectory("jets", "/data/rschoefbeck/JetMET/localJEC/jetTrees/", treeName = "jets")

maxN = -1
#qcd_AllChGood_Flat0to50 = Sample.fromCMGOutput("qcd_AllChGood_Flat0to50", "/data/rschoefbeck/cmgTuples/MC25ns_v2_0l/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_schoef-crab_QCD_Pt-15to7000_AllChlgoodAsymptFlat0to50bx25_AllChannelsGood_v0-v2_RunIISpring15MiniAODv2-74X", treeFilename='tree.root', maxN=maxN)
#qcd_Flat0to50 = Sample.fromCMGOutput("qcd_Flat0to50", "/data/rschoefbeck/cmgTuples/MC25ns_v2_0l/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_schoef-crab_QCD_Pt-15to7000_AsymptFlat0to50bx25Reco_MCRUN2_74_V9-v3_RunIISpring15MiniAODv2-74X_3", treeFilename='tree.root', maxN=maxN)
qcd_AllChGood_noPU = Sample.fromCMGOutput("qcd_AllChGood_noPU", "/data/rschoefbeck/cmgTuples/MC25ns_v2_0l/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_schoef-crab_QCD_Pt-15to7000_AllChlgoodAsymptNoPUbx25_AllChannelsGood_v0-v2_RunIISpring15MiniAODv2-74X", treeFilename='tree.root', maxN=maxN)
qcd_noPU = Sample.fromCMGOutput("qcd_noPU", "/data/rschoefbeck/cmgTuples/MC25ns_v2_0l/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_schoef-crab_QCD_Pt-15to7000_AsymptNoPUbx25Reco_MCRUN2_74_V9-v3_RunIISpring15MiniAODv2-74X_3", treeFilename='tree.root', maxN=maxN)

#Flag_HBHENoiseIsoFilter
#Flag_EcalDeadCellTriggerPrimitiveFilter
#Flag_trkPOG_manystripclus53X
#Flag_ecalLaserCorrFilter
#Flag_trkPOG_toomanystripclus53X
#Flag_hcalLaserEventFilter
#Flag_trkPOG_logErrorTooManyClusters
#Flag_trkPOGFilters
#Flag_trackingFailureFilter
#Flag_CSCTightHaloFilter
#Flag_HBHENoiseFilter
#Flag_goodVertices
#Flag_METFilters
#Flag_eeBadScFilter

#qcd_noPU.chain.Draw("met_rawPt>>h(100,0,2000)")
#qcd_AllChGood_noPU.chain.Draw("met_rawPt","","same")
#qcd_noPU.chain.Draw("met_genPt","","same")
#qcd_AllChGood_noPU.chain.Draw("met_genPt","","same")

#qcd_noPU.chain.Draw("met_rawPt>>h(100,0,2000)", "met_genPt<50&&Flag_EcalDeadCellTriggerPrimitiveFilter")
#qcd_AllChGood_noPU.chain.Draw("met_rawPt","met_genPt<50&&Flag_EcalDeadCellTriggerPrimitiveFilter","same")
#qcd_noPU.chain.Draw("met_genPt","met_genPt<50&&Flag_EcalDeadCellTriggerPrimitiveFilter","same")
#qcd_AllChGood_noPU.chain.Draw("met_genPt","met_genPt<50&&Flag_EcalDeadCellTriggerPrimitiveFilter","same")

