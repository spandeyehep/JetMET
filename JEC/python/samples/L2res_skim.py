from RootTools.core.standard import *

from JetMET.tools.user import skim_ntuple_directory

JetHT_Run2016B    = Sample.fromDirectory("JetHT_Run2016B", skim_ntuple_directory+"/JetHT_Run2016B")
JetHT_Run2016C    = Sample.fromDirectory("JetHT_Run2016C", skim_ntuple_directory+"/JetHT_Run2016C")
JetHT_Run2016D    = Sample.fromDirectory("JetHT_Run2016D", skim_ntuple_directory+"/JetHT_Run2016D")
JetHT_Run2016E    = Sample.fromDirectory("JetHT_Run2016E", skim_ntuple_directory+"/JetHT_Run2016E")
JetHT_Run2016F    = Sample.fromDirectory("JetHT_Run2016F", skim_ntuple_directory+"/JetHT_Run2016F")
JetHT_Run2016G    = Sample.fromDirectory("JetHT_Run2016G", skim_ntuple_directory+"/JetHT_Run2016G")
JetHT_Run2016H_v2 = Sample.fromDirectory("JetHT_Run2016H_v2", skim_ntuple_directory+"/JetHT_Run2016H_v2")
JetHT_Run2016H_v3 = Sample.fromDirectory("JetHT_Run2016H_v3", skim_ntuple_directory+"/JetHT_Run2016H_v3")

JetHT_Run2016BCD     = Sample.combine("JetHT_Run2016BCD", [JetHT_Run2016B, JetHT_Run2016C, JetHT_Run2016D])
JetHT_Run2016EFearly = Sample.combine("JetHT_Run2016EFearly", [JetHT_Run2016E, JetHT_Run2016F])
JetHT_Run2016EFearly.setSelectionString("run<=278801")
JetHT_Run2016FlateG  = Sample.combine("JetHT_Run2016FLateG", [JetHT_Run2016F, JetHT_Run2016G])
JetHT_Run2016FlateG.setSelectionString("run>=278802")
JetHT_Run2016H       = Sample.combine("JetHT_Run2016H", [JetHT_Run2016H_v2, JetHT_Run2016H_v3])

JetHT_Run2016        = Sample.combine("JetHT_Run2016", [JetHT_Run2016B, JetHT_Run2016C, JetHT_Run2016D, JetHT_Run2016E, JetHT_Run2016F, JetHT_Run2016G, JetHT_Run2016H_v2, JetHT_Run2016H_v3] )

QCD_Pt_50to80     = Sample.fromDirectory("QCD_Pt_50to80", skim_ntuple_directory+"/QCD_Pt_50to80")
QCD_Pt_80to120    = Sample.fromDirectory("QCD_Pt_80to120", skim_ntuple_directory+"/QCD_Pt_80to120")
QCD_Pt_120to170   = Sample.fromDirectory("QCD_Pt_120to170", skim_ntuple_directory+"/QCD_Pt_120to170")
QCD_Pt_170to300   = Sample.fromDirectory("QCD_Pt_170to300", skim_ntuple_directory+"/QCD_Pt_170to300")
QCD_Pt_300to470   = Sample.fromDirectory("QCD_Pt_300to470", skim_ntuple_directory+"/QCD_Pt_300to470")
QCD_Pt_470to600   = Sample.fromDirectory("QCD_Pt_470to600", skim_ntuple_directory+"/QCD_Pt_470to600")
QCD_Pt_600to800   = Sample.fromDirectory("QCD_Pt_600to800", skim_ntuple_directory+"/QCD_Pt_600to800")
QCD_Pt_800to1000  = Sample.fromDirectory("QCD_Pt_800to1000", skim_ntuple_directory+"/QCD_Pt_800to1000")
QCD_Pt_1000to1400 = Sample.fromDirectory("QCD_Pt_1000to1400", skim_ntuple_directory+"/QCD_Pt_1000to1400")
QCD_Pt_1400to1800 = Sample.fromDirectory("QCD_Pt_1400to1800", skim_ntuple_directory+"/QCD_Pt_1400to1800")
QCD_Pt_1800to2400 = Sample.fromDirectory("QCD_Pt_1800to2400", skim_ntuple_directory+"/QCD_Pt_1800to2400")
QCD_Pt_2400to3200 = Sample.fromDirectory("QCD_Pt_2400to3200", skim_ntuple_directory+"/QCD_Pt_2400to3200")
QCD_Pt_3200toInf  = Sample.fromDirectory("QCD_Pt_3200toInf", skim_ntuple_directory+"/QCD_Pt_3200toInf")

QCD_Pt = Sample.combine( "QCD_Pt", [ QCD_Pt_50to80, QCD_Pt_80to120, QCD_Pt_120to170, QCD_Pt_170to300, QCD_Pt_300to470, QCD_Pt_470to600, QCD_Pt_600to800, QCD_Pt_800to1000, QCD_Pt_1000to1400, QCD_Pt_1400to1800, QCD_Pt_1800to2400, QCD_Pt_2400to3200, QCD_Pt_3200toInf] )
