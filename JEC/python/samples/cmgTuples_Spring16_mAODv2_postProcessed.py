import copy, os, sys
from RootTools.core.Sample import Sample
import ROOT

from StopsDilepton.samples.color import color

# Data directory
try:
    data_directory = sys.modules['__main__'].data_directory
except:
    from StopsDilepton.tools.user import data_directory as user_data_directory
    data_directory = user_data_directory 

# Take post processing directory if defined in main module
try:
  import sys
  postProcessing_directory = sys.modules['__main__'].postProcessing_directory
except:
  postProcessing_directory = "postProcessed_80X_v23/dilepTiny/"

DY_M5to50_HT = ["DYJetsToLL_M5to50_LO_lheHT100", 
    "DYJetsToLL_M5to50_HT100to200_comb", "DYJetsToLL_M5to50_HT200to400_comb", "DYJetsToLL_M5to50_HT400to600", "DYJetsToLL_M5to50_HT600toInf_comb"] 
DY_M50_HT = ["DYJetsToLL_M50_LO_lheHT100", 
    "DYJetsToLL_M50_HT100to200_comb", "DYJetsToLL_M50_HT200to400_comb", "DYJetsToLL_M50_HT400to600_ext", "DYJetsToLL_M50_HT600toInf_comb"] 


dirs = {}
dirs['DY']               = ["DYJetsToLL_M50", "DYJetsToLL_M10to50" ]
#dirs['DY_LO']            = ["DYJetsToLL_M10to50_LO", "DYJetsToLL_M50_LO"]
#dirs['DY_HT_LO']         =  DY_M50_HT + DY_M5to50_HT
dirs['TTJets']           = ["TTJets"]
dirs['TTLep_pow']        = ["TTLep_pow_ext"]
dirs['TT_pow']           = ["TT_pow_ext3_comb"]

#dirs['TTJets_Lep']       = ["TTJets_DiLepton_comb", "TTJets_SingleLeptonFromTbar_comb", "TTJets_SingleLeptonFromT_comb"]
#dirs['TTJets_Dilep']      = ["TTJets_DiLepton_comb"]
#dirs['TTJets_Singlelep']  = ["TTJets_SingleLeptonFromTbar_comb", "TTJets_SingleLeptonFromT_comb"]
dirs['singleTop']        = ["TBar_tWch", "T_tWch", "TToLeptons_tch_powheg", "TBarToLeptons_tch_powheg"]
dirs['singleTop_tch']    = ["TToLeptons_tch_powheg", "TBarToLeptons_tch_powheg"]
dirs['Top_amc']          = dirs['singleTop'] + dirs['TTJets']
#dirs['Top']              = dirs['singleTop'] + dirs['TTJets_Lep']
dirs['Top_pow']          = dirs['TTLep_pow'] + dirs['singleTop']
dirs['Top_pow_incl']     = dirs['TT_pow'] + dirs['singleTop']
dirs['TZQ']              = ["tZq_ll", "tZq_nunu_reHLT"]
dirs['TWZ']              = ["tWll", "tWnunu"]
dirs['TTW']              = ["TTWToLNu", "TTWToQQ"]
dirs['TTH']              = [ \
        "TTHbb_ext3", 
        "TTHnobb_mWCutfix_ch0"]
dirs['TTZtoLLNuNu']      = ["TTZToLLNuNu"]
dirs['TTZtoQQ']          = ["TTZToQQ"]
dirs['TTZ']              = ["TTZToLLNuNu", "TTZToQQ"]
dirs['TTZ_LO']           = ["TTZ_LO"]
dirs['TTXNoZ']           = dirs['TTH'] + dirs['TTW'] + dirs['TZQ'] + dirs['TWZ']
dirs['TTX']              = dirs['TTXNoZ'] + dirs['TTZ_LO']
dirs['WJetsToLNu']       = ["WJetsToLNu"]
dirs['WJetsToLNu_HT']    = ["WJetsToLNu_HT100to200_comb", "WJetsToLNu_HT200to400_comb", "WJetsToLNu_HT400to600", "WJetsToLNu_HT600to800", "WJetsToLNu_HT800to1200", "WJetsToLNu_HT1200to2500", "WJetsToLNu_HT2500toInf"]
dirs['diBosonInclusive'] = ["WW", "WZ", "ZZ"]
dirs['WW']               = ["WWTo1L1Nu2Q", "WWToLNuQQ_comb"]
dirs['WW_']              = ["WWTo1L1Nu2Q", "WWToLNuQQ_comb","WWTo2L2Nu"]
dirs['WWTo2L2Nu']        = ["WWTo2L2Nu"]
dirs['VVTo2L2Nu']        = ["VVTo2L2Nu"]
dirs['WZ']               = ["WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2L2Q", "WZTo3LNu"]
dirs['ZZ']               = ["ZZTo2L2Q", "ZZTo2Q2Nu"]
dirs['ZZTo2L2Nu']        = ["ZZTo2L2Q"]
dirs['ZZ_']              = ["ZZTo2L2Q", "ZZTo2Q2Nu","ZZTo2L2Nu"]
dirs['diBoson']          = dirs['WW'] + dirs['WZ'] + dirs['ZZ'] + dirs['VVTo2L2Nu']
dirs['diBoson_']         = dirs['WW_'] + dirs['WZ'] + dirs['ZZ_']
dirs['triBoson']         = ["WWZ","WZZ","ZZZ"] 
dirs['multiBoson']       = dirs['diBoson'] + dirs['triBoson']
dirs['EWK']              = dirs['diBoson'] + dirs['triBoson'] + dirs['TTX']
#dirs['QCD_HT']           = ["QCD_HT100to200", "QCD_HT200to300", "QCD_HT300to500", "QCD_HT500to700", "QCD_HT700to1000", "QCD_HT1000to1500", "QCD_HT1500to2000", "QCD_HT2000toInf"]
#dirs['QCD_Mu5']          = ["QCD_Pt20to30_Mu5", "QCD_Pt50to80_Mu5", "QCD_Pt80to120_Mu5", "QCD_Pt120to170_Mu5", "QCD_Pt170to300_Mu5", "QCD_Pt300to470_Mu5", "QCD_Pt470to600_Mu5", "QCD_Pt600to800_Mu5", "QCD_Pt800to1000_Mu5", "QCD_Pt1000toInf_Mu5"]
#dirs['QCD_EM+bcToE']     = ["QCD_Pt_20to30_bcToE", "QCD_Pt_30to80_bcToE", "QCD_Pt_80to170_bcToE", "QCD_Pt_170to250_bcToE", "QCD_Pt_250toInf_bcToE", "QCD_Pt15to20_EMEnriched", "QCD_Pt20to30_EMEnriched", "QCD_Pt50to80_EMEnriched", "QCD_Pt80to120_EMEnriched", "QCD_Pt120to170_EMEnriched", "QCD_Pt170to300_EMEnriched"]
#dirs['QCD_EM+bcToE']    = ["QCD_Pt_15to20_bcToE", "QCD_Pt_20to30_bcToE", "QCD_Pt_30to80_bcToE", "QCD_Pt_80to170_bcToE", "QCD_Pt_170to250_bcToE", "QCD_Pt_250toInf_bcToE", "QCD_Pt15to20_EMEnriched", "QCD_Pt20to30_EMEnriched", "QCD_Pt30to50_EMEnriched", "QCD_Pt50to80_EMEnriched", "QCD_Pt80to120_EMEnriched", "QCD_Pt120to170_EMEnriched", "QCD_Pt170to300_EMEnriched", "QCD_Pt300toInf_EMEnriched"]
#dirs['QCD_Mu5+EM+bcToE'] = ["QCD_Pt20to30_Mu5", "QCD_Pt50to80_Mu5", "QCD_Pt80to120_Mu5", "QCD_Pt120to170_Mu5", "QCD_Pt170to300_Mu5", "QCD_Pt300to470_Mu5", "QCD_Pt470to600_Mu5", "QCD_Pt600to800_Mu5", "QCD_Pt800to1000_Mu5", "QCD_Pt1000toInf_Mu5", "QCD_Pt_20to30_bcToE", "QCD_Pt_30to80_bcToE", "QCD_Pt_80to170_bcToE", "QCD_Pt_170to250_bcToE", "QCD_Pt_250toInf_bcToE", "QCD_Pt15to20_EMEnriched", "QCD_Pt20to30_EMEnriched", "QCD_Pt50to80_EMEnriched", "QCD_Pt80to120_EMEnriched", "QCD_Pt120to170_EMEnriched", "QCD_Pt170to300_EMEnriched"]
#dirs['QCD_Mu5+EM+bcToE']= ["QCD_Pt20to30_Mu5", "QCD_Pt50to80_Mu5", "QCD_Pt80to120_Mu5", "QCD_Pt120to170_Mu5", "QCD_Pt170to300_Mu5", "QCD_Pt300to470_Mu5", "QCD_Pt470to600_Mu5", "QCD_Pt600to800_Mu5", "QCD_Pt800to1000_Mu5", "QCD_Pt1000toInf_Mu5", "QCD_Pt_15to20_bcToE", "QCD_Pt_20to30_bcToE", "QCD_Pt_30to80_bcToE", "QCD_Pt_80to170_bcToE", "QCD_Pt_170to250_bcToE", "QCD_Pt_250toInf_bcToE", "QCD_Pt15to20_EMEnriched", "QCD_Pt20to30_EMEnriched", "QCD_Pt30to50_EMEnriched", "QCD_Pt50to80_EMEnriched", "QCD_Pt80to120_EMEnriched", "QCD_Pt120to170_EMEnriched", "QCD_Pt170to300_EMEnriched", "QCD_Pt300toInf_EMEnriched"]
#dirs['QCD']              = ["QCD_Pt30to50", "QCD_Pt50to80", "QCD_Pt80to120", "QCD_Pt120to170", "QCD_Pt170to300", "QCD_Pt300to470", "QCD_Pt470to600", "QCD_Pt600to800", "QCD_Pt800to1000", "QCD_Pt1000to1400", "QCD_Pt1400to1800", "QCD_Pt1800to2400", "QCD_Pt2400to3200"]
#dirs['QCD']             = ["QCD_Pt10to15", "QCD_Pt15to30", "QCD_Pt30to50", "QCD_Pt50to80", "QCD_Pt80to120", "QCD_Pt120to170", "QCD_Pt170to300", "QCD_Pt300to470", "QCD_Pt470to600", "QCD_Pt600to800", "QCD_Pt800to1000", "QCD_Pt1000to1400", "QCD_Pt1400to1800", "QCD_Pt1800to2400", "QCD_Pt2400to3200", "QCD_Pt3200"]

#dirs['WGToLNuG']     = ["WGToLNuG"]
dirs['ZGTo2LG']      = ["ZGTo2LG"]
#dirs['WGJets']       = ["WGJets"]
dirs['ZGJets']       = ["ZGJets"]
dirs['ZG']           = dirs['ZGTo2LG'] + dirs['ZGJets']

directories = { key : [ os.path.join( data_directory, postProcessing_directory, dir) for dir in dirs[key]] for key in dirs.keys()}

DY              = Sample.fromDirectory(name="DY",               treeName="Events", isData=False, color=color.DY,              texName="DY",                                directory=directories['DY'])
#DY_LO           = Sample.fromDirectory(name="DY_LO",            treeName="Events", isData=False, color=color.DY,              texName="DY (LO)",                           directory=directories['DY_LO'])
#DY_HT_LO        = Sample.fromDirectory(name="DY_HT_LO",         treeName="Events", isData=False, color=color.DY,              texName="DY (LO,HT)",                        directory=directories['DY_HT_LO'])
TTJets          = Sample.fromDirectory(name="TTJets",           treeName="Events", isData=False, color=color.TTJets,          texName="t#bar{t}",                          directory=directories['TTJets'])
#TTJets_Lep      = Sample.fromDirectory(name="TTJets_Lep",       treeName="Events", isData=False, color=color.TTJets,          texName="t#bar{t}(lep)",                     directory=directories['TTJets_Lep'])
#TTJets_Singlelep= Sample.fromDirectory(name="TTJets_Singlelep", treeName="Events", isData=False, color=color.TTJets,          texName="t#bar{t}",                          directory=directories['TTJets_Singlelep'])
#TTJets_Singlelep_singleTop = Sample.fromDirectory(name="TTJets_Singlelep_singleTop", treeName="Events", isData=False, color=color.TTJets, texName="t#bar{t}/t (1l)",       directory=directories['TTJets_Singlelep']+directories['singleTop'])
#TTJets_Dilep    = Sample.fromDirectory(name="TTJets_Dilep",     treeName="Events", isData=False, color=color.TTJets,          texName="t#bar{t}",                          directory=directories['TTJets_Dilep'])
#TTJets_LO       = Sample.fromDirectory(name="TTJets_LO",        treeName="Events", isData=False, color=color.TTJets,          texName="t#bar{t} + Jets (LO)",              directory=directories['TTJets_LO'])
TTLep_pow       = Sample.fromDirectory(name="TTLep_pow",        treeName="Events", isData=False, color=color.TTJets,          texName="t#bar{t} + Jets (lep,pow)",         directory=directories['TTLep_pow'])
TT_pow          = Sample.fromDirectory(name="TT_pow",        treeName="Events", isData=False, color=color.TTJets,          texName="t#bar{t} + Jets (powheg)",         directory=directories['TT_pow'])
#Top             = Sample.fromDirectory(name="Top",              treeName="Events", isData=False, color=color.TTJets,          texName="t#bar{t}/single-t",                 directory=directories['Top'])
Top_pow         = Sample.fromDirectory(name="Top_pow",          treeName="Events", isData=False, color=color.TTJets,          texName="t#bar{t}/single-t(powHeg)",         directory=directories['Top_pow'])
Top_pow_incl    = Sample.fromDirectory(name="Top_pow_incl",     treeName="Events", isData=False, color=color.TTJets,          texName="t#bar{t}/single-t(pow incl)",       directory=directories['Top_pow_incl'])
Top_amc         = Sample.fromDirectory(name="Top_amc",          treeName="Events", isData=False, color=color.TTJets,          texName="t#bar{t}/single-t(amc@NLO)",        directory=directories['Top_amc'])


singleTop      = Sample.fromDirectory(name="singleTop",        treeName="Events", isData=False, color=color.singleTop,       texName="single top",                        directory=directories['singleTop'])
singleTop_tch  = Sample.fromDirectory(name="singleTop_tch",    treeName="Events", isData=False, color=color.singleTop,       texName="single top tch",                    directory=directories['singleTop_tch'])
TTX            = Sample.fromDirectory(name="TTX",              treeName="Events", isData=False, color=color.TTX,             texName="t#bar{t}H/W/Z, tZq",                directory=directories['TTX'])
TTXNoZ         = Sample.fromDirectory(name="TTXNoZ",           treeName="Events", isData=False, color=color.TTXNoZ,          texName="t#bar{t}H/W, tZq, tWZ",             directory=directories['TTXNoZ'])
TTH            = Sample.fromDirectory(name="TTH",              treeName="Events", isData=False, color=color.TTH,             texName="t#bar{t}H",                         directory=directories['TTH'])
TTW            = Sample.fromDirectory(name="TTW",              treeName="Events", isData=False, color=color.TTW,             texName="t#bar{t}W",                         directory=directories['TTW'])
TTZ            = Sample.fromDirectory(name="TTZ",              treeName="Events", isData=False, color=color.TTZ,             texName="t#bar{t}Z",                         directory=directories['TTZ'])
TTZ_LO        = Sample.fromDirectory(name="TTZ_LO",              treeName="Events", isData=False, color=color.TTZ,             texName="t#bar{t}Z",                         directory=directories['TTZ_LO'])
TTZtoLLNuNu    = Sample.fromDirectory(name="TTZtoNuNu",        treeName="Events", isData=False, color=color.TTZtoLLNuNu,     texName="t#bar{t}Z (l#bar{l}/#nu#bar{#nu})", directory=directories['TTZtoLLNuNu'])
TTZtoQQ        = Sample.fromDirectory(name="TTZtoQQ",          treeName="Events", isData=False, color=color.TTZtoQQ,         texName="t#bar{t}Z (q#bar{q})",              directory=directories['TTZtoQQ'])
TZQ            = Sample.fromDirectory(name="TZQ",              treeName="Events", isData=False, color=color.TZQ,             texName="tZq",                               directory=directories['TZQ'])
TWZ            = Sample.fromDirectory(name="TWZ",              treeName="Events", isData=False, color=color.TWZ,             texName="tWZ",                               directory=directories['TWZ'])
WJetsToLNu     = Sample.fromDirectory(name="WJetsToLNu",       treeName="Events", isData=False, color=color.WJetsToLNu,      texName="W(l,#nu) + Jets",                   directory=directories['WJetsToLNu'])
diBoson        = Sample.fromDirectory(name="diBoson",          treeName="Events", isData=False, color=color.diBoson,         texName="VV (excl.)",                        directory=directories['diBoson'])
diBoson_       = Sample.fromDirectory(name="diBoson",          treeName="Events", isData=False, color=color.diBoson,         texName="VV (excl.)",                        directory=directories['diBoson_'])
diBosonInclusive = Sample.fromDirectory(name="diBosonInclusive",treeName="Events", isData=False, color=color.diBoson,        texName="VV (incl.)",                        directory=directories['diBosonInclusive'])
ZZ             = Sample.fromDirectory(name="ZZ",               treeName="Events", isData=False, color=color.ZZ,              texName="ZZ",                                directory=directories['ZZ_'])
ZZNo2L2Nu      = Sample.fromDirectory(name="ZZNo2L2Nu",        treeName="Events", isData=False, color=color.ZZ,              texName="ZZ (no 2L2Nu)",                     directory=directories['ZZ'])
ZZTo2L2Nu      = Sample.fromDirectory(name="ZZTo2L2Nu",        treeName="Events", isData=False, color=color.ZZ,              texName="ZZTo2l2Nu",                         directory=directories['ZZTo2L2Nu'])
WZ             = Sample.fromDirectory(name="WZ",               treeName="Events", isData=False, color=color.WZ,              texName="WZ",                                directory=directories['WZ'])
WW             = Sample.fromDirectory(name="WW",               treeName="Events", isData=False, color=color.WW,              texName="WW",                                directory=directories['WW_'])
WWNo2L2Nu      = Sample.fromDirectory(name="WWNo2L2Nu",        treeName="Events", isData=False, color=color.WW,              texName="WW (no 2L2Nu)",                     directory=directories['WW'])
WWTo2L2Nu      = Sample.fromDirectory(name="WWTo2L2Nu",        treeName="Events", isData=False, color=color.WW,              texName="WWTo2L2Nu",                         directory=directories['WWTo2L2Nu'])
VVTo2L2Nu      = Sample.fromDirectory(name="VVTo2L2Nu",               treeName="Events", isData=False, color=color.VV,              texName="VV to ll#nu#nu",             directory=directories['VVTo2L2Nu'])
triBoson       = Sample.fromDirectory(name="triBoson",         treeName="Events", isData=False, color=color.triBoson,        texName="WWZ,WZZ,ZZZ",                       directory=directories['triBoson'])
multiBoson     = Sample.fromDirectory(name="multiBoson",       treeName="Events", isData=False, color=color.diBoson,         texName="multi boson",                       directory=directories['multiBoson'])
#QCD_HT         = Sample.fromDirectory(name="QCD_HT",           treeName="Events", isData=False, color=color.QCD,             texName="QCD (HT)",                          directory=directories['QCD_HT'])
#QCD_Mu5        = Sample.fromDirectory(name="QCD_Mu5",          treeName="Events", isData=False, color=color.QCD,             texName="QCD (Mu5)",                         directory=directories['QCD_Mu5'])
#QCD_EMbcToE    = Sample.fromDirectory(name="QCD_EM+bcToE",     treeName="Events", isData=False, color=color.QCD,             texName="QCD (Em+bcToE)",                    directory=directories['QCD_EM+bcToE'])
#QCD_Mu5EMbcToE = Sample.fromDirectory(name="QCD_Mu5+EM+bcToE", treeName="Events", isData=False, color=color.QCD,             texName="QCD (Mu5+Em+bcToE)",                directory=directories['QCD_Mu5+EM+bcToE'])
#QCD_Pt         = Sample.fromDirectory(name="QCD",              treeName="Events", isData=False, color=color.QCD,             texName="QCD",                               directory=directories['QCD'])

#WGToLNuG = Sample.fromDirectory(name="WGToLNuG",       treeName="Events", isData=False, color=color.diBoson,       texName="WGToLNuG",                       directory=directories['WGToLNuG'])
#ZGTo2LG  = Sample.fromDirectory(name="ZGTo2LG",       treeName="Events", isData=False, color=color.diBoson,        texName="ZGTo2LG",                       directory=directories['ZGTo2LG'] )
#WGJets   = Sample.fromDirectory(name="WGJets",       treeName="Events", isData=False, color=color.diBoson,         texName="WGJets",                       directory=directories['WGJets']  )
#ZGJets   = Sample.fromDirectory(name="ZGJets",       treeName="Events", isData=False, color=color.diBoson,         texName="ZGJets",                       directory=directories['ZGJets']  )
ZG        = Sample.fromDirectory(name="ZG",            treeName="Events", isData=False, color=color.QCD,             texName="ZG",                           directory=directories['ZG']  )
EWK        = Sample.fromDirectory(name="EWK",            treeName="Events", isData=False, color=color.QCD,             texName="EWK",                           directory=directories['EWK']  )

