from StopsDilepton.tools.helpers import mZ, getVarValue, getObjDict
from math import *
import numbers

jetVars = ['eta','pt','phi','btagCSV','id','area','rawPt']

def getJets(c, jetVars=jetVars, jetColl="Jet"):
    return [getObjDict(c, jetColl+'_', jetVars, i) for i in range(int(getVarValue(c, 'n'+jetColl)))]

def jetId(j, ptCut=30, absEtaCut=2.4, ptVar='pt'):
  return j[ptVar]>ptCut and abs(j['eta'])<absEtaCut and j['id']

def getGoodJets(c, ptCut=30, absEtaCut=2.4, jetVars=jetVars, jetColl="Jet"):
    return filter(lambda j:jetId(j, ptCut=ptCut, absEtaCut=absEtaCut), getJets(c, jetVars, jetColl=jetColl))

def isBJet(j):
    return j['btagCSV']>0.800

def getGoodBJets(c):
    return filter(lambda j:isBJet(j), getGoodJets(c))

def getGenLeps(c):
    return [getObjDict(c, 'genLep_', ['eta','pt','phi','charge', 'pdgId', 'sourceId'], i) for i in range(int(getVarValue(c, 'ngenLep')))]

def getGenParts(c):
    return [getObjDict(c, 'GenPart_', ['eta','pt','phi','charge', 'pdgId', 'motherId', 'grandmotherId'], i) for i in range(int(getVarValue(c, 'nGenPart')))]

genVars = ['eta','pt','phi','mass','charge', 'status', 'pdgId', 'motherId', 'grandmotherId','nDaughters','daughterIndex1','daughterIndex2','nMothers','motherIndex1','motherIndex2','isPromptHard'] 
def getGenPartsAll(c):
    return [getObjDict(c, 'genPartAll_', genVars, i) for i in range(int(getVarValue(c, 'ngenPartAll')))]

def alwaysTrue(*args, **kwargs):
  return True
def alwaysFalse(*args, **kwargs):
  return False

def get_index_str( index ):
    if isinstance(index, int):
        index_str = "["+str(index)+"]"
    elif type(index)==type(""):
        if index.startswith('[') and index.endswith(']'):
            index_str = index
        else:
            index_str = '['+index+']'
    elif index is None:
        index_str=""
    else:
        raise ValueError( "Don't know what to do with index %r" % index )
    return index_str

def miniIsoSelectorString( iso, index = None):
    ''' Cut string for mini Iso'''
    index_str = get_index_str( index )
    if not isinstance(iso, numbers.Number):
        raise ValueError( "Don't know what to do with miniIso %r" % iso )
    return  "LepGood_miniRelIso"+index_str+"<%s"%iso

def relIso03SelectorString( iso, index = None):
    ''' Cut string for relIso03'''
    index_str = get_index_str( index )
    if not isinstance(iso, numbers.Number):
        raise ValueError( "Don't know what to do with relIso03 %r" % iso )
    return  "LepGood_relIso03"+index_str+"<%s"%iso

def isoSelectorString( relIso03, miniIso = None, index = None):
    ''' Cut string for all isos'''

    # isoSelector( x ) defaults to relIso03 selector
    if isinstance(relIso03, numbers.Number): iso_ = relIso03SelectorString( relIso03, index = index )
    # similar for miniIso
    elif isinstance(miniIso, numbers.Number):  iso_ = miniIsoSelectorString( miniIso, index = index )
    else:    raise ValueError( "Don't know what to do with iso args %r %r"%(relIso03, miniIso) )
    return iso_

def miniIsoSelector( miniRelIso ):
    assert isinstance(miniRelIso, numbers.Number), "Don't know what to do with miniRelIso %r"%miniRelIso
    def func(l):
        return  l["miniRelIso"] < miniRelIso
    return func

def relIso03Selector(iso):
    if not isinstance(iso, numbers.Number):
        raise ValueError( "Don't know what to do with relIso03 %r" % iso )
    def func(l):
        return l["relIso03"]<iso 
    return func

def isoSelector( relIso03, miniIso = None):

    # always true if no arguments
    if relIso03 is None and miniIso is None: iso_ = alwaysTrue
    # isoSelector( x ) defaults to relIso03 selector
    elif isinstance(relIso03, numbers.Number): iso_ = relIso03Selector( relIso03 )
    # similar for miniIso
    elif isinstance(miniIso, numbers.Number):  iso_ = miniIsoSelector( miniIso )
    else:    raise ValueError( "Don't know what to do with iso args %r %r"%(relIso03, miniIso) )
    return iso_

# MUONS
def muonSelector(relIso03 = 0.2, miniIso = None, absEtaCut = 2.4, dxy = 0.05, dz = 0.1, loose=False):

    iso_ = isoSelector( relIso03 = relIso03, miniIso = miniIso)

    def func(l, ptCut = 20):
        return \
            l["pt"]>=ptCut\
            and abs(l["pdgId"])==13\
            and abs(l["eta"])<absEtaCut\
            and l["mediumMuonId"]>=1 \
            and iso_(l) \
            and l["sip3d"]<4.0\
            and abs(l["dxy"])<dxy\
            and abs(l["dz"])<dz

    def funcLoose(l, ptCut = 20):
        return \
            l["pt"]>=ptCut\
            and abs(l["pdgId"])==13\
            and abs(l["eta"])<absEtaCut\
            and iso_(l) \
            and abs(l["dxy"])<dxy\
            and abs(l["dz"])<dz

    return func if not loose else funcLoose

default_muon_selector = muonSelector( relIso03 = 0.2, absEtaCut = 2.4 )

def muonSelectorString(relIso03 = 0.2, miniIso = None, ptCut = 20, absEtaCut = 2.4, dxy = 0.05, dz = 0.1, index = "Sum"):
    idx = None if (index is None) or (type(index)==type("") and index.lower()=="sum") else index
    index_str = get_index_str( index  = idx)
    string = [\
                   "LepGood_pt"+index_str+">=%s"%ptCut ,
                   "abs(LepGood_pdgId"+index_str+")==13" ,
                   "abs(LepGood_eta"+index_str+")<%s" % absEtaCut ,
                   "LepGood_mediumMuonId"+index_str+">=1" ,
                   "LepGood_sip3d"+index_str+"<4.0" ,
                   "abs(LepGood_dxy"+index_str+")<%s" % dxy ,
                   "abs(LepGood_dz"+index_str+")<%s" % dz ,
                   isoSelectorString( relIso03 = relIso03, miniIso = miniIso, index = idx),
             ]
    if type(index)==type("") and index.lower()=='sum':
        return 'Sum$('+'&&'.join(string)+')'
    else:
        return '&&'.join(string)

# ELECTRONS

#ele_MVAID =  {'VL': {(0,0.8):-0.16 , (0.8, 1.479):-0.65, (1.57, 999): -0.74},
#              'T':  {(0,0.8):0.87 , (0.8, 1.479):0.60, (1.57, 999):  0.17}
#}
#
#def eleMVAIDSelector( eleId ):
#    ele_mva_WP = ele_MVAID[eleId]
#    def func(l):
#        abs_ele_eta = abs(l["eta"])
#        for abs_ele_bin, mva_threshold in ele_mva_WP.iteritems():
#            if abs_ele_eta>=abs_ele_bin[0] and abs_ele_eta<abs_ele_bin[1] and l["mvaIdSpring15"] > mva_threshold: return True
#        return False
#    return func

def eleCutIDSelector( ele_cut_Id = 4):
    def func(l):
        return l["eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz"]>=ele_cut_Id 
    return func

def eleSelector(relIso03 = 0.2, miniIso = None, eleId = 4, absEtaCut = 2.4, dxy = 0.05, dz = 0.1, noMissingHits=True, loose=False):

    iso_ = isoSelector( relIso03 = relIso03, miniIso = miniIso)

    if isinstance(eleId, numbers.Number): id_ = eleCutIDSelector( eleId )
    # elif type(eleId)==type(""):           id_ = eleMVAIDSelector( eleId )
    else:                                 raise ValueError( "Don't know what to do with eleId %r" % eleId )

    def func(l, ptCut = 20):
        return \
            l["pt"]>=ptCut\
            and abs(l["eta"])<absEtaCut\
            and abs(l["pdgId"])==11\
            and id_(l)\
            and iso_(l)\
            and l["convVeto"]\
            and (l["lostHits"]==0 or not noMissingHits)\
            and l["sip3d"] < 4.0\
            and abs(l["dxy"]) < dxy\
            and abs(l["dz"]) < dz

    def funcLoose(l, ptCut = 20):
        return \
            l["pt"]>=ptCut\
            and abs(l["eta"])<absEtaCut\
            and abs(l["pdgId"])==11\
            and id_(l)\
            and iso_(l)\
            and abs(l["dxy"]) < dxy\
            and abs(l["dz"]) < dz

    return func if not loose else funcLoose

default_ele_selector = eleSelector( relIso03 = 0.2, eleId = 4, absEtaCut = 2.4 )

def eleIDSelectorString( eleId, index = None ):
    index_str = get_index_str( index )

    if isinstance(eleId, numbers.Number):
        return "LepGood_eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz"+index_str+">=%i" % eleId 
    #elif type(eleId)==type(""):
    #    return eleMVAString( eleId, index = index) 
    else:
        raise ValueError( "Don't know what to do with eleId %r" % eleId )

#def eleMVAString( eleId, index = None):
#    index_str = get_index_str( index )
#    ele_mva_WP = ele_MVAID[eleId]
#    abs_ele_eta = "abs(LepGood_eta"+index_str+")"
#    strings = []
#    for abs_ele_bin, mva_threshold in ele_mva_WP.iteritems():
#        strings.append("({abs_ele_eta}>={low_abs_ele_eta}&&{abs_ele_eta}<{high_abs_ele_eta}&&LepGood_mvaIdSpring15{idx_str}>{mva_threshold})".format(\
#            abs_ele_eta=abs_ele_eta, 
#            low_abs_ele_eta = abs_ele_bin[0],
#            high_abs_ele_eta=abs_ele_bin[1], 
#            idx_str = index_str,
#            mva_threshold = mva_threshold))
#
#    return "("+'||'.join(strings)+')'

def eleSelectorString(relIso03 = 0.2, miniIso = None, eleId = 4, ptCut = 20, absEtaCut = 2.4, dxy = 0.05, dz = 0.1, index = "Sum", noMissingHits=True):
    idx = None if (index is None) or (type(index)==type("") and index.lower()=="sum") else index
    index_str = get_index_str( index  = idx)
    string = [\
                   "LepGood_pt"+index_str+">=%s" % ptCut ,
                   "abs(LepGood_eta"+index_str+")<%s" % absEtaCut ,
                   "abs(LepGood_pdgId"+index_str+")==11" ,
                   "LepGood_convVeto"+index_str+"",
                   "LepGood_lostHits"+index_str+"==0" if noMissingHits else "(1)",
                   "LepGood_sip3d"+index_str+"<4.0" ,
                   "abs(LepGood_dxy"+index_str+")<%s" % dxy ,
                   "abs(LepGood_dz"+index_str+")<%s" % dz ,
                   isoSelectorString( relIso03 = relIso03, miniIso = miniIso, index = idx),
                   eleIDSelectorString( eleId, index = idx),
             ]

    if type(index)==type("") and index.lower()=='sum':
        return 'Sum$('+'&&'.join(string)+')'
    else:
        return '&&'.join(string)


leptonVars_data = ['eta','etaSc', 'pt','phi','dxy', 'dz','tightId', 'pdgId', 'mediumMuonId', 'miniRelIso', 'relIso03', 'sip3d', 'mvaIdSpring15', 'convVeto', 'lostHits', 'jetPtRelv2', 'jetPtRatiov2', 'eleCutId_Spring2016_25ns_v1_ConvVetoDxyDz']
leptonVars = leptonVars_data + ['mcMatchId','mcMatchAny']

def getLeptons(c, collVars=leptonVars):
    return [getObjDict(c, 'LepGood_', collVars, i) for i in range(int(getVarValue(c, 'nLepGood')))]
def getOtherLeptons(c, collVars=leptonVars):
    return [getObjDict(c, 'LepOther_', collVars, i) for i in range(int(getVarValue(c, 'nLepOther')))]
def getMuons(c, collVars=leptonVars):
    return [getObjDict(c, 'LepGood_', collVars, i) for i in range(int(getVarValue(c, 'nLepGood'))) if abs(getVarValue(c,"LepGood_pdgId",i))==13]
def getElectrons(c, collVars=leptonVars):
    return [getObjDict(c, 'LepGood_', collVars, i) for i in range(int(getVarValue(c, 'nLepGood'))) if abs(getVarValue(c,"LepGood_pdgId",i))==11]

def getGoodMuons(c, ptCut = 20, collVars=leptonVars, mu_selector = default_muon_selector):
    return [l for l in getMuons(c, collVars) if mu_selector(l, ptCut = ptCut)]
def getGoodElectrons(c, ptCut = 20, collVars=leptonVars, ele_selector = default_ele_selector):
    return [l for l in getElectrons(c, collVars) if ele_selector(l, ptCut = ptCut)]
def getGoodLeptons(c, ptCut=20, collVars=leptonVars, mu_selector = default_muon_selector, ele_selector = default_ele_selector):
    return [l for l in getLeptons(c, collVars) if (abs(l["pdgId"])==11 and ele_selector(l, ptCut = ptCut)) or (abs(l["pdgId"])==13 and mu_selector(l, ptCut = ptCut))]

def getGoodAndOtherLeptons(c, ptCut=20, collVars=leptonVars, mu_selector = default_muon_selector, ele_selector = default_ele_selector):
    good_lep = getLeptons(c, collVars)
    other_lep = getOtherLeptons(c, collVars)
    for l in other_lep: #dirty trick to find back the full lepton if it was in the 'other' collection
        l['index']+=1000
    res = [l for l in good_lep+other_lep if (abs(l["pdgId"])==11 and ele_selector(l, ptCut = ptCut)) or (abs(l["pdgId"])==13 and mu_selector(l, ptCut = ptCut))]
    res.sort( key = lambda l:-l['pt'] )
    return res

tauVars=['eta','pt','phi','pdgId','charge', 'dxy', 'dz', 'idDecayModeNewDMs', 'idCI3hit', 'idAntiMu','idAntiE','mcMatchId']

def getTaus(c, collVars=tauVars):
    return [getObjDict(c, 'TauGood_', collVars, i) for i in range(int(getVarValue(c, 'nTauGood')))]

def looseTauID(l, ptCut=20, absEtaCut=2.4):
    return \
        l["pt"]>=ptCut\
        and abs(l["eta"])<absEtaCut\
        and l["idDecayModeNewDMs"]>=1\
        and l["idCI3hit"]>=1\
        and l["idAntiMu"]>=1\
        and l["idAntiE"]>=1\

def getGoodTaus(c, collVars=tauVars):
    return [l for l in getTaus(c,collVars=collVars) if looseTauID(l)]

idCutBased={'loose':1 ,'medium':2, 'tight':3}
photonVars=['eta','pt','phi','mass','idCutBased','pdgId']
photonVarsMC = photonVars + ['mcPt']
def getPhotons(c, collVars=photonVars, idLevel='loose'):
    return [getObjDict(c, 'gamma_', collVars, i) for i in range(int(getVarValue(c, 'ngamma')))]
def getGoodPhotons(c, ptCut=50, idLevel="loose", isData=True, collVars=None):
    if collVars is None: collVars = photonVars if isData else photonVarsMC
    return [p for p in getPhotons(c, collVars) if p['idCutBased'] >= idCutBased[idLevel] and p['pt'] > ptCut and p['pdgId']==22]

def getFilterCut(isData=False, isFastSim = False, badMuonFilters = "Summer2016"):
    if isFastSim:
        filterCut            = "Flag_goodVertices"
    else:
        filterCut            = "Flag_goodVertices&&Flag_HBHENoiseIsoFilter&&Flag_HBHENoiseFilter&&Flag_globalTightHalo2016Filter&&Flag_eeBadScFilter&&Flag_EcalDeadCellTriggerPrimitiveFilter"
        if badMuonFilters == "Summer2016":
            filterCut += "&&Flag_badChargedHadronSummer2016&&Flag_badMuonSummer2016"
        elif badMuonFilters == "Summer2016_pt20":
            filterCut += "&&Flag_badChargedHadronSummer2016&&Flag_badMuonSummer2016_pt20"
        elif badMuonFilters is None or badMuonFilters == "None":
            pass
    if isData: filterCut += "&&weight>0"
    return filterCut
