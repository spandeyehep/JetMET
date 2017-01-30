#Standard imports
import ROOT
from math import pi, sqrt, cos, sin, sinh, log, cosh
from array import array
import itertools

# Logging
import logging
logger = logging.getLogger(__name__)

#scripts
ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/StopsDilepton/tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()
mZ=91.1876

def deltaPhi(phi1, phi2):
    dphi = phi2-phi1
    if  dphi > pi:
        dphi -= 2.0*pi
    if dphi <= -pi:
        dphi += 2.0*pi
    return abs(dphi)

def deltaR2(l1, l2):
    return deltaPhi(l1['phi'], l2['phi'])**2 + (l1['eta'] - l2['eta'])**2

def deltaR(l1, l2):
    return sqrt(deltaR2(l1,l2))

def bestDRMatchInCollection(l, coll, deltaR = 0.2, deltaRelPt = 0.5 ):
    lst = []
    for l2 in coll:
        dr2 = deltaR2(l, l2)
        if  ( dr2 < deltaR**2 ) and (abs( -1 + l2['pt']/l['pt']) < deltaRelPt or deltaRelPt < 0 ):
            lst.append((dr2, l2))
    lst.sort()
    if len(lst)>0:
        return lst[0][1]
    else:
        return None

def getFileList(dir, histname='histo', maxN=-1):
    import os
    filelist = os.listdir(os.path.expanduser(dir))
    filelist = [dir+'/'+f for f in filelist if histname in f and f.endswith(".root")]
    if maxN>=0:
        filelist = filelist[:maxN]
    return filelist

def checkRootFile(f, checkForObjects=[]):
    rf = ROOT.TFile.Open(f)
    if not rf: return False
    try:
        good = (not rf.IsZombie()) and (not rf.TestBit(ROOT.TFile.kRecovered))
    except:
        if rf: rf.Close()
        return False
    for o in checkForObjects:
        if not rf.GetListOfKeys().Contains(o):
            print "[checkRootFile] Failed to find object %s in file %s"%(o, f)
            rf.Close()
            return False
#    print "Keys recoveredd %i zombie %i tb %i"%(rf.Recover(), rf.IsZombie(), rf.TestBit(ROOT.TFile.kRecovered))
    rf.Close()
    return good

def getObjFromFile(fname, hname):
    f = ROOT.TFile(fname)
    assert not f.IsZombie()
    f.cd()
    htmp = f.Get(hname)
    if not htmp:  return htmp
    ROOT.gDirectory.cd('PyROOT:/')
    res = htmp.Clone()
    f.Close()
    return res

def writeObjToFile(fname, obj):
    gDir = ROOT.gDirectory.GetName()
    f = ROOT.TFile(fname, 'recreate')
    objw = obj.Clone()
    objw.Write()
    f.Close()
    ROOT.gDirectory.cd(gDir+':/')
    return

def getObjDict(c, prefix, variables, i):
    res={var: getVarValue(c, prefix+var, i) for var in variables}
    res['index']=i
    return res

def getCollection(c, prefix, variables, counter_variable):
    return [getObjDict(c, prefix+'_', variables, i) for i in range(int(getVarValue(c, counter_variable)))]

def read_from_subprocess(arglist):
    ''' Read line by line from subprocess
    '''
    import subprocess

    proc = subprocess.Popen(arglist,stdout=subprocess.PIPE)
    res = []
    while True:
        l = proc.stdout.readline()
        if l !=  '':
            res.append( l.rstrip() )
        else:
            break
    return res

def renew_proxy( filename = None, rfc = False, request_time = 192, min_time = 0):
    import os, subprocess

    proxy = None
    timeleft = 0

    # Make voms-proxy-info look for a specific proxy
    if filename is not None:
        os.environ["X509_USER_PROXY"] = filename

    # Check proxy path
    try:
        proxy     = read_from_subprocess( 'voms-proxy-info --path'.split() )[0]
    except IndexError:
        pass

    try:
        tl = read_from_subprocess( 'voms-proxy-info --timeleft'.split() )
        timeleft = int(float( tl[0] ))
    except IndexError:
        pass
    except ValueError:
        print tl
        pass

    # Return existing proxy from $X509_USER_PROXY, the default location or filename
    if proxy is not None and os.path.exists( proxy ):
        if filename is None or os.path.abspath( filename ) == proxy:
            logger.info( "Found proxy %s with lifetime %i hours", proxy, timeleft/3600)
            if timeleft > 0 and timeleft >= min_time*3600 :
                os.environ["X509_USER_PROXY"] = proxy
                return proxy
            else:
                logger.info( "Lifetime %i not sufficient (require %i, will request %i hours).", timeleft/3600, min_time, request_time )

    
    arg_list = ['voms-proxy-init', '-voms', 'cms']

    if filename is not None:
        arg_list += [ '-out', filename ]

    arg_list += ['--valid', "%i:0"%request_time ]

    if rfc:
        arg_list += ['-rfc']

    # make proxy
    p = subprocess.call( arg_list )

    # read path
    new_proxy = None
    try:
        new_proxy     = read_from_subprocess( 'voms-proxy-info --path'.split() )[0]
    except IndexError:
        pass

    if new_proxy is not None and os.path.exists( new_proxy ):
        os.environ["X509_USER_PROXY"] = new_proxy
        logger.info( "Successfully created new proxy %s", new_proxy )
        return new_proxy
    else:
        raise RuntimeError( "Failed to make proxy %s" % new_proxy )

