''' FWLiteReader example: Loop over a sample and write some data to a histogram.
'''
# Standard imports
import os
import logging
import ROOT
import array

#RootTools
from RootTools.core.standard import *

#Helper
import JetMET.tools.helpers as helpers
from math import pi

# argParser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel', 
      action='store',
      nargs='?',
      choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],
      default='INFO',
      help="Log level for logging"
)

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

max_events  = -1
max_files   = -1

#Kenichi Feb. 2016

#sample_prefix = "kenichi_private_qcd_"
#/afs/cern.ch/user/d/deguio/public/ForKen/Plan0_10024.0_TTbar_13.txt
#/afs/cern.ch/user/d/deguio/public/ForKen/Plan0_10043.0_QCDForPF_14TeV.txt
#/afs/cern.ch/user/d/deguio/public/ForKen/Plan1_10024.0_TTbar_13.txt
#/afs/cern.ch/user/d/deguio/public/ForKen/Plan1_10043.0_QCDForPF_14TeV.txt
#files_ttbar_plan0   = [ 'root://eoscms.cern.ch/%s'%s.rstrip() for s in open('/afs/cern.ch/user/d/deguio/public/ForKen/Plan0_10024.0_TTbar_13.txt').readlines() if os.path.split(s)[-1].startswith('step3_') ]
#files_ttbar_plan1   = [ 'root://eoscms.cern.ch/%s'%s.rstrip() for s in open('/afs/cern.ch/user/d/deguio/public/ForKen/Plan1_10024.0_TTbar_13.txt').readlines() if os.path.split(s)[-1].startswith('step3_') ]
#files_qcd_plan0     = [ 'root://eoscms.cern.ch/%s'%s.rstrip() for s in open('/afs/cern.ch/user/d/deguio/public/ForKen/Plan0_10043.0_QCDForPF_14TeV.txt').readlines() if os.path.split(s)[-1].startswith('step3_') ]
#files_qcd_plan1     = [ 'root://eoscms.cern.ch/%s'%s.rstrip() for s in open('/afs/cern.ch/user/d/deguio/public/ForKen/Plan1_10043.0_QCDForPF_14TeV.txt').readlines() if os.path.split(s)[-1].startswith('step3_') ]

veto = ['step3_800.root']

# Federico March 3rd
sample_prefix = "federico_private_qcd_caloJets_"
dir_qcd_plan0 = "/eos/cms/store/group/dpg_hcal/comm_hcal/deguio/Plan0/10043.0_QCDForPF_14TeV+QCDForPF_14TeV_TuneCUETP8M1_2017_GenSimFull+DigiFull_2017+RecoFull_2017+ALCAFull_2017+HARVESTFull_2017/submit_20170303_005403/"
dir_qcd_plan1 = "/eos/cms/store/group/dpg_hcal/comm_hcal/deguio/Plan1/10043.0_QCDForPF_14TeV+QCDForPF_14TeV_TuneCUETP8M1_2017_GenSimFull+DigiFull_2017+RecoFull_2017+ALCAFull_2017+HARVESTFull_2017/submit_20170303_005702/"

files_qcd_plan0 = []
for f in os.listdir(dir_qcd_plan0):
    if os.path.isfile(os.path.join(dir_qcd_plan0, f)) and f.startswith('step3') and not any(v in f for v in veto):
        files_qcd_plan0.append( os.path.join(dir_qcd_plan0, f) )

files_qcd_plan1 = []
for f in os.listdir(dir_qcd_plan1):
    if os.path.isfile(os.path.join(dir_qcd_plan1, f)) and f.startswith('step3') and not any(v in f for v in veto):
        files_qcd_plan1.append( os.path.join(dir_qcd_plan1, f) )

plan0 = FWLiteSample.fromFiles("plan0", files = files_qcd_plan0, maxN = max_files)
plan1 = FWLiteSample.fromFiles("plan1", files = files_qcd_plan1, maxN = max_files)

pt_threshold = 50
preprefix = "%s_pt%i" % ( sample_prefix, pt_threshold )

#plan0    = FWLiteSample.fromFiles("plan0", files = files_qcd_plan0, maxN = max_files)
#plan1    = FWLiteSample.fromFiles("plan1", files = files_qcd_plan1, maxN = max_files)
#preprefix = "kenichi_private_ttbar_pt%i" % pt_threshold
#plan0    = FWLiteSample.fromFiles("plan0", files = files_ttbar_plan0, maxN = max_files)
#plan1    = FWLiteSample.fromFiles("plan1", files = files_ttbar_plan1, maxN = max_files)


# define TProfiles
pt_thresholds = [ 10**(x/10.) for x in range(11,36) ] 


phi_low  = 310/180.*pi - 0.2 - 2*pi
phi_high = 330/180.*pi + 0.2 - 2*pi

resp = {}
resp_eta = {}
resp_eta_HEP17 = {}
resp_eta_nonHEP17 = {}
resp_eta_phi = {}

# 1D Profile inclusive
resp_eta = ROOT.TProfile("resp_eta", "resp_eta", 26, -5.2, 5.2 )
resp_eta.style = styles.lineStyle( ROOT.kBlack )
resp_eta.legendText = "all" 

resp_eta_HEP17 = ROOT.TProfile("resp_eta_HEP17", "resp_eta_HEP17", 26, -5.2, 5.2 )
resp_eta_HEP17.style = styles.lineStyle( ROOT.kRed )
resp_eta_HEP17.legendText = "%3.2f #leq #phi #leq %3.2f"%( phi_low, phi_high)

resp_eta_nonHEP17 = ROOT.TProfile("resp_eta_nonHEP17", "resp_eta_nonHEP17", 26, -5.2, 5.2 )
resp_eta_nonHEP17.style = styles.lineStyle( ROOT.kBlue )
resp_eta_nonHEP17.legendText = "!(%3.2f #leq #phi #leq %3.2f)"%( phi_low, phi_high)

resp_eta_phi = ROOT.TProfile2D("resp_eta_phi", "resp_eta", 26, -5.2, 5.2, 20,-pi,pi )

resp={}
eta_thresholds = [(0,1.1), (1.3,3), (3,5)]
eta_color = {(0,1.1):ROOT.kBlack, (1.3,3):ROOT.kRed, (3,5):ROOT.kBlue}
for i_eta_th, eta_th in enumerate( eta_thresholds ):
    resp[eta_th] = ROOT.TProfile("response", "response", 36, -pi, pi )
    resp[eta_th].style = styles.lineStyle(eta_color[eta_th] )
    resp[eta_th].legendText = "%2.1f #leq #eta < %2.1f"%eta_th

thresholds = [10**(x/10.) for x in range(11,36)] 
resp_pt = ROOT.TProfile("response_pt", "response_pt", len(thresholds)-1, array.array('d', thresholds) )
resp_pt.style = styles.lineStyle( ROOT.kRed )
resp_pt.legendText = "HEP17"
resp_pt_ref = ROOT.TProfile("response_pt_ref", "response_pt_ref", len(thresholds)-1, array.array('d', thresholds) )
resp_pt_ref.style = styles.lineStyle( ROOT.kBlack )
resp_pt_ref.legendText = "other"

products = {
#    'jets':      {'type': 'vector<reco::PFJet>', 'label':"ak4PFJetsCHS"},
    'jets':     {'type':'vector<reco::CaloJet>', 'label': "ak4CaloJets" },

    }

r1 = plan1.fwliteReader( products = products )
r2 = plan0.fwliteReader( products = products )

r1.start()
runs_1 = set()
position_r1 = {}
count=0
while r1.run( readProducts = False ):
    file_n1 =  int(os.path.split(r1.sample.events._filenames[r1.sample.events.fileIndex()])[-1].split('.')[0].split('_')[1])
    position_r1[(r1.evt, file_n1)] = r1.position-1
    count+=1
    if max_events is not None and max_events>0 and count>=max_events:break

r2.start()
runs_2 = set()
position_r2 = {}
count=0
while r2.run( readProducts = False ):
    file_n2 =  int(os.path.split(r2.sample.events._filenames[r2.sample.events.fileIndex()])[-1].split('.')[0].split('_')[1])
    position_r2[(r2.evt, file_n2)] = r2.position-1
    count+=1
    if max_events is not None and max_events>0 and count>=max_events:break

logger.info( "Have %i events in first samle and %i in second", len(position_r1), len(position_r2) )

# Fast intersect
intersec = set(position_r1.keys()).intersection(set(position_r2.keys()))
positions = [(position_r1[i], position_r2[i]) for i in intersec]

# Without sorting, there is a jump between files with almost every event -> extremly slow
positions.sort()
logger.info("Have %i events in common.", len(intersec))

#Looping over common events
for i, p in enumerate(positions):
    p1,p2 = p
    r1.goToPosition(p1)
    r2.goToPosition(p2)
    if i%10000==0: logger.info("At %i/%i of common events.", i, len(positions))
    jets1_ = [ j for j in r1.products['jets'] ] #if helpers.jetID( j )]
    jets2_ = [ j for j in r2.products['jets'] ] #if helpers.jetID( j )]
    jets1 = [{'pt':j.pt(), 'eta':j.eta(), 'phi':j.phi(), 'j':j } for j in jets1_]
    jets2 = [{'pt':j.pt(), 'eta':j.eta(), 'phi':j.phi(), 'j':j } for j in jets2_]
    for c in zip(jets1, jets2):
        if helpers.deltaR2(*c)<0.2**2:

            r = c[0]['pt']/c[1]['pt']

            if c[0]['eta']>1.1 and c[0]['eta']<3.2 and c[0]['phi']>-1.07 and c[0]['phi']<-0.32:
                resp_pt.Fill( c[0]['pt'], r )
            else:
                resp_pt_ref.Fill( c[0]['pt'], r )

            if not c[0]['pt']>pt_threshold: continue
            #if not ( helpers.jetID(c[0]['j']) and helpers.jetID(c[1]['j']) ): continue

            #if not c[1]['nHEF']>0.1: continue
            #r = c[0]['nHEF']/c[1]['nHEF']

            # 2D eta, phi
            resp_eta_phi.Fill( c[0]['eta'], c[0]['phi'], r )

            # inclusive eta
            resp_eta.Fill( c[0]['eta'], r )

            isHEP17 = ( c[0]['phi'] > phi_low and c[0]['phi']<phi_high )
            if isHEP17:
                resp_eta_HEP17.Fill( c[0]['eta'], r )
            else:
                resp_eta_nonHEP17.Fill( c[0]['eta'], r )

            for eta_th in eta_thresholds:
                if c[0]['eta']>=eta_th[0] and c[0]['eta']<eta_th[1]:
                    resp[eta_th].Fill( c[0]['phi'], r )
                    break

# Make plot
profiles = [resp[t] for t in eta_thresholds]
#profiles = [ jetResponse_NJC]
prefix=preprefix + "_" + plan1.name + ("_max_events_%s_"%max_events if max_events is not None and max_events>0 else "" )
histos = [ [h.ProjectionX()] for h in profiles ]
for i, h in enumerate(histos):
    h[0].__dict__.update(profiles[i].__dict__)

jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponseRatio_relval", histos = histos, texX = "plan1 Jet #phi" , texY = "response ratio plan1/plan0" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = False, yRange=(0.7,1.2))

profiles = [resp_eta, resp_eta_HEP17, resp_eta_nonHEP17 ]
#profiles = [ jetResponse_NJC]
prefix=preprefix + "_" + plan1.name + ("_max_events_%s_"%max_events if max_events is not None and max_events>0 else "" )
histos = [ [h.ProjectionX()] for h in profiles ]
for i, h in enumerate(histos):
    h[0].__dict__.update(profiles[i].__dict__)

jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponseRatio_relVal_HEP_comparison", histos = histos, texX = "plan1 Jet #eta" , texY = "response ratio plan1/plan0" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = False, yRange=(0.7,1.2))

h2d = resp_eta_phi.ProjectionXY()
plotting.draw2D(
    plot = Plot2D.fromHisto(name = prefix+"jetResponseRatio_2D", histos = [[h2d]], texX = "plan-1 jet #eta", texY = "plan-1 jet #phi"),
    plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/",
    logX = False, logY = False, logZ = False, zRange = (0.95,1.05),
)

profiles = [resp_pt, resp_pt_ref]
prefix=preprefix + "_" + plan1.name + ("_max_events_%s_"%max_events if max_events is not None and max_events>0 else "" )
histos = [ [h.ProjectionX()] for h in profiles ]
for i, h in enumerate(histos):
    h[0].__dict__.update(profiles[i].__dict__)

jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponseRatio_pt", histos = histos, texX = "plan1 Jet p_{T}" , texY = "response ratio plan1/plan0" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = True, yRange=(0.7,1.2))
