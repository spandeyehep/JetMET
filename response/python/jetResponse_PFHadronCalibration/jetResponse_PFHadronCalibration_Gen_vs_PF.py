''' FWLiteReader example: Loop over a sample and write some data to a histogram.
'''
# Standard imports
import os
import sys
import logging
import ROOT
import array

#RootTools
from RootTools.core.standard import *

#Helper
import JetMET.tools.helpers as helpers
import JetMET.tools.user as user

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

max_events = -1
max_files = 1

dir_ref = [
    "/eos/cms/store/group/phys_jetmet/spandey/QCD/step3/CRAB_UserFiles/crab_QCD_step3_RECO_902_May_29/170529_095645/0000/",
    "/eos/cms/store/group/phys_jetmet/spandey/QCD/step3/CRAB_UserFiles/crab_QCD_step3_RECO_902_May_29/170529_095645/0001/",
    
    ]

dir_new = dir_ref
# dir_new = [
#     "/eos/cms/store/group/phys_jetmet/spandey/SingleNeutrino/step3/CRAB_UserFiles/crab_SingleNeutrino_step3_RECO_902_May_25/170525_153554/0000/",
#     "/eos/cms/store/group/phys_jetmet/spandey/SingleNeutrino/step3/CRAB_UserFiles/crab_SingleNeutrino_step3_RECO_902_May_25/170525_153554/0001/",
# ]

new = FWLiteSample.fromDirectory("new", directory = dir_new, maxN = max_files)
ref = FWLiteSample.fromDirectory("ref", directory = dir_ref, maxN = max_files)

preprefix = "QCD_genJets_1_file_TEST" 

# define TProfiles
pt_thresholds = [ 10**(x/10.) for x in range(11,36) ] 

eta_thresholds = [0, 1.3, 2.5, 3.2]
color={0:ROOT.kBlack, 1.3:ROOT.kBlue, 2.5:ROOT.kRed, 3.2:ROOT.kMagenta, "18Apr_vs_03Feb":ROOT.kBlue}

resp_eta = ROOT.TProfile("resp_eta", "resp_eta", 26, 0, 5.2 )
resp_eta.style = styles.lineStyle(ROOT.kBlack )
#resp_eta.legendText = comp 

resp={}
for i_eta_th, eta_th in enumerate( eta_thresholds ):
    resp[eta_th] = ROOT.TProfile("response", "response", len(pt_thresholds)-1, array.array('d', pt_thresholds) )
    resp[eta_th].style = styles.lineStyle(color[eta_th] )
    resp[eta_th].legendText = "%2.1f<=#eta"%eta_th
    if eta_th!=eta_thresholds[-1]: resp[eta_th].legendText += "<%2.1f"%eta_thresholds[i_eta_th+1]

products_ref = {
    'PFjets':      {'type': 'vector<reco:PFJet>', 'label':"ak4PFJetsCHS"},
    }

products_new = {
    'Genjets':      {'type': 'vector<reco:GenJet>', 'label':"ak4GenJets"},
    }

# products = {
#     'Genjets':      {'type': 'vector<reco:GenJet>', 'label':"ak4GenJets"},
#     'PFjets':      {'type': 'vector<reco:PFJet>', 'label':"ak4PFJetsCHS"},
#     }


r1 = new.fwliteReader( products = products_new )
r2 = ref.fwliteReader( products = products_ref )


r1.start()
evnt=0
while r1.run( readProducts = False ):
    evnt+=1
    if max_events is not None and max_events>0 and evnt>=max_events:break

logger.info( "Have %i events in first sample", evnt)

MatchedJets = 0
dR_list = list()
pt_ratio = list()
ii = 0
matched_jets_Gen_PF = list()
for count in range(evnt):
        
        r1.goToPosition(count)
        r2.goToPosition(count)
        jets1_ = [ j for j in r1.products['Genjets']]
        jets2_ = [ j for j in r2.products['PFjets'] if helpers.jetID( j )]
        
        jets1 = [{'pt':j.pt(), 'eta':j.eta(), 'phi':j.phi(), 'j':j} for j in jets1_]
        jets2 = [{'pt':j.pt(), 'eta':j.eta(), 'phi':j.phi(), 'j':j} for j in jets2_]
        # print "# of GenJet1: ",jets1.__len__()
        # print "# of PFJet1: ",jets2.__len__()
        # print "**********"
        for PFJ in jets2:
            #print "NEW PF"
            for GenJ in jets1:
                #c = [(PFJ,GenJ)]
                dR_list.append((helpers.deltaR2(PFJ,GenJ))**(0.5))
                pt_ratio.append(PFJ['pt']/GenJ['pt'])
                # print "deltaR(PF,Gen):", helpers.deltaR2(PFJ,GenJ)
                # print "Ratio of pT PF/Gen:", PFJ['pt']/GenJ['pt']
                #if (helpers.deltaR2(PFJ,GenJ) < 0.2**2): 
            if (min(dR_list) < 0.2):
                MatchedJets += 1
                #print "**** : index of Min dR: ",dR_list.index(min(dR_list))
                matched_jets_Gen_PF.append((jets1[dR_list.index(min(dR_list))], PFJ))
                # print "PFJet: ", PFJ
                # print "GenJet: ", jets1[dR_list.index(min(dR_list))]
                # print "dR_cross_check:", ((helpers.deltaR2(PFJ,(jets1[dR_list.index(min(dR_list))])))**(0.5))/min(dR_list)
                
            #print "dR:", dR_list
            #print "**** : Min dR: ",min(dR_list), " , pt Ratio:", pt_ratio[dR_list.index(min(dR_list))]
            
            del pt_ratio[:]
            del dR_list[:]

        # print "Event #", count
        # print "# of GenJet1: ",jets1.__len__()
        # print "# of PFJet1: ",jets2.__len__()
        #print "Jet2: ",jets2
        
        # print "Total # of Matched jets:", MatchedJets
        # print "Size of matched_jets_Gen_PF:", matched_jets_Gen_PF.__len__()
        if count%10000==0: logger.info("At %i/%i of total events.", count, evnt)
        for c in matched_jets_Gen_PF:
            resp_eta.Fill( abs(c[0]['eta']), c[0]['pt']/c[1]['pt'] )
            for eta_th in reversed(eta_thresholds):
                if abs(c[0]['eta'])>eta_th:
                    resp[eta_th].Fill( c[0]['pt'], c[0]['pt']/c[1]['pt'] )
                    #break

        # print "**********"
        count+=1
        MatchedJets = 0
        del matched_jets_Gen_PF[:]
        # if count == 720: break
        # if max_events is not None and max_events>0 and count>=max_events:break






# Make plot
profiles = [resp[t] for t in eta_thresholds]
#profiles = [ jetResponse_NJC]
prefix=preprefix + "_" + new.name + ("_max_events_%s_"%max_events if max_events is not None and max_events>0 else "" )
histos = [ [h.ProjectionX()] for h in profiles ]
for i, h in enumerate(histos):
    h[0].__dict__.update(profiles[i].__dict__)

jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponseRatio_relval", histos = histos, texX = "new Jet p_{T}" , texY = "response ratio Gen/PF" )
#plotting.draw(jetResponsePlot, plot_directory = user.plot_directory, ratio = None, logY = False, logX = True, yRange=(0.7,1.2))
plotting.draw(jetResponsePlot, plot_directory = user.plot_directory, ratio = None, logY = False, logX = True, yRange=(0.5,1.5))

# Make eta plot
profiles = [resp_eta  ]
prefix=preprefix + "_eta_" + ("_max_events_%s_"%max_events if max_events is not None and max_events>0 else "" )
histos = [ [h.ProjectionX()] for h in profiles ]
for i, h in enumerate(histos):
    h[0].__dict__.update(profiles[i].__dict__)

jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponseRatio_eta_relval", histos = histos, texX = "new Jet #eta" , texY = "response ratio Gen/PF" )
plotting.draw(jetResponsePlot, plot_directory = user.plot_directory, ratio = None, logY = False, logX = False, yRange=(0.5,1.5))
