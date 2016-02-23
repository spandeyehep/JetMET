#Standard imports
import sys
import logging
import ROOT
import array
from math import log

#RootTools
from RootTools.core.Sample import Sample
from RootTools.core.Variable import Variable, ScalarType, VectorType
from RootTools.core.TreeReader import TreeReader
from RootTools.core.logger import get_logger
from RootTools.plot.Plot import Plot
import RootTools.plot.plotting as plotting

#Helper
from helpers import deltaR2

# argParser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',
      action='store',
      nargs='?',
      choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],
      default='DEBUG',
      help="Log level for logging"
)

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

#make samples. Samples are statisticall compatible.
maxN = -1

QCD_Pt_15to3000_M2_5_500 = Sample.fromCMGOutput("QCD_Pt_15to3000_M2_5_500", baseDirectory = "/data/rschoefbeck/cmgTuples/QCD_HCal/QCD_Pt_15to3000_M2_5_500", treeFilename='tree.root', maxN=maxN)

