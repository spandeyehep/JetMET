import os

if os.environ['USER'] in ['schoef', 'rschoefbeck', 'schoefbeck']:
    # master ntuple
    master_ntuple_directory = "/data/rschoefbeck/cmgTuples/"
    # Where you store cmg output
    cmg_directory      = "/scratch/rschoefbeck/cmgTuples/80X_1l_9"
    # Where postprocessed data goes 
    data_output_directory      = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
    #data_output_directory      = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
    # Where you store the data
    data_directory      = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
    # Where the plots go
    plot_directory      = "/afs/hephy.at/user/r/rschoefbeck/www/"
    #plot_directory      = "/afs/cern.ch/work/s/schoef/www/"
    # Analysis result files
    analysis_results    = '/afs/hephy.at/data/rschoefbeck01/JetMET/results/' #Path to analysis results
    # 715 release for limit calculation 
    combineReleaseLocation = '/afs/hephy.at/work/r/rschoefbeck/CMS/tmp/CMSSW_7_1_5/src'
    runOnGentT2 = False
