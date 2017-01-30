''' JEC on the fly
'''
# Standard imports
import os
import urllib
import tarfile
import ROOT

# Logging

import logging
logger = logging.getLogger(__name__)

# config
correction_levels_data = [ 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' ]
Summer16_23Sep2016_DATA = \
((1, 'Summer16_23Sep2016BCDV3_DATA'),
 (276831, 'Summer16_23Sep2016EFV3_DATA'),
 (278802, 'Summer16_23Sep2016GV3_DATA'),
 (280919, 'Summer16_23Sep2016HV3_DATA'))

# config
correction_levels_mc = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]
Summer16_23Sep2016_MC = ((1, 'Summer16_23Sep2016V3_MC'),)

#Summer16_23Sep2016GV3_DATA_L2L3Residual_AK8PFchs.txt

def wget( source, target ):
    ''' Download source to target and make directories
    '''
    if not os.path.exists( os.path.dirname( target ) ):
        os.makedirs( os.path.dirname( target ) )

    urllib.urlretrieve( source, target )

class JetCorrector:

    def __init__( self, 
            iovs, 
            correctionLevels, 
            baseurl     = "https://github.com/cms-jet/JECDatabase/raw/master/tarballs/",
            jetFlavour  = 'AK4PFchs',
            extension   = "tar.gz",
            directory   = "$CMSSW_BASE/src/JetMET/JetCorrector/data/",
            overwrite   = False,
            ):

        for runnumber, filename in iovs:
            # Download file
            source = os.path.join( baseurl, "%s.%s"%(filename, extension) )
            target = os.path.join( os.path.expandvars( directory ), "%s.%s"%(filename, extension) )
            if not os.path.exists( target ) or overwrite:
                logger.info( "Downloading %s to %s.", source, target )
                wget( source, target )
                # Extract txt files
                with tarfile.open(target, 'r:gz') as tar:
                    for member in tar.getmembers():
                        if not member.name.endswith( jetFlavour+'.txt' ): continue
                        logger.debug( "Found file %s in %s", member.name, target )
                        with file( os.path.join( os.path.expandvars(directory), member.name), 'w') as f_out:
                            f_out.writelines( tar.extractfile( member ).readlines() )
            else:
                logger.info( "Found %s.", target )

        # Make corrector
        self.vPar = ROOT.vector(ROOT.JetCorrectorParameters)()
        for level in correctionLevels:
            txtfile = os.path.join( os.path.expandvars(directory), "%s_%s_%s.txt"%( filename, level, jetFlavour)  )
            self.vPar.push_back( ROOT.JetCorrectorParameters( txtfile, "" ) )

        self.JetCorrector = ROOT.FactorizedJetCorrector(self.vPar) 

if __name__ == "__main__":

    # Logging
    import JetMET.tools.logger as logger
    logger  = logger.get_logger('INFO', logFile = None)

    jetCorrector_data = JetCorrector( Summer16_23Sep2016_DATA, correctionLevels = correction_levels_data ) 
    jetCorrector_mc   = JetCorrector( Summer16_23Sep2016_MC, correctionLevels = correction_levels_mc ) 
