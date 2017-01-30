''' JEC on the fly
'''
# Standard imports
import os

# Logging

import logging
logger = logging.getLogger(__name__)

# config
Summer16_23Sep2016_DATA = \
((1, 'Summer16_23Sep2016BCDV3_DATA'),
 (276831, 'Summer16_23Sep2016EFV3_DATA'),
 (278802, 'Summer16_23Sep2016GV3_DATA'),
 (280919, 'Summer16_23Sep2016HV3_DATA'))

# config
Summer16_23Sep2016_MC = ((1, 'Summer16_23Sep2016BCDV3_MC'))

def wget( source, target ):
    ''' Download source to target and make directories
    '''
    if not os.path.exists( os.path.dirname( target ) ):
        os.makedirs( os.path.dirname( target ) )

    import urllib
    urllib.urlretrieve (source, target)

class JetCorrector:

    def __init__( self, 
            iovs, 
            baseurl     = "https://github.com/cms-jet/JECDatabase/raw/master/tarballs/",
            extension = "tar.gz",
            directory   = "$CMSSW_BASE/src/JetMET/JetCorrector/data/",
            overwrite = False,
            ):

        for runnumber, filename in iovs:

            # Download file
            source = os.path.join( baseurl, "%s.%s"%(filename, extension) )
            target = os.path.join( os.path.expandvars( directory ), "%s.%s"%(filename, extension) )
            if not os.path.exists( target ) or overwrite:
                logger.info( "Downloading %s to %s.", source, target )
                wget( source, target )
            else:
                logger.info( "Found %s. Do nothing.", source )


   
if __name__ == "__main__":

    # Logging
    import JetMET.tools.logger as logger
    logger  = logger.get_logger('INFO', logFile = None)

    jetCorrector = JetCorrector( Summer16_23Sep2016_DATA ) 
