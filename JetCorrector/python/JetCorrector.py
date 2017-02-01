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

def wget( source, target ):
    ''' Download source to target and make directories
    '''
    if not os.path.exists( os.path.dirname( target ) ):
        os.makedirs( os.path.dirname( target ) )

    urllib.urlretrieve( source, target )

correction_levels_data = [ 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' ]
correction_levels_mc = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]

class JetCorrector:

    def __init__( self, 
            iovs, 
            correctionLevels, 
            baseurl     = "https://github.com/cms-jet/JECDatabase/raw/master/tarballs/",
            jetFlavour  = 'AK4PFchs',
            directory   = "$CMSSW_BASE/src/JetMET/JetCorrector/data/",
            ):

        extension   = "tar.gz"
        self.jetCorrectors = []
        for runnumber, filename in iovs:
            # Download file
            source = os.path.join( baseurl, "%s.%s"%(filename, extension) )
            target = os.path.join( os.path.expandvars( directory ), "%s.%s"%(filename, extension) )

            params = ROOT.vector(ROOT.JetCorrectorParameters)()
            
            for level in correctionLevels:
                txtfile = os.path.join( os.path.expandvars(directory), "%s_%s_%s.txt"%( filename, level, jetFlavour)  )

                # Do we have the txt file?
                if not os.path.exists( txtfile ):
                    logger.info( "txt file %s not found.", txtfile )

                    # Do we actually have the tar.gz?
                    if not os.path.exists( target ):
                        logger.info( "%s not found. Downloading from %s.", target, source )
                        wget( source, target )

                    # Extract txt files that match the bill
                    logger.info( "Extracting %s", target )
                    with tarfile.open(target, 'r:gz') as tar:
                        for member in tar.getmembers():
                            if not member.name.endswith( jetFlavour+'.txt' ): continue
                            logger.debug( "Found file %s in %s", member.name, target )
                            with file( os.path.join( os.path.expandvars(directory), member.name), 'w') as f_out:
                                f_out.writelines( tar.extractfile( member ).readlines() )

                params.push_back( ROOT.JetCorrectorParameters( txtfile, "" ) )

            # Make corrector
            self.jetCorrectors.append( ( runnumber,  ROOT.FactorizedJetCorrector( params )) )

        # Sort wrt IOVs
        self.jetCorrectors.sort( key = lambda p: p[0] )

    def correction(self, rawPt, eta, area, rho, run = 1):

        for runnumber, corrector in reversed( self.jetCorrectors ):
            if run >= runnumber:
                corrector.setJetPt( rawPt )
                corrector.setJetEta( eta )
                corrector.setJetA( area )
                corrector.setRho( rho )
                return corrector.getCorrection()
