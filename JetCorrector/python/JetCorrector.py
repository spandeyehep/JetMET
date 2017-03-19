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

correction_levels_data  = [ 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' ]
correction_levels_mc    = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]

class JetCorrector:

    data_directory   = "$CMSSW_BASE/src/JetMET/JetCorrector/data/"
    extension        = "tar.gz"


    def __init__( self, iovs ):

        self.iovs = iovs
        self.makeNewCorrectors()
    
    def makeNewCorrectors( self, require = None ):
        self.jetCorrectors = []
        for runnumber, txtfiles in self.iovs:
            params = ROOT.vector(ROOT.JetCorrectorParameters)()
            for txtfile in txtfiles:

                # Skip the file if not any of elements of require is found
                if require is not None:
                    if not any( r in txtfile for r in require):
                        logger.debug( "Could not find any of %s in txt file %s. Skip.", ",".join(require), txtfile)
                        continue 

                logger.debug( "Including %s in corrector.", txtfile )

                params.push_back( ROOT.JetCorrectorParameters( txtfile, "" ) )

            self.jetCorrectors.append( ( runnumber,  ROOT.FactorizedJetCorrector( params )) )

        # Sort wrt IOVs
        self.jetCorrectors.sort( key = lambda p: p[0] )

    @classmethod        
    def fromTarBalls( cls, 
            iovs, 
            correctionLevels, 
            baseurl     = "https://github.com/cms-jet/JECDatabase/raw/master/tarballs/",
            jetflavour  = 'AK4PFchs',
            ):

        _iovs = []
        for runnumber, filename in iovs:

            txtfiles = []

            # Download file
            source = os.path.join( baseurl, "%s.%s"%(filename, JetCorrector.extension) )
            target = os.path.join( os.path.expandvars( JetCorrector.data_directory ), "%s.%s"%(filename, JetCorrector.extension) )

            for level in correctionLevels:
                txtfile = os.path.join( os.path.expandvars( JetCorrector.data_directory ), "%s_%s_%s.txt"%( filename, level, jetflavour)  )

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
                            if not member.name.endswith( jetflavour+'.txt' ): continue
                            member_filename = os.path.basename( member.name )
                            logger.debug( "Found file %s in %s", member_filename, target )
                            with file( os.path.join( os.path.expandvars( JetCorrector.data_directory ), member_filename), 'w') as f_out:
                                f_out.writelines( tar.extractfile( member ).readlines() )

                logger.debug( "Adding txtfile %s", txtfile )
                txtfiles.append( txtfile )

            _iovs.append( ( runnumber, txtfiles ) )

        return cls( _iovs )

    def reduceLevels( self, correctionLevels ):
        self.makeNewCorrectors( require = correctionLevels )
        return self 

    def correction(self, rawPt, eta, area, rho, run ):

        for runnumber, corrector in reversed( self.jetCorrectors ):
            if run >= runnumber:
                corrector.setJetPt( rawPt )
                corrector.setJetEta( eta )
                corrector.setJetA( area )
                corrector.setRho( rho )
                return corrector.getCorrection()
