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

    def __init__( self, 
            iovs, 
            correctionLevels, 
            baseurl     = "https://github.com/cms-jet/JECDatabase/raw/master/tarballs/",
            jetflavour  = 'AK4PFchs',
            directory   = "$CMSSW_BASE/src/JetMET/JetCorrector/data/",
            ):

        self.extension   = "tar.gz"
        self.iovs        = iovs
        self.correctionLevels = correctionLevels
        self.baseurl     = baseurl
        self.jetflavour  = jetflavour
        self.directory   = directory

        self.makeCorrectors()

    def makeCorrectors( self ):
        self.jetCorrectors = []
        for runnumber, filename in self.iovs:
            # Download file
            source = os.path.join( self.baseurl, "%s.%s"%(filename, self.extension) )
            target = os.path.join( os.path.expandvars( self.directory ), "%s.%s"%(filename, self.extension) )

            params = ROOT.vector(ROOT.JetCorrectorParameters)()
            for level in self.correctionLevels:
                txtfile = os.path.join( os.path.expandvars( self.directory ), "%s_%s_%s.txt"%( filename, level, self.jetflavour)  )

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
                            if not member.name.endswith( self.jetflavour+'.txt' ): continue
                            member_filename = os.path.basename( member.name )
                            logger.debug( "Found file %s in %s", member_filename, target )
                            with file( os.path.join( os.path.expandvars( self.directory ), member_filename), 'w') as f_out:
                                f_out.writelines( tar.extractfile( member ).readlines() )

                logger.debug( "Adding runnumber %i txtfile %s", runnumber, txtfile )

                params.push_back( ROOT.JetCorrectorParameters( txtfile, "" ) )

            # Make corrector
            self.jetCorrectors.append( ( runnumber,  ROOT.FactorizedJetCorrector( params )) )

        # Sort wrt IOVs
        self.jetCorrectors.sort( key = lambda p: p[0] )

    
    def fromLevels( self, correctionLevels ):
        return JetCorrector( iovs = self.iovs, correctionLevels = correctionLevels, baseurl = self.baseurl, jetflavour = self.jetflavour, directory = self.directory )

    def correction(self, rawPt, eta, area, rho, run ):

        for runnumber, corrector in reversed( self.jetCorrectors ):
            if run >= runnumber:
                corrector.setJetPt( rawPt )
                corrector.setJetEta( eta )
                corrector.setJetA( area )
                corrector.setRho( rho )
                return corrector.getCorrection()
