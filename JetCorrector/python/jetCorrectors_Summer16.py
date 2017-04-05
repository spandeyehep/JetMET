from JetMET.JetCorrector.JetCorrector import JetCorrector 
from JetMET.JetCorrector.JetCorrector import correction_levels_data
from JetMET.JetCorrector.JetCorrector import correction_levels_mc

# config
Summer16_23Sep2016_DATA = \
[(1, 'Summer16_23Sep2016BCDV3_DATA'),
 (276831, 'Summer16_23Sep2016EFV3_DATA'),
 (278802, 'Summer16_23Sep2016GV3_DATA'),
 (280919, 'Summer16_23Sep2016HV3_DATA')]

Summer16_23Sep2016_MC = [(1, 'Summer16_23Sep2016V3_MC') ]

if __name__ == "__main__":

    # Logging
    import JetMET.tools.logger as logger
    logger  = logger.get_logger('DEBUG', logFile = None)

    
jetCorrector_data = JetCorrector.fromTarBalls( Summer16_23Sep2016_DATA, correctionLevels = correction_levels_data ) 
jetCorrector_mc   = JetCorrector.fromTarBalls( Summer16_23Sep2016_MC,   correctionLevels = correction_levels_mc )
