from JetMET.JetCorrector.JetCorrector import JetCorrector 
from JetMET.JetCorrector.JetCorrector import correction_levels_data
from JetMET.JetCorrector.JetCorrector import correction_levels_mc

# config
Spring16_23Sep2016_DATA = \
((1, 'Spring16_23Sep2016BCDV2_DATA'),
 (276831, 'Spring16_23Sep2016EFV2_DATA'),
 (278802, 'Spring16_23Sep2016GV2_DATA'),
 (280919, 'Spring16_23Sep2016HV2_DATA'))

Spring16_23Sep2016_MC = ((1, 'Spring16_23Sep2016V2_MC'),)

if __name__ == "__main__":

    # Logging
    import JetMET.tools.logger as logger
    logger  = logger.get_logger('INFO', logFile = None)

    
jetCorrector_data = JetCorrector( Spring16_23Sep2016_DATA, correctionLevels = correction_levels_data ) 
jetCorrector_mc   = JetCorrector( Spring16_23Sep2016_MC, correctionLevels = correction_levels_mc ) 

