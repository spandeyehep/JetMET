# aequidistant in log 
log_ptll_thresholds     = [10**(x/10.) for x in range(11,36)]

# L3res
default_ptll_thresholds = [30, 40, 50, 60, 85, 105, 130, 175, 230, 300, 400, 500, 700, 1000, 1500]

ptll_thresholds = default_ptll_thresholds
ptll_thresholds.sort()

ptll_bins    = [ (ptll_thresholds[0], -1) ] + [(ptll_thresholds[i], ptll_thresholds[i+1]) for i in xrange(len(ptll_thresholds)-1)] + [ (ptll_thresholds[-1], -1) ]
coarse_ptll_bins    = [ (30, -1), (100, -1), (30, 100), (100,200), (200,500), (500, 1000), (1000, -1) ]

#abs_eta_bins = [ (0, 5.2), (0, 0.8), (0.8, 1.3), (1.3, 1.9), (1.9, 2.5), (2.5, 3), (3, 3.2), (3.2, 5.2) ]
abs_eta_bins = [ (0, 5.2), (0, 0.8), (0, 1.3), (0.8, 1.3), (1.3, 1.9), (1.9, 2.5), (2.5, 3), (3, 3.2), (3.2, 5.2) ]
coarse_abs_eta_bins = [ (0, 1.3), ( 1.3, 2.5 ), (2.5, 3), (3, 3.2), (3.2, 5.2) ]

# L2res
L2res_abs_eta_thresholds = [0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191]
L2res_abs_eta_bins       = [(L2res_abs_eta_thresholds[i], L2res_abs_eta_thresholds[i+1]) for i in range( len( L2res_abs_eta_thresholds ) -1 ) ]
