log_ptll_thresholds     = [10**(x/10.) for x in range(11,36)]
default_ptll_thresholds = [30, 40, 50, 60, 85, 105, 130, 175, 230, 300, 400, 500, 700, 1000, 1500]

ptll_thresholds = default_ptll_thresholds
ptll_thresholds.sort()

ptll_bins    = [ (ptll_thresholds[0], -1) ] + [(ptll_thresholds[i], ptll_thresholds[i+1]) for i in xrange(len(ptll_thresholds)-1)] + [ (ptll_thresholds[-1], -1) ]
abs_eta_bins = [ (0, 5.2), (0, 0.8), (0.8, 1.3), (1.3, 1.9), (1.9, 2.5), (2.5, 3), (3, 3.2), (3.2, 5.2) ]
