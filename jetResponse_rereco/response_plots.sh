#!/bin/sh
python jetResponseEraPlot.py $1 &
python jetResponsePlot.py --era Run2016B $1 &
python jetResponsePlot.py --era Run2016C $1 &
python jetResponsePlot.py --era Run2016D $1 &
python jetResponsePlot.py --era Run2016E $1 &
python jetResponsePlot.py --era Run2016F $1 &
python jetResponsePlot.py --era Run2016F_early $1 &
python jetResponsePlot.py --era Run2016F_late $1 &
python jetResponsePlot.py --era Run2016G $1 &
