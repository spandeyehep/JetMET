# JetMET/JERC analysis code and CMG setup in 80X

```
cmsrel CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src
cmsenv
git cms-init
# Analysis code
git clone https://github.com/schoef/JetMET
# CMG + Heppy setup for ntuple production
./JetMET/setup/setup.sh
```
