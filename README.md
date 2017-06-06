# JetMET/JERC analysis code and CMG setup in 80X

```
cmsrel CMSSW_9_0_2
cd CMSSW_9_0_2/src
cmsenv
git clone -b PF_jets_comp git@github.com:spandeyehep/JetMET.git
git clone https://github.com/schoef/RootTools.git
<edit output directory https://github.com/spandeyehep/JetMET/blob/PF_jets_comp/tools/python/user.py#L11-L13>
scram b
cmsRun JetMET/response/python/jetResponse_PFHadronCalibration/jetResponse_PFHadronCalibration.py
```
