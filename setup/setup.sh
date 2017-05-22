# 
# Recipe to continue the setup of our 80X analysis after the checkout of StopsDilepton package
#
eval `scram runtime -sh`
cd $CMSSW_BASE/src

git clone git@github.com:schoef/RootTools

#
# Setting up CMG
#
git remote add cmg-central https://github.com/CERN-PH-CMG/cmg-cmssw.git -f -t heppy_80X
cp $CMSSW_BASE/src/JetMET/setup/heppy_sparse-checkout $CMSSW_BASE/src/.git/info/sparse-checkout
git checkout -b heppy_80X cmg-central/heppy_80X

# add your mirror, and push the 80X branch to it
git remote add origin git@github.com:schoef/cmg-cmssw.git
git push -u origin heppy_80X

# now get the CMGTools subsystem from the cmgtools-lite repository
git clone -o cmg-central https://github.com/CERN-PH-CMG/cmgtools-lite.git -b 80X CMGTools
cd CMGTools 
# only take what we need
cp $CMSSW_BASE/src/JetMET/setup/cmgtools_sparse-checkout $CMSSW_BASE/src/CMGTools/.git/info/sparse-checkout
git config core.sparsecheckout True
git read-tree -mu HEAD

# add your fork, and push the 80X branch to it
git remote add origin git@github.com:schoef/cmgtools-lite.git 
git push -u origin 80X_JetMET

cd $CMSSW_BASE/src
git fetch origin
#git checkout -b 80X_JetMET origin/80X_JetMET

#compile
cd $CMSSW_BASE/src && scram b -j 8 
