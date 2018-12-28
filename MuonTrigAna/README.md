- login to lxplus

- download code from gitHub:

- setup standard ATLAS sofwtare release:

{noformat}
cd MuonTrigAna
mkdir build
mkdir run
cd build/
asetup 21.2.38,AnalysisBase,here
cmake ../source/
make
source x86_64-*/setup.sh
{noformat}

- NOTE: after first build, next times only needs:

{noformat}
cd build/
setup --restore
{noformat}

- run the code: 

-- with root:

{noformat}
cd run/
rm -rf pippo
root -b -q '$ROOTCOREDIR/scripts/load_packages.C' '../source/TARun.cxx("pippo")'
{noformat}

-- with python: 

{noformat}
cd run/
rm -rf pippo
../source/TARun.py --submission-dir='pippo'
{noformat}
