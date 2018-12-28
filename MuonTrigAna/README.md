- login to lxplus

- download code from gitHub:

- setup standard ATLAS sofwtare release:

cd MuonTrigAna
mkdir build
mkdir run

cd build/
asetup 21.2.38,AnalysisBase,here
cmake ../source/
make
source x86_64-*/setup.sh

- NOTE: after first build, next times only needs:

cd build/
setup --restore

- run the code: 

-- with root:

cd run/
rm -rf pippo
root -b -q '$ROOTCOREDIR/scripts/load_packages.C' '../source/TARun.cxx("pippo")'

-- with python: 

cd run/
rm -rf pippo
../source/TARun.py --submission-dir='pippo'
