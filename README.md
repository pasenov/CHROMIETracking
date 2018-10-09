# CHROMIETracking


# to install and compile the code
#on lxplus

cmsrel CMSSW_10_1_2

cd CMSSW_10_1_2/src/

cmsenv

git clone https://github.com/jandrea/CHROMIETracking.git

mv CHROMIETracking/* .

scramv1 b -j4


# to run the DQM plot producer

cd Geometry/PixelTelescope/test

cmsRun PixelTelescope_BeamData_RECO_jeremy.py

