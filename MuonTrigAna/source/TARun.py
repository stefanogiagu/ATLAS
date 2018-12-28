#!/usr/bin/env python

# Read the submission directory as a command line argument. You can
# extend the list of arguments with your private ones later on.
import optparse
parser = optparse.OptionParser()
parser.add_option( '-s', '--submission-dir', dest = 'submission_dir',
                   action = 'store', type = 'string', default = 'submitDir',
                   help = 'Submission directory for EventLoop' )
( options, args ) = parser.parse_args()

# Set up (Py)ROOT.
import ROOT
ROOT.xAOD.Init().ignore()

# Set up the sample handler object. See comments from the C++ macro
# for the details about these lines.
import os
sh = ROOT.SH.SampleHandler()
sh.setMetaString( 'nc_tree', 'CollectionTree' )
inputFilePath = '/afs/cern.ch/work/g/giagu/public/mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.AOD.e3601_s3126_r10114_r10148'
ROOT.SH.ScanDir().filePattern( 'AOD.12631155._000001.pool.root.1' ).scan( sh, inputFilePath )
sh.printContent()

# Create an EventLoop job.
job = ROOT.EL.Job()
job.sampleHandler( sh )
job.options().setDouble( ROOT.EL.Job.optMaxEvents, 500 )

# Create the algorithm's configuration. Note that we'll be able to add
# algorithm property settings here later on.
from AnaAlgorithm.AnaAlgorithmConfig import AnaAlgorithmConfig
#std trigger to be validated:
#L1: L1_MU20, L1_MU21
#HLT single iso: HLT_mu26_ivarmedium, HLT_mu28_ivarmedium, HLT_mu26_ivartight, HLT_mu28_ivartight
#HLT single: HLT_mu50, HLT_mu60, HLT_mu80, HLT_mu60_msonly_3layerEC, HLT_mu80_msonly_3layerEC
#HLT multi: HLT_mu22_mu8noL1, HLT_mu20_2mu4noL1, HLT_2mu14, HLT_3mu6_msonly
config = AnaAlgorithmConfig( 'TrigxAODAnalysis/TrigAnalysisAlg',
                             ChainName=['HLT_mu24_ivarmedium','HLT_mu50','HLT_2mu14','HLT_3mu6_msonly'],
                             MuonQuality=1,
                             Verbose=2)

job.algsAdd( config )

job.outputAdd (ROOT.EL.OutputStream ('ANALYSIS'))

# Run the job using the direct driver.
driver = ROOT.EL.DirectDriver()
driver.submit( job, options.submission_dir )
