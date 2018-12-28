void TARun (const std::string& submitDir)
{
  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  // create a new sample handler to describe the data files we use
  SH::SampleHandler sh;

  // scan for datasets in the given directory
  // this works if you are on lxplus, otherwise you'd want to copy over files
  // to your local machine and use a local path.  if you do so, make sure
  // that you copy all subdirectories and point this to the directory
  // containing all the files, not the subdirectories.

  // use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
  const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/g/giagu/public/mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.AOD.e3601_s3126_r10114_r10148");
  SH::ScanDir().filePattern("AOD.12631155._000001.pool.root.1").scan(sh,inputFilePath);


  // set the name of the tree in our files
  // in the xAOD the TTree containing the EDM containers is "CollectionTree"
  sh.setMetaString ("nc_tree", "CollectionTree");

  // further sample handler configuration may go here

  // print out the samples we found
  sh.print ();

  // this is the basic description of our job
  EL::Job job;
  job.sampleHandler (sh); // use SampleHandler in this job
  job.options()->setDouble (EL::Job::optMaxEvents, 500); // for testing purposes, limit to run over the first 500 events only!

  // add our algorithm to the job
  EL::AnaAlgorithmConfig config;
  config.setType ("TrigxAODAnalysis");

  // set the name of the algorithm (this is the name use with
  // messages)
  config.setName ("TrigAnalysisAlg");

  // configuration options for algorithm go here
  std::vector<std::string> chains;

  // std trigger to be validated:
  // L1: L1_MU20, L1_MU21
  // HLT single iso: HLT_mu26_ivarmedium, HLT_mu28_ivarmedium, HLT_mu26_ivartight, HLT_mu28_ivartight
  // HLT single: HLT_mu50, HLT_mu60, HLT_mu80, HLT_mu60_msonly_3layerEC, HLT_mu80_msonly_3layerEC
  // HLT multi: HLT_mu22_mu8noL1, HLT_mu20_2mu4noL1, HLT_2mu14, HLT_3mu6_msonly
  chains.push_back("HLT_mu24_ivarmedium");
  chains.push_back("HLT_mu50");
  chains.push_back("HLT_2mu14");
  chains.push_back("HLT_3mu6_msonly");
  config.setProperty( "ChainName", chains ).ignore();
  config.setProperty( "MuonQuality", 1 ).ignore();
  config.setProperty( "Verbose", 2 ).ignore();

  job.algsAdd (config);

  // output ntuple
  job.outputAdd (EL::OutputStream ("ANALYSIS"));

  // make the driver we want to use:
  // this one works by running the algorithm directly:
  EL::DirectDriver driver;
  // we can use other drivers to run things on the Grid, with PROOF, etc.

  // process the job using the driver
  driver.submit (job, submitDir);
}
