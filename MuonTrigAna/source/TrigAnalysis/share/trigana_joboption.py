import AthenaRootComps.ReadAthenaxAODHybrid

theApp.EvtMax = 500
testFile = '/afs/cern.ch/work/g/giagu/public/mc16_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.AOD.e3601_s3126_r10114_r10148/AOD.12631155._000001.pool.root.1'
svcMgr.EventSelector.InputCollections = [testFile]

algSeq = CfgMgr.AthSequencer("AthAlgSeq")

# create our algorithm with teh given name
alg = CfgMgr.TrigxAODAnalysis()

# later on we'll add some configuration options for our algorithm that go here

# // std trigger to be validated:
# // L1: L1_MU20, L1_MU21
# // HLT single iso: HLT_mu26_ivarmedium, HLT_mu28_ivarmedium, HLT_mu26_ivartight, HLT_mu28_ivartight
# // HLT single: HLT_mu50, HLT_mu60, HLT_mu80, HLT_mu60_msonly_3layerEC, HLT_mu80_msonly_3layerEC
# // HLT multi: HLT_mu22_mu8noL1, HLT_mu20_2mu4noL1, HLT_2mu14, HLT_3mu6_msonly
alg.ChainName = ['HLT_mu24_ivarmedium','HLT_mu50','HLT_2mu14','HLT_3mu6_msonly']
alg.MuonQuality=1
alg.Verbose=2

algSeq += alg

ServiceMgr += CfgMgr.THistSvc()
ServiceMgr.THistSvc.Output += [
    "ANALYSIS DATAFILE='TrigxAODAnalysis.outputs.root' OPT='RECREATE'"
        ]

# optional include for reducing printout from athena
include("AthAnalysisBaseComps/SuppressLogging.py")
