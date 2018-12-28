#include <AsgTools/MessageCheck.h>
#include <TrigAnalysis/TrigxAODAnalysis.h>

#include <xAODEventInfo/EventInfo.h>
#include <xAODTruth/TruthEventContainer.h>
#include <xAODTruth/TruthParticleContainer.h>
#include <xAODTruth/TruthVertexContainer.h>
#include <xAODMuon/MuonContainer.h>
#include <xAODMuon/Muon.h>

#include <xAODTrigger/MuonRoI.h>
#include <xAODTrigger/MuonRoIContainer.h>

#include <xAODTrigMuon/L2StandAloneMuon.h>
#include <xAODTrigMuon/L2StandAloneMuonContainer.h>
#include <xAODTrigMuon/L2CombinedMuon.h>
#include <xAODTrigMuon/L2CombinedMuonContainer.h>

#include <TrigDecisionTool/Feature.h>



TrigxAODAnalysis :: TrigxAODAnalysis (const std::string& name, ISvcLocator *pSvcLocator) : 
                                      EL::AnaAlgorithm (name, pSvcLocator), 
                                      m_muonSelection ("CP::MuonSelectionTool", this),
                                      m_trigConfigTool("TrigConf::xAODConfigTool/xAODConfigTool"),
                                      m_trigDecisionTool ("Trig::TrigDecisionTool/TrigDecisionTool")
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0. Note that things like resetting
  // statistics variables should rather go into the initialize() function.
  declareProperty( "ChainName", m_chainName, "Muon trigger chains to be analysed" );
  declareProperty( "MuonQuality", m_muonQuality = 1,
                    "Offline Muon Quality" );
  declareProperty( "Verbose", m_verbose = 2, ">0 to activate verbose printout" );

  declareProperty("DeltaR_L1_truth", m_dr_L1_truth = 0.15, "DeltaR match L1-truth" );
  declareProperty("DeltaR_L2SA_truth", m_dr_L2SA_truth = 0.1, "DeltaR match L2SA-truth" );
  declareProperty("DeltaR_L2CB_truth", m_dr_L2CB_truth = 0.1, "DeltaR match L2CB-truth" );
  declareProperty("DeltaR_EF_truth", m_dr_EF_truth = 0.1, "DeltaR match EF-truth" );
  declareProperty("DeltaR_L1_offl", m_dr_L1_offl = 0.15, "DeltaR match L1-offline" );
  declareProperty("DeltaR_L2SA_offl", m_dr_L2SA_offl = 0.1, "DeltaR match L2SA-offline" );
  declareProperty("DeltaR_L2CB_offl", m_dr_L2CB_offl = 0.1, "DeltaR match L2CB-offline" );
  declareProperty("DeltaR_EF_offl", m_dr_EF_offl = 0.05, "DeltaR match EF-offline" );
  declareProperty("DeltaR_EF_L1", m_dr_EF_L1 = 0.15, "DeltaR match EF-offline" );
  declareProperty("DeltaR_EF_L2SA", m_dr_EF_L2SA = 0.1, "DeltaR match EF-L2SA" );
  declareProperty("DeltaR_EF_L2CB", m_dr_EF_L2CB = 0.1, "DeltaR match EF-L2CB" );

}

TrigxAODAnalysis::~TrigxAODAnalysis() {

   if (m_muonTruthEta) delete m_muonTruthEta;
   if (m_muonTruthPhi) delete m_muonTruthPhi;
   if (m_muonTruthPt)  delete m_muonTruthPt;
   if (m_muonTruthE)   delete m_muonTruthE;
   if (m_muonTruthQ)   delete m_muonTruthQ;

   if (m_ZTruthEta)    delete   m_ZTruthEta;
   if (m_ZTruthPhi)    delete   m_ZTruthPhi;
   if (m_ZTruthPt)     delete   m_ZTruthPt;
   if (m_ZTruthE)      delete   m_ZTruthE;
   if (m_JpsiTruthEta) delete   m_JpsiTruthEta;
   if (m_JpsiTruthPhi) delete   m_JpsiTruthPhi;
   if (m_JpsiTruthPt)  delete   m_JpsiTruthPt;
   if (m_JpsiTruthE)   delete   m_JpsiTruthE;

   if (m_muonEta) delete m_muonEta;
   if (m_muonPhi) delete m_muonPhi;
   if (m_muonPt)  delete m_muonPt;
   if (m_muonE)   delete m_muonE;
   if (m_muonQ)   delete m_muonQ;
   if (m_muonIso) delete m_muonIso;

   if (m_Z_id1) delete m_Z_id1;
   if (m_Z_id2) delete m_Z_id2;
   if (m_Z_m) delete m_Z_m;
   if (m_Jpsi_id1) delete m_Jpsi_id1;
   if (m_Jpsi_id2) delete m_Jpsi_id2;
   if (m_Jpsi_m) delete m_Jpsi_m;

   if (m_l1roiEta)   delete m_l1roiEta;
   if (m_l1roiPhi)   delete m_l1roiPhi;
   if (m_l1roiPtThr) delete m_l1roiPtThr;
   if (m_l1roiRoI)   delete m_l1roiRoI;
   if (m_l1roiMatchTruth)  delete m_l1roiMatchTruth;
   if (m_l1roiMatchOffl)   delete m_l1roiMatchOffl;

   if (m_chainNameTrg) delete m_chainNameTrg;
   if (m_chainPassTrg) delete m_chainPassTrg;
   if (m_chainL1PassTrg) delete m_chainL1PassTrg;
   if (m_chainL1PrescaledTrg)   delete m_chainL1PrescaledTrg;
   if (m_chainPSTrg)   delete m_chainPSTrg;

   if (m_l2saEta)    delete m_l2saEta;
   if (m_l2saPhi)    delete m_l2saPhi;
   if (m_l2saPt)     delete m_l2saPt;
   if (m_l2saQ)      delete m_l2saQ;
   if (m_l2saRoI)    delete m_l2saRoI;
   if (m_l2saRoId)   delete m_l2saRoId;
   if (m_l2saMatchTruth)   delete m_l2saMatchTruth;
   if (m_l2saMatchOffl)   delete m_l2saMatchOffl;
   if (m_l2saMatchL1)   delete m_l2saMatchL1;

   if (m_l2cbEta)    delete m_l2cbEta;
   if (m_l2cbPhi)    delete m_l2cbPhi;
   if (m_l2cbPt)     delete m_l2cbPt;
   if (m_l2cbQ)      delete m_l2cbQ;
   if (m_l2cbRoI)    delete m_l2cbRoI;
   if (m_l2cbRoId)   delete m_l2cbRoId;
   if (m_l2cbMatchTruth)  delete m_l2cbMatchTruth;
   if (m_l2cbMatchOffl)   delete m_l2cbMatchOffl;
   if (m_l2cbMatchL1)     delete m_l2cbMatchL1;
   if (m_l2cbMatchL2SA)   delete m_l2cbMatchL2SA;

   if (m_efmuEta) delete m_efmuEta;
   if (m_efmuPhi) delete m_efmuPhi;
   if (m_efmuPt)  delete m_efmuPt;
   if (m_efmuQ)   delete m_efmuQ;
   if (m_efmuMatchTruth)  delete m_efmuMatchTruth;
   if (m_efmuMatchOffl)   delete m_efmuMatchOffl;
   if (m_efmuMatchL1)     delete m_efmuMatchL1;
   if (m_efmuMatchL2SA)   delete m_efmuMatchL2SA;
   if (m_efmuMatchL2CB)   delete m_efmuMatchL2CB;

}

int TrigxAODAnalysis :: match_dr(std::vector<float> eta0, std::vector<float> phi0, float eta, float phi, float dr) 
{
  // utility function for DeltaR match
  int n_entries = eta0.size();
  int idx = -1;
  float mindr = 99999999999999.;
  for (int j=0; j<n_entries; ++j) {
    float deta = (eta0.at(j)-eta);
    float dphi = (phi0.at(j)-phi);
    if(dphi > 3.1415927) dphi = 2*3.1415927 - dphi;
    float DeltaR = sqrt(pow(deta,2)+pow(dphi,2));
    ANA_MSG_INFO ("Match parameters: j,eta0,phi0,eta,phi,dr: " << j << " / " << eta0.at(j) << " / " << phi0.at(j) << " / " << eta << " / " << phi << " / " << DeltaR);
    if (DeltaR < dr) {
      if (DeltaR < mindr) {
        mindr = DeltaR;
        idx = j;
      }
    }
  }
  ANA_MSG_INFO ("Match idx: " << idx);
  return idx;
}


// J/psi and Z search
void TrigxAODAnalysis :: match_zjpsi(std::vector<float> eta, std::vector<float> phi, std::vector<float> pt, std::vector<float> Q, std::vector<float> iso, std::vector<unsigned int> & idZ1, std::vector<unsigned int> & idZ2, std::vector<float> & mZ, std::vector<unsigned int> & idJpsi1, std::vector<unsigned int> & idJpsi2, std::vector<float> & mJpsi ) 
{
  float PI = 3.1415927;
  float muMass = 0.1056583745;
  float ZMass = 91.1876;
  float DeltaZMass = 10.0;
  float isoZmu = 0.15;
  float ptmuZ = 20.0;
  float dphiZ = PI*(1. - 10./180);
  float dRJpsi_min = 0.2;
  float dRJpsi_max = 3.5;
  float ptmuJpsi = 4.0;
  float JpsiMass = 3.096900;
  float DeltaJpsiMass = 0.3;
  for (unsigned int i=0; i<eta.size(); ++i) {
     if (pt.at(i) < ptmuJpsi) continue;
     TLorentzVector m1;
     m1.SetPtEtaPhiM(pt.at(i),eta.at(i),phi.at(i), muMass);
     for (unsigned int j=i+1; j<eta.size(); ++j) {
        if (pt.at(j) < ptmuJpsi) continue;
        if (Q.at(i)*Q.at(j) > 0) continue;
        TLorentzVector m2;
        m2.SetPtEtaPhiM(pt.at(j),eta.at(j),phi.at(j), muMass);
        if (fabs((m1+m2).M() - ZMass) < DeltaZMass) { // Z Match
          if (m1.DeltaPhi(m2) > dphiZ) { // back to back
            if (pt.at(i) > ptmuZ && pt.at(j) > ptmuZ) { // high pt
               if (iso.at(i) < isoZmu) { // one muon isolated
                  // found a Z
                  idZ1.push_back(i);
                  idZ2.push_back(j);
                  mZ.push_back((m1+m2).M());
               } else if (iso.at(j) < isoZmu) {
                  idZ1.push_back(j);
                  idZ2.push_back(i);
                  mZ.push_back((m1+m2).M());
               }
            }
          }
        }
        if (fabs((m1+m2).M() - JpsiMass) < DeltaJpsiMass) { // J/psi Match
           if (m1.DeltaR(m2) > dRJpsi_min) { // DeltaR > 0.2 (to avoid overlaps in same ROI)
              if (m1.DeltaR(m2) < dRJpsi_max) { // DeltaR < 3.5
                 // found a J/psi
                 idJpsi1.push_back(i);
                 idJpsi2.push_back(j);
                 mJpsi.push_back((m1+m2).M());
              }
           }
        }
     } 
  }

}

StatusCode TrigxAODAnalysis :: initialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  ANA_MSG_INFO ("in initialize");

  ANA_CHECK (book (TTree ("muontrig", "Muon Trigger analysis ntuple")));

  TTree* mytree = tree ("muontrig");
  mytree->Branch ("RunNumber", &m_runNumber);
  mytree->Branch ("EventNumber", &m_eventNumber);

  m_muonTruthEta = new std::vector<float>();
  mytree->Branch ("MuonTruthEta", &m_muonTruthEta);
  m_muonTruthPhi = new std::vector<float>();
  mytree->Branch ("MuonTruthPhi", &m_muonTruthPhi);
  m_muonTruthPt = new std::vector<float>();
  mytree->Branch ("MuonTruthPt", &m_muonTruthPt);
  m_muonTruthE = new std::vector<float>();
  mytree->Branch ("MuonTruthE", &m_muonTruthE);
  m_muonTruthQ = new std::vector<float>();
  mytree->Branch ("MuonTruthQ", &m_muonTruthQ);

  m_ZTruthEta = new std::vector<float>();
  mytree->Branch ("ZTruthEta", &m_ZTruthEta);
  m_ZTruthPhi = new std::vector<float>();
  mytree->Branch ("ZTruthPhi", &m_ZTruthPhi);
  m_ZTruthPt = new std::vector<float>();
  mytree->Branch ("ZTruthPt", &m_ZTruthPt);
  m_ZTruthE = new std::vector<float>();
  mytree->Branch ("ZTruthE", &m_ZTruthE);
  m_JpsiTruthEta = new std::vector<float>();
  mytree->Branch ("JpsiTruthEta", &m_JpsiTruthEta);
  m_JpsiTruthPhi = new std::vector<float>();
  mytree->Branch ("JpsiTruthPhi", &m_JpsiTruthPhi);
  m_JpsiTruthPt = new std::vector<float>();
  mytree->Branch ("JpsiTruthPt", &m_JpsiTruthPt);
  m_JpsiTruthE = new std::vector<float>();
  mytree->Branch ("JpsiTruthE", &m_JpsiTruthE);

  m_muonEta = new std::vector<float>();
  mytree->Branch ("MuonEta", &m_muonEta);
  m_muonPhi = new std::vector<float>();
  mytree->Branch ("MuonPhi", &m_muonPhi);
  m_muonPt = new std::vector<float>();
  mytree->Branch ("MuonPt", &m_muonPt);
  m_muonE = new std::vector<float>();
  mytree->Branch ("MuonE", &m_muonE);
  m_muonQ = new std::vector<float>();
  mytree->Branch ("MuonQ", &m_muonQ);
  m_muonIso = new std::vector<float>();
  mytree->Branch ("MuonIso", &m_muonIso);

  
  m_Z_id1 = new std::vector<unsigned int>();
  mytree->Branch ("Z_id1", &m_Z_id1);
  m_Z_id2 = new std::vector<unsigned int>();
  mytree->Branch ("Z_id2", &m_Z_id2);
  m_Z_m = new std::vector<float>();
  mytree->Branch ("Z_m", &m_Z_m);
  m_Jpsi_id1 = new std::vector<unsigned int>();
  mytree->Branch ("Jpsi_id1", &m_Jpsi_id1);
  m_Jpsi_id2 = new std::vector<unsigned int>();
  mytree->Branch ("Jpsi_id2", &m_Jpsi_id2);
  m_Jpsi_m = new std::vector<float>();
  mytree->Branch ("Jpsi_m", &m_Jpsi_m);

  m_l1roiEta = new std::vector<float>();
  mytree->Branch ("L1RoIEta", &m_l1roiEta);
  m_l1roiPhi = new std::vector<float>();
  mytree->Branch ("L1RoIPhi", &m_l1roiPhi);
  m_l1roiPtThr = new std::vector<float>();
  mytree->Branch ("L1RoIPtThr", &m_l1roiPtThr);
  m_l1roiRoI = new std::vector<unsigned int>();
  mytree->Branch ("L1RoI", &m_l1roiRoI);
  m_l1roiMatchTruth = new std::vector<unsigned int>();
  mytree->Branch ("L1RoIMatchTruth", &m_l1roiMatchTruth);
  m_l1roiMatchOffl = new std::vector<unsigned int>();
  mytree->Branch ("L1RoIMatchOffl", &m_l1roiMatchOffl);

  m_chainNameTrg = new std::vector<std::string>();
  mytree->Branch ("ChainNameTrg", &m_chainNameTrg);
  m_chainPassTrg = new std::vector<int>(); 
  mytree->Branch ("ChainPassTrg", &m_chainPassTrg);
  m_chainL1PassTrg = new std::vector<int>(); 
  mytree->Branch ("ChainL1PassTrg", &m_chainL1PassTrg);
  m_chainPSTrg = new std::vector<float>(); 
  mytree->Branch ("ChainPSTrg", &m_chainPSTrg);
  m_chainL1PrescaledTrg = new std::vector<int>(); 
  mytree->Branch ("ChainL1PrescaledTrg", &m_chainL1PrescaledTrg);

  m_l2saEta = new std::vector<std::vector<float>>();
  mytree->Branch ("L2SAEta", &m_l2saEta);
  m_l2saPhi = new std::vector<std::vector<float>>();
  mytree->Branch ("L2SAPhi", &m_l2saPhi);
  m_l2saPt = new std::vector<std::vector<float>>();
  mytree->Branch ("L2SAPt", &m_l2saPt);
  m_l2saQ = new std::vector<std::vector<float>>();
  mytree->Branch ("L2SAQ", &m_l2saQ);
  m_l2saRoI = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("L2SARoI", &m_l2saRoI);
  m_l2saRoId = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("L2SARoId", &m_l2saRoId);
  m_l2saMatchTruth = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("L2SAMatchTruth", &m_l2saMatchTruth);
  m_l2saMatchOffl = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("L2SAMatchOffl", &m_l2saMatchOffl);
  m_l2saMatchL1 = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("L2SAMatchL1", &m_l2saMatchL1);
  
  m_l2cbEta = new std::vector<std::vector<float>>();
  mytree->Branch ("L2CBEta", &m_l2cbEta);
  m_l2cbPhi = new std::vector<std::vector<float>>();
  mytree->Branch ("L2CBPhi", &m_l2cbPhi);
  m_l2cbPt = new std::vector<std::vector<float>>();
  mytree->Branch ("L2CBPt", &m_l2cbPt);
  m_l2cbQ = new std::vector<std::vector<float>>();
  mytree->Branch ("L2CBQ", &m_l2cbQ);
  m_l2cbRoI = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("L2CBRoI", &m_l2cbRoI);
  m_l2cbRoId = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("L2CBRoId", &m_l2cbRoId);
  m_l2cbMatchTruth = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("L2CBMatchTruth", &m_l2cbMatchTruth);
  m_l2cbMatchOffl = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("L2CBMatchOffl", &m_l2cbMatchOffl);
  m_l2cbMatchL1 = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("L2CBMatchL1", &m_l2cbMatchL1);
  m_l2cbMatchL2SA = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("L2CBMatchL2SA", &m_l2cbMatchL2SA);

  m_efmuEta = new std::vector<std::vector<float>>();
  mytree->Branch ("EFmuEta", &m_efmuEta);
  m_efmuPhi = new std::vector<std::vector<float>>();
  mytree->Branch ("EFmuPhi", &m_efmuPhi);
  m_efmuPt = new std::vector<std::vector<float>>();
  mytree->Branch ("EFmuPt", &m_efmuPt);
  m_efmuQ = new std::vector<std::vector<float>>();
  mytree->Branch ("EFmuQ", &m_efmuQ);
  m_efmuMatchTruth = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("EFmuMatchTruth", &m_efmuMatchTruth);
  m_efmuMatchOffl = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("EFmuMatchOffl", &m_efmuMatchOffl);
  m_efmuMatchL1 = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("EFmuMatchL1", &m_efmuMatchL1);
  m_efmuMatchL2SA = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("EFmuMatchL2SA", &m_efmuMatchL2SA);
  m_efmuMatchL2CB = new std::vector<std::vector<unsigned int>>();
  mytree->Branch ("EFmuMatchL2CB", &m_efmuMatchL2CB);

  // Properties
  ANA_MSG_INFO( "Muon trigger chains analised" );
  for (auto &chName : m_chainName) {
     ANA_MSG_INFO( "ChainName = " << chName );
  }
  ANA_MSG_INFO( "MuonQuality  = " << m_muonQuality );
  ANA_MSG_INFO( "Verbose Output level  = " << m_verbose );

  ANA_MSG_INFO( "DeltaR_L1_truth = " << m_dr_L1_truth);
  ANA_MSG_INFO( "DeltaR_L2SA_truth = " << m_dr_L2SA_truth);
  ANA_MSG_INFO( "DeltaR_L2CB_truth = " << m_dr_L2CB_truth);
  ANA_MSG_INFO( "DeltaR_EF_truth = " << m_dr_EF_truth);
  ANA_MSG_INFO( "DeltaR_L1_offl = " << m_dr_L1_offl);
  ANA_MSG_INFO( "DeltaR_L2SA_offl = " << m_dr_L2SA_offl);
  ANA_MSG_INFO( "DeltaR_L2CB_offl = " << m_dr_L2CB_offl);
  ANA_MSG_INFO( "DeltaR_EF_offl = " << m_dr_EF_offl);
  ANA_MSG_INFO( "DeltaR_EF_L1 = " << m_dr_EF_L1);
  ANA_MSG_INFO( "DeltaR_EF_L2SA = " << m_dr_EF_L2SA);
  ANA_MSG_INFO( "DeltaR_EF_L2CB = " << m_dr_EF_L2CB);

  // set up and initialize the muon selection tool in initialize()
  // MuQuality = {0,1,2,3,4,5} === {Tight, Medium, Loose, VeryLoose, HighPt and LowPt} Efficiency muon qualities
  ANA_CHECK (m_muonSelection.setProperty("MaxEta", 2.5)); 
  ANA_CHECK (m_muonSelection.setProperty("MuQuality", m_muonQuality)); //medium
  ANA_CHECK (m_muonSelection.initialize());

  // Initialize and configure trigger tools
  ANA_CHECK (m_trigConfigTool.initialize());
  ANA_CHECK (m_trigDecisionTool.setProperty ("ConfigTool", m_trigConfigTool.getHandle())); // connect the TrigDecisionTool to the ConfigTool
  ANA_CHECK (m_trigDecisionTool.setProperty ("TrigDecisionKey", "xTrigDecision"));
  ANA_CHECK (m_trigDecisionTool.initialize());


 

  return StatusCode::SUCCESS;
}



StatusCode TrigxAODAnalysis :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  ANA_MSG_INFO ("in execute");

  // retrieve the eventInfo object from the event store
  const xAOD::EventInfo *eventInfo = nullptr;
  ANA_CHECK (evtStore()->retrieve (eventInfo, "EventInfo"));

  // check if the event is data or MC
  bool isMC = false;
  // check if the event is MC
  if (eventInfo->eventType (xAOD::EventInfo::IS_SIMULATION)) {
    isMC = true;
  }

  // print out run and event number from retrieved object
  ANA_MSG_INFO ("in execute, runNumber = " << eventInfo->runNumber() << ", eventNumber = " << eventInfo->eventNumber());
  m_runNumber = eventInfo->runNumber ();
  m_eventNumber = eventInfo->eventNumber ();

  // Read/fill truth muons variables:
  m_muonTruthEta->clear();
  m_muonTruthPhi->clear();
  m_muonTruthPt->clear();
  m_muonTruthE->clear();
  m_muonTruthQ->clear();

  m_ZTruthEta->clear();
  m_ZTruthPhi->clear();
  m_ZTruthPt->clear();
  m_ZTruthE->clear();
  m_JpsiTruthEta->clear();
  m_JpsiTruthPhi->clear();
  m_JpsiTruthPt->clear();
  m_JpsiTruthE->clear();

  if (isMC) {
     const xAOD::TruthEventContainer* truthEvts = nullptr;
     ANA_CHECK (evtStore()->retrieve (truthEvts, "TruthEvents"));
     int n_truth_muons = 0;
     int n_Z_truth = 0;
     int n_Jpsi_truth = 0;
     for (auto truth : *truthEvts) {
        unsigned int nPart = truth->nTruthParticles();
        if (m_verbose > 0) ANA_MSG_INFO ("Number of truth particles: " << nPart);
        for (unsigned int iPart=0; iPart<nPart; ++iPart) {
           const xAOD::TruthParticle *particle = (truth)->truthParticle(iPart);
           if (particle) {
              if (particle->isZ()) {
                 if (particle->status() == 22) {
                    n_Z_truth++;
                    m_ZTruthEta->push_back (particle->eta ());
                    m_ZTruthPhi->push_back (particle->phi ());
                    m_ZTruthPt-> push_back (particle->pt () / 1000.);
                    m_ZTruthE->  push_back (particle->e ());
                 }
                 if (m_verbose > 1) {
                    ANA_MSG_INFO ("particle: " << iPart << " PDGId: " << particle->pdgId());
                    ANA_MSG_INFO ("particle: " << iPart << " Status: " << particle->status());
                    ANA_MSG_INFO ("particle: " << iPart << " hasDecayVtx: " << particle->hasDecayVtx());
                 }
              }
              if (fabs(particle->pdgId()) == 443) {
                 if (particle->status() == 22) {
                    n_Jpsi_truth++;
                    m_JpsiTruthEta->push_back (particle->eta ());
                    m_JpsiTruthPhi->push_back (particle->phi ());
                    m_JpsiTruthPt-> push_back (particle->pt () / 1000.);
                    m_JpsiTruthE->  push_back (particle->e ());
                 }
                 if (m_verbose > 1) {
                    ANA_MSG_INFO ("particle: " << iPart << " PDGId: " << particle->pdgId());
                    ANA_MSG_INFO ("particle: " << iPart << " Status: " << particle->status());
                    ANA_MSG_INFO ("particle: " << iPart << " hasDecayVtx: " << particle->hasDecayVtx());
                 }
              }
              if (particle->status() != 1) continue;
              if (particle->hasDecayVtx()) continue;
              if (particle->isMuon()) {
                 if (m_verbose > 1) {
                    ANA_MSG_INFO ("particle: " << iPart << " PDGId: " << particle->pdgId());
                    ANA_MSG_INFO ("particle: " << iPart << " Status: " << particle->status());
                    ANA_MSG_INFO ("particle: " << iPart << " hasDecayVtx: " << particle->hasDecayVtx());
                 }
                 n_truth_muons++;
                 m_muonTruthEta->push_back (particle->eta ());
                 m_muonTruthPhi->push_back (particle->phi ());
                 m_muonTruthPt-> push_back (particle->pt () / 1000.);
                 m_muonTruthE->  push_back (particle->e ());
                 m_muonTruthQ->  push_back (particle->charge ());
              }
           }
        }
     }
     if (m_verbose > 0) ANA_MSG_INFO ("Number of truth muons found: " << n_truth_muons);
     if (m_verbose > 0) ANA_MSG_INFO ("Number of truth Z found: " << n_Z_truth);
     if (m_verbose > 0) ANA_MSG_INFO ("Number of truth J/psi found: " << n_Jpsi_truth);
  }

  // Read/fill offline muons variables:
  const xAOD::MuonContainer* muons = nullptr;
  ANA_CHECK (evtStore()->retrieve (muons, "Muons"));
  m_muonEta->clear();
  m_muonPhi->clear();
  m_muonPt->clear();
  m_muonE->clear();
  m_muonQ->clear();
  m_muonIso->clear();
  int n_muons = 0;
  for (auto muon : *muons) {
     if ( !m_muonSelection->accept(*muon) ) continue;
     n_muons++;
     m_muonEta->push_back (muon->eta ());
     m_muonPhi->push_back (muon->phi ());
     m_muonPt-> push_back (muon->pt ()/1000.);
     m_muonE->  push_back (muon->e ());
     m_muonQ->  push_back (muon->charge ());
     m_muonIso->  push_back (muon->auxdata< float >("ptvarcone30")/muon->pt());
  }
  if (m_verbose > 0) ANA_MSG_INFO ("Number of offline muons found: " << n_muons);

  // J/psi and Z search
  m_Z_id1->clear();
  m_Z_id2->clear();
  m_Z_m->clear();
  m_Jpsi_id1->clear();
  m_Jpsi_id2->clear();
  m_Jpsi_m->clear();

  match_zjpsi(*m_muonEta, *m_muonPhi, *m_muonPt, *m_muonQ, *m_muonIso, *m_Z_id1, *m_Z_id2, *m_Z_m, *m_Jpsi_id1, *m_Jpsi_id2, *m_Jpsi_m);

  if (m_verbose > 0) ANA_MSG_INFO ("Number of Z found: " << m_Z_m->size());
  if (m_verbose > 0) ANA_MSG_INFO ("Number of J/psi found: " << m_Jpsi_m->size());

  // L1 Muons (common to all chains as ancestor acess doesn't work with standalone code)
  // LVL1MuonRoI
  m_l1roiEta->clear();
  m_l1roiPhi->clear();
  m_l1roiPtThr->clear();
  m_l1roiRoI->clear();
  m_l1roiMatchTruth->clear();
  m_l1roiMatchOffl->clear();
  const xAOD::MuonRoIContainer* l1roic = nullptr;
  ANA_CHECK (evtStore()->retrieve (l1roic, "LVL1MuonRoIs"));
  int n_l1roi = 0;
  for (auto l1 : *l1roic) {
     n_l1roi++;
     if (m_verbose > 1) ANA_MSG_INFO ("execute(): L1 RoI: RoI number: " << l1->getRoI() << ", pt/eta/phi: " << l1->thrValue()/1000. << " / " << l1->eta() << " / " << l1->phi());
     m_l1roiEta->push_back(l1->eta());
     m_l1roiPhi->push_back(l1->phi());
     m_l1roiPtThr->push_back(l1->thrValue()/1000.);
     m_l1roiRoI->push_back(l1->getRoI());
     // match with truth and offline
     int idx_truth = match_dr(*m_muonTruthEta,*m_muonTruthPhi,l1->eta(),l1->phi(), m_dr_L1_truth);
     int idx_offl = match_dr(*m_muonEta,*m_muonPhi,l1->eta(),l1->phi(), m_dr_L1_offl);
     m_l1roiMatchTruth->push_back(idx_truth);
     m_l1roiMatchOffl->push_back(idx_offl);
  }
  if (m_verbose > 0) ANA_MSG_INFO ("Number of L1 muons found: " << n_l1roi);

  // HLT trigger infos

  // check if the event passed selected trigger chains (save result and L1*HLT prescale)
  m_chainNameTrg->clear();
  m_chainPassTrg->clear();
  m_chainL1PassTrg->clear();
  m_chainL1PrescaledTrg->clear();
  m_chainPSTrg->clear();
  for (auto &chName : m_chainName) {
     auto cg = m_trigDecisionTool->getChainGroup(chName);
     int passed = cg->isPassed(); 
     float prescale = cg->getPrescale();

     const unsigned int bits = cg->isPassedBits();
     // L1?
     bool tbp = bits&TrigDefs::L1_isPassedBeforePrescale;
     bool tap = bits&TrigDefs::L1_isPassedAfterPrescale;
     bool tav = bits&TrigDefs::L1_isPassedAfterVeto;
     // L1 is prescaled if tbp and not tap
     bool l1prescale=tbp && !tap;

     m_chainNameTrg->push_back (chName);
     m_chainPassTrg->push_back (passed);
     if (tap) m_chainL1PassTrg->push_back (1);
     else     m_chainL1PassTrg->push_back (0);
     if (l1prescale) m_chainL1PrescaledTrg->push_back (1);
     else            m_chainL1PrescaledTrg->push_back (0);
     m_chainPSTrg->push_back (prescale);
     if (m_verbose > 1) ANA_MSG_INFO ("execute(): " << chName << ", chain passed(1)/failed(0) = " << passed << ", total chain prescale (L1*HLT) = " << prescale);
     if (m_verbose > 1) ANA_MSG_INFO ("execute(): " << chName << ", chain L1 passed(1)/failed(0) = " << tap << ", L1 prescaled? = " << l1prescale);
  }

  // HLT Muons (for each analysed chain)

  m_l2saEta->clear();
  m_l2saPhi->clear();
  m_l2saPt->clear();
  m_l2saQ->clear();
  m_l2saRoI->clear();
  m_l2saRoId->clear();
  m_l2saMatchTruth->clear();
  m_l2saMatchOffl->clear();
  m_l2saMatchL1->clear();

  m_l2cbEta->clear();
  m_l2cbPhi->clear();
  m_l2cbPt->clear();
  m_l2cbQ->clear();
  m_l2cbRoI->clear();
  m_l2cbRoId->clear();
  m_l2cbMatchTruth->clear();
  m_l2cbMatchOffl->clear();
  m_l2cbMatchL1->clear();
  m_l2cbMatchL2SA->clear();

  m_efmuEta->clear();
  m_efmuPhi->clear();
  m_efmuPt->clear();
  m_efmuQ->clear();
  m_efmuMatchTruth->clear();
  m_efmuMatchOffl->clear();
  m_efmuMatchL1->clear();
  m_efmuMatchL2CB->clear();
  m_efmuMatchL2SA->clear();

  for (auto &chName : m_chainName) {

     int n_l2sa = 0;
     int n_l2cb = 0;
     auto cg = m_trigDecisionTool->getChainGroup(chName);
     auto fc = cg->features(TrigDefs::alsoDeactivateTEs);
     auto l2sac = fc.containerFeature<xAOD::L2StandAloneMuonContainer>();
     auto l2cbc = fc.containerFeature<xAOD::L2CombinedMuonContainer>();
     if (m_verbose > 0) ANA_MSG_INFO("Looking at L2 features for the event for chain: " << chName);

     std::vector<float> x_l2saEta;
     std::vector<float> x_l2saPhi;
     std::vector<float> x_l2saPt;
     std::vector<float> x_l2saQ;
     std::vector<unsigned int> x_l2saRoI;
     std::vector<unsigned int> x_l2saRoId;
     std::vector<unsigned int> x_l2saMatchTruth;
     std::vector<unsigned int> x_l2saMatchOffl;
     std::vector<unsigned int> x_l2saMatchL1;

     for (auto &l2samuc : l2sac) {
        if (m_verbose > 0) ANA_MSG_INFO(" -> xAOD::L2StandaloneMuonContainer: " << l2samuc.label());
        for (auto l2samu : *l2samuc.cptr()) {
           if (m_verbose > 1) ANA_MSG_INFO("    -> muon pt = " << l2samu->pt() << " [GeV]");
           if (m_verbose > 1) ANA_MSG_INFO("    -> muon roiId = " << l2samu->roiId());
           if (m_verbose > 1) ANA_MSG_INFO("    -> muon roiN = " << l2samu->roiNumber());
           n_l2sa++;
           x_l2saEta.push_back(l2samu->eta());
           x_l2saPhi.push_back(l2samu->phi());
           x_l2saPt.push_back(fabs(l2samu->pt()));
           float q = (l2samu->pt() > 0 ? 1.0 : -1.0);
           x_l2saQ.push_back(q);
           x_l2saRoI.push_back(l2samu->roiNumber());
           x_l2saRoId.push_back(l2samu->roiId());
           int idx_truth = match_dr(*m_muonTruthEta,*m_muonTruthPhi,l2samu->eta(),l2samu->phi(), m_dr_L2SA_truth);
           int idx_offl = match_dr(*m_muonEta,*m_muonPhi,l2samu->eta(),l2samu->phi(), m_dr_L2SA_offl);
           x_l2saMatchTruth.push_back(idx_truth);
           x_l2saMatchOffl.push_back(idx_offl);
           int idx_L1 = -1;
           for (unsigned int j=0; j<m_l1roiRoI->size(); ++j) {
              if (m_l1roiRoI->at(j) == (l2samu->roiNumber())) { idx_L1 = j; }
           }
           x_l2saMatchL1.push_back(idx_L1);
        }
     }
     m_l2saEta->push_back(x_l2saEta);
     m_l2saPhi->push_back(x_l2saPhi);
     m_l2saPt->push_back(x_l2saPt);
     m_l2saQ->push_back(x_l2saQ);
     m_l2saRoI->push_back(x_l2saRoI);
     m_l2saRoId->push_back(x_l2saRoId);
     m_l2saMatchTruth->push_back(x_l2saMatchTruth);
     m_l2saMatchOffl->push_back(x_l2saMatchOffl);
     m_l2saMatchL1->push_back(x_l2saMatchL1);
     if (m_verbose > 0) ANA_MSG_INFO ("Number of L2SA muons found for chain: " << chName << " : " << n_l2sa);

     std::vector<float> x_l2cbEta;
     std::vector<float> x_l2cbPhi;
     std::vector<float> x_l2cbPt;
     std::vector<float> x_l2cbQ;
     std::vector<unsigned int> x_l2cbRoI;
     std::vector<unsigned int> x_l2cbRoId;
     std::vector<unsigned int> x_l2cbMatchTruth;
     std::vector<unsigned int> x_l2cbMatchOffl;
     std::vector<unsigned int> x_l2cbMatchL1;
     std::vector<unsigned int> x_l2cbMatchL2SA;

     for (auto &l2cbmuc : l2cbc) {
        if (m_verbose > 0) ANA_MSG_INFO(" -> xAOD::L2CombinedMuonContainer: " << l2cbmuc.label());
        for (auto l2cbmu : *l2cbmuc.cptr()) {
           if (m_verbose > 1) ANA_MSG_INFO("    -> muon pt = " << l2cbmu->pt()/1000.0 << " [GeV]");
           if (m_verbose > 1) ANA_MSG_INFO("    -> muon roiId = " << (l2cbmu->muSATrack())->roiId());
           n_l2cb++;
           x_l2cbEta.push_back(l2cbmu->eta());
           x_l2cbPhi.push_back(l2cbmu->phi());
           x_l2cbPt.push_back(l2cbmu->pt()/1000.);
           x_l2cbQ.push_back(l2cbmu->charge());
           x_l2cbRoI.push_back((l2cbmu->muSATrack())->roiNumber());
           x_l2cbRoId.push_back((l2cbmu->muSATrack())->roiId());
           int idx_truth = match_dr(*m_muonTruthEta,*m_muonTruthPhi,l2cbmu->eta(),l2cbmu->phi(), m_dr_L2CB_truth);
           int idx_offl = match_dr(*m_muonEta,*m_muonPhi,l2cbmu->eta(),l2cbmu->phi(), m_dr_L2CB_offl);
           int idx_L1 = -1;
           for (unsigned int j=0; j<m_l1roiRoI->size(); ++j) {
              if (m_l1roiRoI->at(j) == ((l2cbmu->muSATrack())->roiNumber())) { idx_L1 = j; }
           }
           x_l2cbMatchTruth.push_back(idx_truth);
           x_l2cbMatchOffl.push_back(idx_offl);
           x_l2cbMatchL1.push_back(idx_L1);
           x_l2cbMatchL2SA.push_back((l2cbmu->muSATrack())->roiId());
        }
     }
     m_l2cbEta->push_back(x_l2cbEta);
     m_l2cbPhi->push_back(x_l2cbPhi);
     m_l2cbPt->push_back(x_l2cbPt);
     m_l2cbQ->push_back(x_l2cbQ);
     m_l2cbRoI->push_back(x_l2cbRoI);
     m_l2cbMatchTruth->push_back(x_l2cbMatchTruth);
     m_l2cbMatchOffl->push_back(x_l2cbMatchOffl);
     m_l2cbMatchL1->push_back(x_l2cbMatchL1);
     m_l2cbMatchL2SA->push_back(x_l2cbMatchL2SA);
     if (m_verbose > 0) ANA_MSG_INFO ("Number of L2CB muons found for chain: " << chName << " : " << n_l2cb);
  }


  // EF
  // HLT_xAOD__MuonContainer_MuonEFInfo

  int n_chan = 0;
  for (auto &chName : m_chainName) {

     int n_efmu = 0;
     auto cg = m_trigDecisionTool->getChainGroup(chName);
     //auto fc = cg->features();
     auto fc = cg->features(TrigDefs::alsoDeactivateTEs);
     auto efc = fc.containerFeature<xAOD::MuonContainer>();

     if (m_verbose > 0) ANA_MSG_INFO("Looking at EF features for the event for chain: " << chName);

     std::vector<float> x_efmuEta;
     std::vector<float> x_efmuPhi;
     std::vector<float> x_efmuPt;
     std::vector<float> x_efmuQ;
     std::vector<unsigned int> x_efmuMatchTruth;
     std::vector<unsigned int> x_efmuMatchOffl;
     std::vector<unsigned int> x_efmuMatchL1;
     std::vector<unsigned int> x_efmuMatchL2SA;
     std::vector<unsigned int> x_efmuMatchL2CB;

     for (auto &efmuc : efc) {
        if (m_verbose > 0) ANA_MSG_INFO(" -> xAOD::MuonContainer: " << efmuc.label());
        for (auto efmu : *efmuc.cptr()) {
           if (m_verbose > 1) ANA_MSG_INFO("    -> muon pt = " << efmu->pt()/1000.0 << " [GeV]");
           n_efmu++;
           x_efmuEta.push_back(efmu->eta());
           x_efmuPhi.push_back(efmu->phi());
           x_efmuPt.push_back(efmu->pt()/1000.);
           x_efmuQ.push_back(efmu->charge());
           int idx_truth = match_dr(*m_muonTruthEta,*m_muonTruthPhi,efmu->eta(),efmu->phi(), m_dr_EF_truth);
           int idx_offl = match_dr(*m_muonEta,*m_muonPhi,efmu->eta(),efmu->phi(), m_dr_EF_offl);
           int idx_L1 = match_dr(*m_l1roiEta,*m_l1roiPhi,efmu->eta(),efmu->phi(), m_dr_EF_L1);
           int idx_L2SA = match_dr((m_l2saEta->at(n_chan)), (m_l2saPhi->at(n_chan)),efmu->eta(),efmu->phi(), m_dr_EF_L2SA);
           int idx_L2CB = match_dr((m_l2cbEta->at(n_chan)), (m_l2cbPhi->at(n_chan)),efmu->eta(),efmu->phi(), m_dr_EF_L2CB);
           x_efmuMatchTruth.push_back(idx_truth);
           x_efmuMatchOffl.push_back(idx_offl);
           x_efmuMatchL1.push_back(idx_L1);
           x_efmuMatchL2SA.push_back(idx_L2SA);
           x_efmuMatchL2CB.push_back(idx_L2CB);
        }
     }
     m_efmuEta->push_back(x_efmuEta);
     m_efmuPhi->push_back(x_efmuPhi);
     m_efmuPt->push_back(x_efmuPt);
     m_efmuQ->push_back(x_efmuQ);
     m_efmuMatchTruth->push_back(x_efmuMatchTruth);
     m_efmuMatchOffl->push_back(x_efmuMatchOffl);
     m_efmuMatchL1->push_back(x_efmuMatchL1);
     m_efmuMatchL2SA->push_back(x_efmuMatchL2SA);
     m_efmuMatchL2CB->push_back(x_efmuMatchL2CB);
     if (m_verbose > 0) ANA_MSG_INFO ("Number of EF muons found for chain: " << chName << " : " << n_efmu);
     n_chan++;
  }


  // Fill the event into the tree:
  tree ("muontrig")->Fill ();
  
  return StatusCode::SUCCESS;
}



StatusCode TrigxAODAnalysis :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.
  return StatusCode::SUCCESS;
}
