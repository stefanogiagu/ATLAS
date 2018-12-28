#ifndef TrigAnalysis_TrigxAODAnalysis_H
#define TrigAnalysis_TrigxAODAnalysis_H

#include <AnaAlgorithm/AnaAlgorithm.h>

#include <TTree.h>
#include <vector>

#include <MuonAnalysisInterfaces/IMuonSelectionTool.h>
#include <AsgTools/AnaToolHandle.h>

#include <TrigConfInterfaces/ITrigConfigTool.h>
#include <TrigDecisionTool/TrigDecisionTool.h>

class TrigxAODAnalysis : public EL::AnaAlgorithm
{
public:
  // this is a standard algorithm constructor
  TrigxAODAnalysis (const std::string& name, ISvcLocator* pSvcLocator);

  ~TrigxAODAnalysis();

  // these are the functions inherited from Algorithm
  virtual StatusCode initialize () override;
  virtual StatusCode execute () override;
  virtual StatusCode finalize () override;

  /// utility functions
  int match_dr(std::vector<float>, std::vector<float>, float, float, float);
  void match_zjpsi(std::vector<float>,  std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<unsigned int> &, std::vector<unsigned int> &, std::vector<float> &, std::vector<unsigned int> &, std::vector<unsigned int> &, std::vector<float> &);

  /// MuonSelectionTool
  asg::AnaToolHandle<CP::IMuonSelectionTool> m_muonSelection; //! 

  /// trigger tools member variables
  asg::AnaToolHandle<Trig::TrigDecisionTool> m_trigDecisionTool; //!
  asg::AnaToolHandle<TrigConf::ITrigConfigTool> m_trigConfigTool; //!

private:
  // Configuration, and any other types of variables go here.

  /// trigger chains to analyse
  std::vector<std::string> m_chainName;

  /// offline muon quality
  int m_muonQuality;

  /// printout
  int m_verbose;

  /// DeltaR match pars
  float m_dr_L1_truth;
  float m_dr_L2SA_truth;
  float m_dr_L2CB_truth;
  float m_dr_EF_truth;
  float m_dr_L1_offl;
  float m_dr_L2SA_offl;
  float m_dr_L2CB_offl;
  float m_dr_EF_offl;
  float m_dr_EF_L1;
  float m_dr_EF_L2SA;
  float m_dr_EF_L2CB;


  /// Global quantities 
  unsigned int m_runNumber = 0; ///< Run number for the current event
  unsigned long long m_eventNumber = 0; ///< Event number

  /// Truth Muon 4-momentum variables
  std::vector<float> *m_muonTruthEta = nullptr;
  std::vector<float> *m_muonTruthPhi = nullptr;
  std::vector<float> *m_muonTruthPt = nullptr;
  std::vector<float> *m_muonTruthE = nullptr;
  std::vector<float> *m_muonTruthQ = nullptr;

  /// Z and J/psi truth
  std::vector<float> *m_ZTruthEta = nullptr;
  std::vector<float> *m_ZTruthPhi = nullptr;
  std::vector<float> *m_ZTruthPt = nullptr;
  std::vector<float> *m_ZTruthE = nullptr;
  std::vector<float> *m_JpsiTruthEta = nullptr;
  std::vector<float> *m_JpsiTruthPhi = nullptr;
  std::vector<float> *m_JpsiTruthPt = nullptr;
  std::vector<float> *m_JpsiTruthE = nullptr;

  /// Offline Muon 4-momentum variables
  std::vector<float> *m_muonEta = nullptr;
  std::vector<float> *m_muonPhi = nullptr;
  std::vector<float> *m_muonPt = nullptr;
  std::vector<float> *m_muonE = nullptr;
  std::vector<float> *m_muonQ = nullptr;
  std::vector<float> *m_muonIso = nullptr;

  /// Z and J/psi tagged
  std::vector<unsigned int> *m_Z_id1 = nullptr;
  std::vector<unsigned int> *m_Z_id2 = nullptr;
  std::vector<float> *m_Z_m = nullptr;
  std::vector<unsigned int> *m_Jpsi_id1 = nullptr;
  std::vector<unsigned int> *m_Jpsi_id2 = nullptr;
  std::vector<float> *m_Jpsi_m = nullptr;

  /// L1RoI Muon variables
  std::vector<float>        *m_l1roiEta   = nullptr;
  std::vector<float>        *m_l1roiPhi   = nullptr;
  std::vector<float>        *m_l1roiPtThr = nullptr;
  std::vector<unsigned int> *m_l1roiRoI   = nullptr;
  std::vector<unsigned int> *m_l1roiMatchTruth = nullptr;
  std::vector<unsigned int> *m_l1roiMatchOffl  = nullptr;

  /// Global trigger info
  std::vector<std::string> *m_chainNameTrg = nullptr;
  std::vector<int>         *m_chainPassTrg = nullptr;
  std::vector<int>         *m_chainL1PassTrg = nullptr;
  std::vector<float>       *m_chainPSTrg = nullptr;
  std::vector<int>         *m_chainL1PrescaledTrg = nullptr;

  /// L2SA Muon variables
  std::vector<std::vector<float>>        *m_l2saEta = nullptr;
  std::vector<std::vector<float>>        *m_l2saPhi = nullptr;
  std::vector<std::vector<float>>        *m_l2saPt = nullptr;
  std::vector<std::vector<float>>        *m_l2saQ = nullptr;
  std::vector<std::vector<unsigned int>> *m_l2saRoI = nullptr;
  std::vector<std::vector<unsigned int>> *m_l2saRoId = nullptr;
  std::vector<std::vector<unsigned int>> *m_l2saMatchTruth = nullptr;
  std::vector<std::vector<unsigned int>> *m_l2saMatchOffl  = nullptr;
  std::vector<std::vector<unsigned int>> *m_l2saMatchL1  = nullptr;

  /// L2CB Muon variables
  std::vector<std::vector<float>>        *m_l2cbEta = nullptr;
  std::vector<std::vector<float>>        *m_l2cbPhi = nullptr;
  std::vector<std::vector<float>>        *m_l2cbPt = nullptr;
  std::vector<std::vector<float>>        *m_l2cbQ = nullptr;
  std::vector<std::vector<unsigned int>> *m_l2cbRoI = nullptr;
  std::vector<std::vector<unsigned int>> *m_l2cbRoId = nullptr;
  std::vector<std::vector<unsigned int>> *m_l2cbMatchTruth = nullptr;
  std::vector<std::vector<unsigned int>> *m_l2cbMatchOffl  = nullptr;
  std::vector<std::vector<unsigned int>> *m_l2cbMatchL1  = nullptr;
  std::vector<std::vector<unsigned int>> *m_l2cbMatchL2SA= nullptr;

  /// EF Muon variables
  std::vector<std::vector<float>>       *m_efmuEta = nullptr;
  std::vector<std::vector<float>>       *m_efmuPhi = nullptr;
  std::vector<std::vector<float>>       *m_efmuPt = nullptr;
  std::vector<std::vector<float>>       *m_efmuQ = nullptr;
  std::vector<std::vector<unsigned int>> *m_efmuMatchTruth = nullptr;
  std::vector<std::vector<unsigned int>> *m_efmuMatchOffl  = nullptr;
  std::vector<std::vector<unsigned int>> *m_efmuMatchL1  = nullptr;
  std::vector<std::vector<unsigned int>> *m_efmuMatchL2SA= nullptr;
  std::vector<std::vector<unsigned int>> *m_efmuMatchL2CB= nullptr;

};

#endif
