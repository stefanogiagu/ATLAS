# The name of the package:
atlas_subdir (TrigAnalysis)

# Add the shared library:
atlas_add_library (TrigAnalysisLib
  TrigAnalysis/*.h Root/*.cxx
  PUBLIC_HEADERS TrigAnalysis
  LINK_LIBRARIES AnaAlgorithmLib xAODEventInfo AsgAnalysisInterfaces xAODTruth xAODMuon xAODTrigger xAODTrigEgamma xAODTrigMuon MuonAnalysisInterfacesLib TrigDecisionToolLib TrigConfInterfaces)

# Add the dictionary:
atlas_add_dictionary (TrigAnalysisDict
  TrigAnalysis/TrigAnalysisDict.h
  TrigAnalysis/selection.xml
  LINK_LIBRARIES TrigAnalysisLib)

if (NOT XAOD_STANDALONE)
  # Add a component library for AthAnalysis only:
  atlas_add_component (TrigAnalysis
    src/components/*.cxx
    LINK_LIBRARIES TrigAnalysisLib)
endif ()

# Install files from the package:
atlas_install_joboptions( share/*.py )
