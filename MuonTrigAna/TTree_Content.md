******************************************************************************
*Tree    :muontrig  : Muon Trigger analysis ntuple                           *
******************************************************************************
* 
* Run Number and Event Number
*
*Br    0 :RunNumber : RunNumber/i                                            *
*Br    1 :EventNumber : EventNumber/l                                        *
*............................................................................*
* 
* For each truth muon the 4mom + charge info
*
*Br    2 :MuonTruthEta : vector<float>                                       *
*Br    3 :MuonTruthPhi : vector<float>                                       *
*Br    4 :MuonTruthPt : vector<float>                                        *
*Br    5 :MuonTruthE : vector<float>                                         *
*Br    6 :MuonTruthQ : vector<float>                                         *
*............................................................................*
* 
* For each truth Z the 4mom info
*
*Br    7 :ZTruthEta : vector<float>                                          *
*Br    8 :ZTruthPhi : vector<float>                                          *
*Br    9 :ZTruthPt  : vector<float>                                          *
*Br   10 :ZTruthE   : vector<float>                                          *
*............................................................................*
* 
* For each truth J/psi the 4mom info
*
*Br   11 :JpsiTruthEta : vector<float>                                       *
*Br   12 :JpsiTruthPhi : vector<float>                                       *
*Br   13 :JpsiTruthPt : vector<float>                                        *
*Br   14 :JpsiTruthE : vector<float>                                         *
*............................................................................*
* 
* For each selected offline muon the 4mom info + charge + iso
*
*Br   15 :MuonEta   : vector<float>                                          *
*Br   16 :MuonPhi   : vector<float>                                          *
*Br   17 :MuonPt    : vector<float>                                          *
*Br   18 :MuonE     : vector<float>                                          *
*Br   19 :MuonQ     : vector<float>                                          *
*Br   20 :MuonIso   : vector<float>   (def as ptvarcone3/ptmu)               *
*............................................................................*
* 
* For each reconstructed (from offline muons) Z: 
*
*Br   21 :Z_id1     : vector<unsigned int> ((isolated) muon entry index in the Muon block)
*Br   22 :Z_id2     : vector<unsigned int> (second muon entry in index the Muon block)
*Br   23 :Z_m       : vector<float>        (reconstructed Z inv. mass)       *
*............................................................................*
* 
* For each reconstructed (from offline muons) J/psi:
*
*Br   24 :Jpsi_id1  : vector<unsigned int>                                   *
*Br   25 :Jpsi_id2  : vector<unsigned int>                                   *
*Br   26 :Jpsi_m    : vector<float>                                          *
*............................................................................*
* 
* For each Muon L1 RoI
*
*Br   27 :L1RoIEta  : vector<float>                                          *
*Br   28 :L1RoIPhi  : vector<float>                                          *
*Br   29 :L1RoIPtThr : vector<float>  (pt threshold value)                   *
*Br   30 :L1RoI     : vector<unsigned int>   (RoI number)                    *
*Br   31 :L1RoIMatchTruth : vector<unsigned int> (matched entry index with a truth muon in the MuonTruth block (-1 not matched))
*Br   32 :L1RoIMatchOffl : vector<unsigned int>  (matched entry index with an offline muon in the Muon block (-1 not matched))
*............................................................................*
* 
* For each analysed HLT muon chain (selectable via property)
*
*Br   33 :ChainNameTrg : vector<string> (name of the chain)                  *
*Br   34 :ChainPassTrg : vector<int>    (1 passed 0 not-passed)              *
*Br   XX :ChainL1PassTrg : vector<int>  L1 passed (1 passed 0 not-passed)    *
*Br   35 :ChainPSTrg : vector<float>    (total L1*HLT prescale)              *
*Br   XX :ChainL1PrescaleTrg : vector<int>   (1 if prescaled at L1)          *
*............................................................................*
* 
* L2 StandAlone muons:
* For each analysed HLT muon chain (same order as in the previous block) a std:vector with the info for each feature)
*
*Br   36 :L2SAEta   : vector<vector<float> >                                 *
*Br   37 :L2SAPhi   : vector<vector<float> >                                 *
*Br   38 :L2SAPt    : vector<vector<float> >                                 *
*Br   39 :L2SAQ     : vector<vector<float> >                                 *
*Br   40 :L2SARoI   : vector<vector<unsigned int> >                          *
*Br   41 :L2SARoId  : vector<vector<unsigned int> >                          *
*Br   42 :L2SAMatchTruth : vector<vector<unsigned int> > (matched entry index with TruthMuon)
*Br   43 :L2SAMatchOffl : vector<vector<unsigned int> >  (matched entry index with offline Muon)
*Br   44 :L2SAMatchL1 : vector<vector<unsigned int> >    (matched entry index with L1 RoI)
*............................................................................*
*
* L2 Combined muons:
* For each analysed HLT muon chain (same order as in the previous block) a std:vector with the info for each feature)
*
*Br   45 :L2CBEta   : vector<vector<float> >                                 *
*Br   46 :L2CBPhi   : vector<vector<float> >                                 *
*Br   47 :L2CBPt    : vector<vector<float> >                                 *
*Br   48 :L2CBQ     : vector<vector<float> >                                 *
*Br   49 :L2CBRoI   : vector<vector<unsigned int> >                          *
*Br   50 :L2CBRoId  : vector<vector<unsigned int> >                          *
*Br   51 :L2CBMatchTruth : vector<vector<unsigned int> >                     *
*Br   52 :L2CBMatchOffl : vector<vector<unsigned int> >                      *
*Br   53 :L2CBMatchL1 : vector<vector<unsigned int> >                        *
*Br   54 :L2CBMatchL2SA : vector<vector<unsigned int> > (matched entry index with L2 SA)
*............................................................................*
*
* EF Muons:
* For each analysed HLT muon chain (same order as in the previous block) a std:vector with the info for each feature)
*
*Br   55 :EFmuEta   : vector<vector<float> >                                 *
*Br   56 :EFmuPhi   : vector<vector<float> >                                 *
*Br   57 :EFmuPt    : vector<vector<float> >                                 *
*Br   58 :EFmuQ     : vector<vector<float> >                                 *
*Br   59 :EFmuMatchTruth : vector<vector<unsigned int> >                     *
*Br   60 :EFmuMatchOffl : vector<vector<unsigned int> >                      *
*Br   61 :EFmuMatchL1 : vector<vector<unsigned int> >                        *
*Br   62 :EFmuMatchL2SA : vector<vector<unsigned int> >                      *
*Br   63 :EFmuMatchL2CB : vector<vector<unsigned int> > (matched entry index with L2 CB)
*............................................................................*
