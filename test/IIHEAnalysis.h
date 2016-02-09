//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan  4 15:45:36 2016 by ROOT version 5.34/32
// from TTree IIHEAnalysis/IIHEAnalysis
// found on file: outfile.root
//////////////////////////////////////////////////////////

#ifndef IIHEAnalysis_h
#define IIHEAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class IIHEAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          ev_event;
   UInt_t          ev_run;
   UInt_t          ev_luminosityBlock;
   UInt_t          ev_time;
   UInt_t          ev_time_unixTime;
   UInt_t          ev_time_microsecondOffset;
   Float_t         ev_fixedGridRhoAll;
   Float_t         ev_rho_kt6PFJetsForIsolation;
   UInt_t          pv_n;
   vector<float>   *pv_x;
   vector<float>   *pv_y;
   vector<float>   *pv_z;
   vector<int>     *pv_isValid;
   vector<float>   *pv_normalizedChi2;
   vector<int>     *pv_ndof;
   vector<int>     *pv_nTracks;
   vector<int>     *pv_totTrackSize;
   UInt_t          sc_n;
   vector<float>   *sc_energy;
   vector<float>   *sc_eta;
   vector<float>   *sc_etacorr;
   vector<float>   *sc_theta;
   vector<float>   *sc_thetacorr;
   vector<float>   *sc_et;
   vector<float>   *sc_phi;
   vector<float>   *sc_px;
   vector<float>   *sc_py;
   vector<float>   *sc_pz;
   vector<float>   *sc_x;
   vector<float>   *sc_y;
   vector<float>   *sc_z;
   UInt_t          ph_n;
   vector<float>   *ph_px;
   vector<float>   *ph_pt;
   vector<float>   *ph_eta;
   vector<float>   *ph_theta;
   vector<float>   *ph_phi;
   vector<float>   *ph_energy;
   vector<float>   *ph_mass;
   vector<int>     *ph_isPFlowPhoton;
   vector<int>     *ph_isStandardPhoton;
   vector<int>     *ph_hasConversionTracks;
   vector<int>     *ph_hasPixelSeed;
   vector<int>     *ph_isEB;
   vector<int>     *ph_isEE;
   vector<int>     *ph_isEBGap;
   vector<int>     *ph_isEBEtaGap;
   vector<int>     *ph_isEBPhiGap;
   vector<int>     *ph_isEEGap;
   vector<int>     *ph_isEERingGap;
   vector<int>     *ph_isEEDeeGap;
   vector<int>     *ph_isEBEEGap;
   vector<float>   *ph_hadronicOverEm;
   vector<float>   *ph_hadronicDepth1OverEm;
   vector<float>   *ph_hadronicDepth2OverEm;
   vector<float>   *ph_hadTowOverEm;
   vector<float>   *ph_hadTowDepth1OverEm;
   vector<float>   *ph_hadTowDepth2OverEm;
   vector<float>   *ph_e1x5;
   vector<float>   *ph_e2x5;
   vector<float>   *ph_e3x3;
   vector<float>   *ph_e5x5;
   vector<float>   *ph_maxEnergyXtal;
   vector<float>   *ph_sigmaEtaEta;
   vector<float>   *ph_sigmaIetaIeta;
   vector<float>   *ph_r1x5;
   vector<float>   *ph_r2x5;
   vector<float>   *ph_r9;
   vector<float>   *ph_mipChi2;
   vector<float>   *ph_mipTotEnergy;
   vector<float>   *ph_mipSlope;
   vector<float>   *ph_mipIntercept;
   vector<int>     *ph_mipNhitCone;
   vector<int>     *ph_mipIsHalo;
   vector<float>   *ph_ecalRecHitSumEtConeDR04;
   vector<float>   *ph_hcalTowerSumEtConeDR04;
   vector<float>   *ph_hcalDepth1TowerSumEtConeDR04;
   vector<float>   *ph_hcalDepth2TowerSumEtConeDR04;
   vector<float>   *ph_hcalTowerSumEtBcConeDR04;
   vector<float>   *ph_hcalDepth1TowerSumEtBcConeDR04;
   vector<float>   *ph_hcalDepth2TowerSumEtBcConeDR04;
   vector<float>   *ph_trkSumPtSolidConeDR04;
   vector<float>   *ph_trkSumPtHollowConeDR04;
   vector<int>     *ph_nTrkSolidConeDR04;
   vector<int>     *ph_nTrkHollowConeDR04;
   vector<float>   *ph_ecalRecHitSumEtConeDR03;
   vector<float>   *ph_hcalTowerSumEtConeDR03;
   vector<float>   *ph_hcalDepth1TowerSumEtConeDR03;
   vector<float>   *ph_hcalDepth2TowerSumEtConeDR03;
   vector<float>   *ph_hcalTowerSumEtBcConeDR03;
   vector<float>   *ph_hcalDepth1TowerSumEtBcConeDR03;
   vector<float>   *ph_hcalDepth2TowerSumEtBcConeDR03;
   vector<float>   *ph_trkSumPtSolidConeDR03;
   vector<float>   *ph_trkSumPtHollowConeDR03;
   vector<int>     *ph_nTrkSolidConeDR03;
   vector<int>     *ph_nTrkHollowConeDR03;
   vector<float>   *ph_chargedHadronIso;
   vector<float>   *ph_neutralHadronIso;
   vector<float>   *ph_photonIso;
   vector<float>   *ph_chargedHadronIsoWrongVtx;
   vector<float>   *ph_sumChargedParticlePt;
   vector<float>   *ph_sumNeutralHadronEtHighThreshold;
   vector<float>   *ph_sumPhotonEtHighThreshold;
   vector<float>   *ph_sumPUPt;
   vector<int>     *ph_nClusterOutsideMustache;
   vector<float>   *ph_etOutsideMustache;
   vector<float>   *ph_pfMVA;
   vector<float>   *ph_mc_bestDR;
   vector<int>     *ph_mc_index;
   vector<float>   *ph_mc_ERatio;
   UInt_t          gsf_n;
   vector<int>     *gsf_classification;
   vector<float>   *gsf_energy;
   vector<float>   *gsf_p;
   vector<float>   *gsf_pt;
   vector<float>   *gsf_scE1x5;
   vector<float>   *gsf_scE5x5;
   vector<float>   *gsf_scE2x5Max;
   vector<float>   *gsf_eta;
   vector<float>   *gsf_phi;
   vector<float>   *gsf_theta;
   vector<float>   *gsf_px;
   vector<float>   *gsf_py;
   vector<float>   *gsf_pz;
   vector<float>   *gsf_superClusterEta;
   vector<float>   *gsf_superClusterEnergy;
   vector<float>   *gsf_caloEnergy;
   vector<float>   *gsf_deltaEtaSuperClusterTrackAtVtx;
   vector<float>   *gsf_deltaPhiSuperClusterTrackAtVtx;
   vector<float>   *gsf_hadronicOverEm;
   vector<float>   *gsf_hcalDepth1OverEcal;
   vector<float>   *gsf_hcalDepth2OverEcal;
   vector<float>   *gsf_dr03TkSumPt;
   vector<float>   *gsf_dr03EcalRecHitSumEt;
   vector<float>   *gsf_dr03HcalDepth1TowerSumEt;
   vector<float>   *gsf_dr03HcalDepth2TowerSumEt;
   vector<int>     *gsf_charge;
   vector<float>   *gsf_sigmaIetaIeta;
   vector<int>     *gsf_ecaldrivenSeed;
   vector<int>     *gsf_trackerdrivenSeed;
   vector<int>     *gsf_isEB;
   vector<int>     *gsf_isEE;
   vector<float>   *gsf_deltaEtaSeedClusterTrackAtCalo;
   vector<float>   *gsf_deltaPhiSeedClusterTrackAtCalo;
   vector<float>   *gsf_ecalEnergy;
   vector<float>   *gsf_eSuperClusterOverP;
   vector<float>   *gsf_dxy;
   vector<float>   *gsf_dxy_beamSpot;
   vector<float>   *gsf_dxy_firstPVtx;
   vector<float>   *gsf_dxyError;
   vector<float>   *gsf_dz;
   vector<float>   *gsf_dz_beamSpot;
   vector<float>   *gsf_dz_firstPVtx;
   vector<float>   *gsf_dzError;
   vector<float>   *gsf_vz;
   vector<int>     *gsf_numberOfValidHits;
   vector<int>     *gsf_nLostInnerHits;
   vector<int>     *gsf_nLostOuterHits;
   vector<int>     *gsf_convFlags;
   vector<float>   *gsf_convDist;
   vector<float>   *gsf_convDcot;
   vector<float>   *gsf_convRadius;
   vector<float>   *gsf_fBrem;
   vector<float>   *gsf_e1x5;
   vector<float>   *gsf_e2x5Max;
   vector<float>   *gsf_e5x5;
   vector<float>   *gsf_r9;
   vector<vector<int> > *gsf_hitsinfo;
   vector<float>   *gsf_pixelMatch_dPhi1;
   vector<float>   *gsf_pixelMatch_dPhi2;
   vector<float>   *gsf_pixelMatch_dRz1;
   vector<float>   *gsf_pixelMatch_dRz2;
   vector<int>     *gsf_pixelMatch_subDetector1;
   vector<int>     *gsf_pixelMatch_subDetector2;
   vector<float>   *gsf_mc_bestDR;
   vector<int>     *gsf_mc_index;
   vector<float>   *gsf_mc_ERatio;
   UInt_t          mu_n;
   UInt_t          mu_n_saved;
   UInt_t          mu_gt_n;
   UInt_t          mu_ot_n;
   UInt_t          mu_it_n;
   vector<float>   *mu_gt_qoverp;
   vector<int>     *mu_gt_charge;
   vector<float>   *mu_gt_pt;
   vector<float>   *mu_gt_eta;
   vector<float>   *mu_gt_phi;
   vector<float>   *mu_gt_p;
   vector<float>   *mu_gt_px;
   vector<float>   *mu_gt_py;
   vector<float>   *mu_gt_pz;
   vector<float>   *mu_gt_theta;
   vector<float>   *mu_gt_lambda;
   vector<float>   *mu_gt_d0;
   vector<float>   *mu_gt_dz;
   vector<float>   *mu_gt_dz_beamspot;
   vector<float>   *mu_gt_dz_firstPVtx;
   vector<float>   *mu_gt_dxy;
   vector<float>   *mu_gt_dxy_beamspot;
   vector<float>   *mu_gt_dxy_firstPVtx;
   vector<float>   *mu_gt_dsz;
   vector<float>   *mu_gt_vx;
   vector<float>   *mu_gt_vy;
   vector<float>   *mu_gt_vz;
   vector<float>   *mu_gt_qoverpError;
   vector<float>   *mu_gt_ptError;
   vector<float>   *mu_gt_thetaError;
   vector<float>   *mu_gt_lambdaError;
   vector<float>   *mu_gt_phiError;
   vector<float>   *mu_gt_dxyError;
   vector<float>   *mu_gt_d0Error;
   vector<float>   *mu_gt_dszError;
   vector<float>   *mu_gt_dzError;
   vector<float>   *mu_gt_etaError;
   vector<float>   *mu_gt_chi2;
   vector<float>   *mu_gt_ndof;
   vector<float>   *mu_gt_normalizedChi2;
   vector<float>   *mu_ot_qoverp;
   vector<int>     *mu_ot_charge;
   vector<float>   *mu_ot_pt;
   vector<float>   *mu_ot_eta;
   vector<float>   *mu_ot_phi;
   vector<float>   *mu_ot_p;
   vector<float>   *mu_ot_px;
   vector<float>   *mu_ot_py;
   vector<float>   *mu_ot_pz;
   vector<float>   *mu_ot_theta;
   vector<float>   *mu_ot_lambda;
   vector<float>   *mu_ot_d0;
   vector<float>   *mu_ot_dz;
   vector<float>   *mu_ot_dz_beamspot;
   vector<float>   *mu_ot_dz_firstPVtx;
   vector<float>   *mu_ot_dxy;
   vector<float>   *mu_ot_dxy_beamspot;
   vector<float>   *mu_ot_dxy_firstPVtx;
   vector<float>   *mu_ot_dsz;
   vector<float>   *mu_ot_vx;
   vector<float>   *mu_ot_vy;
   vector<float>   *mu_ot_vz;
   vector<float>   *mu_ot_qoverpError;
   vector<float>   *mu_ot_ptError;
   vector<float>   *mu_ot_thetaError;
   vector<float>   *mu_ot_lambdaError;
   vector<float>   *mu_ot_phiError;
   vector<float>   *mu_ot_dxyError;
   vector<float>   *mu_ot_d0Error;
   vector<float>   *mu_ot_dszError;
   vector<float>   *mu_ot_dzError;
   vector<float>   *mu_ot_etaError;
   vector<float>   *mu_ot_chi2;
   vector<float>   *mu_ot_ndof;
   vector<float>   *mu_ot_normalizedChi2;
   vector<float>   *mu_it_qoverp;
   vector<int>     *mu_it_charge;
   vector<float>   *mu_it_pt;
   vector<float>   *mu_it_eta;
   vector<float>   *mu_it_phi;
   vector<float>   *mu_it_p;
   vector<float>   *mu_it_px;
   vector<float>   *mu_it_py;
   vector<float>   *mu_it_pz;
   vector<float>   *mu_it_theta;
   vector<float>   *mu_it_lambda;
   vector<float>   *mu_it_d0;
   vector<float>   *mu_it_dz;
   vector<float>   *mu_it_dz_beamspot;
   vector<float>   *mu_it_dz_firstPVtx;
   vector<float>   *mu_it_dxy;
   vector<float>   *mu_it_dxy_beamspot;
   vector<float>   *mu_it_dxy_firstPVtx;
   vector<float>   *mu_it_dsz;
   vector<float>   *mu_it_vx;
   vector<float>   *mu_it_vy;
   vector<float>   *mu_it_vz;
   vector<float>   *mu_it_qoverpError;
   vector<float>   *mu_it_ptError;
   vector<float>   *mu_it_thetaError;
   vector<float>   *mu_it_lambdaError;
   vector<float>   *mu_it_phiError;
   vector<float>   *mu_it_dxyError;
   vector<float>   *mu_it_d0Error;
   vector<float>   *mu_it_dszError;
   vector<float>   *mu_it_dzError;
   vector<float>   *mu_it_etaError;
   vector<float>   *mu_it_chi2;
   vector<float>   *mu_it_ndof;
   vector<float>   *mu_it_normalizedChi2;
   vector<int>     *mu_isGlobalMuon;
   vector<int>     *mu_isStandAloneMuon;
   vector<int>     *mu_isTrackerMuon;
   vector<int>     *mu_isPFMuon;
   vector<int>     *mu_isPFIsolationValid;
   vector<int>     *mu_numberOfMatchedStations;
   vector<int>     *mu_numberOfValidPixelHits;
   vector<int>     *mu_numberOfValidTrackerHits;
   vector<int>     *mu_numberOfValidMuonHits;
   vector<int>     *mu_tevOptimized_charge;
   vector<float>   *mu_tevOptimized_pt;
   vector<float>   *mu_tevOptimized_eta;
   vector<float>   *mu_tevOptimized_phi;
   vector<float>   *mu_tevOptimized_theta;
   vector<float>   *mu_tevOptimized_px;
   vector<float>   *mu_tevOptimized_py;
   vector<float>   *mu_tevOptimized_pz;
   vector<float>   *mu_tevOptimized_d0;
   vector<float>   *mu_tevOptimized_dz;
   vector<float>   *mu_tevOptimized_dz_beamSpot;
   vector<float>   *mu_tevOptimized_dz_firstPVtx;
   vector<float>   *mu_tevOptimized_dxy;
   vector<float>   *mu_tevOptimized_dxy_beamSpot;
   vector<float>   *mu_tevOptimized_dxy_firstPVtx;
   vector<float>   *mu_tevOptimized_ptError;
   vector<float>   *mu_tevOptimized_etaError;
   vector<float>   *mu_tevOptimized_phiError;
   vector<float>   *mu_tevOptimized_thetaError;
   vector<float>   *mu_tevOptimized_d0Error;
   vector<float>   *mu_tevOptimized_dzError;
   vector<float>   *mu_tevOptimized_dxyError;
   vector<float>   *mu_isolationR03_sumPt;
   vector<float>   *mu_isolationR03_trackerVetoPt;
   vector<float>   *mu_isolationR03_emEt;
   vector<float>   *mu_isolationR03_emVetoEt;
   vector<float>   *mu_isolationR03_hadEt;
   vector<float>   *mu_isolationR03_hadVetoEt;
   vector<float>   *mu_isolationR05_sumPt;
   vector<float>   *mu_isolationR05_trackerVetoPt;
   vector<float>   *mu_isolationR05_emEt;
   vector<float>   *mu_isolationR05_emVetoEt;
   vector<float>   *mu_isolationR05_hadEt;
   vector<float>   *mu_isolationR05_hadVetoEt;
   vector<float>   *mu_pfIsolationR03_sumChargedHadronPt;
   vector<float>   *mu_pfIsolationR03_sumChargedParticlePt;
   vector<float>   *mu_pfIsolationR03_sumPhotonEt;
   vector<float>   *mu_pfIsolationR03_sumNeutralHadronEtHighThreshold;
   vector<float>   *mu_pfIsolationR03_sumPhotonEtHighThreshold;
   vector<float>   *mu_pfIsolationR03_sumPUPt;
   vector<float>   *mu_pfIsolationR04_sumChargedHadronPt;
   vector<float>   *mu_pfIsolationR04_sumChargedParticlePt;
   vector<float>   *mu_pfIsolationR04_sumPhotonEt;
   vector<float>   *mu_pfIsolationR04_sumNeutralHadronEtHighThreshold;
   vector<float>   *mu_pfIsolationR04_sumPhotonEtHighThreshold;
   vector<float>   *mu_pfIsolationR04_sumPUPt;
   vector<float>   *mu_pfMeanDRIsoProfileR03_sumChargedHadronPt;
   vector<float>   *mu_pfMeanDRIsoProfileR03_sumChargedParticlePt;
   vector<float>   *mu_pfMeanDRIsoProfileR03_sumPhotonEt;
   vector<float>   *mu_pfMeanDRIsoProfileR03_sumNeutralHadronEtHighThreshold;
   vector<float>   *mu_pfMeanDRIsoProfileR03_sumPhotonEtHighThreshold;
   vector<float>   *mu_pfMeanDRIsoProfileR03_sumPUPt;
   vector<float>   *mu_pfMeanDRIsoProfileR04_sumChargedHadronPt;
   vector<float>   *mu_pfMeanDRIsoProfileR04_sumChargedParticlePt;
   vector<float>   *mu_pfMeanDRIsoProfileR04_sumPhotonEt;
   vector<float>   *mu_pfMeanDRIsoProfileR04_sumNeutralHadronEtHighThreshold;
   vector<float>   *mu_pfMeanDRIsoProfileR04_sumPhotonEtHighThreshold;
   vector<float>   *mu_pfMeanDRIsoProfileR04_sumPUPt;
   vector<float>   *mu_mc_bestDR;
   vector<int>     *mu_mc_index;
   vector<float>   *mu_mc_ERatio;
   Float_t         MET_caloMet_et;
   Float_t         MET_caloMet_phi;
   Float_t         MET_pfMet_et;
   Float_t         MET_pfMet_phi;
   vector<float>   *HEEP_eseffsixix;
   vector<float>   *HEEP_eseffsiyiy;
   vector<float>   *HEEP_eseffsirir;
   vector<float>   *HEEP_preshowerEnergy;
   vector<float>   *HEEP_e1x3;
   vector<float>   *HEEP_eMax;
   vector<float>   *HEEP_e5x5;
   vector<float>   *HEEP_e2x5Right;
   vector<float>   *HEEP_e2x5Left;
   vector<float>   *HEEP_e2x5Top;
   vector<float>   *HEEP_e2x5Bottom;
   vector<float>   *HEEP_eRight;
   vector<float>   *HEEP_eLeft;
   vector<float>   *HEEP_eTop;
   vector<float>   *HEEP_eBottom;
   vector<float>   *HEEP_basicClusterSeedTime;
   vector<int>     *EBHits_rawId;
   vector<int>     *EBHits_iRechit;
   vector<float>   *EBHits_energy;
   vector<int>     *EBHits_ieta;
   vector<int>     *EBHits_iphi;
   vector<int>     *EBHits_RecoFlag;
   vector<int>     *EBHits_kSaturated;
   vector<int>     *EBHits_kLeadingEdgeRecovered;
   vector<int>     *EBHits_kNeighboursRecovered;
   vector<int>     *EBHits_kWeird;
   vector<int>     *EEHits_rawId;
   vector<int>     *EEHits_iRechit;
   vector<float>   *EEHits_energy;
   vector<int>     *EEHits_ieta;
   vector<int>     *EEHits_iphi;
   vector<int>     *EEHits_RecoFlag;
   vector<int>     *EEHits_kSaturated;
   vector<int>     *EEHits_kLeadingEdgeRecovered;
   vector<int>     *EEHits_kNeighboursRecovered;
   vector<int>     *EEHits_kWeird;
   vector<vector<float> > *HEEP_crystal_energy;
   vector<vector<float> > *HEEP_crystal_eta;
   vector<vector<float> > *HEEP_eshitsixix;
   vector<vector<float> > *HEEP_eshitsiyiy;
   vector<vector<int> > *HEEP_crystal_ietaorix;
   vector<vector<int> > *HEEP_crystal_iphioriy;
   vector<int>     *HEEP_cutflow41_Et;
   Int_t           HEEP_cutflow41_Et_n;
   Int_t           HEEP_cutflow41_Et_nCumulative;
   vector<float>   *HEEP_cutflow41_Et_value;
   vector<int>     *HEEP_cutflow41_eta;
   Int_t           HEEP_cutflow41_eta_n;
   Int_t           HEEP_cutflow41_eta_nCumulative;
   vector<float>   *HEEP_cutflow41_eta_value;
   vector<int>     *HEEP_cutflow41_acceptance;
   Int_t           HEEP_cutflow41_acceptance_n;
   vector<int>     *HEEP_cutflow41_EcalDriven;
   Int_t           HEEP_cutflow41_EcalDriven_n;
   Int_t           HEEP_cutflow41_EcalDriven_nCumulative;
   vector<float>   *HEEP_cutflow41_EcalDriven_value;
   vector<int>     *HEEP_cutflow41_dEtaIn;
   Int_t           HEEP_cutflow41_dEtaIn_n;
   Int_t           HEEP_cutflow41_dEtaIn_nCumulative;
   vector<float>   *HEEP_cutflow41_dEtaIn_value;
   vector<int>     *HEEP_cutflow41_dPhiIn;
   Int_t           HEEP_cutflow41_dPhiIn_n;
   Int_t           HEEP_cutflow41_dPhiIn_nCumulative;
   vector<float>   *HEEP_cutflow41_dPhiIn_value;
   vector<int>     *HEEP_cutflow41_HOverE;
   Int_t           HEEP_cutflow41_HOverE_n;
   Int_t           HEEP_cutflow41_HOverE_nCumulative;
   vector<float>   *HEEP_cutflow41_HOverE_value;
   vector<int>     *HEEP_cutflow41_SigmaIetaIeta;
   Int_t           HEEP_cutflow41_SigmaIetaIeta_n;
   Int_t           HEEP_cutflow41_SigmaIetaIeta_nCumulative;
   vector<float>   *HEEP_cutflow41_SigmaIetaIeta_value;
   vector<int>     *HEEP_cutflow41_E1x5OverE5x5;
   Int_t           HEEP_cutflow41_E1x5OverE5x5_n;
   Int_t           HEEP_cutflow41_E1x5OverE5x5_nCumulative;
   vector<float>   *HEEP_cutflow41_E1x5OverE5x5_value;
   vector<int>     *HEEP_cutflow41_E2x5OverE5x5;
   Int_t           HEEP_cutflow41_E2x5OverE5x5_n;
   Int_t           HEEP_cutflow41_E2x5OverE5x5_nCumulative;
   vector<float>   *HEEP_cutflow41_E2x5OverE5x5_value;
   vector<int>     *HEEP_cutflow41_missingHits;
   Int_t           HEEP_cutflow41_missingHits_n;
   Int_t           HEEP_cutflow41_missingHits_nCumulative;
   vector<float>   *HEEP_cutflow41_missingHits_value;
   vector<int>     *HEEP_cutflow41_dxyFirstPV;
   Int_t           HEEP_cutflow41_dxyFirstPV_n;
   Int_t           HEEP_cutflow41_dxyFirstPV_nCumulative;
   vector<float>   *HEEP_cutflow41_dxyFirstPV_value;
   vector<int>     *HEEP_cutflow41_ID;
   Int_t           HEEP_cutflow41_ID_n;
   vector<int>     *HEEP_cutflow41_isolEMHadDepth1;
   Int_t           HEEP_cutflow41_isolEMHadDepth1_n;
   Int_t           HEEP_cutflow41_isolEMHadDepth1_nCumulative;
   vector<float>   *HEEP_cutflow41_isolEMHadDepth1_value;
   vector<int>     *HEEP_cutflow41_IsolPtTrks;
   Int_t           HEEP_cutflow41_IsolPtTrks_n;
   Int_t           HEEP_cutflow41_IsolPtTrks_nCumulative;
   vector<float>   *HEEP_cutflow41_IsolPtTrks_value;
   vector<int>     *HEEP_cutflow41_isolation;
   Int_t           HEEP_cutflow41_isolation_n;
   vector<int>     *HEEP_cutflow41_total;
   Int_t           HEEP_cutflow41_total_n;
   vector<int>     *HEEP_cutflow50_25ns_Et;
   Int_t           HEEP_cutflow50_25ns_Et_n;
   Int_t           HEEP_cutflow50_25ns_Et_nCumulative;
   vector<float>   *HEEP_cutflow50_25ns_Et_value;
   vector<int>     *HEEP_cutflow50_25ns_eta;
   Int_t           HEEP_cutflow50_25ns_eta_n;
   Int_t           HEEP_cutflow50_25ns_eta_nCumulative;
   vector<float>   *HEEP_cutflow50_25ns_eta_value;
   vector<int>     *HEEP_cutflow50_25ns_acceptance;
   Int_t           HEEP_cutflow50_25ns_acceptance_n;
   vector<int>     *HEEP_cutflow50_25ns_EcalDriven;
   Int_t           HEEP_cutflow50_25ns_EcalDriven_n;
   Int_t           HEEP_cutflow50_25ns_EcalDriven_nCumulative;
   vector<float>   *HEEP_cutflow50_25ns_EcalDriven_value;
   vector<int>     *HEEP_cutflow50_25ns_dEtaIn;
   Int_t           HEEP_cutflow50_25ns_dEtaIn_n;
   Int_t           HEEP_cutflow50_25ns_dEtaIn_nCumulative;
   vector<float>   *HEEP_cutflow50_25ns_dEtaIn_value;
   vector<int>     *HEEP_cutflow50_25ns_dPhiIn;
   Int_t           HEEP_cutflow50_25ns_dPhiIn_n;
   Int_t           HEEP_cutflow50_25ns_dPhiIn_nCumulative;
   vector<float>   *HEEP_cutflow50_25ns_dPhiIn_value;
   vector<int>     *HEEP_cutflow50_25ns_HOverE;
   Int_t           HEEP_cutflow50_25ns_HOverE_n;
   Int_t           HEEP_cutflow50_25ns_HOverE_nCumulative;
   vector<float>   *HEEP_cutflow50_25ns_HOverE_value;
   vector<int>     *HEEP_cutflow50_25ns_SigmaIetaIeta;
   Int_t           HEEP_cutflow50_25ns_SigmaIetaIeta_n;
   Int_t           HEEP_cutflow50_25ns_SigmaIetaIeta_nCumulative;
   vector<float>   *HEEP_cutflow50_25ns_SigmaIetaIeta_value;
   vector<int>     *HEEP_cutflow50_25ns_E1x5OverE5x5;
   Int_t           HEEP_cutflow50_25ns_E1x5OverE5x5_n;
   Int_t           HEEP_cutflow50_25ns_E1x5OverE5x5_nCumulative;
   vector<float>   *HEEP_cutflow50_25ns_E1x5OverE5x5_value;
   vector<int>     *HEEP_cutflow50_25ns_E2x5OverE5x5;
   Int_t           HEEP_cutflow50_25ns_E2x5OverE5x5_n;
   Int_t           HEEP_cutflow50_25ns_E2x5OverE5x5_nCumulative;
   vector<float>   *HEEP_cutflow50_25ns_E2x5OverE5x5_value;
   vector<int>     *HEEP_cutflow50_25ns_missingHits;
   Int_t           HEEP_cutflow50_25ns_missingHits_n;
   Int_t           HEEP_cutflow50_25ns_missingHits_nCumulative;
   vector<float>   *HEEP_cutflow50_25ns_missingHits_value;
   vector<int>     *HEEP_cutflow50_25ns_dxyFirstPV;
   Int_t           HEEP_cutflow50_25ns_dxyFirstPV_n;
   Int_t           HEEP_cutflow50_25ns_dxyFirstPV_nCumulative;
   vector<float>   *HEEP_cutflow50_25ns_dxyFirstPV_value;
   vector<int>     *HEEP_cutflow50_25ns_ID;
   Int_t           HEEP_cutflow50_25ns_ID_n;
   vector<int>     *HEEP_cutflow50_25ns_isolEMHadDepth1;
   Int_t           HEEP_cutflow50_25ns_isolEMHadDepth1_n;
   Int_t           HEEP_cutflow50_25ns_isolEMHadDepth1_nCumulative;
   vector<float>   *HEEP_cutflow50_25ns_isolEMHadDepth1_value;
   vector<int>     *HEEP_cutflow50_25ns_IsolPtTrks;
   Int_t           HEEP_cutflow50_25ns_IsolPtTrks_n;
   Int_t           HEEP_cutflow50_25ns_IsolPtTrks_nCumulative;
   vector<float>   *HEEP_cutflow50_25ns_IsolPtTrks_value;
   vector<int>     *HEEP_cutflow50_25ns_isolation;
   Int_t           HEEP_cutflow50_25ns_isolation_n;
   vector<int>     *HEEP_cutflow50_25ns_total;
   Int_t           HEEP_cutflow50_25ns_total_n;
   vector<int>     *HEEP_cutflow50_50ns_Et;
   Int_t           HEEP_cutflow50_50ns_Et_n;
   Int_t           HEEP_cutflow50_50ns_Et_nCumulative;
   vector<float>   *HEEP_cutflow50_50ns_Et_value;
   vector<int>     *HEEP_cutflow50_50ns_eta;
   Int_t           HEEP_cutflow50_50ns_eta_n;
   Int_t           HEEP_cutflow50_50ns_eta_nCumulative;
   vector<float>   *HEEP_cutflow50_50ns_eta_value;
   vector<int>     *HEEP_cutflow50_50ns_acceptance;
   Int_t           HEEP_cutflow50_50ns_acceptance_n;
   vector<int>     *HEEP_cutflow50_50ns_EcalDriven;
   Int_t           HEEP_cutflow50_50ns_EcalDriven_n;
   Int_t           HEEP_cutflow50_50ns_EcalDriven_nCumulative;
   vector<float>   *HEEP_cutflow50_50ns_EcalDriven_value;
   vector<int>     *HEEP_cutflow50_50ns_dEtaIn;
   Int_t           HEEP_cutflow50_50ns_dEtaIn_n;
   Int_t           HEEP_cutflow50_50ns_dEtaIn_nCumulative;
   vector<float>   *HEEP_cutflow50_50ns_dEtaIn_value;
   vector<int>     *HEEP_cutflow50_50ns_dPhiIn;
   Int_t           HEEP_cutflow50_50ns_dPhiIn_n;
   Int_t           HEEP_cutflow50_50ns_dPhiIn_nCumulative;
   vector<float>   *HEEP_cutflow50_50ns_dPhiIn_value;
   vector<int>     *HEEP_cutflow50_50ns_HOverE;
   Int_t           HEEP_cutflow50_50ns_HOverE_n;
   Int_t           HEEP_cutflow50_50ns_HOverE_nCumulative;
   vector<float>   *HEEP_cutflow50_50ns_HOverE_value;
   vector<int>     *HEEP_cutflow50_50ns_SigmaIetaIeta;
   Int_t           HEEP_cutflow50_50ns_SigmaIetaIeta_n;
   Int_t           HEEP_cutflow50_50ns_SigmaIetaIeta_nCumulative;
   vector<float>   *HEEP_cutflow50_50ns_SigmaIetaIeta_value;
   vector<int>     *HEEP_cutflow50_50ns_E1x5OverE5x5;
   Int_t           HEEP_cutflow50_50ns_E1x5OverE5x5_n;
   Int_t           HEEP_cutflow50_50ns_E1x5OverE5x5_nCumulative;
   vector<float>   *HEEP_cutflow50_50ns_E1x5OverE5x5_value;
   vector<int>     *HEEP_cutflow50_50ns_E2x5OverE5x5;
   Int_t           HEEP_cutflow50_50ns_E2x5OverE5x5_n;
   Int_t           HEEP_cutflow50_50ns_E2x5OverE5x5_nCumulative;
   vector<float>   *HEEP_cutflow50_50ns_E2x5OverE5x5_value;
   vector<int>     *HEEP_cutflow50_50ns_missingHits;
   Int_t           HEEP_cutflow50_50ns_missingHits_n;
   Int_t           HEEP_cutflow50_50ns_missingHits_nCumulative;
   vector<float>   *HEEP_cutflow50_50ns_missingHits_value;
   vector<int>     *HEEP_cutflow50_50ns_dxyFirstPV;
   Int_t           HEEP_cutflow50_50ns_dxyFirstPV_n;
   Int_t           HEEP_cutflow50_50ns_dxyFirstPV_nCumulative;
   vector<float>   *HEEP_cutflow50_50ns_dxyFirstPV_value;
   vector<int>     *HEEP_cutflow50_50ns_ID;
   Int_t           HEEP_cutflow50_50ns_ID_n;
   vector<int>     *HEEP_cutflow50_50ns_isolEMHadDepth1;
   Int_t           HEEP_cutflow50_50ns_isolEMHadDepth1_n;
   Int_t           HEEP_cutflow50_50ns_isolEMHadDepth1_nCumulative;
   vector<float>   *HEEP_cutflow50_50ns_isolEMHadDepth1_value;
   vector<int>     *HEEP_cutflow50_50ns_IsolPtTrks;
   Int_t           HEEP_cutflow50_50ns_IsolPtTrks_n;
   Int_t           HEEP_cutflow50_50ns_IsolPtTrks_nCumulative;
   vector<float>   *HEEP_cutflow50_50ns_IsolPtTrks_value;
   vector<int>     *HEEP_cutflow50_50ns_isolation;
   Int_t           HEEP_cutflow50_50ns_isolation_n;
   vector<int>     *HEEP_cutflow50_50ns_total;
   Int_t           HEEP_cutflow50_50ns_total_n;
   vector<int>     *HEEP_cutflow50_Et;
   Int_t           HEEP_cutflow50_Et_n;
   Int_t           HEEP_cutflow50_Et_nCumulative;
   vector<float>   *HEEP_cutflow50_Et_value;
   vector<int>     *HEEP_cutflow50_eta;
   Int_t           HEEP_cutflow50_eta_n;
   Int_t           HEEP_cutflow50_eta_nCumulative;
   vector<float>   *HEEP_cutflow50_eta_value;
   vector<int>     *HEEP_cutflow50_acceptance;
   Int_t           HEEP_cutflow50_acceptance_n;
   vector<int>     *HEEP_cutflow50_EcalDriven;
   Int_t           HEEP_cutflow50_EcalDriven_n;
   Int_t           HEEP_cutflow50_EcalDriven_nCumulative;
   vector<float>   *HEEP_cutflow50_EcalDriven_value;
   vector<int>     *HEEP_cutflow50_dEtaIn;
   Int_t           HEEP_cutflow50_dEtaIn_n;
   Int_t           HEEP_cutflow50_dEtaIn_nCumulative;
   vector<float>   *HEEP_cutflow50_dEtaIn_value;
   vector<int>     *HEEP_cutflow50_dPhiIn;
   Int_t           HEEP_cutflow50_dPhiIn_n;
   Int_t           HEEP_cutflow50_dPhiIn_nCumulative;
   vector<float>   *HEEP_cutflow50_dPhiIn_value;
   vector<int>     *HEEP_cutflow50_HOverE;
   Int_t           HEEP_cutflow50_HOverE_n;
   Int_t           HEEP_cutflow50_HOverE_nCumulative;
   vector<float>   *HEEP_cutflow50_HOverE_value;
   vector<int>     *HEEP_cutflow50_SigmaIetaIeta;
   Int_t           HEEP_cutflow50_SigmaIetaIeta_n;
   Int_t           HEEP_cutflow50_SigmaIetaIeta_nCumulative;
   vector<float>   *HEEP_cutflow50_SigmaIetaIeta_value;
   vector<int>     *HEEP_cutflow50_E1x5OverE5x5;
   Int_t           HEEP_cutflow50_E1x5OverE5x5_n;
   Int_t           HEEP_cutflow50_E1x5OverE5x5_nCumulative;
   vector<float>   *HEEP_cutflow50_E1x5OverE5x5_value;
   vector<int>     *HEEP_cutflow50_E2x5OverE5x5;
   Int_t           HEEP_cutflow50_E2x5OverE5x5_n;
   Int_t           HEEP_cutflow50_E2x5OverE5x5_nCumulative;
   vector<float>   *HEEP_cutflow50_E2x5OverE5x5_value;
   vector<int>     *HEEP_cutflow50_missingHits;
   Int_t           HEEP_cutflow50_missingHits_n;
   Int_t           HEEP_cutflow50_missingHits_nCumulative;
   vector<float>   *HEEP_cutflow50_missingHits_value;
   vector<int>     *HEEP_cutflow50_dxyFirstPV;
   Int_t           HEEP_cutflow50_dxyFirstPV_n;
   Int_t           HEEP_cutflow50_dxyFirstPV_nCumulative;
   vector<float>   *HEEP_cutflow50_dxyFirstPV_value;
   vector<int>     *HEEP_cutflow50_ID;
   Int_t           HEEP_cutflow50_ID_n;
   vector<int>     *HEEP_cutflow50_isolEMHadDepth1;
   Int_t           HEEP_cutflow50_isolEMHadDepth1_n;
   Int_t           HEEP_cutflow50_isolEMHadDepth1_nCumulative;
   vector<float>   *HEEP_cutflow50_isolEMHadDepth1_value;
   vector<int>     *HEEP_cutflow50_IsolPtTrks;
   Int_t           HEEP_cutflow50_IsolPtTrks_n;
   Int_t           HEEP_cutflow50_IsolPtTrks_nCumulative;
   vector<float>   *HEEP_cutflow50_IsolPtTrks_value;
   vector<int>     *HEEP_cutflow50_isolation;
   Int_t           HEEP_cutflow50_isolation_n;
   vector<int>     *HEEP_cutflow50_total;
   Int_t           HEEP_cutflow50_total_n;
   vector<int>     *HEEP_cutflow51_Et;
   Int_t           HEEP_cutflow51_Et_n;
   Int_t           HEEP_cutflow51_Et_nCumulative;
   vector<float>   *HEEP_cutflow51_Et_value;
   vector<int>     *HEEP_cutflow51_eta;
   Int_t           HEEP_cutflow51_eta_n;
   Int_t           HEEP_cutflow51_eta_nCumulative;
   vector<float>   *HEEP_cutflow51_eta_value;
   vector<int>     *HEEP_cutflow51_acceptance;
   Int_t           HEEP_cutflow51_acceptance_n;
   vector<int>     *HEEP_cutflow51_EcalDriven;
   Int_t           HEEP_cutflow51_EcalDriven_n;
   Int_t           HEEP_cutflow51_EcalDriven_nCumulative;
   vector<float>   *HEEP_cutflow51_EcalDriven_value;
   vector<int>     *HEEP_cutflow51_dEtaIn;
   Int_t           HEEP_cutflow51_dEtaIn_n;
   Int_t           HEEP_cutflow51_dEtaIn_nCumulative;
   vector<float>   *HEEP_cutflow51_dEtaIn_value;
   vector<int>     *HEEP_cutflow51_dPhiIn;
   Int_t           HEEP_cutflow51_dPhiIn_n;
   Int_t           HEEP_cutflow51_dPhiIn_nCumulative;
   vector<float>   *HEEP_cutflow51_dPhiIn_value;
   vector<int>     *HEEP_cutflow51_HOverE;
   Int_t           HEEP_cutflow51_HOverE_n;
   Int_t           HEEP_cutflow51_HOverE_nCumulative;
   vector<float>   *HEEP_cutflow51_HOverE_value;
   vector<int>     *HEEP_cutflow51_SigmaIetaIeta;
   Int_t           HEEP_cutflow51_SigmaIetaIeta_n;
   Int_t           HEEP_cutflow51_SigmaIetaIeta_nCumulative;
   vector<float>   *HEEP_cutflow51_SigmaIetaIeta_value;
   vector<int>     *HEEP_cutflow51_E1x5OverE5x5;
   Int_t           HEEP_cutflow51_E1x5OverE5x5_n;
   Int_t           HEEP_cutflow51_E1x5OverE5x5_nCumulative;
   vector<float>   *HEEP_cutflow51_E1x5OverE5x5_value;
   vector<int>     *HEEP_cutflow51_E2x5OverE5x5;
   Int_t           HEEP_cutflow51_E2x5OverE5x5_n;
   Int_t           HEEP_cutflow51_E2x5OverE5x5_nCumulative;
   vector<float>   *HEEP_cutflow51_E2x5OverE5x5_value;
   vector<int>     *HEEP_cutflow51_missingHits;
   Int_t           HEEP_cutflow51_missingHits_n;
   Int_t           HEEP_cutflow51_missingHits_nCumulative;
   vector<float>   *HEEP_cutflow51_missingHits_value;
   vector<int>     *HEEP_cutflow51_dxyFirstPV;
   Int_t           HEEP_cutflow51_dxyFirstPV_n;
   Int_t           HEEP_cutflow51_dxyFirstPV_nCumulative;
   vector<float>   *HEEP_cutflow51_dxyFirstPV_value;
   vector<int>     *HEEP_cutflow51_ID;
   Int_t           HEEP_cutflow51_ID_n;
   vector<int>     *HEEP_cutflow51_isolEMHadDepth1;
   Int_t           HEEP_cutflow51_isolEMHadDepth1_n;
   Int_t           HEEP_cutflow51_isolEMHadDepth1_nCumulative;
   vector<float>   *HEEP_cutflow51_isolEMHadDepth1_value;
   vector<int>     *HEEP_cutflow51_IsolPtTrks;
   Int_t           HEEP_cutflow51_IsolPtTrks_n;
   Int_t           HEEP_cutflow51_IsolPtTrks_nCumulative;
   vector<float>   *HEEP_cutflow51_IsolPtTrks_value;
   vector<int>     *HEEP_cutflow51_isolation;
   Int_t           HEEP_cutflow51_isolation_n;
   vector<int>     *HEEP_cutflow51_total;
   Int_t           HEEP_cutflow51_total_n;
   vector<int>     *HEEP_cutflow60_Et;
   Int_t           HEEP_cutflow60_Et_n;
   Int_t           HEEP_cutflow60_Et_nCumulative;
   vector<float>   *HEEP_cutflow60_Et_value;
   vector<int>     *HEEP_cutflow60_eta;
   Int_t           HEEP_cutflow60_eta_n;
   Int_t           HEEP_cutflow60_eta_nCumulative;
   vector<float>   *HEEP_cutflow60_eta_value;
   vector<int>     *HEEP_cutflow60_acceptance;
   Int_t           HEEP_cutflow60_acceptance_n;
   vector<int>     *HEEP_cutflow60_EcalDriven;
   Int_t           HEEP_cutflow60_EcalDriven_n;
   Int_t           HEEP_cutflow60_EcalDriven_nCumulative;
   vector<float>   *HEEP_cutflow60_EcalDriven_value;
   vector<int>     *HEEP_cutflow60_dEtaIn;
   Int_t           HEEP_cutflow60_dEtaIn_n;
   Int_t           HEEP_cutflow60_dEtaIn_nCumulative;
   vector<float>   *HEEP_cutflow60_dEtaIn_value;
   vector<int>     *HEEP_cutflow60_dPhiIn;
   Int_t           HEEP_cutflow60_dPhiIn_n;
   Int_t           HEEP_cutflow60_dPhiIn_nCumulative;
   vector<float>   *HEEP_cutflow60_dPhiIn_value;
   vector<int>     *HEEP_cutflow60_HOverE;
   Int_t           HEEP_cutflow60_HOverE_n;
   Int_t           HEEP_cutflow60_HOverE_nCumulative;
   vector<float>   *HEEP_cutflow60_HOverE_value;
   vector<int>     *HEEP_cutflow60_SigmaIetaIeta;
   Int_t           HEEP_cutflow60_SigmaIetaIeta_n;
   Int_t           HEEP_cutflow60_SigmaIetaIeta_nCumulative;
   vector<float>   *HEEP_cutflow60_SigmaIetaIeta_value;
   vector<int>     *HEEP_cutflow60_E1x5OverE5x5;
   Int_t           HEEP_cutflow60_E1x5OverE5x5_n;
   Int_t           HEEP_cutflow60_E1x5OverE5x5_nCumulative;
   vector<float>   *HEEP_cutflow60_E1x5OverE5x5_value;
   vector<int>     *HEEP_cutflow60_E2x5OverE5x5;
   Int_t           HEEP_cutflow60_E2x5OverE5x5_n;
   Int_t           HEEP_cutflow60_E2x5OverE5x5_nCumulative;
   vector<float>   *HEEP_cutflow60_E2x5OverE5x5_value;
   vector<int>     *HEEP_cutflow60_missingHits;
   Int_t           HEEP_cutflow60_missingHits_n;
   Int_t           HEEP_cutflow60_missingHits_nCumulative;
   vector<float>   *HEEP_cutflow60_missingHits_value;
   vector<int>     *HEEP_cutflow60_dxyFirstPV;
   Int_t           HEEP_cutflow60_dxyFirstPV_n;
   Int_t           HEEP_cutflow60_dxyFirstPV_nCumulative;
   vector<float>   *HEEP_cutflow60_dxyFirstPV_value;
   vector<int>     *HEEP_cutflow60_ID;
   Int_t           HEEP_cutflow60_ID_n;
   vector<int>     *HEEP_cutflow60_isolEMHadDepth1;
   Int_t           HEEP_cutflow60_isolEMHadDepth1_n;
   Int_t           HEEP_cutflow60_isolEMHadDepth1_nCumulative;
   vector<float>   *HEEP_cutflow60_isolEMHadDepth1_value;
   vector<int>     *HEEP_cutflow60_IsolPtTrks;
   Int_t           HEEP_cutflow60_IsolPtTrks_n;
   Int_t           HEEP_cutflow60_IsolPtTrks_nCumulative;
   vector<float>   *HEEP_cutflow60_IsolPtTrks_value;
   vector<int>     *HEEP_cutflow60_isolation;
   Int_t           HEEP_cutflow60_isolation_n;
   vector<int>     *HEEP_cutflow60_total;
   Int_t           HEEP_cutflow60_total_n;
   Int_t           Zee_n;
   vector<float>   *Zee_mass;
   vector<float>   *Zee_mass_HEEP;
   vector<int>     *Zee_i1;
   vector<int>     *Zee_i2;
   Int_t           Zee_highestMass;
   Int_t           Zmm_n;
   vector<float>   *Zmm_mass;
   vector<int>     *Zmm_i1;
   vector<int>     *Zmm_i2;
   Int_t           Zmm_highestMass;
   Int_t           Zem_n;
   vector<float>   *Zem_mass;
   vector<float>   *Zem_mass_HEEP;
   vector<int>     *Zem_i1;
   vector<int>     *Zem_i2;
   Int_t           Zem_highestMass;
   UInt_t          trk_n;
   vector<float>   *trk_pt;
   vector<float>   *trk_eta;
   vector<float>   *trk_phi;
   Int_t           trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_accept;
   Int_t           trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_prescale;
   vector<float>   *trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltL1sL1DoubleEG2210_eta;
   vector<float>   *trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltL1sL1DoubleEG2210_phi;
   vector<float>   *trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg1TrackIsoFilter_eta;
   vector<float>   *trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg1TrackIsoFilter_phi;
   vector<float>   *trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg2TrackIsoFilter_eta;
   vector<float>   *trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg2TrackIsoFilter_phi;
   Int_t           trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_accept;
   Int_t           trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_prescale;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLPixelMatchFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLPixelMatchFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLNewPixelMatchFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLNewPixelMatchFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEG33EtUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEG33EtUnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLNewPixelMatchUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLNewPixelMatchUnseededFilter_phi;
   Int_t           trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_accept;
   Int_t           trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_prescale;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLPixelMatchFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLPixelMatchFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEG33EtUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEG33EtUnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_eta;
   vector<float>   *trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_phi;
   Int_t           trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_accept;
   Int_t           trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_prescale;
   vector<float>   *trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_eta;
   vector<float>   *trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_phi;
   vector<float>   *trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltSingleEle22WPLooseGsfTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltSingleEle22WPLooseGsfTrackIsoFilter_phi;
   Int_t           trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_accept;
   Int_t           trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_prescale;
   vector<float>   *trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_eta;
   vector<float>   *trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_phi;
   vector<float>   *trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltSingleEle22WPTightGsfTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltSingleEle22WPTightGsfTrackIsoFilter_phi;
   Int_t           trig_HLT_Ele30WP60_Ele8_Mass55_v2_accept;
   Int_t           trig_HLT_Ele30WP60_Ele8_Mass55_v2_prescale;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltL1sL1SingleEG25OrSingleIsoEG30er_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltL1sL1SingleEG25OrSingleIsoEG30er_phi;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtFilter_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtFilter_phi;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8ClusterShapeFilter_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8ClusterShapeFilter_phi;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HEFilter_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HEFilter_phi;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EcalIsoFilter_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EcalIsoFilter_phi;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HcalIsoFilter_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HcalIsoFilter_phi;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchFilter_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchFilter_phi;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8OneOEMinusOneOPFilter_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8OneOEMinusOneOPFilter_phi;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DetaFilter_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DetaFilter_phi;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DphiFilter_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DphiFilter_phi;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8TrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8TrackIsoFilter_phi;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtUnseededFilter_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtUnseededFilter_phi;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchUnseededFilter_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchUnseededFilter_phi;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8Mass55Filter_eta;
   vector<float>   *trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8Mass55Filter_phi;
   Int_t           trig_HLT_Ele23_WPLoose_Gsf_v2_accept;
   Int_t           trig_HLT_Ele23_WPLoose_Gsf_v2_prescale;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_v2_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_eta;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_v2_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_phi;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_v2_hltEle23WPLooseGsfTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_v2_hltEle23WPLooseGsfTrackIsoFilter_phi;
   Int_t           trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_accept;
   Int_t           trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_prescale;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_eta;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_phi;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltEle23WPLooseGsfTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltEle23WPLooseGsfTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltJetFilterEle23WPLoose_eta;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltJetFilterEle23WPLoose_phi;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltHCand80NoEle23WPLoose_eta;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltHCand80NoEle23WPLoose_phi;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMhtFilter0_eta;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMhtFilter0_phi;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMetFilter0_eta;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMetFilter0_phi;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand80NoEle23WPLooseMET_eta;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand80NoEle23WPLooseMET_phi;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand70NoEle23WPLooseMHTIDTight_eta;
   vector<float>   *trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand70NoEle23WPLooseMHTIDTight_phi;
   Int_t           trig_HLT_Ele27_WPLoose_Gsf_v1_accept;
   Int_t           trig_HLT_Ele27_WPLoose_Gsf_v1_prescale;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_v1_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_eta;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_v1_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_phi;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_v1_hltEle27noerWPLooseGsfTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_v1_hltEle27noerWPLooseGsfTrackIsoFilter_phi;
   Int_t           trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_accept;
   Int_t           trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_prescale;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_eta;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_phi;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltEle27noerWPLooseGsfTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltEle27noerWPLooseGsfTrackIsoFilter_phi;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltJetFilterEle27WPLoose_eta;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltJetFilterEle27WPLoose_phi;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltHCand80NoEle27WPLoose_eta;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltHCand80NoEle27WPLoose_phi;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMhtFilter0_eta;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMhtFilter0_phi;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMetFilter0_eta;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMetFilter0_phi;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand80NoEle27WPLooseMET_eta;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand80NoEle27WPLooseMET_phi;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand70NoEle27WPLooseMHTIDTight_eta;
   vector<float>   *trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand70NoEle27WPLooseMHTIDTight_phi;
   Int_t           trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_accept;
   Int_t           trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_prescale;
   vector<float>   *trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_eta;
   vector<float>   *trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_phi;
   vector<float>   *trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltEle27WPLooseGsfTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltEle27WPLooseGsfTrackIsoFilter_phi;
   Int_t           trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_accept;
   Int_t           trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_prescale;
   vector<float>   *trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_eta;
   vector<float>   *trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_phi;
   vector<float>   *trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltEle27WPTightGsfTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltEle27WPTightGsfTrackIsoFilter_phi;
   Int_t           trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_accept;
   Int_t           trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_prescale;
   vector<float>   *trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG30erOrSingleEG40_eta;
   vector<float>   *trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG30erOrSingleEG40_phi;
   vector<float>   *trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltEle32WPTightGsfTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltEle32WPTightGsfTrackIsoFilter_phi;
   Int_t           trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_accept;
   Int_t           trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_prescale;
   vector<float>   *trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_eta;
   vector<float>   *trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_phi;
   vector<float>   *trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter_eta;
   vector<float>   *trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter_phi;
   Int_t           trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_accept;
   Int_t           trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_prescale;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG2210_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG2210_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_phi;
   Int_t           trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_accept;
   Int_t           trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_prescale;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG1510_eta;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG1510_phi;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter_eta;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter_phi;
   Int_t           trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept;
   Int_t           trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG2210_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG2210_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi;
   Int_t           trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept;
   Int_t           trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG1510_eta;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG1510_phi;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta;
   vector<float>   *trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi;
   Int_t           trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_accept;
   Int_t           trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_prescale;
   vector<float>   *trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG20_eta;
   vector<float>   *trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG20_phi;
   vector<float>   *trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter_phi;
   Int_t           trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_accept;
   Int_t           trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_prescale;
   vector<float>   *trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltL1sL1SingleEG15_eta;
   vector<float>   *trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltL1sL1SingleEG15_phi;
   vector<float>   *trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter_phi;
   Int_t           trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept;
   Int_t           trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale;
   vector<float>   *trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG15ORSingleEG10_eta;
   vector<float>   *trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG15ORSingleEG10_phi;
   vector<float>   *trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter_eta;
   vector<float>   *trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter_phi;
   Int_t           trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_accept;
   Int_t           trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_prescale;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_eta;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_phi;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_eta;
   vector<float>   *trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_phi;
   Int_t           trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_accept;
   Int_t           trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_prescale;
   vector<float>   *trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltL1sL1SingleIsoEG22erOrSingleIsoEG20OrSingleEG25OrSingleEG20_eta;
   vector<float>   *trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltL1sL1SingleIsoEG22erOrSingleIsoEG20OrSingleEG25OrSingleEG20_phi;
   vector<float>   *trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltEle23WPLooseGsfTrackIsoFilterForEleStream_eta;
   vector<float>   *trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltEle23WPLooseGsfTrackIsoFilterForEleStream_phi;

   // List of branches
   TBranch        *b_ev_event;   //!
   TBranch        *b_ev_run;   //!
   TBranch        *b_ev_luminosityBlock;   //!
   TBranch        *b_ev_time;   //!
   TBranch        *b_ev_time_unixTime;   //!
   TBranch        *b_ev_time_microsecondOffset;   //!
   TBranch        *b_ev_fixedGridRhoAll;   //!
   TBranch        *b_ev_rho_kt6PFJetsForIsolation;   //!
   TBranch        *b_pv_n;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_isValid;   //!
   TBranch        *b_pv_normalizedChi2;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_pv_nTracks;   //!
   TBranch        *b_pv_totTrackSize;   //!
   TBranch        *b_sc_n;   //!
   TBranch        *b_sc_energy;   //!
   TBranch        *b_sc_eta;   //!
   TBranch        *b_sc_etacorr;   //!
   TBranch        *b_sc_theta;   //!
   TBranch        *b_sc_thetacorr;   //!
   TBranch        *b_sc_et;   //!
   TBranch        *b_sc_phi;   //!
   TBranch        *b_sc_px;   //!
   TBranch        *b_sc_py;   //!
   TBranch        *b_sc_pz;   //!
   TBranch        *b_sc_x;   //!
   TBranch        *b_sc_y;   //!
   TBranch        *b_sc_z;   //!
   TBranch        *b_ph_n;   //!
   TBranch        *b_ph_px;   //!
   TBranch        *b_ph_pt;   //!
   TBranch        *b_ph_eta;   //!
   TBranch        *b_ph_theta;   //!
   TBranch        *b_ph_phi;   //!
   TBranch        *b_ph_energy;   //!
   TBranch        *b_ph_mass;   //!
   TBranch        *b_ph_isPFlowPhoton;   //!
   TBranch        *b_ph_isStandardPhoton;   //!
   TBranch        *b_ph_hasConversionTracks;   //!
   TBranch        *b_ph_hasPixelSeed;   //!
   TBranch        *b_ph_isEB;   //!
   TBranch        *b_ph_isEE;   //!
   TBranch        *b_ph_isEBGap;   //!
   TBranch        *b_ph_isEBEtaGap;   //!
   TBranch        *b_ph_isEBPhiGap;   //!
   TBranch        *b_ph_isEEGap;   //!
   TBranch        *b_ph_isEERingGap;   //!
   TBranch        *b_ph_isEEDeeGap;   //!
   TBranch        *b_ph_isEBEEGap;   //!
   TBranch        *b_ph_hadronicOverEm;   //!
   TBranch        *b_ph_hadronicDepth1OverEm;   //!
   TBranch        *b_ph_hadronicDepth2OverEm;   //!
   TBranch        *b_ph_hadTowOverEm;   //!
   TBranch        *b_ph_hadTowDepth1OverEm;   //!
   TBranch        *b_ph_hadTowDepth2OverEm;   //!
   TBranch        *b_ph_e1x5;   //!
   TBranch        *b_ph_e2x5;   //!
   TBranch        *b_ph_e3x3;   //!
   TBranch        *b_ph_e5x5;   //!
   TBranch        *b_ph_maxEnergyXtal;   //!
   TBranch        *b_ph_sigmaEtaEta;   //!
   TBranch        *b_ph_sigmaIetaIeta;   //!
   TBranch        *b_ph_r1x5;   //!
   TBranch        *b_ph_r2x5;   //!
   TBranch        *b_ph_r9;   //!
   TBranch        *b_ph_mipChi2;   //!
   TBranch        *b_ph_mipTotEnergy;   //!
   TBranch        *b_ph_mipSlope;   //!
   TBranch        *b_ph_mipIntercept;   //!
   TBranch        *b_ph_mipNhitCone;   //!
   TBranch        *b_ph_mipIsHalo;   //!
   TBranch        *b_ph_ecalRecHitSumEtConeDR04;   //!
   TBranch        *b_ph_hcalTowerSumEtConeDR04;   //!
   TBranch        *b_ph_hcalDepth1TowerSumEtConeDR04;   //!
   TBranch        *b_ph_hcalDepth2TowerSumEtConeDR04;   //!
   TBranch        *b_ph_hcalTowerSumEtBcConeDR04;   //!
   TBranch        *b_ph_hcalDepth1TowerSumEtBcConeDR04;   //!
   TBranch        *b_ph_hcalDepth2TowerSumEtBcConeDR04;   //!
   TBranch        *b_ph_trkSumPtSolidConeDR04;   //!
   TBranch        *b_ph_trkSumPtHollowConeDR04;   //!
   TBranch        *b_ph_nTrkSolidConeDR04;   //!
   TBranch        *b_ph_nTrkHollowConeDR04;   //!
   TBranch        *b_ph_ecalRecHitSumEtConeDR03;   //!
   TBranch        *b_ph_hcalTowerSumEtConeDR03;   //!
   TBranch        *b_ph_hcalDepth1TowerSumEtConeDR03;   //!
   TBranch        *b_ph_hcalDepth2TowerSumEtConeDR03;   //!
   TBranch        *b_ph_hcalTowerSumEtBcConeDR03;   //!
   TBranch        *b_ph_hcalDepth1TowerSumEtBcConeDR03;   //!
   TBranch        *b_ph_hcalDepth2TowerSumEtBcConeDR03;   //!
   TBranch        *b_ph_trkSumPtSolidConeDR03;   //!
   TBranch        *b_ph_trkSumPtHollowConeDR03;   //!
   TBranch        *b_ph_nTrkSolidConeDR03;   //!
   TBranch        *b_ph_nTrkHollowConeDR03;   //!
   TBranch        *b_ph_chargedHadronIso;   //!
   TBranch        *b_ph_neutralHadronIso;   //!
   TBranch        *b_ph_photonIso;   //!
   TBranch        *b_ph_chargedHadronIsoWrongVtx;   //!
   TBranch        *b_ph_sumChargedParticlePt;   //!
   TBranch        *b_ph_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_ph_sumPhotonEtHighThreshold;   //!
   TBranch        *b_ph_sumPUPt;   //!
   TBranch        *b_ph_nClusterOutsideMustache;   //!
   TBranch        *b_ph_etOutsideMustache;   //!
   TBranch        *b_ph_pfMVA;   //!
   TBranch        *b_ph_mc_bestDR;   //!
   TBranch        *b_ph_mc_index;   //!
   TBranch        *b_ph_mc_ERatio;   //!
   TBranch        *b_gsf_n;   //!
   TBranch        *b_gsf_classification;   //!
   TBranch        *b_gsf_energy;   //!
   TBranch        *b_gsf_p;   //!
   TBranch        *b_gsf_pt;   //!
   TBranch        *b_gsf_scE1x5;   //!
   TBranch        *b_gsf_scE5x5;   //!
   TBranch        *b_gsf_scE2x5Max;   //!
   TBranch        *b_gsf_eta;   //!
   TBranch        *b_gsf_phi;   //!
   TBranch        *b_gsf_theta;   //!
   TBranch        *b_gsf_px;   //!
   TBranch        *b_gsf_py;   //!
   TBranch        *b_gsf_pz;   //!
   TBranch        *b_gsf_superClusterEta;   //!
   TBranch        *b_gsf_superClusterEnergy;   //!
   TBranch        *b_gsf_caloEnergy;   //!
   TBranch        *b_gsf_deltaEtaSuperClusterTrackAtVtx;   //!
   TBranch        *b_gsf_deltaPhiSuperClusterTrackAtVtx;   //!
   TBranch        *b_gsf_hadronicOverEm;   //!
   TBranch        *b_gsf_hcalDepth1OverEcal;   //!
   TBranch        *b_gsf_hcalDepth2OverEcal;   //!
   TBranch        *b_gsf_dr03TkSumPt;   //!
   TBranch        *b_gsf_dr03EcalRecHitSumEt;   //!
   TBranch        *b_gsf_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_gsf_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_gsf_charge;   //!
   TBranch        *b_gsf_sigmaIetaIeta;   //!
   TBranch        *b_gsf_ecaldrivenSeed;   //!
   TBranch        *b_gsf_trackerdrivenSeed;   //!
   TBranch        *b_gsf_isEB;   //!
   TBranch        *b_gsf_isEE;   //!
   TBranch        *b_gsf_deltaEtaSeedClusterTrackAtCalo;   //!
   TBranch        *b_gsf_deltaPhiSeedClusterTrackAtCalo;   //!
   TBranch        *b_gsf_ecalEnergy;   //!
   TBranch        *b_gsf_eSuperClusterOverP;   //!
   TBranch        *b_gsf_dxy;   //!
   TBranch        *b_gsf_dxy_beamSpot;   //!
   TBranch        *b_gsf_dxy_firstPVtx;   //!
   TBranch        *b_gsf_dxyError;   //!
   TBranch        *b_gsf_dz;   //!
   TBranch        *b_gsf_dz_beamSpot;   //!
   TBranch        *b_gsf_dz_firstPVtx;   //!
   TBranch        *b_gsf_dzError;   //!
   TBranch        *b_gsf_vz;   //!
   TBranch        *b_gsf_numberOfValidHits;   //!
   TBranch        *b_gsf_nLostInnerHits;   //!
   TBranch        *b_gsf_nLostOuterHits;   //!
   TBranch        *b_gsf_convFlags;   //!
   TBranch        *b_gsf_convDist;   //!
   TBranch        *b_gsf_convDcot;   //!
   TBranch        *b_gsf_convRadius;   //!
   TBranch        *b_gsf_fBrem;   //!
   TBranch        *b_gsf_e1x5;   //!
   TBranch        *b_gsf_e2x5Max;   //!
   TBranch        *b_gsf_e5x5;   //!
   TBranch        *b_gsf_r9;   //!
   TBranch        *b_gsf_hitsinfo;   //!
   TBranch        *b_gsf_pixelMatch_dPhi1;   //!
   TBranch        *b_gsf_pixelMatch_dPhi2;   //!
   TBranch        *b_gsf_pixelMatch_dRz1;   //!
   TBranch        *b_gsf_pixelMatch_dRz2;   //!
   TBranch        *b_gsf_pixelMatch_subDetector1;   //!
   TBranch        *b_gsf_pixelMatch_subDetector2;   //!
   TBranch        *b_gsf_mc_bestDR;   //!
   TBranch        *b_gsf_mc_index;   //!
   TBranch        *b_gsf_mc_ERatio;   //!
   TBranch        *b_mu_n;   //!
   TBranch        *b_mu_n_saved;   //!
   TBranch        *b_mu_gt_n;   //!
   TBranch        *b_mu_ot_n;   //!
   TBranch        *b_mu_it_n;   //!
   TBranch        *b_mu_gt_qoverp;   //!
   TBranch        *b_mu_gt_charge;   //!
   TBranch        *b_mu_gt_pt;   //!
   TBranch        *b_mu_gt_eta;   //!
   TBranch        *b_mu_gt_phi;   //!
   TBranch        *b_mu_gt_p;   //!
   TBranch        *b_mu_gt_px;   //!
   TBranch        *b_mu_gt_py;   //!
   TBranch        *b_mu_gt_pz;   //!
   TBranch        *b_mu_gt_theta;   //!
   TBranch        *b_mu_gt_lambda;   //!
   TBranch        *b_mu_gt_d0;   //!
   TBranch        *b_mu_gt_dz;   //!
   TBranch        *b_mu_gt_dz_beamspot;   //!
   TBranch        *b_mu_gt_dz_firstPVtx;   //!
   TBranch        *b_mu_gt_dxy;   //!
   TBranch        *b_mu_gt_dxy_beamspot;   //!
   TBranch        *b_mu_gt_dxy_firstPVtx;   //!
   TBranch        *b_mu_gt_dsz;   //!
   TBranch        *b_mu_gt_vx;   //!
   TBranch        *b_mu_gt_vy;   //!
   TBranch        *b_mu_gt_vz;   //!
   TBranch        *b_mu_gt_qoverpError;   //!
   TBranch        *b_mu_gt_ptError;   //!
   TBranch        *b_mu_gt_thetaError;   //!
   TBranch        *b_mu_gt_lambdaError;   //!
   TBranch        *b_mu_gt_phiError;   //!
   TBranch        *b_mu_gt_dxyError;   //!
   TBranch        *b_mu_gt_d0Error;   //!
   TBranch        *b_mu_gt_dszError;   //!
   TBranch        *b_mu_gt_dzError;   //!
   TBranch        *b_mu_gt_etaError;   //!
   TBranch        *b_mu_gt_chi2;   //!
   TBranch        *b_mu_gt_ndof;   //!
   TBranch        *b_mu_gt_normalizedChi2;   //!
   TBranch        *b_mu_ot_qoverp;   //!
   TBranch        *b_mu_ot_charge;   //!
   TBranch        *b_mu_ot_pt;   //!
   TBranch        *b_mu_ot_eta;   //!
   TBranch        *b_mu_ot_phi;   //!
   TBranch        *b_mu_ot_p;   //!
   TBranch        *b_mu_ot_px;   //!
   TBranch        *b_mu_ot_py;   //!
   TBranch        *b_mu_ot_pz;   //!
   TBranch        *b_mu_ot_theta;   //!
   TBranch        *b_mu_ot_lambda;   //!
   TBranch        *b_mu_ot_d0;   //!
   TBranch        *b_mu_ot_dz;   //!
   TBranch        *b_mu_ot_dz_beamspot;   //!
   TBranch        *b_mu_ot_dz_firstPVtx;   //!
   TBranch        *b_mu_ot_dxy;   //!
   TBranch        *b_mu_ot_dxy_beamspot;   //!
   TBranch        *b_mu_ot_dxy_firstPVtx;   //!
   TBranch        *b_mu_ot_dsz;   //!
   TBranch        *b_mu_ot_vx;   //!
   TBranch        *b_mu_ot_vy;   //!
   TBranch        *b_mu_ot_vz;   //!
   TBranch        *b_mu_ot_qoverpError;   //!
   TBranch        *b_mu_ot_ptError;   //!
   TBranch        *b_mu_ot_thetaError;   //!
   TBranch        *b_mu_ot_lambdaError;   //!
   TBranch        *b_mu_ot_phiError;   //!
   TBranch        *b_mu_ot_dxyError;   //!
   TBranch        *b_mu_ot_d0Error;   //!
   TBranch        *b_mu_ot_dszError;   //!
   TBranch        *b_mu_ot_dzError;   //!
   TBranch        *b_mu_ot_etaError;   //!
   TBranch        *b_mu_ot_chi2;   //!
   TBranch        *b_mu_ot_ndof;   //!
   TBranch        *b_mu_ot_normalizedChi2;   //!
   TBranch        *b_mu_it_qoverp;   //!
   TBranch        *b_mu_it_charge;   //!
   TBranch        *b_mu_it_pt;   //!
   TBranch        *b_mu_it_eta;   //!
   TBranch        *b_mu_it_phi;   //!
   TBranch        *b_mu_it_p;   //!
   TBranch        *b_mu_it_px;   //!
   TBranch        *b_mu_it_py;   //!
   TBranch        *b_mu_it_pz;   //!
   TBranch        *b_mu_it_theta;   //!
   TBranch        *b_mu_it_lambda;   //!
   TBranch        *b_mu_it_d0;   //!
   TBranch        *b_mu_it_dz;   //!
   TBranch        *b_mu_it_dz_beamspot;   //!
   TBranch        *b_mu_it_dz_firstPVtx;   //!
   TBranch        *b_mu_it_dxy;   //!
   TBranch        *b_mu_it_dxy_beamspot;   //!
   TBranch        *b_mu_it_dxy_firstPVtx;   //!
   TBranch        *b_mu_it_dsz;   //!
   TBranch        *b_mu_it_vx;   //!
   TBranch        *b_mu_it_vy;   //!
   TBranch        *b_mu_it_vz;   //!
   TBranch        *b_mu_it_qoverpError;   //!
   TBranch        *b_mu_it_ptError;   //!
   TBranch        *b_mu_it_thetaError;   //!
   TBranch        *b_mu_it_lambdaError;   //!
   TBranch        *b_mu_it_phiError;   //!
   TBranch        *b_mu_it_dxyError;   //!
   TBranch        *b_mu_it_d0Error;   //!
   TBranch        *b_mu_it_dszError;   //!
   TBranch        *b_mu_it_dzError;   //!
   TBranch        *b_mu_it_etaError;   //!
   TBranch        *b_mu_it_chi2;   //!
   TBranch        *b_mu_it_ndof;   //!
   TBranch        *b_mu_it_normalizedChi2;   //!
   TBranch        *b_mu_isGlobalMuon;   //!
   TBranch        *b_mu_isStandAloneMuon;   //!
   TBranch        *b_mu_isTrackerMuon;   //!
   TBranch        *b_mu_isPFMuon;   //!
   TBranch        *b_mu_isPFIsolationValid;   //!
   TBranch        *b_mu_numberOfMatchedStations;   //!
   TBranch        *b_mu_numberOfValidPixelHits;   //!
   TBranch        *b_mu_numberOfValidTrackerHits;   //!
   TBranch        *b_mu_numberOfValidMuonHits;   //!
   TBranch        *b_mu_tevOptimized_charge;   //!
   TBranch        *b_mu_tevOptimized_pt;   //!
   TBranch        *b_mu_tevOptimized_eta;   //!
   TBranch        *b_mu_tevOptimized_phi;   //!
   TBranch        *b_mu_tevOptimized_theta;   //!
   TBranch        *b_mu_tevOptimized_px;   //!
   TBranch        *b_mu_tevOptimized_py;   //!
   TBranch        *b_mu_tevOptimized_pz;   //!
   TBranch        *b_mu_tevOptimized_d0;   //!
   TBranch        *b_mu_tevOptimized_dz;   //!
   TBranch        *b_mu_tevOptimized_dz_beamSpot;   //!
   TBranch        *b_mu_tevOptimized_dz_firstPVtx;   //!
   TBranch        *b_mu_tevOptimized_dxy;   //!
   TBranch        *b_mu_tevOptimized_dxy_beamSpot;   //!
   TBranch        *b_mu_tevOptimized_dxy_firstPVtx;   //!
   TBranch        *b_mu_tevOptimized_ptError;   //!
   TBranch        *b_mu_tevOptimized_etaError;   //!
   TBranch        *b_mu_tevOptimized_phiError;   //!
   TBranch        *b_mu_tevOptimized_thetaError;   //!
   TBranch        *b_mu_tevOptimized_d0Error;   //!
   TBranch        *b_mu_tevOptimized_dzError;   //!
   TBranch        *b_mu_tevOptimized_dxyError;   //!
   TBranch        *b_mu_isolationR03_sumPt;   //!
   TBranch        *b_mu_isolationR03_trackerVetoPt;   //!
   TBranch        *b_mu_isolationR03_emEt;   //!
   TBranch        *b_mu_isolationR03_emVetoEt;   //!
   TBranch        *b_mu_isolationR03_hadEt;   //!
   TBranch        *b_mu_isolationR03_hadVetoEt;   //!
   TBranch        *b_mu_isolationR05_sumPt;   //!
   TBranch        *b_mu_isolationR05_trackerVetoPt;   //!
   TBranch        *b_mu_isolationR05_emEt;   //!
   TBranch        *b_mu_isolationR05_emVetoEt;   //!
   TBranch        *b_mu_isolationR05_hadEt;   //!
   TBranch        *b_mu_isolationR05_hadVetoEt;   //!
   TBranch        *b_mu_pfIsolationR03_sumChargedHadronPt;   //!
   TBranch        *b_mu_pfIsolationR03_sumChargedParticlePt;   //!
   TBranch        *b_mu_pfIsolationR03_sumPhotonEt;   //!
   TBranch        *b_mu_pfIsolationR03_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mu_pfIsolationR03_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mu_pfIsolationR03_sumPUPt;   //!
   TBranch        *b_mu_pfIsolationR04_sumChargedHadronPt;   //!
   TBranch        *b_mu_pfIsolationR04_sumChargedParticlePt;   //!
   TBranch        *b_mu_pfIsolationR04_sumPhotonEt;   //!
   TBranch        *b_mu_pfIsolationR04_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mu_pfIsolationR04_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mu_pfIsolationR04_sumPUPt;   //!
   TBranch        *b_mu_pfMeanDRIsoProfileR03_sumChargedHadronPt;   //!
   TBranch        *b_mu_pfMeanDRIsoProfileR03_sumChargedParticlePt;   //!
   TBranch        *b_mu_pfMeanDRIsoProfileR03_sumPhotonEt;   //!
   TBranch        *b_mu_pfMeanDRIsoProfileR03_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mu_pfMeanDRIsoProfileR03_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mu_pfMeanDRIsoProfileR03_sumPUPt;   //!
   TBranch        *b_mu_pfMeanDRIsoProfileR04_sumChargedHadronPt;   //!
   TBranch        *b_mu_pfMeanDRIsoProfileR04_sumChargedParticlePt;   //!
   TBranch        *b_mu_pfMeanDRIsoProfileR04_sumPhotonEt;   //!
   TBranch        *b_mu_pfMeanDRIsoProfileR04_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mu_pfMeanDRIsoProfileR04_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mu_pfMeanDRIsoProfileR04_sumPUPt;   //!
   TBranch        *b_mu_mc_bestDR;   //!
   TBranch        *b_mu_mc_index;   //!
   TBranch        *b_mu_mc_ERatio;   //!
   TBranch        *b_MET_caloMet_et;   //!
   TBranch        *b_MET_caloMet_phi;   //!
   TBranch        *b_MET_pfMet_et;   //!
   TBranch        *b_MET_pfMet_phi;   //!
   TBranch        *b_HEEP_eseffsixix;   //!
   TBranch        *b_HEEP_eseffsiyiy;   //!
   TBranch        *b_HEEP_eseffsirir;   //!
   TBranch        *b_HEEP_preshowerEnergy;   //!
   TBranch        *b_HEEP_e1x3;   //!
   TBranch        *b_HEEP_eMax;   //!
   TBranch        *b_HEEP_e5x5;   //!
   TBranch        *b_HEEP_e2x5Right;   //!
   TBranch        *b_HEEP_e2x5Left;   //!
   TBranch        *b_HEEP_e2x5Top;   //!
   TBranch        *b_HEEP_e2x5Bottom;   //!
   TBranch        *b_HEEP_eRight;   //!
   TBranch        *b_HEEP_eLeft;   //!
   TBranch        *b_HEEP_eTop;   //!
   TBranch        *b_HEEP_eBottom;   //!
   TBranch        *b_HEEP_basicClusterSeedTime;   //!
   TBranch        *b_EBHits_rawId;   //!
   TBranch        *b_EBHits_iRechit;   //!
   TBranch        *b_EBHits_energy;   //!
   TBranch        *b_EBHits_ieta;   //!
   TBranch        *b_EBHits_iphi;   //!
   TBranch        *b_EBHits_RecoFlag;   //!
   TBranch        *b_EBHits_kSaturated;   //!
   TBranch        *b_EBHits_kLeadingEdgeRecovered;   //!
   TBranch        *b_EBHits_kNeighboursRecovered;   //!
   TBranch        *b_EBHits_kWeird;   //!
   TBranch        *b_EEHits_rawId;   //!
   TBranch        *b_EEHits_iRechit;   //!
   TBranch        *b_EEHits_energy;   //!
   TBranch        *b_EEHits_ieta;   //!
   TBranch        *b_EEHits_iphi;   //!
   TBranch        *b_EEHits_RecoFlag;   //!
   TBranch        *b_EEHits_kSaturated;   //!
   TBranch        *b_EEHits_kLeadingEdgeRecovered;   //!
   TBranch        *b_EEHits_kNeighboursRecovered;   //!
   TBranch        *b_EEHits_kWeird;   //!
   TBranch        *b_HEEP_crystal_energy;   //!
   TBranch        *b_HEEP_crystal_eta;   //!
   TBranch        *b_HEEP_eshitsixix;   //!
   TBranch        *b_HEEP_eshitsiyiy;   //!
   TBranch        *b_HEEP_crystal_ietaorix;   //!
   TBranch        *b_HEEP_crystal_iphioriy;   //!
   TBranch        *b_HEEP_cutflow41_Et;   //!
   TBranch        *b_HEEP_cutflow41_Et_n;   //!
   TBranch        *b_HEEP_cutflow41_Et_nCumulative;   //!
   TBranch        *b_HEEP_cutflow41_Et_value;   //!
   TBranch        *b_HEEP_cutflow41_eta;   //!
   TBranch        *b_HEEP_cutflow41_eta_n;   //!
   TBranch        *b_HEEP_cutflow41_eta_nCumulative;   //!
   TBranch        *b_HEEP_cutflow41_eta_value;   //!
   TBranch        *b_HEEP_cutflow41_acceptance;   //!
   TBranch        *b_HEEP_cutflow41_acceptance_n;   //!
   TBranch        *b_HEEP_cutflow41_EcalDriven;   //!
   TBranch        *b_HEEP_cutflow41_EcalDriven_n;   //!
   TBranch        *b_HEEP_cutflow41_EcalDriven_nCumulative;   //!
   TBranch        *b_HEEP_cutflow41_EcalDriven_value;   //!
   TBranch        *b_HEEP_cutflow41_dEtaIn;   //!
   TBranch        *b_HEEP_cutflow41_dEtaIn_n;   //!
   TBranch        *b_HEEP_cutflow41_dEtaIn_nCumulative;   //!
   TBranch        *b_HEEP_cutflow41_dEtaIn_value;   //!
   TBranch        *b_HEEP_cutflow41_dPhiIn;   //!
   TBranch        *b_HEEP_cutflow41_dPhiIn_n;   //!
   TBranch        *b_HEEP_cutflow41_dPhiIn_nCumulative;   //!
   TBranch        *b_HEEP_cutflow41_dPhiIn_value;   //!
   TBranch        *b_HEEP_cutflow41_HOverE;   //!
   TBranch        *b_HEEP_cutflow41_HOverE_n;   //!
   TBranch        *b_HEEP_cutflow41_HOverE_nCumulative;   //!
   TBranch        *b_HEEP_cutflow41_HOverE_value;   //!
   TBranch        *b_HEEP_cutflow41_SigmaIetaIeta;   //!
   TBranch        *b_HEEP_cutflow41_SigmaIetaIeta_n;   //!
   TBranch        *b_HEEP_cutflow41_SigmaIetaIeta_nCumulative;   //!
   TBranch        *b_HEEP_cutflow41_SigmaIetaIeta_value;   //!
   TBranch        *b_HEEP_cutflow41_E1x5OverE5x5;   //!
   TBranch        *b_HEEP_cutflow41_E1x5OverE5x5_n;   //!
   TBranch        *b_HEEP_cutflow41_E1x5OverE5x5_nCumulative;   //!
   TBranch        *b_HEEP_cutflow41_E1x5OverE5x5_value;   //!
   TBranch        *b_HEEP_cutflow41_E2x5OverE5x5;   //!
   TBranch        *b_HEEP_cutflow41_E2x5OverE5x5_n;   //!
   TBranch        *b_HEEP_cutflow41_E2x5OverE5x5_nCumulative;   //!
   TBranch        *b_HEEP_cutflow41_E2x5OverE5x5_value;   //!
   TBranch        *b_HEEP_cutflow41_missingHits;   //!
   TBranch        *b_HEEP_cutflow41_missingHits_n;   //!
   TBranch        *b_HEEP_cutflow41_missingHits_nCumulative;   //!
   TBranch        *b_HEEP_cutflow41_missingHits_value;   //!
   TBranch        *b_HEEP_cutflow41_dxyFirstPV;   //!
   TBranch        *b_HEEP_cutflow41_dxyFirstPV_n;   //!
   TBranch        *b_HEEP_cutflow41_dxyFirstPV_nCumulative;   //!
   TBranch        *b_HEEP_cutflow41_dxyFirstPV_value;   //!
   TBranch        *b_HEEP_cutflow41_ID;   //!
   TBranch        *b_HEEP_cutflow41_ID_n;   //!
   TBranch        *b_HEEP_cutflow41_isolEMHadDepth1;   //!
   TBranch        *b_HEEP_cutflow41_isolEMHadDepth1_n;   //!
   TBranch        *b_HEEP_cutflow41_isolEMHadDepth1_nCumulative;   //!
   TBranch        *b_HEEP_cutflow41_isolEMHadDepth1_value;   //!
   TBranch        *b_HEEP_cutflow41_IsolPtTrks;   //!
   TBranch        *b_HEEP_cutflow41_IsolPtTrks_n;   //!
   TBranch        *b_HEEP_cutflow41_IsolPtTrks_nCumulative;   //!
   TBranch        *b_HEEP_cutflow41_IsolPtTrks_value;   //!
   TBranch        *b_HEEP_cutflow41_isolation;   //!
   TBranch        *b_HEEP_cutflow41_isolation_n;   //!
   TBranch        *b_HEEP_cutflow41_total;   //!
   TBranch        *b_HEEP_cutflow41_total_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_Et;   //!
   TBranch        *b_HEEP_cutflow50_25ns_Et_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_Et_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_25ns_Et_value;   //!
   TBranch        *b_HEEP_cutflow50_25ns_eta;   //!
   TBranch        *b_HEEP_cutflow50_25ns_eta_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_eta_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_25ns_eta_value;   //!
   TBranch        *b_HEEP_cutflow50_25ns_acceptance;   //!
   TBranch        *b_HEEP_cutflow50_25ns_acceptance_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_EcalDriven;   //!
   TBranch        *b_HEEP_cutflow50_25ns_EcalDriven_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_EcalDriven_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_25ns_EcalDriven_value;   //!
   TBranch        *b_HEEP_cutflow50_25ns_dEtaIn;   //!
   TBranch        *b_HEEP_cutflow50_25ns_dEtaIn_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_dEtaIn_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_25ns_dEtaIn_value;   //!
   TBranch        *b_HEEP_cutflow50_25ns_dPhiIn;   //!
   TBranch        *b_HEEP_cutflow50_25ns_dPhiIn_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_dPhiIn_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_25ns_dPhiIn_value;   //!
   TBranch        *b_HEEP_cutflow50_25ns_HOverE;   //!
   TBranch        *b_HEEP_cutflow50_25ns_HOverE_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_HOverE_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_25ns_HOverE_value;   //!
   TBranch        *b_HEEP_cutflow50_25ns_SigmaIetaIeta;   //!
   TBranch        *b_HEEP_cutflow50_25ns_SigmaIetaIeta_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_SigmaIetaIeta_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_25ns_SigmaIetaIeta_value;   //!
   TBranch        *b_HEEP_cutflow50_25ns_E1x5OverE5x5;   //!
   TBranch        *b_HEEP_cutflow50_25ns_E1x5OverE5x5_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_E1x5OverE5x5_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_25ns_E1x5OverE5x5_value;   //!
   TBranch        *b_HEEP_cutflow50_25ns_E2x5OverE5x5;   //!
   TBranch        *b_HEEP_cutflow50_25ns_E2x5OverE5x5_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_E2x5OverE5x5_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_25ns_E2x5OverE5x5_value;   //!
   TBranch        *b_HEEP_cutflow50_25ns_missingHits;   //!
   TBranch        *b_HEEP_cutflow50_25ns_missingHits_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_missingHits_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_25ns_missingHits_value;   //!
   TBranch        *b_HEEP_cutflow50_25ns_dxyFirstPV;   //!
   TBranch        *b_HEEP_cutflow50_25ns_dxyFirstPV_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_dxyFirstPV_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_25ns_dxyFirstPV_value;   //!
   TBranch        *b_HEEP_cutflow50_25ns_ID;   //!
   TBranch        *b_HEEP_cutflow50_25ns_ID_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_isolEMHadDepth1;   //!
   TBranch        *b_HEEP_cutflow50_25ns_isolEMHadDepth1_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_isolEMHadDepth1_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_25ns_isolEMHadDepth1_value;   //!
   TBranch        *b_HEEP_cutflow50_25ns_IsolPtTrks;   //!
   TBranch        *b_HEEP_cutflow50_25ns_IsolPtTrks_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_IsolPtTrks_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_25ns_IsolPtTrks_value;   //!
   TBranch        *b_HEEP_cutflow50_25ns_isolation;   //!
   TBranch        *b_HEEP_cutflow50_25ns_isolation_n;   //!
   TBranch        *b_HEEP_cutflow50_25ns_total;   //!
   TBranch        *b_HEEP_cutflow50_25ns_total_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_Et;   //!
   TBranch        *b_HEEP_cutflow50_50ns_Et_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_Et_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_50ns_Et_value;   //!
   TBranch        *b_HEEP_cutflow50_50ns_eta;   //!
   TBranch        *b_HEEP_cutflow50_50ns_eta_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_eta_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_50ns_eta_value;   //!
   TBranch        *b_HEEP_cutflow50_50ns_acceptance;   //!
   TBranch        *b_HEEP_cutflow50_50ns_acceptance_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_EcalDriven;   //!
   TBranch        *b_HEEP_cutflow50_50ns_EcalDriven_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_EcalDriven_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_50ns_EcalDriven_value;   //!
   TBranch        *b_HEEP_cutflow50_50ns_dEtaIn;   //!
   TBranch        *b_HEEP_cutflow50_50ns_dEtaIn_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_dEtaIn_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_50ns_dEtaIn_value;   //!
   TBranch        *b_HEEP_cutflow50_50ns_dPhiIn;   //!
   TBranch        *b_HEEP_cutflow50_50ns_dPhiIn_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_dPhiIn_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_50ns_dPhiIn_value;   //!
   TBranch        *b_HEEP_cutflow50_50ns_HOverE;   //!
   TBranch        *b_HEEP_cutflow50_50ns_HOverE_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_HOverE_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_50ns_HOverE_value;   //!
   TBranch        *b_HEEP_cutflow50_50ns_SigmaIetaIeta;   //!
   TBranch        *b_HEEP_cutflow50_50ns_SigmaIetaIeta_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_SigmaIetaIeta_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_50ns_SigmaIetaIeta_value;   //!
   TBranch        *b_HEEP_cutflow50_50ns_E1x5OverE5x5;   //!
   TBranch        *b_HEEP_cutflow50_50ns_E1x5OverE5x5_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_E1x5OverE5x5_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_50ns_E1x5OverE5x5_value;   //!
   TBranch        *b_HEEP_cutflow50_50ns_E2x5OverE5x5;   //!
   TBranch        *b_HEEP_cutflow50_50ns_E2x5OverE5x5_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_E2x5OverE5x5_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_50ns_E2x5OverE5x5_value;   //!
   TBranch        *b_HEEP_cutflow50_50ns_missingHits;   //!
   TBranch        *b_HEEP_cutflow50_50ns_missingHits_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_missingHits_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_50ns_missingHits_value;   //!
   TBranch        *b_HEEP_cutflow50_50ns_dxyFirstPV;   //!
   TBranch        *b_HEEP_cutflow50_50ns_dxyFirstPV_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_dxyFirstPV_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_50ns_dxyFirstPV_value;   //!
   TBranch        *b_HEEP_cutflow50_50ns_ID;   //!
   TBranch        *b_HEEP_cutflow50_50ns_ID_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_isolEMHadDepth1;   //!
   TBranch        *b_HEEP_cutflow50_50ns_isolEMHadDepth1_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_isolEMHadDepth1_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_50ns_isolEMHadDepth1_value;   //!
   TBranch        *b_HEEP_cutflow50_50ns_IsolPtTrks;   //!
   TBranch        *b_HEEP_cutflow50_50ns_IsolPtTrks_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_IsolPtTrks_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_50ns_IsolPtTrks_value;   //!
   TBranch        *b_HEEP_cutflow50_50ns_isolation;   //!
   TBranch        *b_HEEP_cutflow50_50ns_isolation_n;   //!
   TBranch        *b_HEEP_cutflow50_50ns_total;   //!
   TBranch        *b_HEEP_cutflow50_50ns_total_n;   //!
   TBranch        *b_HEEP_cutflow50_Et;   //!
   TBranch        *b_HEEP_cutflow50_Et_n;   //!
   TBranch        *b_HEEP_cutflow50_Et_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_Et_value;   //!
   TBranch        *b_HEEP_cutflow50_eta;   //!
   TBranch        *b_HEEP_cutflow50_eta_n;   //!
   TBranch        *b_HEEP_cutflow50_eta_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_eta_value;   //!
   TBranch        *b_HEEP_cutflow50_acceptance;   //!
   TBranch        *b_HEEP_cutflow50_acceptance_n;   //!
   TBranch        *b_HEEP_cutflow50_EcalDriven;   //!
   TBranch        *b_HEEP_cutflow50_EcalDriven_n;   //!
   TBranch        *b_HEEP_cutflow50_EcalDriven_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_EcalDriven_value;   //!
   TBranch        *b_HEEP_cutflow50_dEtaIn;   //!
   TBranch        *b_HEEP_cutflow50_dEtaIn_n;   //!
   TBranch        *b_HEEP_cutflow50_dEtaIn_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_dEtaIn_value;   //!
   TBranch        *b_HEEP_cutflow50_dPhiIn;   //!
   TBranch        *b_HEEP_cutflow50_dPhiIn_n;   //!
   TBranch        *b_HEEP_cutflow50_dPhiIn_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_dPhiIn_value;   //!
   TBranch        *b_HEEP_cutflow50_HOverE;   //!
   TBranch        *b_HEEP_cutflow50_HOverE_n;   //!
   TBranch        *b_HEEP_cutflow50_HOverE_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_HOverE_value;   //!
   TBranch        *b_HEEP_cutflow50_SigmaIetaIeta;   //!
   TBranch        *b_HEEP_cutflow50_SigmaIetaIeta_n;   //!
   TBranch        *b_HEEP_cutflow50_SigmaIetaIeta_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_SigmaIetaIeta_value;   //!
   TBranch        *b_HEEP_cutflow50_E1x5OverE5x5;   //!
   TBranch        *b_HEEP_cutflow50_E1x5OverE5x5_n;   //!
   TBranch        *b_HEEP_cutflow50_E1x5OverE5x5_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_E1x5OverE5x5_value;   //!
   TBranch        *b_HEEP_cutflow50_E2x5OverE5x5;   //!
   TBranch        *b_HEEP_cutflow50_E2x5OverE5x5_n;   //!
   TBranch        *b_HEEP_cutflow50_E2x5OverE5x5_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_E2x5OverE5x5_value;   //!
   TBranch        *b_HEEP_cutflow50_missingHits;   //!
   TBranch        *b_HEEP_cutflow50_missingHits_n;   //!
   TBranch        *b_HEEP_cutflow50_missingHits_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_missingHits_value;   //!
   TBranch        *b_HEEP_cutflow50_dxyFirstPV;   //!
   TBranch        *b_HEEP_cutflow50_dxyFirstPV_n;   //!
   TBranch        *b_HEEP_cutflow50_dxyFirstPV_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_dxyFirstPV_value;   //!
   TBranch        *b_HEEP_cutflow50_ID;   //!
   TBranch        *b_HEEP_cutflow50_ID_n;   //!
   TBranch        *b_HEEP_cutflow50_isolEMHadDepth1;   //!
   TBranch        *b_HEEP_cutflow50_isolEMHadDepth1_n;   //!
   TBranch        *b_HEEP_cutflow50_isolEMHadDepth1_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_isolEMHadDepth1_value;   //!
   TBranch        *b_HEEP_cutflow50_IsolPtTrks;   //!
   TBranch        *b_HEEP_cutflow50_IsolPtTrks_n;   //!
   TBranch        *b_HEEP_cutflow50_IsolPtTrks_nCumulative;   //!
   TBranch        *b_HEEP_cutflow50_IsolPtTrks_value;   //!
   TBranch        *b_HEEP_cutflow50_isolation;   //!
   TBranch        *b_HEEP_cutflow50_isolation_n;   //!
   TBranch        *b_HEEP_cutflow50_total;   //!
   TBranch        *b_HEEP_cutflow50_total_n;   //!
   TBranch        *b_HEEP_cutflow51_Et;   //!
   TBranch        *b_HEEP_cutflow51_Et_n;   //!
   TBranch        *b_HEEP_cutflow51_Et_nCumulative;   //!
   TBranch        *b_HEEP_cutflow51_Et_value;   //!
   TBranch        *b_HEEP_cutflow51_eta;   //!
   TBranch        *b_HEEP_cutflow51_eta_n;   //!
   TBranch        *b_HEEP_cutflow51_eta_nCumulative;   //!
   TBranch        *b_HEEP_cutflow51_eta_value;   //!
   TBranch        *b_HEEP_cutflow51_acceptance;   //!
   TBranch        *b_HEEP_cutflow51_acceptance_n;   //!
   TBranch        *b_HEEP_cutflow51_EcalDriven;   //!
   TBranch        *b_HEEP_cutflow51_EcalDriven_n;   //!
   TBranch        *b_HEEP_cutflow51_EcalDriven_nCumulative;   //!
   TBranch        *b_HEEP_cutflow51_EcalDriven_value;   //!
   TBranch        *b_HEEP_cutflow51_dEtaIn;   //!
   TBranch        *b_HEEP_cutflow51_dEtaIn_n;   //!
   TBranch        *b_HEEP_cutflow51_dEtaIn_nCumulative;   //!
   TBranch        *b_HEEP_cutflow51_dEtaIn_value;   //!
   TBranch        *b_HEEP_cutflow51_dPhiIn;   //!
   TBranch        *b_HEEP_cutflow51_dPhiIn_n;   //!
   TBranch        *b_HEEP_cutflow51_dPhiIn_nCumulative;   //!
   TBranch        *b_HEEP_cutflow51_dPhiIn_value;   //!
   TBranch        *b_HEEP_cutflow51_HOverE;   //!
   TBranch        *b_HEEP_cutflow51_HOverE_n;   //!
   TBranch        *b_HEEP_cutflow51_HOverE_nCumulative;   //!
   TBranch        *b_HEEP_cutflow51_HOverE_value;   //!
   TBranch        *b_HEEP_cutflow51_SigmaIetaIeta;   //!
   TBranch        *b_HEEP_cutflow51_SigmaIetaIeta_n;   //!
   TBranch        *b_HEEP_cutflow51_SigmaIetaIeta_nCumulative;   //!
   TBranch        *b_HEEP_cutflow51_SigmaIetaIeta_value;   //!
   TBranch        *b_HEEP_cutflow51_E1x5OverE5x5;   //!
   TBranch        *b_HEEP_cutflow51_E1x5OverE5x5_n;   //!
   TBranch        *b_HEEP_cutflow51_E1x5OverE5x5_nCumulative;   //!
   TBranch        *b_HEEP_cutflow51_E1x5OverE5x5_value;   //!
   TBranch        *b_HEEP_cutflow51_E2x5OverE5x5;   //!
   TBranch        *b_HEEP_cutflow51_E2x5OverE5x5_n;   //!
   TBranch        *b_HEEP_cutflow51_E2x5OverE5x5_nCumulative;   //!
   TBranch        *b_HEEP_cutflow51_E2x5OverE5x5_value;   //!
   TBranch        *b_HEEP_cutflow51_missingHits;   //!
   TBranch        *b_HEEP_cutflow51_missingHits_n;   //!
   TBranch        *b_HEEP_cutflow51_missingHits_nCumulative;   //!
   TBranch        *b_HEEP_cutflow51_missingHits_value;   //!
   TBranch        *b_HEEP_cutflow51_dxyFirstPV;   //!
   TBranch        *b_HEEP_cutflow51_dxyFirstPV_n;   //!
   TBranch        *b_HEEP_cutflow51_dxyFirstPV_nCumulative;   //!
   TBranch        *b_HEEP_cutflow51_dxyFirstPV_value;   //!
   TBranch        *b_HEEP_cutflow51_ID;   //!
   TBranch        *b_HEEP_cutflow51_ID_n;   //!
   TBranch        *b_HEEP_cutflow51_isolEMHadDepth1;   //!
   TBranch        *b_HEEP_cutflow51_isolEMHadDepth1_n;   //!
   TBranch        *b_HEEP_cutflow51_isolEMHadDepth1_nCumulative;   //!
   TBranch        *b_HEEP_cutflow51_isolEMHadDepth1_value;   //!
   TBranch        *b_HEEP_cutflow51_IsolPtTrks;   //!
   TBranch        *b_HEEP_cutflow51_IsolPtTrks_n;   //!
   TBranch        *b_HEEP_cutflow51_IsolPtTrks_nCumulative;   //!
   TBranch        *b_HEEP_cutflow51_IsolPtTrks_value;   //!
   TBranch        *b_HEEP_cutflow51_isolation;   //!
   TBranch        *b_HEEP_cutflow51_isolation_n;   //!
   TBranch        *b_HEEP_cutflow51_total;   //!
   TBranch        *b_HEEP_cutflow51_total_n;   //!
   TBranch        *b_HEEP_cutflow60_Et;   //!
   TBranch        *b_HEEP_cutflow60_Et_n;   //!
   TBranch        *b_HEEP_cutflow60_Et_nCumulative;   //!
   TBranch        *b_HEEP_cutflow60_Et_value;   //!
   TBranch        *b_HEEP_cutflow60_eta;   //!
   TBranch        *b_HEEP_cutflow60_eta_n;   //!
   TBranch        *b_HEEP_cutflow60_eta_nCumulative;   //!
   TBranch        *b_HEEP_cutflow60_eta_value;   //!
   TBranch        *b_HEEP_cutflow60_acceptance;   //!
   TBranch        *b_HEEP_cutflow60_acceptance_n;   //!
   TBranch        *b_HEEP_cutflow60_EcalDriven;   //!
   TBranch        *b_HEEP_cutflow60_EcalDriven_n;   //!
   TBranch        *b_HEEP_cutflow60_EcalDriven_nCumulative;   //!
   TBranch        *b_HEEP_cutflow60_EcalDriven_value;   //!
   TBranch        *b_HEEP_cutflow60_dEtaIn;   //!
   TBranch        *b_HEEP_cutflow60_dEtaIn_n;   //!
   TBranch        *b_HEEP_cutflow60_dEtaIn_nCumulative;   //!
   TBranch        *b_HEEP_cutflow60_dEtaIn_value;   //!
   TBranch        *b_HEEP_cutflow60_dPhiIn;   //!
   TBranch        *b_HEEP_cutflow60_dPhiIn_n;   //!
   TBranch        *b_HEEP_cutflow60_dPhiIn_nCumulative;   //!
   TBranch        *b_HEEP_cutflow60_dPhiIn_value;   //!
   TBranch        *b_HEEP_cutflow60_HOverE;   //!
   TBranch        *b_HEEP_cutflow60_HOverE_n;   //!
   TBranch        *b_HEEP_cutflow60_HOverE_nCumulative;   //!
   TBranch        *b_HEEP_cutflow60_HOverE_value;   //!
   TBranch        *b_HEEP_cutflow60_SigmaIetaIeta;   //!
   TBranch        *b_HEEP_cutflow60_SigmaIetaIeta_n;   //!
   TBranch        *b_HEEP_cutflow60_SigmaIetaIeta_nCumulative;   //!
   TBranch        *b_HEEP_cutflow60_SigmaIetaIeta_value;   //!
   TBranch        *b_HEEP_cutflow60_E1x5OverE5x5;   //!
   TBranch        *b_HEEP_cutflow60_E1x5OverE5x5_n;   //!
   TBranch        *b_HEEP_cutflow60_E1x5OverE5x5_nCumulative;   //!
   TBranch        *b_HEEP_cutflow60_E1x5OverE5x5_value;   //!
   TBranch        *b_HEEP_cutflow60_E2x5OverE5x5;   //!
   TBranch        *b_HEEP_cutflow60_E2x5OverE5x5_n;   //!
   TBranch        *b_HEEP_cutflow60_E2x5OverE5x5_nCumulative;   //!
   TBranch        *b_HEEP_cutflow60_E2x5OverE5x5_value;   //!
   TBranch        *b_HEEP_cutflow60_missingHits;   //!
   TBranch        *b_HEEP_cutflow60_missingHits_n;   //!
   TBranch        *b_HEEP_cutflow60_missingHits_nCumulative;   //!
   TBranch        *b_HEEP_cutflow60_missingHits_value;   //!
   TBranch        *b_HEEP_cutflow60_dxyFirstPV;   //!
   TBranch        *b_HEEP_cutflow60_dxyFirstPV_n;   //!
   TBranch        *b_HEEP_cutflow60_dxyFirstPV_nCumulative;   //!
   TBranch        *b_HEEP_cutflow60_dxyFirstPV_value;   //!
   TBranch        *b_HEEP_cutflow60_ID;   //!
   TBranch        *b_HEEP_cutflow60_ID_n;   //!
   TBranch        *b_HEEP_cutflow60_isolEMHadDepth1;   //!
   TBranch        *b_HEEP_cutflow60_isolEMHadDepth1_n;   //!
   TBranch        *b_HEEP_cutflow60_isolEMHadDepth1_nCumulative;   //!
   TBranch        *b_HEEP_cutflow60_isolEMHadDepth1_value;   //!
   TBranch        *b_HEEP_cutflow60_IsolPtTrks;   //!
   TBranch        *b_HEEP_cutflow60_IsolPtTrks_n;   //!
   TBranch        *b_HEEP_cutflow60_IsolPtTrks_nCumulative;   //!
   TBranch        *b_HEEP_cutflow60_IsolPtTrks_value;   //!
   TBranch        *b_HEEP_cutflow60_isolation;   //!
   TBranch        *b_HEEP_cutflow60_isolation_n;   //!
   TBranch        *b_HEEP_cutflow60_total;   //!
   TBranch        *b_HEEP_cutflow60_total_n;   //!
   TBranch        *b_Zee_n;   //!
   TBranch        *b_Zee_mass;   //!
   TBranch        *b_Zee_mass_HEEP;   //!
   TBranch        *b_Zee_i1;   //!
   TBranch        *b_Zee_i2;   //!
   TBranch        *b_Zee_highestMass;   //!
   TBranch        *b_Zmm_n;   //!
   TBranch        *b_Zmm_mass;   //!
   TBranch        *b_Zmm_i1;   //!
   TBranch        *b_Zmm_i2;   //!
   TBranch        *b_Zmm_highestMass;   //!
   TBranch        *b_Zem_n;   //!
   TBranch        *b_Zem_mass;   //!
   TBranch        *b_Zem_mass_HEEP;   //!
   TBranch        *b_Zem_i1;   //!
   TBranch        *b_Zem_i2;   //!
   TBranch        *b_Zem_highestMass;   //!
   TBranch        *b_trk_n;   //!
   TBranch        *b_trk_pt;   //!
   TBranch        *b_trk_eta;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_prescale;   //!
   TBranch        *b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltL1sL1DoubleEG2210_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltL1sL1DoubleEG2210_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg1TrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg1TrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg2TrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg2TrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_prescale;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLPixelMatchFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLPixelMatchFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLNewPixelMatchFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLNewPixelMatchFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEG33EtUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEG33EtUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLNewPixelMatchUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLNewPixelMatchUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_prescale;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLPixelMatchFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLPixelMatchFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEG33EtUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEG33EtUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_accept;   //!
   TBranch        *b_trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_prescale;   //!
   TBranch        *b_trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_eta;   //!
   TBranch        *b_trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_phi;   //!
   TBranch        *b_trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltSingleEle22WPLooseGsfTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltSingleEle22WPLooseGsfTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_accept;   //!
   TBranch        *b_trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_prescale;   //!
   TBranch        *b_trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_eta;   //!
   TBranch        *b_trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_phi;   //!
   TBranch        *b_trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltSingleEle22WPTightGsfTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltSingleEle22WPTightGsfTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_accept;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_prescale;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltL1sL1SingleEG25OrSingleIsoEG30er_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltL1sL1SingleEG25OrSingleIsoEG30er_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8ClusterShapeFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8ClusterShapeFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HEFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HEFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HcalIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HcalIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8OneOEMinusOneOPFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8OneOEMinusOneOPFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DetaFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DetaFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DphiFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DphiFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8TrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8TrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchUnseededFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchUnseededFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8Mass55Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8Mass55Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_v2_accept;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_v2_prescale;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_v2_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_eta;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_v2_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_phi;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_v2_hltEle23WPLooseGsfTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_v2_hltEle23WPLooseGsfTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_accept;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_prescale;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_eta;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_phi;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltEle23WPLooseGsfTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltEle23WPLooseGsfTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltJetFilterEle23WPLoose_eta;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltJetFilterEle23WPLoose_phi;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltHCand80NoEle23WPLoose_eta;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltHCand80NoEle23WPLoose_phi;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMhtFilter0_eta;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMhtFilter0_phi;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMetFilter0_eta;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMetFilter0_phi;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand80NoEle23WPLooseMET_eta;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand80NoEle23WPLooseMET_phi;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand70NoEle23WPLooseMHTIDTight_eta;   //!
   TBranch        *b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand70NoEle23WPLooseMHTIDTight_phi;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_v1_accept;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_v1_prescale;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_v1_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_eta;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_v1_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_phi;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_v1_hltEle27noerWPLooseGsfTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_v1_hltEle27noerWPLooseGsfTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_accept;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_prescale;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_eta;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_phi;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltEle27noerWPLooseGsfTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltEle27noerWPLooseGsfTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltJetFilterEle27WPLoose_eta;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltJetFilterEle27WPLoose_phi;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltHCand80NoEle27WPLoose_eta;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltHCand80NoEle27WPLoose_phi;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMhtFilter0_eta;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMhtFilter0_phi;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMetFilter0_eta;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMetFilter0_phi;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand80NoEle27WPLooseMET_eta;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand80NoEle27WPLooseMET_phi;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand70NoEle27WPLooseMHTIDTight_eta;   //!
   TBranch        *b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand70NoEle27WPLooseMHTIDTight_phi;   //!
   TBranch        *b_trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_accept;   //!
   TBranch        *b_trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_prescale;   //!
   TBranch        *b_trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_eta;   //!
   TBranch        *b_trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_phi;   //!
   TBranch        *b_trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltEle27WPLooseGsfTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltEle27WPLooseGsfTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_accept;   //!
   TBranch        *b_trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_prescale;   //!
   TBranch        *b_trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_eta;   //!
   TBranch        *b_trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_phi;   //!
   TBranch        *b_trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltEle27WPTightGsfTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltEle27WPTightGsfTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_accept;   //!
   TBranch        *b_trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_prescale;   //!
   TBranch        *b_trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG30erOrSingleEG40_eta;   //!
   TBranch        *b_trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG30erOrSingleEG40_phi;   //!
   TBranch        *b_trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltEle32WPTightGsfTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltEle32WPTightGsfTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_accept;   //!
   TBranch        *b_trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_prescale;   //!
   TBranch        *b_trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_eta;   //!
   TBranch        *b_trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_phi;   //!
   TBranch        *b_trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_accept;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_prescale;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG2210_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG2210_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_accept;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_prescale;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG1510_eta;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG1510_phi;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG2210_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG2210_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG1510_eta;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG1510_phi;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta;   //!
   TBranch        *b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi;   //!
   TBranch        *b_trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_accept;   //!
   TBranch        *b_trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_prescale;   //!
   TBranch        *b_trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG20_eta;   //!
   TBranch        *b_trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG20_phi;   //!
   TBranch        *b_trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_accept;   //!
   TBranch        *b_trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_prescale;   //!
   TBranch        *b_trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltL1sL1SingleEG15_eta;   //!
   TBranch        *b_trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltL1sL1SingleEG15_phi;   //!
   TBranch        *b_trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept;   //!
   TBranch        *b_trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale;   //!
   TBranch        *b_trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG15ORSingleEG10_eta;   //!
   TBranch        *b_trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG15ORSingleEG10_phi;   //!
   TBranch        *b_trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter_phi;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_accept;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_prescale;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_eta;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_phi;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_eta;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_phi;   //!
   TBranch        *b_trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_accept;   //!
   TBranch        *b_trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_prescale;   //!
   TBranch        *b_trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltL1sL1SingleIsoEG22erOrSingleIsoEG20OrSingleEG25OrSingleEG20_eta;   //!
   TBranch        *b_trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltL1sL1SingleIsoEG22erOrSingleIsoEG20OrSingleEG25OrSingleEG20_phi;   //!
   TBranch        *b_trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltEle23WPLooseGsfTrackIsoFilterForEleStream_eta;   //!
   TBranch        *b_trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltEle23WPLooseGsfTrackIsoFilterForEleStream_phi;   //!

   IIHEAnalysis(TTree *tree=0);
   virtual ~IIHEAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef IIHEAnalysis_cxx
IIHEAnalysis::IIHEAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("outfile.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("outfile.root");
      }
      f->GetObject("IIHEAnalysis",tree);

   }
   Init(tree);
}

IIHEAnalysis::~IIHEAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t IIHEAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t IIHEAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void IIHEAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pv_x = 0;
   pv_y = 0;
   pv_z = 0;
   pv_isValid = 0;
   pv_normalizedChi2 = 0;
   pv_ndof = 0;
   pv_nTracks = 0;
   pv_totTrackSize = 0;
   sc_energy = 0;
   sc_eta = 0;
   sc_etacorr = 0;
   sc_theta = 0;
   sc_thetacorr = 0;
   sc_et = 0;
   sc_phi = 0;
   sc_px = 0;
   sc_py = 0;
   sc_pz = 0;
   sc_x = 0;
   sc_y = 0;
   sc_z = 0;
   ph_px = 0;
   ph_pt = 0;
   ph_eta = 0;
   ph_theta = 0;
   ph_phi = 0;
   ph_energy = 0;
   ph_mass = 0;
   ph_isPFlowPhoton = 0;
   ph_isStandardPhoton = 0;
   ph_hasConversionTracks = 0;
   ph_hasPixelSeed = 0;
   ph_isEB = 0;
   ph_isEE = 0;
   ph_isEBGap = 0;
   ph_isEBEtaGap = 0;
   ph_isEBPhiGap = 0;
   ph_isEEGap = 0;
   ph_isEERingGap = 0;
   ph_isEEDeeGap = 0;
   ph_isEBEEGap = 0;
   ph_hadronicOverEm = 0;
   ph_hadronicDepth1OverEm = 0;
   ph_hadronicDepth2OverEm = 0;
   ph_hadTowOverEm = 0;
   ph_hadTowDepth1OverEm = 0;
   ph_hadTowDepth2OverEm = 0;
   ph_e1x5 = 0;
   ph_e2x5 = 0;
   ph_e3x3 = 0;
   ph_e5x5 = 0;
   ph_maxEnergyXtal = 0;
   ph_sigmaEtaEta = 0;
   ph_sigmaIetaIeta = 0;
   ph_r1x5 = 0;
   ph_r2x5 = 0;
   ph_r9 = 0;
   ph_mipChi2 = 0;
   ph_mipTotEnergy = 0;
   ph_mipSlope = 0;
   ph_mipIntercept = 0;
   ph_mipNhitCone = 0;
   ph_mipIsHalo = 0;
   ph_ecalRecHitSumEtConeDR04 = 0;
   ph_hcalTowerSumEtConeDR04 = 0;
   ph_hcalDepth1TowerSumEtConeDR04 = 0;
   ph_hcalDepth2TowerSumEtConeDR04 = 0;
   ph_hcalTowerSumEtBcConeDR04 = 0;
   ph_hcalDepth1TowerSumEtBcConeDR04 = 0;
   ph_hcalDepth2TowerSumEtBcConeDR04 = 0;
   ph_trkSumPtSolidConeDR04 = 0;
   ph_trkSumPtHollowConeDR04 = 0;
   ph_nTrkSolidConeDR04 = 0;
   ph_nTrkHollowConeDR04 = 0;
   ph_ecalRecHitSumEtConeDR03 = 0;
   ph_hcalTowerSumEtConeDR03 = 0;
   ph_hcalDepth1TowerSumEtConeDR03 = 0;
   ph_hcalDepth2TowerSumEtConeDR03 = 0;
   ph_hcalTowerSumEtBcConeDR03 = 0;
   ph_hcalDepth1TowerSumEtBcConeDR03 = 0;
   ph_hcalDepth2TowerSumEtBcConeDR03 = 0;
   ph_trkSumPtSolidConeDR03 = 0;
   ph_trkSumPtHollowConeDR03 = 0;
   ph_nTrkSolidConeDR03 = 0;
   ph_nTrkHollowConeDR03 = 0;
   ph_chargedHadronIso = 0;
   ph_neutralHadronIso = 0;
   ph_photonIso = 0;
   ph_chargedHadronIsoWrongVtx = 0;
   ph_sumChargedParticlePt = 0;
   ph_sumNeutralHadronEtHighThreshold = 0;
   ph_sumPhotonEtHighThreshold = 0;
   ph_sumPUPt = 0;
   ph_nClusterOutsideMustache = 0;
   ph_etOutsideMustache = 0;
   ph_pfMVA = 0;
   ph_mc_bestDR = 0;
   ph_mc_index = 0;
   ph_mc_ERatio = 0;
   gsf_classification = 0;
   gsf_energy = 0;
   gsf_p = 0;
   gsf_pt = 0;
   gsf_scE1x5 = 0;
   gsf_scE5x5 = 0;
   gsf_scE2x5Max = 0;
   gsf_eta = 0;
   gsf_phi = 0;
   gsf_theta = 0;
   gsf_px = 0;
   gsf_py = 0;
   gsf_pz = 0;
   gsf_superClusterEta = 0;
   gsf_superClusterEnergy = 0;
   gsf_caloEnergy = 0;
   gsf_deltaEtaSuperClusterTrackAtVtx = 0;
   gsf_deltaPhiSuperClusterTrackAtVtx = 0;
   gsf_hadronicOverEm = 0;
   gsf_hcalDepth1OverEcal = 0;
   gsf_hcalDepth2OverEcal = 0;
   gsf_dr03TkSumPt = 0;
   gsf_dr03EcalRecHitSumEt = 0;
   gsf_dr03HcalDepth1TowerSumEt = 0;
   gsf_dr03HcalDepth2TowerSumEt = 0;
   gsf_charge = 0;
   gsf_sigmaIetaIeta = 0;
   gsf_ecaldrivenSeed = 0;
   gsf_trackerdrivenSeed = 0;
   gsf_isEB = 0;
   gsf_isEE = 0;
   gsf_deltaEtaSeedClusterTrackAtCalo = 0;
   gsf_deltaPhiSeedClusterTrackAtCalo = 0;
   gsf_ecalEnergy = 0;
   gsf_eSuperClusterOverP = 0;
   gsf_dxy = 0;
   gsf_dxy_beamSpot = 0;
   gsf_dxy_firstPVtx = 0;
   gsf_dxyError = 0;
   gsf_dz = 0;
   gsf_dz_beamSpot = 0;
   gsf_dz_firstPVtx = 0;
   gsf_dzError = 0;
   gsf_vz = 0;
   gsf_numberOfValidHits = 0;
   gsf_nLostInnerHits = 0;
   gsf_nLostOuterHits = 0;
   gsf_convFlags = 0;
   gsf_convDist = 0;
   gsf_convDcot = 0;
   gsf_convRadius = 0;
   gsf_fBrem = 0;
   gsf_e1x5 = 0;
   gsf_e2x5Max = 0;
   gsf_e5x5 = 0;
   gsf_r9 = 0;
   gsf_hitsinfo = 0;
   gsf_pixelMatch_dPhi1 = 0;
   gsf_pixelMatch_dPhi2 = 0;
   gsf_pixelMatch_dRz1 = 0;
   gsf_pixelMatch_dRz2 = 0;
   gsf_pixelMatch_subDetector1 = 0;
   gsf_pixelMatch_subDetector2 = 0;
   gsf_mc_bestDR = 0;
   gsf_mc_index = 0;
   gsf_mc_ERatio = 0;
   mu_gt_qoverp = 0;
   mu_gt_charge = 0;
   mu_gt_pt = 0;
   mu_gt_eta = 0;
   mu_gt_phi = 0;
   mu_gt_p = 0;
   mu_gt_px = 0;
   mu_gt_py = 0;
   mu_gt_pz = 0;
   mu_gt_theta = 0;
   mu_gt_lambda = 0;
   mu_gt_d0 = 0;
   mu_gt_dz = 0;
   mu_gt_dz_beamspot = 0;
   mu_gt_dz_firstPVtx = 0;
   mu_gt_dxy = 0;
   mu_gt_dxy_beamspot = 0;
   mu_gt_dxy_firstPVtx = 0;
   mu_gt_dsz = 0;
   mu_gt_vx = 0;
   mu_gt_vy = 0;
   mu_gt_vz = 0;
   mu_gt_qoverpError = 0;
   mu_gt_ptError = 0;
   mu_gt_thetaError = 0;
   mu_gt_lambdaError = 0;
   mu_gt_phiError = 0;
   mu_gt_dxyError = 0;
   mu_gt_d0Error = 0;
   mu_gt_dszError = 0;
   mu_gt_dzError = 0;
   mu_gt_etaError = 0;
   mu_gt_chi2 = 0;
   mu_gt_ndof = 0;
   mu_gt_normalizedChi2 = 0;
   mu_ot_qoverp = 0;
   mu_ot_charge = 0;
   mu_ot_pt = 0;
   mu_ot_eta = 0;
   mu_ot_phi = 0;
   mu_ot_p = 0;
   mu_ot_px = 0;
   mu_ot_py = 0;
   mu_ot_pz = 0;
   mu_ot_theta = 0;
   mu_ot_lambda = 0;
   mu_ot_d0 = 0;
   mu_ot_dz = 0;
   mu_ot_dz_beamspot = 0;
   mu_ot_dz_firstPVtx = 0;
   mu_ot_dxy = 0;
   mu_ot_dxy_beamspot = 0;
   mu_ot_dxy_firstPVtx = 0;
   mu_ot_dsz = 0;
   mu_ot_vx = 0;
   mu_ot_vy = 0;
   mu_ot_vz = 0;
   mu_ot_qoverpError = 0;
   mu_ot_ptError = 0;
   mu_ot_thetaError = 0;
   mu_ot_lambdaError = 0;
   mu_ot_phiError = 0;
   mu_ot_dxyError = 0;
   mu_ot_d0Error = 0;
   mu_ot_dszError = 0;
   mu_ot_dzError = 0;
   mu_ot_etaError = 0;
   mu_ot_chi2 = 0;
   mu_ot_ndof = 0;
   mu_ot_normalizedChi2 = 0;
   mu_it_qoverp = 0;
   mu_it_charge = 0;
   mu_it_pt = 0;
   mu_it_eta = 0;
   mu_it_phi = 0;
   mu_it_p = 0;
   mu_it_px = 0;
   mu_it_py = 0;
   mu_it_pz = 0;
   mu_it_theta = 0;
   mu_it_lambda = 0;
   mu_it_d0 = 0;
   mu_it_dz = 0;
   mu_it_dz_beamspot = 0;
   mu_it_dz_firstPVtx = 0;
   mu_it_dxy = 0;
   mu_it_dxy_beamspot = 0;
   mu_it_dxy_firstPVtx = 0;
   mu_it_dsz = 0;
   mu_it_vx = 0;
   mu_it_vy = 0;
   mu_it_vz = 0;
   mu_it_qoverpError = 0;
   mu_it_ptError = 0;
   mu_it_thetaError = 0;
   mu_it_lambdaError = 0;
   mu_it_phiError = 0;
   mu_it_dxyError = 0;
   mu_it_d0Error = 0;
   mu_it_dszError = 0;
   mu_it_dzError = 0;
   mu_it_etaError = 0;
   mu_it_chi2 = 0;
   mu_it_ndof = 0;
   mu_it_normalizedChi2 = 0;
   mu_isGlobalMuon = 0;
   mu_isStandAloneMuon = 0;
   mu_isTrackerMuon = 0;
   mu_isPFMuon = 0;
   mu_isPFIsolationValid = 0;
   mu_numberOfMatchedStations = 0;
   mu_numberOfValidPixelHits = 0;
   mu_numberOfValidTrackerHits = 0;
   mu_numberOfValidMuonHits = 0;
   mu_tevOptimized_charge = 0;
   mu_tevOptimized_pt = 0;
   mu_tevOptimized_eta = 0;
   mu_tevOptimized_phi = 0;
   mu_tevOptimized_theta = 0;
   mu_tevOptimized_px = 0;
   mu_tevOptimized_py = 0;
   mu_tevOptimized_pz = 0;
   mu_tevOptimized_d0 = 0;
   mu_tevOptimized_dz = 0;
   mu_tevOptimized_dz_beamSpot = 0;
   mu_tevOptimized_dz_firstPVtx = 0;
   mu_tevOptimized_dxy = 0;
   mu_tevOptimized_dxy_beamSpot = 0;
   mu_tevOptimized_dxy_firstPVtx = 0;
   mu_tevOptimized_ptError = 0;
   mu_tevOptimized_etaError = 0;
   mu_tevOptimized_phiError = 0;
   mu_tevOptimized_thetaError = 0;
   mu_tevOptimized_d0Error = 0;
   mu_tevOptimized_dzError = 0;
   mu_tevOptimized_dxyError = 0;
   mu_isolationR03_sumPt = 0;
   mu_isolationR03_trackerVetoPt = 0;
   mu_isolationR03_emEt = 0;
   mu_isolationR03_emVetoEt = 0;
   mu_isolationR03_hadEt = 0;
   mu_isolationR03_hadVetoEt = 0;
   mu_isolationR05_sumPt = 0;
   mu_isolationR05_trackerVetoPt = 0;
   mu_isolationR05_emEt = 0;
   mu_isolationR05_emVetoEt = 0;
   mu_isolationR05_hadEt = 0;
   mu_isolationR05_hadVetoEt = 0;
   mu_pfIsolationR03_sumChargedHadronPt = 0;
   mu_pfIsolationR03_sumChargedParticlePt = 0;
   mu_pfIsolationR03_sumPhotonEt = 0;
   mu_pfIsolationR03_sumNeutralHadronEtHighThreshold = 0;
   mu_pfIsolationR03_sumPhotonEtHighThreshold = 0;
   mu_pfIsolationR03_sumPUPt = 0;
   mu_pfIsolationR04_sumChargedHadronPt = 0;
   mu_pfIsolationR04_sumChargedParticlePt = 0;
   mu_pfIsolationR04_sumPhotonEt = 0;
   mu_pfIsolationR04_sumNeutralHadronEtHighThreshold = 0;
   mu_pfIsolationR04_sumPhotonEtHighThreshold = 0;
   mu_pfIsolationR04_sumPUPt = 0;
   mu_pfMeanDRIsoProfileR03_sumChargedHadronPt = 0;
   mu_pfMeanDRIsoProfileR03_sumChargedParticlePt = 0;
   mu_pfMeanDRIsoProfileR03_sumPhotonEt = 0;
   mu_pfMeanDRIsoProfileR03_sumNeutralHadronEtHighThreshold = 0;
   mu_pfMeanDRIsoProfileR03_sumPhotonEtHighThreshold = 0;
   mu_pfMeanDRIsoProfileR03_sumPUPt = 0;
   mu_pfMeanDRIsoProfileR04_sumChargedHadronPt = 0;
   mu_pfMeanDRIsoProfileR04_sumChargedParticlePt = 0;
   mu_pfMeanDRIsoProfileR04_sumPhotonEt = 0;
   mu_pfMeanDRIsoProfileR04_sumNeutralHadronEtHighThreshold = 0;
   mu_pfMeanDRIsoProfileR04_sumPhotonEtHighThreshold = 0;
   mu_pfMeanDRIsoProfileR04_sumPUPt = 0;
   mu_mc_bestDR = 0;
   mu_mc_index = 0;
   mu_mc_ERatio = 0;
   HEEP_eseffsixix = 0;
   HEEP_eseffsiyiy = 0;
   HEEP_eseffsirir = 0;
   HEEP_preshowerEnergy = 0;
   HEEP_e1x3 = 0;
   HEEP_eMax = 0;
   HEEP_e5x5 = 0;
   HEEP_e2x5Right = 0;
   HEEP_e2x5Left = 0;
   HEEP_e2x5Top = 0;
   HEEP_e2x5Bottom = 0;
   HEEP_eRight = 0;
   HEEP_eLeft = 0;
   HEEP_eTop = 0;
   HEEP_eBottom = 0;
   HEEP_basicClusterSeedTime = 0;
   EBHits_rawId = 0;
   EBHits_iRechit = 0;
   EBHits_energy = 0;
   EBHits_ieta = 0;
   EBHits_iphi = 0;
   EBHits_RecoFlag = 0;
   EBHits_kSaturated = 0;
   EBHits_kLeadingEdgeRecovered = 0;
   EBHits_kNeighboursRecovered = 0;
   EBHits_kWeird = 0;
   EEHits_rawId = 0;
   EEHits_iRechit = 0;
   EEHits_energy = 0;
   EEHits_ieta = 0;
   EEHits_iphi = 0;
   EEHits_RecoFlag = 0;
   EEHits_kSaturated = 0;
   EEHits_kLeadingEdgeRecovered = 0;
   EEHits_kNeighboursRecovered = 0;
   EEHits_kWeird = 0;
   HEEP_crystal_energy = 0;
   HEEP_crystal_eta = 0;
   HEEP_eshitsixix = 0;
   HEEP_eshitsiyiy = 0;
   HEEP_crystal_ietaorix = 0;
   HEEP_crystal_iphioriy = 0;
   HEEP_cutflow41_Et = 0;
   HEEP_cutflow41_Et_value = 0;
   HEEP_cutflow41_eta = 0;
   HEEP_cutflow41_eta_value = 0;
   HEEP_cutflow41_acceptance = 0;
   HEEP_cutflow41_EcalDriven = 0;
   HEEP_cutflow41_EcalDriven_value = 0;
   HEEP_cutflow41_dEtaIn = 0;
   HEEP_cutflow41_dEtaIn_value = 0;
   HEEP_cutflow41_dPhiIn = 0;
   HEEP_cutflow41_dPhiIn_value = 0;
   HEEP_cutflow41_HOverE = 0;
   HEEP_cutflow41_HOverE_value = 0;
   HEEP_cutflow41_SigmaIetaIeta = 0;
   HEEP_cutflow41_SigmaIetaIeta_value = 0;
   HEEP_cutflow41_E1x5OverE5x5 = 0;
   HEEP_cutflow41_E1x5OverE5x5_value = 0;
   HEEP_cutflow41_E2x5OverE5x5 = 0;
   HEEP_cutflow41_E2x5OverE5x5_value = 0;
   HEEP_cutflow41_missingHits = 0;
   HEEP_cutflow41_missingHits_value = 0;
   HEEP_cutflow41_dxyFirstPV = 0;
   HEEP_cutflow41_dxyFirstPV_value = 0;
   HEEP_cutflow41_ID = 0;
   HEEP_cutflow41_isolEMHadDepth1 = 0;
   HEEP_cutflow41_isolEMHadDepth1_value = 0;
   HEEP_cutflow41_IsolPtTrks = 0;
   HEEP_cutflow41_IsolPtTrks_value = 0;
   HEEP_cutflow41_isolation = 0;
   HEEP_cutflow41_total = 0;
   HEEP_cutflow50_25ns_Et = 0;
   HEEP_cutflow50_25ns_Et_value = 0;
   HEEP_cutflow50_25ns_eta = 0;
   HEEP_cutflow50_25ns_eta_value = 0;
   HEEP_cutflow50_25ns_acceptance = 0;
   HEEP_cutflow50_25ns_EcalDriven = 0;
   HEEP_cutflow50_25ns_EcalDriven_value = 0;
   HEEP_cutflow50_25ns_dEtaIn = 0;
   HEEP_cutflow50_25ns_dEtaIn_value = 0;
   HEEP_cutflow50_25ns_dPhiIn = 0;
   HEEP_cutflow50_25ns_dPhiIn_value = 0;
   HEEP_cutflow50_25ns_HOverE = 0;
   HEEP_cutflow50_25ns_HOverE_value = 0;
   HEEP_cutflow50_25ns_SigmaIetaIeta = 0;
   HEEP_cutflow50_25ns_SigmaIetaIeta_value = 0;
   HEEP_cutflow50_25ns_E1x5OverE5x5 = 0;
   HEEP_cutflow50_25ns_E1x5OverE5x5_value = 0;
   HEEP_cutflow50_25ns_E2x5OverE5x5 = 0;
   HEEP_cutflow50_25ns_E2x5OverE5x5_value = 0;
   HEEP_cutflow50_25ns_missingHits = 0;
   HEEP_cutflow50_25ns_missingHits_value = 0;
   HEEP_cutflow50_25ns_dxyFirstPV = 0;
   HEEP_cutflow50_25ns_dxyFirstPV_value = 0;
   HEEP_cutflow50_25ns_ID = 0;
   HEEP_cutflow50_25ns_isolEMHadDepth1 = 0;
   HEEP_cutflow50_25ns_isolEMHadDepth1_value = 0;
   HEEP_cutflow50_25ns_IsolPtTrks = 0;
   HEEP_cutflow50_25ns_IsolPtTrks_value = 0;
   HEEP_cutflow50_25ns_isolation = 0;
   HEEP_cutflow50_25ns_total = 0;
   HEEP_cutflow50_50ns_Et = 0;
   HEEP_cutflow50_50ns_Et_value = 0;
   HEEP_cutflow50_50ns_eta = 0;
   HEEP_cutflow50_50ns_eta_value = 0;
   HEEP_cutflow50_50ns_acceptance = 0;
   HEEP_cutflow50_50ns_EcalDriven = 0;
   HEEP_cutflow50_50ns_EcalDriven_value = 0;
   HEEP_cutflow50_50ns_dEtaIn = 0;
   HEEP_cutflow50_50ns_dEtaIn_value = 0;
   HEEP_cutflow50_50ns_dPhiIn = 0;
   HEEP_cutflow50_50ns_dPhiIn_value = 0;
   HEEP_cutflow50_50ns_HOverE = 0;
   HEEP_cutflow50_50ns_HOverE_value = 0;
   HEEP_cutflow50_50ns_SigmaIetaIeta = 0;
   HEEP_cutflow50_50ns_SigmaIetaIeta_value = 0;
   HEEP_cutflow50_50ns_E1x5OverE5x5 = 0;
   HEEP_cutflow50_50ns_E1x5OverE5x5_value = 0;
   HEEP_cutflow50_50ns_E2x5OverE5x5 = 0;
   HEEP_cutflow50_50ns_E2x5OverE5x5_value = 0;
   HEEP_cutflow50_50ns_missingHits = 0;
   HEEP_cutflow50_50ns_missingHits_value = 0;
   HEEP_cutflow50_50ns_dxyFirstPV = 0;
   HEEP_cutflow50_50ns_dxyFirstPV_value = 0;
   HEEP_cutflow50_50ns_ID = 0;
   HEEP_cutflow50_50ns_isolEMHadDepth1 = 0;
   HEEP_cutflow50_50ns_isolEMHadDepth1_value = 0;
   HEEP_cutflow50_50ns_IsolPtTrks = 0;
   HEEP_cutflow50_50ns_IsolPtTrks_value = 0;
   HEEP_cutflow50_50ns_isolation = 0;
   HEEP_cutflow50_50ns_total = 0;
   HEEP_cutflow50_Et = 0;
   HEEP_cutflow50_Et_value = 0;
   HEEP_cutflow50_eta = 0;
   HEEP_cutflow50_eta_value = 0;
   HEEP_cutflow50_acceptance = 0;
   HEEP_cutflow50_EcalDriven = 0;
   HEEP_cutflow50_EcalDriven_value = 0;
   HEEP_cutflow50_dEtaIn = 0;
   HEEP_cutflow50_dEtaIn_value = 0;
   HEEP_cutflow50_dPhiIn = 0;
   HEEP_cutflow50_dPhiIn_value = 0;
   HEEP_cutflow50_HOverE = 0;
   HEEP_cutflow50_HOverE_value = 0;
   HEEP_cutflow50_SigmaIetaIeta = 0;
   HEEP_cutflow50_SigmaIetaIeta_value = 0;
   HEEP_cutflow50_E1x5OverE5x5 = 0;
   HEEP_cutflow50_E1x5OverE5x5_value = 0;
   HEEP_cutflow50_E2x5OverE5x5 = 0;
   HEEP_cutflow50_E2x5OverE5x5_value = 0;
   HEEP_cutflow50_missingHits = 0;
   HEEP_cutflow50_missingHits_value = 0;
   HEEP_cutflow50_dxyFirstPV = 0;
   HEEP_cutflow50_dxyFirstPV_value = 0;
   HEEP_cutflow50_ID = 0;
   HEEP_cutflow50_isolEMHadDepth1 = 0;
   HEEP_cutflow50_isolEMHadDepth1_value = 0;
   HEEP_cutflow50_IsolPtTrks = 0;
   HEEP_cutflow50_IsolPtTrks_value = 0;
   HEEP_cutflow50_isolation = 0;
   HEEP_cutflow50_total = 0;
   HEEP_cutflow51_Et = 0;
   HEEP_cutflow51_Et_value = 0;
   HEEP_cutflow51_eta = 0;
   HEEP_cutflow51_eta_value = 0;
   HEEP_cutflow51_acceptance = 0;
   HEEP_cutflow51_EcalDriven = 0;
   HEEP_cutflow51_EcalDriven_value = 0;
   HEEP_cutflow51_dEtaIn = 0;
   HEEP_cutflow51_dEtaIn_value = 0;
   HEEP_cutflow51_dPhiIn = 0;
   HEEP_cutflow51_dPhiIn_value = 0;
   HEEP_cutflow51_HOverE = 0;
   HEEP_cutflow51_HOverE_value = 0;
   HEEP_cutflow51_SigmaIetaIeta = 0;
   HEEP_cutflow51_SigmaIetaIeta_value = 0;
   HEEP_cutflow51_E1x5OverE5x5 = 0;
   HEEP_cutflow51_E1x5OverE5x5_value = 0;
   HEEP_cutflow51_E2x5OverE5x5 = 0;
   HEEP_cutflow51_E2x5OverE5x5_value = 0;
   HEEP_cutflow51_missingHits = 0;
   HEEP_cutflow51_missingHits_value = 0;
   HEEP_cutflow51_dxyFirstPV = 0;
   HEEP_cutflow51_dxyFirstPV_value = 0;
   HEEP_cutflow51_ID = 0;
   HEEP_cutflow51_isolEMHadDepth1 = 0;
   HEEP_cutflow51_isolEMHadDepth1_value = 0;
   HEEP_cutflow51_IsolPtTrks = 0;
   HEEP_cutflow51_IsolPtTrks_value = 0;
   HEEP_cutflow51_isolation = 0;
   HEEP_cutflow51_total = 0;
   HEEP_cutflow60_Et = 0;
   HEEP_cutflow60_Et_value = 0;
   HEEP_cutflow60_eta = 0;
   HEEP_cutflow60_eta_value = 0;
   HEEP_cutflow60_acceptance = 0;
   HEEP_cutflow60_EcalDriven = 0;
   HEEP_cutflow60_EcalDriven_value = 0;
   HEEP_cutflow60_dEtaIn = 0;
   HEEP_cutflow60_dEtaIn_value = 0;
   HEEP_cutflow60_dPhiIn = 0;
   HEEP_cutflow60_dPhiIn_value = 0;
   HEEP_cutflow60_HOverE = 0;
   HEEP_cutflow60_HOverE_value = 0;
   HEEP_cutflow60_SigmaIetaIeta = 0;
   HEEP_cutflow60_SigmaIetaIeta_value = 0;
   HEEP_cutflow60_E1x5OverE5x5 = 0;
   HEEP_cutflow60_E1x5OverE5x5_value = 0;
   HEEP_cutflow60_E2x5OverE5x5 = 0;
   HEEP_cutflow60_E2x5OverE5x5_value = 0;
   HEEP_cutflow60_missingHits = 0;
   HEEP_cutflow60_missingHits_value = 0;
   HEEP_cutflow60_dxyFirstPV = 0;
   HEEP_cutflow60_dxyFirstPV_value = 0;
   HEEP_cutflow60_ID = 0;
   HEEP_cutflow60_isolEMHadDepth1 = 0;
   HEEP_cutflow60_isolEMHadDepth1_value = 0;
   HEEP_cutflow60_IsolPtTrks = 0;
   HEEP_cutflow60_IsolPtTrks_value = 0;
   HEEP_cutflow60_isolation = 0;
   HEEP_cutflow60_total = 0;
   Zee_mass = 0;
   Zee_mass_HEEP = 0;
   Zee_i1 = 0;
   Zee_i2 = 0;
   Zmm_mass = 0;
   Zmm_i1 = 0;
   Zmm_i2 = 0;
   Zem_mass = 0;
   Zem_mass_HEEP = 0;
   Zem_i1 = 0;
   Zem_i2 = 0;
   trk_pt = 0;
   trk_eta = 0;
   trk_phi = 0;
   trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltL1sL1DoubleEG2210_eta = 0;
   trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltL1sL1DoubleEG2210_phi = 0;
   trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg1TrackIsoFilter_eta = 0;
   trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg1TrackIsoFilter_phi = 0;
   trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg2TrackIsoFilter_eta = 0;
   trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg2TrackIsoFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLPixelMatchFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLPixelMatchFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLNewPixelMatchFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLNewPixelMatchFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEG33EtUnseededFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEG33EtUnseededFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLNewPixelMatchUnseededFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLNewPixelMatchUnseededFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLPixelMatchFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLPixelMatchFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEG33EtUnseededFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEG33EtUnseededFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_eta = 0;
   trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_phi = 0;
   trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_eta = 0;
   trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_phi = 0;
   trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltSingleEle22WPLooseGsfTrackIsoFilter_eta = 0;
   trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltSingleEle22WPLooseGsfTrackIsoFilter_phi = 0;
   trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_eta = 0;
   trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_phi = 0;
   trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltSingleEle22WPTightGsfTrackIsoFilter_eta = 0;
   trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltSingleEle22WPTightGsfTrackIsoFilter_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltL1sL1SingleEG25OrSingleIsoEG30er_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltL1sL1SingleEG25OrSingleIsoEG30er_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtFilter_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtFilter_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8ClusterShapeFilter_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8ClusterShapeFilter_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HEFilter_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HEFilter_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EcalIsoFilter_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EcalIsoFilter_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HcalIsoFilter_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HcalIsoFilter_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchFilter_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchFilter_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8OneOEMinusOneOPFilter_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8OneOEMinusOneOPFilter_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DetaFilter_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DetaFilter_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DphiFilter_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DphiFilter_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8TrackIsoFilter_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8TrackIsoFilter_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtUnseededFilter_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtUnseededFilter_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchUnseededFilter_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchUnseededFilter_phi = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8Mass55Filter_eta = 0;
   trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8Mass55Filter_phi = 0;
   trig_HLT_Ele23_WPLoose_Gsf_v2_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_eta = 0;
   trig_HLT_Ele23_WPLoose_Gsf_v2_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_phi = 0;
   trig_HLT_Ele23_WPLoose_Gsf_v2_hltEle23WPLooseGsfTrackIsoFilter_eta = 0;
   trig_HLT_Ele23_WPLoose_Gsf_v2_hltEle23WPLooseGsfTrackIsoFilter_phi = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_eta = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_phi = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltEle23WPLooseGsfTrackIsoFilter_eta = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltEle23WPLooseGsfTrackIsoFilter_phi = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltJetFilterEle23WPLoose_eta = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltJetFilterEle23WPLoose_phi = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltHCand80NoEle23WPLoose_eta = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltHCand80NoEle23WPLoose_phi = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMhtFilter0_eta = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMhtFilter0_phi = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMetFilter0_eta = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMetFilter0_phi = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand80NoEle23WPLooseMET_eta = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand80NoEle23WPLooseMET_phi = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand70NoEle23WPLooseMHTIDTight_eta = 0;
   trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand70NoEle23WPLooseMHTIDTight_phi = 0;
   trig_HLT_Ele27_WPLoose_Gsf_v1_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_eta = 0;
   trig_HLT_Ele27_WPLoose_Gsf_v1_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_phi = 0;
   trig_HLT_Ele27_WPLoose_Gsf_v1_hltEle27noerWPLooseGsfTrackIsoFilter_eta = 0;
   trig_HLT_Ele27_WPLoose_Gsf_v1_hltEle27noerWPLooseGsfTrackIsoFilter_phi = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_eta = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_phi = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltEle27noerWPLooseGsfTrackIsoFilter_eta = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltEle27noerWPLooseGsfTrackIsoFilter_phi = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltJetFilterEle27WPLoose_eta = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltJetFilterEle27WPLoose_phi = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltHCand80NoEle27WPLoose_eta = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltHCand80NoEle27WPLoose_phi = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMhtFilter0_eta = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMhtFilter0_phi = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMetFilter0_eta = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMetFilter0_phi = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand80NoEle27WPLooseMET_eta = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand80NoEle27WPLooseMET_phi = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand70NoEle27WPLooseMHTIDTight_eta = 0;
   trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand70NoEle27WPLooseMHTIDTight_phi = 0;
   trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_eta = 0;
   trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_phi = 0;
   trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltEle27WPLooseGsfTrackIsoFilter_eta = 0;
   trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltEle27WPLooseGsfTrackIsoFilter_phi = 0;
   trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_eta = 0;
   trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_phi = 0;
   trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltEle27WPTightGsfTrackIsoFilter_eta = 0;
   trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltEle27WPTightGsfTrackIsoFilter_phi = 0;
   trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG30erOrSingleEG40_eta = 0;
   trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG30erOrSingleEG40_phi = 0;
   trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltEle32WPTightGsfTrackIsoFilter_eta = 0;
   trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltEle32WPTightGsfTrackIsoFilter_phi = 0;
   trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_eta = 0;
   trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_phi = 0;
   trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter_eta = 0;
   trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG2210_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG2210_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_phi = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG1510_eta = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG1510_phi = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter_eta = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG2210_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG2210_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta = 0;
   trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG1510_eta = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG1510_phi = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta = 0;
   trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi = 0;
   trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG20_eta = 0;
   trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG20_phi = 0;
   trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter_eta = 0;
   trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter_phi = 0;
   trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltL1sL1SingleEG15_eta = 0;
   trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltL1sL1SingleEG15_phi = 0;
   trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter_eta = 0;
   trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter_phi = 0;
   trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG15ORSingleEG10_eta = 0;
   trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG15ORSingleEG10_phi = 0;
   trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter_eta = 0;
   trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter_phi = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_eta = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_phi = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_eta = 0;
   trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_phi = 0;
   trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltL1sL1SingleIsoEG22erOrSingleIsoEG20OrSingleEG25OrSingleEG20_eta = 0;
   trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltL1sL1SingleIsoEG22erOrSingleIsoEG20OrSingleEG25OrSingleEG20_phi = 0;
   trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltEle23WPLooseGsfTrackIsoFilterForEleStream_eta = 0;
   trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltEle23WPLooseGsfTrackIsoFilterForEleStream_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev_event", &ev_event, &b_ev_event);
   fChain->SetBranchAddress("ev_run", &ev_run, &b_ev_run);
   fChain->SetBranchAddress("ev_luminosityBlock", &ev_luminosityBlock, &b_ev_luminosityBlock);
   fChain->SetBranchAddress("ev_time", &ev_time, &b_ev_time);
   fChain->SetBranchAddress("ev_time_unixTime", &ev_time_unixTime, &b_ev_time_unixTime);
   fChain->SetBranchAddress("ev_time_microsecondOffset", &ev_time_microsecondOffset, &b_ev_time_microsecondOffset);
   fChain->SetBranchAddress("ev_fixedGridRhoAll", &ev_fixedGridRhoAll, &b_ev_fixedGridRhoAll);
   fChain->SetBranchAddress("ev_rho_kt6PFJetsForIsolation", &ev_rho_kt6PFJetsForIsolation, &b_ev_rho_kt6PFJetsForIsolation);
   fChain->SetBranchAddress("pv_n", &pv_n, &b_pv_n);
   fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_isValid", &pv_isValid, &b_pv_isValid);
   fChain->SetBranchAddress("pv_normalizedChi2", &pv_normalizedChi2, &b_pv_normalizedChi2);
   fChain->SetBranchAddress("pv_ndof", &pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("pv_nTracks", &pv_nTracks, &b_pv_nTracks);
   fChain->SetBranchAddress("pv_totTrackSize", &pv_totTrackSize, &b_pv_totTrackSize);
   fChain->SetBranchAddress("sc_n", &sc_n, &b_sc_n);
   fChain->SetBranchAddress("sc_energy", &sc_energy, &b_sc_energy);
   fChain->SetBranchAddress("sc_eta", &sc_eta, &b_sc_eta);
   fChain->SetBranchAddress("sc_etacorr", &sc_etacorr, &b_sc_etacorr);
   fChain->SetBranchAddress("sc_theta", &sc_theta, &b_sc_theta);
   fChain->SetBranchAddress("sc_thetacorr", &sc_thetacorr, &b_sc_thetacorr);
   fChain->SetBranchAddress("sc_et", &sc_et, &b_sc_et);
   fChain->SetBranchAddress("sc_phi", &sc_phi, &b_sc_phi);
   fChain->SetBranchAddress("sc_px", &sc_px, &b_sc_px);
   fChain->SetBranchAddress("sc_py", &sc_py, &b_sc_py);
   fChain->SetBranchAddress("sc_pz", &sc_pz, &b_sc_pz);
   fChain->SetBranchAddress("sc_x", &sc_x, &b_sc_x);
   fChain->SetBranchAddress("sc_y", &sc_y, &b_sc_y);
   fChain->SetBranchAddress("sc_z", &sc_z, &b_sc_z);
   fChain->SetBranchAddress("ph_n", &ph_n, &b_ph_n);
   fChain->SetBranchAddress("ph_px", &ph_px, &b_ph_px);
   fChain->SetBranchAddress("ph_pt", &ph_pt, &b_ph_pt);
   fChain->SetBranchAddress("ph_eta", &ph_eta, &b_ph_eta);
   fChain->SetBranchAddress("ph_theta", &ph_theta, &b_ph_theta);
   fChain->SetBranchAddress("ph_phi", &ph_phi, &b_ph_phi);
   fChain->SetBranchAddress("ph_energy", &ph_energy, &b_ph_energy);
   fChain->SetBranchAddress("ph_mass", &ph_mass, &b_ph_mass);
   fChain->SetBranchAddress("ph_isPFlowPhoton", &ph_isPFlowPhoton, &b_ph_isPFlowPhoton);
   fChain->SetBranchAddress("ph_isStandardPhoton", &ph_isStandardPhoton, &b_ph_isStandardPhoton);
   fChain->SetBranchAddress("ph_hasConversionTracks", &ph_hasConversionTracks, &b_ph_hasConversionTracks);
   fChain->SetBranchAddress("ph_hasPixelSeed", &ph_hasPixelSeed, &b_ph_hasPixelSeed);
   fChain->SetBranchAddress("ph_isEB", &ph_isEB, &b_ph_isEB);
   fChain->SetBranchAddress("ph_isEE", &ph_isEE, &b_ph_isEE);
   fChain->SetBranchAddress("ph_isEBGap", &ph_isEBGap, &b_ph_isEBGap);
   fChain->SetBranchAddress("ph_isEBEtaGap", &ph_isEBEtaGap, &b_ph_isEBEtaGap);
   fChain->SetBranchAddress("ph_isEBPhiGap", &ph_isEBPhiGap, &b_ph_isEBPhiGap);
   fChain->SetBranchAddress("ph_isEEGap", &ph_isEEGap, &b_ph_isEEGap);
   fChain->SetBranchAddress("ph_isEERingGap", &ph_isEERingGap, &b_ph_isEERingGap);
   fChain->SetBranchAddress("ph_isEEDeeGap", &ph_isEEDeeGap, &b_ph_isEEDeeGap);
   fChain->SetBranchAddress("ph_isEBEEGap", &ph_isEBEEGap, &b_ph_isEBEEGap);
   fChain->SetBranchAddress("ph_hadronicOverEm", &ph_hadronicOverEm, &b_ph_hadronicOverEm);
   fChain->SetBranchAddress("ph_hadronicDepth1OverEm", &ph_hadronicDepth1OverEm, &b_ph_hadronicDepth1OverEm);
   fChain->SetBranchAddress("ph_hadronicDepth2OverEm", &ph_hadronicDepth2OverEm, &b_ph_hadronicDepth2OverEm);
   fChain->SetBranchAddress("ph_hadTowOverEm", &ph_hadTowOverEm, &b_ph_hadTowOverEm);
   fChain->SetBranchAddress("ph_hadTowDepth1OverEm", &ph_hadTowDepth1OverEm, &b_ph_hadTowDepth1OverEm);
   fChain->SetBranchAddress("ph_hadTowDepth2OverEm", &ph_hadTowDepth2OverEm, &b_ph_hadTowDepth2OverEm);
   fChain->SetBranchAddress("ph_e1x5", &ph_e1x5, &b_ph_e1x5);
   fChain->SetBranchAddress("ph_e2x5", &ph_e2x5, &b_ph_e2x5);
   fChain->SetBranchAddress("ph_e3x3", &ph_e3x3, &b_ph_e3x3);
   fChain->SetBranchAddress("ph_e5x5", &ph_e5x5, &b_ph_e5x5);
   fChain->SetBranchAddress("ph_maxEnergyXtal", &ph_maxEnergyXtal, &b_ph_maxEnergyXtal);
   fChain->SetBranchAddress("ph_sigmaEtaEta", &ph_sigmaEtaEta, &b_ph_sigmaEtaEta);
   fChain->SetBranchAddress("ph_sigmaIetaIeta", &ph_sigmaIetaIeta, &b_ph_sigmaIetaIeta);
   fChain->SetBranchAddress("ph_r1x5", &ph_r1x5, &b_ph_r1x5);
   fChain->SetBranchAddress("ph_r2x5", &ph_r2x5, &b_ph_r2x5);
   fChain->SetBranchAddress("ph_r9", &ph_r9, &b_ph_r9);
   fChain->SetBranchAddress("ph_mipChi2", &ph_mipChi2, &b_ph_mipChi2);
   fChain->SetBranchAddress("ph_mipTotEnergy", &ph_mipTotEnergy, &b_ph_mipTotEnergy);
   fChain->SetBranchAddress("ph_mipSlope", &ph_mipSlope, &b_ph_mipSlope);
   fChain->SetBranchAddress("ph_mipIntercept", &ph_mipIntercept, &b_ph_mipIntercept);
   fChain->SetBranchAddress("ph_mipNhitCone", &ph_mipNhitCone, &b_ph_mipNhitCone);
   fChain->SetBranchAddress("ph_mipIsHalo", &ph_mipIsHalo, &b_ph_mipIsHalo);
   fChain->SetBranchAddress("ph_ecalRecHitSumEtConeDR04", &ph_ecalRecHitSumEtConeDR04, &b_ph_ecalRecHitSumEtConeDR04);
   fChain->SetBranchAddress("ph_hcalTowerSumEtConeDR04", &ph_hcalTowerSumEtConeDR04, &b_ph_hcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("ph_hcalDepth1TowerSumEtConeDR04", &ph_hcalDepth1TowerSumEtConeDR04, &b_ph_hcalDepth1TowerSumEtConeDR04);
   fChain->SetBranchAddress("ph_hcalDepth2TowerSumEtConeDR04", &ph_hcalDepth2TowerSumEtConeDR04, &b_ph_hcalDepth2TowerSumEtConeDR04);
   fChain->SetBranchAddress("ph_hcalTowerSumEtBcConeDR04", &ph_hcalTowerSumEtBcConeDR04, &b_ph_hcalTowerSumEtBcConeDR04);
   fChain->SetBranchAddress("ph_hcalDepth1TowerSumEtBcConeDR04", &ph_hcalDepth1TowerSumEtBcConeDR04, &b_ph_hcalDepth1TowerSumEtBcConeDR04);
   fChain->SetBranchAddress("ph_hcalDepth2TowerSumEtBcConeDR04", &ph_hcalDepth2TowerSumEtBcConeDR04, &b_ph_hcalDepth2TowerSumEtBcConeDR04);
   fChain->SetBranchAddress("ph_trkSumPtSolidConeDR04", &ph_trkSumPtSolidConeDR04, &b_ph_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("ph_trkSumPtHollowConeDR04", &ph_trkSumPtHollowConeDR04, &b_ph_trkSumPtHollowConeDR04);
   fChain->SetBranchAddress("ph_nTrkSolidConeDR04", &ph_nTrkSolidConeDR04, &b_ph_nTrkSolidConeDR04);
   fChain->SetBranchAddress("ph_nTrkHollowConeDR04", &ph_nTrkHollowConeDR04, &b_ph_nTrkHollowConeDR04);
   fChain->SetBranchAddress("ph_ecalRecHitSumEtConeDR03", &ph_ecalRecHitSumEtConeDR03, &b_ph_ecalRecHitSumEtConeDR03);
   fChain->SetBranchAddress("ph_hcalTowerSumEtConeDR03", &ph_hcalTowerSumEtConeDR03, &b_ph_hcalTowerSumEtConeDR03);
   fChain->SetBranchAddress("ph_hcalDepth1TowerSumEtConeDR03", &ph_hcalDepth1TowerSumEtConeDR03, &b_ph_hcalDepth1TowerSumEtConeDR03);
   fChain->SetBranchAddress("ph_hcalDepth2TowerSumEtConeDR03", &ph_hcalDepth2TowerSumEtConeDR03, &b_ph_hcalDepth2TowerSumEtConeDR03);
   fChain->SetBranchAddress("ph_hcalTowerSumEtBcConeDR03", &ph_hcalTowerSumEtBcConeDR03, &b_ph_hcalTowerSumEtBcConeDR03);
   fChain->SetBranchAddress("ph_hcalDepth1TowerSumEtBcConeDR03", &ph_hcalDepth1TowerSumEtBcConeDR03, &b_ph_hcalDepth1TowerSumEtBcConeDR03);
   fChain->SetBranchAddress("ph_hcalDepth2TowerSumEtBcConeDR03", &ph_hcalDepth2TowerSumEtBcConeDR03, &b_ph_hcalDepth2TowerSumEtBcConeDR03);
   fChain->SetBranchAddress("ph_trkSumPtSolidConeDR03", &ph_trkSumPtSolidConeDR03, &b_ph_trkSumPtSolidConeDR03);
   fChain->SetBranchAddress("ph_trkSumPtHollowConeDR03", &ph_trkSumPtHollowConeDR03, &b_ph_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("ph_nTrkSolidConeDR03", &ph_nTrkSolidConeDR03, &b_ph_nTrkSolidConeDR03);
   fChain->SetBranchAddress("ph_nTrkHollowConeDR03", &ph_nTrkHollowConeDR03, &b_ph_nTrkHollowConeDR03);
   fChain->SetBranchAddress("ph_chargedHadronIso", &ph_chargedHadronIso, &b_ph_chargedHadronIso);
   fChain->SetBranchAddress("ph_neutralHadronIso", &ph_neutralHadronIso, &b_ph_neutralHadronIso);
   fChain->SetBranchAddress("ph_photonIso", &ph_photonIso, &b_ph_photonIso);
   fChain->SetBranchAddress("ph_chargedHadronIsoWrongVtx", &ph_chargedHadronIsoWrongVtx, &b_ph_chargedHadronIsoWrongVtx);
   fChain->SetBranchAddress("ph_sumChargedParticlePt", &ph_sumChargedParticlePt, &b_ph_sumChargedParticlePt);
   fChain->SetBranchAddress("ph_sumNeutralHadronEtHighThreshold", &ph_sumNeutralHadronEtHighThreshold, &b_ph_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("ph_sumPhotonEtHighThreshold", &ph_sumPhotonEtHighThreshold, &b_ph_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("ph_sumPUPt", &ph_sumPUPt, &b_ph_sumPUPt);
   fChain->SetBranchAddress("ph_nClusterOutsideMustache", &ph_nClusterOutsideMustache, &b_ph_nClusterOutsideMustache);
   fChain->SetBranchAddress("ph_etOutsideMustache", &ph_etOutsideMustache, &b_ph_etOutsideMustache);
   fChain->SetBranchAddress("ph_pfMVA", &ph_pfMVA, &b_ph_pfMVA);
   fChain->SetBranchAddress("ph_mc_bestDR", &ph_mc_bestDR, &b_ph_mc_bestDR);
   fChain->SetBranchAddress("ph_mc_index", &ph_mc_index, &b_ph_mc_index);
   fChain->SetBranchAddress("ph_mc_ERatio", &ph_mc_ERatio, &b_ph_mc_ERatio);
   fChain->SetBranchAddress("gsf_n", &gsf_n, &b_gsf_n);
   fChain->SetBranchAddress("gsf_classification", &gsf_classification, &b_gsf_classification);
   fChain->SetBranchAddress("gsf_energy", &gsf_energy, &b_gsf_energy);
   fChain->SetBranchAddress("gsf_p", &gsf_p, &b_gsf_p);
   fChain->SetBranchAddress("gsf_pt", &gsf_pt, &b_gsf_pt);
   fChain->SetBranchAddress("gsf_scE1x5", &gsf_scE1x5, &b_gsf_scE1x5);
   fChain->SetBranchAddress("gsf_scE5x5", &gsf_scE5x5, &b_gsf_scE5x5);
   fChain->SetBranchAddress("gsf_scE2x5Max", &gsf_scE2x5Max, &b_gsf_scE2x5Max);
   fChain->SetBranchAddress("gsf_eta", &gsf_eta, &b_gsf_eta);
   fChain->SetBranchAddress("gsf_phi", &gsf_phi, &b_gsf_phi);
   fChain->SetBranchAddress("gsf_theta", &gsf_theta, &b_gsf_theta);
   fChain->SetBranchAddress("gsf_px", &gsf_px, &b_gsf_px);
   fChain->SetBranchAddress("gsf_py", &gsf_py, &b_gsf_py);
   fChain->SetBranchAddress("gsf_pz", &gsf_pz, &b_gsf_pz);
   fChain->SetBranchAddress("gsf_superClusterEta", &gsf_superClusterEta, &b_gsf_superClusterEta);
   fChain->SetBranchAddress("gsf_superClusterEnergy", &gsf_superClusterEnergy, &b_gsf_superClusterEnergy);
   fChain->SetBranchAddress("gsf_caloEnergy", &gsf_caloEnergy, &b_gsf_caloEnergy);
   fChain->SetBranchAddress("gsf_deltaEtaSuperClusterTrackAtVtx", &gsf_deltaEtaSuperClusterTrackAtVtx, &b_gsf_deltaEtaSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("gsf_deltaPhiSuperClusterTrackAtVtx", &gsf_deltaPhiSuperClusterTrackAtVtx, &b_gsf_deltaPhiSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("gsf_hadronicOverEm", &gsf_hadronicOverEm, &b_gsf_hadronicOverEm);
   fChain->SetBranchAddress("gsf_hcalDepth1OverEcal", &gsf_hcalDepth1OverEcal, &b_gsf_hcalDepth1OverEcal);
   fChain->SetBranchAddress("gsf_hcalDepth2OverEcal", &gsf_hcalDepth2OverEcal, &b_gsf_hcalDepth2OverEcal);
   fChain->SetBranchAddress("gsf_dr03TkSumPt", &gsf_dr03TkSumPt, &b_gsf_dr03TkSumPt);
   fChain->SetBranchAddress("gsf_dr03EcalRecHitSumEt", &gsf_dr03EcalRecHitSumEt, &b_gsf_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("gsf_dr03HcalDepth1TowerSumEt", &gsf_dr03HcalDepth1TowerSumEt, &b_gsf_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("gsf_dr03HcalDepth2TowerSumEt", &gsf_dr03HcalDepth2TowerSumEt, &b_gsf_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("gsf_charge", &gsf_charge, &b_gsf_charge);
   fChain->SetBranchAddress("gsf_sigmaIetaIeta", &gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
   fChain->SetBranchAddress("gsf_ecaldrivenSeed", &gsf_ecaldrivenSeed, &b_gsf_ecaldrivenSeed);
   fChain->SetBranchAddress("gsf_trackerdrivenSeed", &gsf_trackerdrivenSeed, &b_gsf_trackerdrivenSeed);
   fChain->SetBranchAddress("gsf_isEB", &gsf_isEB, &b_gsf_isEB);
   fChain->SetBranchAddress("gsf_isEE", &gsf_isEE, &b_gsf_isEE);
   fChain->SetBranchAddress("gsf_deltaEtaSeedClusterTrackAtCalo", &gsf_deltaEtaSeedClusterTrackAtCalo, &b_gsf_deltaEtaSeedClusterTrackAtCalo);
   fChain->SetBranchAddress("gsf_deltaPhiSeedClusterTrackAtCalo", &gsf_deltaPhiSeedClusterTrackAtCalo, &b_gsf_deltaPhiSeedClusterTrackAtCalo);
   fChain->SetBranchAddress("gsf_ecalEnergy", &gsf_ecalEnergy, &b_gsf_ecalEnergy);
   fChain->SetBranchAddress("gsf_eSuperClusterOverP", &gsf_eSuperClusterOverP, &b_gsf_eSuperClusterOverP);
   fChain->SetBranchAddress("gsf_dxy", &gsf_dxy, &b_gsf_dxy);
   fChain->SetBranchAddress("gsf_dxy_beamSpot", &gsf_dxy_beamSpot, &b_gsf_dxy_beamSpot);
   fChain->SetBranchAddress("gsf_dxy_firstPVtx", &gsf_dxy_firstPVtx, &b_gsf_dxy_firstPVtx);
   fChain->SetBranchAddress("gsf_dxyError", &gsf_dxyError, &b_gsf_dxyError);
   fChain->SetBranchAddress("gsf_dz", &gsf_dz, &b_gsf_dz);
   fChain->SetBranchAddress("gsf_dz_beamSpot", &gsf_dz_beamSpot, &b_gsf_dz_beamSpot);
   fChain->SetBranchAddress("gsf_dz_firstPVtx", &gsf_dz_firstPVtx, &b_gsf_dz_firstPVtx);
   fChain->SetBranchAddress("gsf_dzError", &gsf_dzError, &b_gsf_dzError);
   fChain->SetBranchAddress("gsf_vz", &gsf_vz, &b_gsf_vz);
   fChain->SetBranchAddress("gsf_numberOfValidHits", &gsf_numberOfValidHits, &b_gsf_numberOfValidHits);
   fChain->SetBranchAddress("gsf_nLostInnerHits", &gsf_nLostInnerHits, &b_gsf_nLostInnerHits);
   fChain->SetBranchAddress("gsf_nLostOuterHits", &gsf_nLostOuterHits, &b_gsf_nLostOuterHits);
   fChain->SetBranchAddress("gsf_convFlags", &gsf_convFlags, &b_gsf_convFlags);
   fChain->SetBranchAddress("gsf_convDist", &gsf_convDist, &b_gsf_convDist);
   fChain->SetBranchAddress("gsf_convDcot", &gsf_convDcot, &b_gsf_convDcot);
   fChain->SetBranchAddress("gsf_convRadius", &gsf_convRadius, &b_gsf_convRadius);
   fChain->SetBranchAddress("gsf_fBrem", &gsf_fBrem, &b_gsf_fBrem);
   fChain->SetBranchAddress("gsf_e1x5", &gsf_e1x5, &b_gsf_e1x5);
   fChain->SetBranchAddress("gsf_e2x5Max", &gsf_e2x5Max, &b_gsf_e2x5Max);
   fChain->SetBranchAddress("gsf_e5x5", &gsf_e5x5, &b_gsf_e5x5);
   fChain->SetBranchAddress("gsf_r9", &gsf_r9, &b_gsf_r9);
   fChain->SetBranchAddress("gsf_hitsinfo", &gsf_hitsinfo, &b_gsf_hitsinfo);
   fChain->SetBranchAddress("gsf_pixelMatch_dPhi1", &gsf_pixelMatch_dPhi1, &b_gsf_pixelMatch_dPhi1);
   fChain->SetBranchAddress("gsf_pixelMatch_dPhi2", &gsf_pixelMatch_dPhi2, &b_gsf_pixelMatch_dPhi2);
   fChain->SetBranchAddress("gsf_pixelMatch_dRz1", &gsf_pixelMatch_dRz1, &b_gsf_pixelMatch_dRz1);
   fChain->SetBranchAddress("gsf_pixelMatch_dRz2", &gsf_pixelMatch_dRz2, &b_gsf_pixelMatch_dRz2);
   fChain->SetBranchAddress("gsf_pixelMatch_subDetector1", &gsf_pixelMatch_subDetector1, &b_gsf_pixelMatch_subDetector1);
   fChain->SetBranchAddress("gsf_pixelMatch_subDetector2", &gsf_pixelMatch_subDetector2, &b_gsf_pixelMatch_subDetector2);
   fChain->SetBranchAddress("gsf_mc_bestDR", &gsf_mc_bestDR, &b_gsf_mc_bestDR);
   fChain->SetBranchAddress("gsf_mc_index", &gsf_mc_index, &b_gsf_mc_index);
   fChain->SetBranchAddress("gsf_mc_ERatio", &gsf_mc_ERatio, &b_gsf_mc_ERatio);
   fChain->SetBranchAddress("mu_n", &mu_n, &b_mu_n);
   fChain->SetBranchAddress("mu_n_saved", &mu_n_saved, &b_mu_n_saved);
   fChain->SetBranchAddress("mu_gt_n", &mu_gt_n, &b_mu_gt_n);
   fChain->SetBranchAddress("mu_ot_n", &mu_ot_n, &b_mu_ot_n);
   fChain->SetBranchAddress("mu_it_n", &mu_it_n, &b_mu_it_n);
   fChain->SetBranchAddress("mu_gt_qoverp", &mu_gt_qoverp, &b_mu_gt_qoverp);
   fChain->SetBranchAddress("mu_gt_charge", &mu_gt_charge, &b_mu_gt_charge);
   fChain->SetBranchAddress("mu_gt_pt", &mu_gt_pt, &b_mu_gt_pt);
   fChain->SetBranchAddress("mu_gt_eta", &mu_gt_eta, &b_mu_gt_eta);
   fChain->SetBranchAddress("mu_gt_phi", &mu_gt_phi, &b_mu_gt_phi);
   fChain->SetBranchAddress("mu_gt_p", &mu_gt_p, &b_mu_gt_p);
   fChain->SetBranchAddress("mu_gt_px", &mu_gt_px, &b_mu_gt_px);
   fChain->SetBranchAddress("mu_gt_py", &mu_gt_py, &b_mu_gt_py);
   fChain->SetBranchAddress("mu_gt_pz", &mu_gt_pz, &b_mu_gt_pz);
   fChain->SetBranchAddress("mu_gt_theta", &mu_gt_theta, &b_mu_gt_theta);
   fChain->SetBranchAddress("mu_gt_lambda", &mu_gt_lambda, &b_mu_gt_lambda);
   fChain->SetBranchAddress("mu_gt_d0", &mu_gt_d0, &b_mu_gt_d0);
   fChain->SetBranchAddress("mu_gt_dz", &mu_gt_dz, &b_mu_gt_dz);
   fChain->SetBranchAddress("mu_gt_dz_beamspot", &mu_gt_dz_beamspot, &b_mu_gt_dz_beamspot);
   fChain->SetBranchAddress("mu_gt_dz_firstPVtx", &mu_gt_dz_firstPVtx, &b_mu_gt_dz_firstPVtx);
   fChain->SetBranchAddress("mu_gt_dxy", &mu_gt_dxy, &b_mu_gt_dxy);
   fChain->SetBranchAddress("mu_gt_dxy_beamspot", &mu_gt_dxy_beamspot, &b_mu_gt_dxy_beamspot);
   fChain->SetBranchAddress("mu_gt_dxy_firstPVtx", &mu_gt_dxy_firstPVtx, &b_mu_gt_dxy_firstPVtx);
   fChain->SetBranchAddress("mu_gt_dsz", &mu_gt_dsz, &b_mu_gt_dsz);
   fChain->SetBranchAddress("mu_gt_vx", &mu_gt_vx, &b_mu_gt_vx);
   fChain->SetBranchAddress("mu_gt_vy", &mu_gt_vy, &b_mu_gt_vy);
   fChain->SetBranchAddress("mu_gt_vz", &mu_gt_vz, &b_mu_gt_vz);
   fChain->SetBranchAddress("mu_gt_qoverpError", &mu_gt_qoverpError, &b_mu_gt_qoverpError);
   fChain->SetBranchAddress("mu_gt_ptError", &mu_gt_ptError, &b_mu_gt_ptError);
   fChain->SetBranchAddress("mu_gt_thetaError", &mu_gt_thetaError, &b_mu_gt_thetaError);
   fChain->SetBranchAddress("mu_gt_lambdaError", &mu_gt_lambdaError, &b_mu_gt_lambdaError);
   fChain->SetBranchAddress("mu_gt_phiError", &mu_gt_phiError, &b_mu_gt_phiError);
   fChain->SetBranchAddress("mu_gt_dxyError", &mu_gt_dxyError, &b_mu_gt_dxyError);
   fChain->SetBranchAddress("mu_gt_d0Error", &mu_gt_d0Error, &b_mu_gt_d0Error);
   fChain->SetBranchAddress("mu_gt_dszError", &mu_gt_dszError, &b_mu_gt_dszError);
   fChain->SetBranchAddress("mu_gt_dzError", &mu_gt_dzError, &b_mu_gt_dzError);
   fChain->SetBranchAddress("mu_gt_etaError", &mu_gt_etaError, &b_mu_gt_etaError);
   fChain->SetBranchAddress("mu_gt_chi2", &mu_gt_chi2, &b_mu_gt_chi2);
   fChain->SetBranchAddress("mu_gt_ndof", &mu_gt_ndof, &b_mu_gt_ndof);
   fChain->SetBranchAddress("mu_gt_normalizedChi2", &mu_gt_normalizedChi2, &b_mu_gt_normalizedChi2);
   fChain->SetBranchAddress("mu_ot_qoverp", &mu_ot_qoverp, &b_mu_ot_qoverp);
   fChain->SetBranchAddress("mu_ot_charge", &mu_ot_charge, &b_mu_ot_charge);
   fChain->SetBranchAddress("mu_ot_pt", &mu_ot_pt, &b_mu_ot_pt);
   fChain->SetBranchAddress("mu_ot_eta", &mu_ot_eta, &b_mu_ot_eta);
   fChain->SetBranchAddress("mu_ot_phi", &mu_ot_phi, &b_mu_ot_phi);
   fChain->SetBranchAddress("mu_ot_p", &mu_ot_p, &b_mu_ot_p);
   fChain->SetBranchAddress("mu_ot_px", &mu_ot_px, &b_mu_ot_px);
   fChain->SetBranchAddress("mu_ot_py", &mu_ot_py, &b_mu_ot_py);
   fChain->SetBranchAddress("mu_ot_pz", &mu_ot_pz, &b_mu_ot_pz);
   fChain->SetBranchAddress("mu_ot_theta", &mu_ot_theta, &b_mu_ot_theta);
   fChain->SetBranchAddress("mu_ot_lambda", &mu_ot_lambda, &b_mu_ot_lambda);
   fChain->SetBranchAddress("mu_ot_d0", &mu_ot_d0, &b_mu_ot_d0);
   fChain->SetBranchAddress("mu_ot_dz", &mu_ot_dz, &b_mu_ot_dz);
   fChain->SetBranchAddress("mu_ot_dz_beamspot", &mu_ot_dz_beamspot, &b_mu_ot_dz_beamspot);
   fChain->SetBranchAddress("mu_ot_dz_firstPVtx", &mu_ot_dz_firstPVtx, &b_mu_ot_dz_firstPVtx);
   fChain->SetBranchAddress("mu_ot_dxy", &mu_ot_dxy, &b_mu_ot_dxy);
   fChain->SetBranchAddress("mu_ot_dxy_beamspot", &mu_ot_dxy_beamspot, &b_mu_ot_dxy_beamspot);
   fChain->SetBranchAddress("mu_ot_dxy_firstPVtx", &mu_ot_dxy_firstPVtx, &b_mu_ot_dxy_firstPVtx);
   fChain->SetBranchAddress("mu_ot_dsz", &mu_ot_dsz, &b_mu_ot_dsz);
   fChain->SetBranchAddress("mu_ot_vx", &mu_ot_vx, &b_mu_ot_vx);
   fChain->SetBranchAddress("mu_ot_vy", &mu_ot_vy, &b_mu_ot_vy);
   fChain->SetBranchAddress("mu_ot_vz", &mu_ot_vz, &b_mu_ot_vz);
   fChain->SetBranchAddress("mu_ot_qoverpError", &mu_ot_qoverpError, &b_mu_ot_qoverpError);
   fChain->SetBranchAddress("mu_ot_ptError", &mu_ot_ptError, &b_mu_ot_ptError);
   fChain->SetBranchAddress("mu_ot_thetaError", &mu_ot_thetaError, &b_mu_ot_thetaError);
   fChain->SetBranchAddress("mu_ot_lambdaError", &mu_ot_lambdaError, &b_mu_ot_lambdaError);
   fChain->SetBranchAddress("mu_ot_phiError", &mu_ot_phiError, &b_mu_ot_phiError);
   fChain->SetBranchAddress("mu_ot_dxyError", &mu_ot_dxyError, &b_mu_ot_dxyError);
   fChain->SetBranchAddress("mu_ot_d0Error", &mu_ot_d0Error, &b_mu_ot_d0Error);
   fChain->SetBranchAddress("mu_ot_dszError", &mu_ot_dszError, &b_mu_ot_dszError);
   fChain->SetBranchAddress("mu_ot_dzError", &mu_ot_dzError, &b_mu_ot_dzError);
   fChain->SetBranchAddress("mu_ot_etaError", &mu_ot_etaError, &b_mu_ot_etaError);
   fChain->SetBranchAddress("mu_ot_chi2", &mu_ot_chi2, &b_mu_ot_chi2);
   fChain->SetBranchAddress("mu_ot_ndof", &mu_ot_ndof, &b_mu_ot_ndof);
   fChain->SetBranchAddress("mu_ot_normalizedChi2", &mu_ot_normalizedChi2, &b_mu_ot_normalizedChi2);
   fChain->SetBranchAddress("mu_it_qoverp", &mu_it_qoverp, &b_mu_it_qoverp);
   fChain->SetBranchAddress("mu_it_charge", &mu_it_charge, &b_mu_it_charge);
   fChain->SetBranchAddress("mu_it_pt", &mu_it_pt, &b_mu_it_pt);
   fChain->SetBranchAddress("mu_it_eta", &mu_it_eta, &b_mu_it_eta);
   fChain->SetBranchAddress("mu_it_phi", &mu_it_phi, &b_mu_it_phi);
   fChain->SetBranchAddress("mu_it_p", &mu_it_p, &b_mu_it_p);
   fChain->SetBranchAddress("mu_it_px", &mu_it_px, &b_mu_it_px);
   fChain->SetBranchAddress("mu_it_py", &mu_it_py, &b_mu_it_py);
   fChain->SetBranchAddress("mu_it_pz", &mu_it_pz, &b_mu_it_pz);
   fChain->SetBranchAddress("mu_it_theta", &mu_it_theta, &b_mu_it_theta);
   fChain->SetBranchAddress("mu_it_lambda", &mu_it_lambda, &b_mu_it_lambda);
   fChain->SetBranchAddress("mu_it_d0", &mu_it_d0, &b_mu_it_d0);
   fChain->SetBranchAddress("mu_it_dz", &mu_it_dz, &b_mu_it_dz);
   fChain->SetBranchAddress("mu_it_dz_beamspot", &mu_it_dz_beamspot, &b_mu_it_dz_beamspot);
   fChain->SetBranchAddress("mu_it_dz_firstPVtx", &mu_it_dz_firstPVtx, &b_mu_it_dz_firstPVtx);
   fChain->SetBranchAddress("mu_it_dxy", &mu_it_dxy, &b_mu_it_dxy);
   fChain->SetBranchAddress("mu_it_dxy_beamspot", &mu_it_dxy_beamspot, &b_mu_it_dxy_beamspot);
   fChain->SetBranchAddress("mu_it_dxy_firstPVtx", &mu_it_dxy_firstPVtx, &b_mu_it_dxy_firstPVtx);
   fChain->SetBranchAddress("mu_it_dsz", &mu_it_dsz, &b_mu_it_dsz);
   fChain->SetBranchAddress("mu_it_vx", &mu_it_vx, &b_mu_it_vx);
   fChain->SetBranchAddress("mu_it_vy", &mu_it_vy, &b_mu_it_vy);
   fChain->SetBranchAddress("mu_it_vz", &mu_it_vz, &b_mu_it_vz);
   fChain->SetBranchAddress("mu_it_qoverpError", &mu_it_qoverpError, &b_mu_it_qoverpError);
   fChain->SetBranchAddress("mu_it_ptError", &mu_it_ptError, &b_mu_it_ptError);
   fChain->SetBranchAddress("mu_it_thetaError", &mu_it_thetaError, &b_mu_it_thetaError);
   fChain->SetBranchAddress("mu_it_lambdaError", &mu_it_lambdaError, &b_mu_it_lambdaError);
   fChain->SetBranchAddress("mu_it_phiError", &mu_it_phiError, &b_mu_it_phiError);
   fChain->SetBranchAddress("mu_it_dxyError", &mu_it_dxyError, &b_mu_it_dxyError);
   fChain->SetBranchAddress("mu_it_d0Error", &mu_it_d0Error, &b_mu_it_d0Error);
   fChain->SetBranchAddress("mu_it_dszError", &mu_it_dszError, &b_mu_it_dszError);
   fChain->SetBranchAddress("mu_it_dzError", &mu_it_dzError, &b_mu_it_dzError);
   fChain->SetBranchAddress("mu_it_etaError", &mu_it_etaError, &b_mu_it_etaError);
   fChain->SetBranchAddress("mu_it_chi2", &mu_it_chi2, &b_mu_it_chi2);
   fChain->SetBranchAddress("mu_it_ndof", &mu_it_ndof, &b_mu_it_ndof);
   fChain->SetBranchAddress("mu_it_normalizedChi2", &mu_it_normalizedChi2, &b_mu_it_normalizedChi2);
   fChain->SetBranchAddress("mu_isGlobalMuon", &mu_isGlobalMuon, &b_mu_isGlobalMuon);
   fChain->SetBranchAddress("mu_isStandAloneMuon", &mu_isStandAloneMuon, &b_mu_isStandAloneMuon);
   fChain->SetBranchAddress("mu_isTrackerMuon", &mu_isTrackerMuon, &b_mu_isTrackerMuon);
   fChain->SetBranchAddress("mu_isPFMuon", &mu_isPFMuon, &b_mu_isPFMuon);
   fChain->SetBranchAddress("mu_isPFIsolationValid", &mu_isPFIsolationValid, &b_mu_isPFIsolationValid);
   fChain->SetBranchAddress("mu_numberOfMatchedStations", &mu_numberOfMatchedStations, &b_mu_numberOfMatchedStations);
   fChain->SetBranchAddress("mu_numberOfValidPixelHits", &mu_numberOfValidPixelHits, &b_mu_numberOfValidPixelHits);
   fChain->SetBranchAddress("mu_numberOfValidTrackerHits", &mu_numberOfValidTrackerHits, &b_mu_numberOfValidTrackerHits);
   fChain->SetBranchAddress("mu_numberOfValidMuonHits", &mu_numberOfValidMuonHits, &b_mu_numberOfValidMuonHits);
   fChain->SetBranchAddress("mu_tevOptimized_charge", &mu_tevOptimized_charge, &b_mu_tevOptimized_charge);
   fChain->SetBranchAddress("mu_tevOptimized_pt", &mu_tevOptimized_pt, &b_mu_tevOptimized_pt);
   fChain->SetBranchAddress("mu_tevOptimized_eta", &mu_tevOptimized_eta, &b_mu_tevOptimized_eta);
   fChain->SetBranchAddress("mu_tevOptimized_phi", &mu_tevOptimized_phi, &b_mu_tevOptimized_phi);
   fChain->SetBranchAddress("mu_tevOptimized_theta", &mu_tevOptimized_theta, &b_mu_tevOptimized_theta);
   fChain->SetBranchAddress("mu_tevOptimized_px", &mu_tevOptimized_px, &b_mu_tevOptimized_px);
   fChain->SetBranchAddress("mu_tevOptimized_py", &mu_tevOptimized_py, &b_mu_tevOptimized_py);
   fChain->SetBranchAddress("mu_tevOptimized_pz", &mu_tevOptimized_pz, &b_mu_tevOptimized_pz);
   fChain->SetBranchAddress("mu_tevOptimized_d0", &mu_tevOptimized_d0, &b_mu_tevOptimized_d0);
   fChain->SetBranchAddress("mu_tevOptimized_dz", &mu_tevOptimized_dz, &b_mu_tevOptimized_dz);
   fChain->SetBranchAddress("mu_tevOptimized_dz_beamSpot", &mu_tevOptimized_dz_beamSpot, &b_mu_tevOptimized_dz_beamSpot);
   fChain->SetBranchAddress("mu_tevOptimized_dz_firstPVtx", &mu_tevOptimized_dz_firstPVtx, &b_mu_tevOptimized_dz_firstPVtx);
   fChain->SetBranchAddress("mu_tevOptimized_dxy", &mu_tevOptimized_dxy, &b_mu_tevOptimized_dxy);
   fChain->SetBranchAddress("mu_tevOptimized_dxy_beamSpot", &mu_tevOptimized_dxy_beamSpot, &b_mu_tevOptimized_dxy_beamSpot);
   fChain->SetBranchAddress("mu_tevOptimized_dxy_firstPVtx", &mu_tevOptimized_dxy_firstPVtx, &b_mu_tevOptimized_dxy_firstPVtx);
   fChain->SetBranchAddress("mu_tevOptimized_ptError", &mu_tevOptimized_ptError, &b_mu_tevOptimized_ptError);
   fChain->SetBranchAddress("mu_tevOptimized_etaError", &mu_tevOptimized_etaError, &b_mu_tevOptimized_etaError);
   fChain->SetBranchAddress("mu_tevOptimized_phiError", &mu_tevOptimized_phiError, &b_mu_tevOptimized_phiError);
   fChain->SetBranchAddress("mu_tevOptimized_thetaError", &mu_tevOptimized_thetaError, &b_mu_tevOptimized_thetaError);
   fChain->SetBranchAddress("mu_tevOptimized_d0Error", &mu_tevOptimized_d0Error, &b_mu_tevOptimized_d0Error);
   fChain->SetBranchAddress("mu_tevOptimized_dzError", &mu_tevOptimized_dzError, &b_mu_tevOptimized_dzError);
   fChain->SetBranchAddress("mu_tevOptimized_dxyError", &mu_tevOptimized_dxyError, &b_mu_tevOptimized_dxyError);
   fChain->SetBranchAddress("mu_isolationR03_sumPt", &mu_isolationR03_sumPt, &b_mu_isolationR03_sumPt);
   fChain->SetBranchAddress("mu_isolationR03_trackerVetoPt", &mu_isolationR03_trackerVetoPt, &b_mu_isolationR03_trackerVetoPt);
   fChain->SetBranchAddress("mu_isolationR03_emEt", &mu_isolationR03_emEt, &b_mu_isolationR03_emEt);
   fChain->SetBranchAddress("mu_isolationR03_emVetoEt", &mu_isolationR03_emVetoEt, &b_mu_isolationR03_emVetoEt);
   fChain->SetBranchAddress("mu_isolationR03_hadEt", &mu_isolationR03_hadEt, &b_mu_isolationR03_hadEt);
   fChain->SetBranchAddress("mu_isolationR03_hadVetoEt", &mu_isolationR03_hadVetoEt, &b_mu_isolationR03_hadVetoEt);
   fChain->SetBranchAddress("mu_isolationR05_sumPt", &mu_isolationR05_sumPt, &b_mu_isolationR05_sumPt);
   fChain->SetBranchAddress("mu_isolationR05_trackerVetoPt", &mu_isolationR05_trackerVetoPt, &b_mu_isolationR05_trackerVetoPt);
   fChain->SetBranchAddress("mu_isolationR05_emEt", &mu_isolationR05_emEt, &b_mu_isolationR05_emEt);
   fChain->SetBranchAddress("mu_isolationR05_emVetoEt", &mu_isolationR05_emVetoEt, &b_mu_isolationR05_emVetoEt);
   fChain->SetBranchAddress("mu_isolationR05_hadEt", &mu_isolationR05_hadEt, &b_mu_isolationR05_hadEt);
   fChain->SetBranchAddress("mu_isolationR05_hadVetoEt", &mu_isolationR05_hadVetoEt, &b_mu_isolationR05_hadVetoEt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumChargedHadronPt", &mu_pfIsolationR03_sumChargedHadronPt, &b_mu_pfIsolationR03_sumChargedHadronPt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumChargedParticlePt", &mu_pfIsolationR03_sumChargedParticlePt, &b_mu_pfIsolationR03_sumChargedParticlePt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumPhotonEt", &mu_pfIsolationR03_sumPhotonEt, &b_mu_pfIsolationR03_sumPhotonEt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumNeutralHadronEtHighThreshold", &mu_pfIsolationR03_sumNeutralHadronEtHighThreshold, &b_mu_pfIsolationR03_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumPhotonEtHighThreshold", &mu_pfIsolationR03_sumPhotonEtHighThreshold, &b_mu_pfIsolationR03_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumPUPt", &mu_pfIsolationR03_sumPUPt, &b_mu_pfIsolationR03_sumPUPt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumChargedHadronPt", &mu_pfIsolationR04_sumChargedHadronPt, &b_mu_pfIsolationR04_sumChargedHadronPt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumChargedParticlePt", &mu_pfIsolationR04_sumChargedParticlePt, &b_mu_pfIsolationR04_sumChargedParticlePt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumPhotonEt", &mu_pfIsolationR04_sumPhotonEt, &b_mu_pfIsolationR04_sumPhotonEt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumNeutralHadronEtHighThreshold", &mu_pfIsolationR04_sumNeutralHadronEtHighThreshold, &b_mu_pfIsolationR04_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumPhotonEtHighThreshold", &mu_pfIsolationR04_sumPhotonEtHighThreshold, &b_mu_pfIsolationR04_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumPUPt", &mu_pfIsolationR04_sumPUPt, &b_mu_pfIsolationR04_sumPUPt);
   fChain->SetBranchAddress("mu_pfMeanDRIsoProfileR03_sumChargedHadronPt", &mu_pfMeanDRIsoProfileR03_sumChargedHadronPt, &b_mu_pfMeanDRIsoProfileR03_sumChargedHadronPt);
   fChain->SetBranchAddress("mu_pfMeanDRIsoProfileR03_sumChargedParticlePt", &mu_pfMeanDRIsoProfileR03_sumChargedParticlePt, &b_mu_pfMeanDRIsoProfileR03_sumChargedParticlePt);
   fChain->SetBranchAddress("mu_pfMeanDRIsoProfileR03_sumPhotonEt", &mu_pfMeanDRIsoProfileR03_sumPhotonEt, &b_mu_pfMeanDRIsoProfileR03_sumPhotonEt);
   fChain->SetBranchAddress("mu_pfMeanDRIsoProfileR03_sumNeutralHadronEtHighThreshold", &mu_pfMeanDRIsoProfileR03_sumNeutralHadronEtHighThreshold, &b_mu_pfMeanDRIsoProfileR03_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("mu_pfMeanDRIsoProfileR03_sumPhotonEtHighThreshold", &mu_pfMeanDRIsoProfileR03_sumPhotonEtHighThreshold, &b_mu_pfMeanDRIsoProfileR03_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("mu_pfMeanDRIsoProfileR03_sumPUPt", &mu_pfMeanDRIsoProfileR03_sumPUPt, &b_mu_pfMeanDRIsoProfileR03_sumPUPt);
   fChain->SetBranchAddress("mu_pfMeanDRIsoProfileR04_sumChargedHadronPt", &mu_pfMeanDRIsoProfileR04_sumChargedHadronPt, &b_mu_pfMeanDRIsoProfileR04_sumChargedHadronPt);
   fChain->SetBranchAddress("mu_pfMeanDRIsoProfileR04_sumChargedParticlePt", &mu_pfMeanDRIsoProfileR04_sumChargedParticlePt, &b_mu_pfMeanDRIsoProfileR04_sumChargedParticlePt);
   fChain->SetBranchAddress("mu_pfMeanDRIsoProfileR04_sumPhotonEt", &mu_pfMeanDRIsoProfileR04_sumPhotonEt, &b_mu_pfMeanDRIsoProfileR04_sumPhotonEt);
   fChain->SetBranchAddress("mu_pfMeanDRIsoProfileR04_sumNeutralHadronEtHighThreshold", &mu_pfMeanDRIsoProfileR04_sumNeutralHadronEtHighThreshold, &b_mu_pfMeanDRIsoProfileR04_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("mu_pfMeanDRIsoProfileR04_sumPhotonEtHighThreshold", &mu_pfMeanDRIsoProfileR04_sumPhotonEtHighThreshold, &b_mu_pfMeanDRIsoProfileR04_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("mu_pfMeanDRIsoProfileR04_sumPUPt", &mu_pfMeanDRIsoProfileR04_sumPUPt, &b_mu_pfMeanDRIsoProfileR04_sumPUPt);
   fChain->SetBranchAddress("mu_mc_bestDR", &mu_mc_bestDR, &b_mu_mc_bestDR);
   fChain->SetBranchAddress("mu_mc_index", &mu_mc_index, &b_mu_mc_index);
   fChain->SetBranchAddress("mu_mc_ERatio", &mu_mc_ERatio, &b_mu_mc_ERatio);
   fChain->SetBranchAddress("MET_caloMet_et", &MET_caloMet_et, &b_MET_caloMet_et);
   fChain->SetBranchAddress("MET_caloMet_phi", &MET_caloMet_phi, &b_MET_caloMet_phi);
   fChain->SetBranchAddress("MET_pfMet_et", &MET_pfMet_et, &b_MET_pfMet_et);
   fChain->SetBranchAddress("MET_pfMet_phi", &MET_pfMet_phi, &b_MET_pfMet_phi);
   fChain->SetBranchAddress("HEEP_eseffsixix", &HEEP_eseffsixix, &b_HEEP_eseffsixix);
   fChain->SetBranchAddress("HEEP_eseffsiyiy", &HEEP_eseffsiyiy, &b_HEEP_eseffsiyiy);
   fChain->SetBranchAddress("HEEP_eseffsirir", &HEEP_eseffsirir, &b_HEEP_eseffsirir);
   fChain->SetBranchAddress("HEEP_preshowerEnergy", &HEEP_preshowerEnergy, &b_HEEP_preshowerEnergy);
   fChain->SetBranchAddress("HEEP_e1x3", &HEEP_e1x3, &b_HEEP_e1x3);
   fChain->SetBranchAddress("HEEP_eMax", &HEEP_eMax, &b_HEEP_eMax);
   fChain->SetBranchAddress("HEEP_e5x5", &HEEP_e5x5, &b_HEEP_e5x5);
   fChain->SetBranchAddress("HEEP_e2x5Right", &HEEP_e2x5Right, &b_HEEP_e2x5Right);
   fChain->SetBranchAddress("HEEP_e2x5Left", &HEEP_e2x5Left, &b_HEEP_e2x5Left);
   fChain->SetBranchAddress("HEEP_e2x5Top", &HEEP_e2x5Top, &b_HEEP_e2x5Top);
   fChain->SetBranchAddress("HEEP_e2x5Bottom", &HEEP_e2x5Bottom, &b_HEEP_e2x5Bottom);
   fChain->SetBranchAddress("HEEP_eRight", &HEEP_eRight, &b_HEEP_eRight);
   fChain->SetBranchAddress("HEEP_eLeft", &HEEP_eLeft, &b_HEEP_eLeft);
   fChain->SetBranchAddress("HEEP_eTop", &HEEP_eTop, &b_HEEP_eTop);
   fChain->SetBranchAddress("HEEP_eBottom", &HEEP_eBottom, &b_HEEP_eBottom);
   fChain->SetBranchAddress("HEEP_basicClusterSeedTime", &HEEP_basicClusterSeedTime, &b_HEEP_basicClusterSeedTime);
   fChain->SetBranchAddress("EBHits_rawId", &EBHits_rawId, &b_EBHits_rawId);
   fChain->SetBranchAddress("EBHits_iRechit", &EBHits_iRechit, &b_EBHits_iRechit);
   fChain->SetBranchAddress("EBHits_energy", &EBHits_energy, &b_EBHits_energy);
   fChain->SetBranchAddress("EBHits_ieta", &EBHits_ieta, &b_EBHits_ieta);
   fChain->SetBranchAddress("EBHits_iphi", &EBHits_iphi, &b_EBHits_iphi);
   fChain->SetBranchAddress("EBHits_RecoFlag", &EBHits_RecoFlag, &b_EBHits_RecoFlag);
   fChain->SetBranchAddress("EBHits_kSaturated", &EBHits_kSaturated, &b_EBHits_kSaturated);
   fChain->SetBranchAddress("EBHits_kLeadingEdgeRecovered", &EBHits_kLeadingEdgeRecovered, &b_EBHits_kLeadingEdgeRecovered);
   fChain->SetBranchAddress("EBHits_kNeighboursRecovered", &EBHits_kNeighboursRecovered, &b_EBHits_kNeighboursRecovered);
   fChain->SetBranchAddress("EBHits_kWeird", &EBHits_kWeird, &b_EBHits_kWeird);
   fChain->SetBranchAddress("EEHits_rawId", &EEHits_rawId, &b_EEHits_rawId);
   fChain->SetBranchAddress("EEHits_iRechit", &EEHits_iRechit, &b_EEHits_iRechit);
   fChain->SetBranchAddress("EEHits_energy", &EEHits_energy, &b_EEHits_energy);
   fChain->SetBranchAddress("EEHits_ieta", &EEHits_ieta, &b_EEHits_ieta);
   fChain->SetBranchAddress("EEHits_iphi", &EEHits_iphi, &b_EEHits_iphi);
   fChain->SetBranchAddress("EEHits_RecoFlag", &EEHits_RecoFlag, &b_EEHits_RecoFlag);
   fChain->SetBranchAddress("EEHits_kSaturated", &EEHits_kSaturated, &b_EEHits_kSaturated);
   fChain->SetBranchAddress("EEHits_kLeadingEdgeRecovered", &EEHits_kLeadingEdgeRecovered, &b_EEHits_kLeadingEdgeRecovered);
   fChain->SetBranchAddress("EEHits_kNeighboursRecovered", &EEHits_kNeighboursRecovered, &b_EEHits_kNeighboursRecovered);
   fChain->SetBranchAddress("EEHits_kWeird", &EEHits_kWeird, &b_EEHits_kWeird);
   fChain->SetBranchAddress("HEEP_crystal_energy", &HEEP_crystal_energy, &b_HEEP_crystal_energy);
   fChain->SetBranchAddress("HEEP_crystal_eta", &HEEP_crystal_eta, &b_HEEP_crystal_eta);
   fChain->SetBranchAddress("HEEP_eshitsixix", &HEEP_eshitsixix, &b_HEEP_eshitsixix);
   fChain->SetBranchAddress("HEEP_eshitsiyiy", &HEEP_eshitsiyiy, &b_HEEP_eshitsiyiy);
   fChain->SetBranchAddress("HEEP_crystal_ietaorix", &HEEP_crystal_ietaorix, &b_HEEP_crystal_ietaorix);
   fChain->SetBranchAddress("HEEP_crystal_iphioriy", &HEEP_crystal_iphioriy, &b_HEEP_crystal_iphioriy);
   fChain->SetBranchAddress("HEEP_cutflow41_Et", &HEEP_cutflow41_Et, &b_HEEP_cutflow41_Et);
   fChain->SetBranchAddress("HEEP_cutflow41_Et_n", &HEEP_cutflow41_Et_n, &b_HEEP_cutflow41_Et_n);
   fChain->SetBranchAddress("HEEP_cutflow41_Et_nCumulative", &HEEP_cutflow41_Et_nCumulative, &b_HEEP_cutflow41_Et_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow41_Et_value", &HEEP_cutflow41_Et_value, &b_HEEP_cutflow41_Et_value);
   fChain->SetBranchAddress("HEEP_cutflow41_eta", &HEEP_cutflow41_eta, &b_HEEP_cutflow41_eta);
   fChain->SetBranchAddress("HEEP_cutflow41_eta_n", &HEEP_cutflow41_eta_n, &b_HEEP_cutflow41_eta_n);
   fChain->SetBranchAddress("HEEP_cutflow41_eta_nCumulative", &HEEP_cutflow41_eta_nCumulative, &b_HEEP_cutflow41_eta_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow41_eta_value", &HEEP_cutflow41_eta_value, &b_HEEP_cutflow41_eta_value);
   fChain->SetBranchAddress("HEEP_cutflow41_acceptance", &HEEP_cutflow41_acceptance, &b_HEEP_cutflow41_acceptance);
   fChain->SetBranchAddress("HEEP_cutflow41_acceptance_n", &HEEP_cutflow41_acceptance_n, &b_HEEP_cutflow41_acceptance_n);
   fChain->SetBranchAddress("HEEP_cutflow41_EcalDriven", &HEEP_cutflow41_EcalDriven, &b_HEEP_cutflow41_EcalDriven);
   fChain->SetBranchAddress("HEEP_cutflow41_EcalDriven_n", &HEEP_cutflow41_EcalDriven_n, &b_HEEP_cutflow41_EcalDriven_n);
   fChain->SetBranchAddress("HEEP_cutflow41_EcalDriven_nCumulative", &HEEP_cutflow41_EcalDriven_nCumulative, &b_HEEP_cutflow41_EcalDriven_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow41_EcalDriven_value", &HEEP_cutflow41_EcalDriven_value, &b_HEEP_cutflow41_EcalDriven_value);
   fChain->SetBranchAddress("HEEP_cutflow41_dEtaIn", &HEEP_cutflow41_dEtaIn, &b_HEEP_cutflow41_dEtaIn);
   fChain->SetBranchAddress("HEEP_cutflow41_dEtaIn_n", &HEEP_cutflow41_dEtaIn_n, &b_HEEP_cutflow41_dEtaIn_n);
   fChain->SetBranchAddress("HEEP_cutflow41_dEtaIn_nCumulative", &HEEP_cutflow41_dEtaIn_nCumulative, &b_HEEP_cutflow41_dEtaIn_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow41_dEtaIn_value", &HEEP_cutflow41_dEtaIn_value, &b_HEEP_cutflow41_dEtaIn_value);
   fChain->SetBranchAddress("HEEP_cutflow41_dPhiIn", &HEEP_cutflow41_dPhiIn, &b_HEEP_cutflow41_dPhiIn);
   fChain->SetBranchAddress("HEEP_cutflow41_dPhiIn_n", &HEEP_cutflow41_dPhiIn_n, &b_HEEP_cutflow41_dPhiIn_n);
   fChain->SetBranchAddress("HEEP_cutflow41_dPhiIn_nCumulative", &HEEP_cutflow41_dPhiIn_nCumulative, &b_HEEP_cutflow41_dPhiIn_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow41_dPhiIn_value", &HEEP_cutflow41_dPhiIn_value, &b_HEEP_cutflow41_dPhiIn_value);
   fChain->SetBranchAddress("HEEP_cutflow41_HOverE", &HEEP_cutflow41_HOverE, &b_HEEP_cutflow41_HOverE);
   fChain->SetBranchAddress("HEEP_cutflow41_HOverE_n", &HEEP_cutflow41_HOverE_n, &b_HEEP_cutflow41_HOverE_n);
   fChain->SetBranchAddress("HEEP_cutflow41_HOverE_nCumulative", &HEEP_cutflow41_HOverE_nCumulative, &b_HEEP_cutflow41_HOverE_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow41_HOverE_value", &HEEP_cutflow41_HOverE_value, &b_HEEP_cutflow41_HOverE_value);
   fChain->SetBranchAddress("HEEP_cutflow41_SigmaIetaIeta", &HEEP_cutflow41_SigmaIetaIeta, &b_HEEP_cutflow41_SigmaIetaIeta);
   fChain->SetBranchAddress("HEEP_cutflow41_SigmaIetaIeta_n", &HEEP_cutflow41_SigmaIetaIeta_n, &b_HEEP_cutflow41_SigmaIetaIeta_n);
   fChain->SetBranchAddress("HEEP_cutflow41_SigmaIetaIeta_nCumulative", &HEEP_cutflow41_SigmaIetaIeta_nCumulative, &b_HEEP_cutflow41_SigmaIetaIeta_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow41_SigmaIetaIeta_value", &HEEP_cutflow41_SigmaIetaIeta_value, &b_HEEP_cutflow41_SigmaIetaIeta_value);
   fChain->SetBranchAddress("HEEP_cutflow41_E1x5OverE5x5", &HEEP_cutflow41_E1x5OverE5x5, &b_HEEP_cutflow41_E1x5OverE5x5);
   fChain->SetBranchAddress("HEEP_cutflow41_E1x5OverE5x5_n", &HEEP_cutflow41_E1x5OverE5x5_n, &b_HEEP_cutflow41_E1x5OverE5x5_n);
   fChain->SetBranchAddress("HEEP_cutflow41_E1x5OverE5x5_nCumulative", &HEEP_cutflow41_E1x5OverE5x5_nCumulative, &b_HEEP_cutflow41_E1x5OverE5x5_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow41_E1x5OverE5x5_value", &HEEP_cutflow41_E1x5OverE5x5_value, &b_HEEP_cutflow41_E1x5OverE5x5_value);
   fChain->SetBranchAddress("HEEP_cutflow41_E2x5OverE5x5", &HEEP_cutflow41_E2x5OverE5x5, &b_HEEP_cutflow41_E2x5OverE5x5);
   fChain->SetBranchAddress("HEEP_cutflow41_E2x5OverE5x5_n", &HEEP_cutflow41_E2x5OverE5x5_n, &b_HEEP_cutflow41_E2x5OverE5x5_n);
   fChain->SetBranchAddress("HEEP_cutflow41_E2x5OverE5x5_nCumulative", &HEEP_cutflow41_E2x5OverE5x5_nCumulative, &b_HEEP_cutflow41_E2x5OverE5x5_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow41_E2x5OverE5x5_value", &HEEP_cutflow41_E2x5OverE5x5_value, &b_HEEP_cutflow41_E2x5OverE5x5_value);
   fChain->SetBranchAddress("HEEP_cutflow41_missingHits", &HEEP_cutflow41_missingHits, &b_HEEP_cutflow41_missingHits);
   fChain->SetBranchAddress("HEEP_cutflow41_missingHits_n", &HEEP_cutflow41_missingHits_n, &b_HEEP_cutflow41_missingHits_n);
   fChain->SetBranchAddress("HEEP_cutflow41_missingHits_nCumulative", &HEEP_cutflow41_missingHits_nCumulative, &b_HEEP_cutflow41_missingHits_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow41_missingHits_value", &HEEP_cutflow41_missingHits_value, &b_HEEP_cutflow41_missingHits_value);
   fChain->SetBranchAddress("HEEP_cutflow41_dxyFirstPV", &HEEP_cutflow41_dxyFirstPV, &b_HEEP_cutflow41_dxyFirstPV);
   fChain->SetBranchAddress("HEEP_cutflow41_dxyFirstPV_n", &HEEP_cutflow41_dxyFirstPV_n, &b_HEEP_cutflow41_dxyFirstPV_n);
   fChain->SetBranchAddress("HEEP_cutflow41_dxyFirstPV_nCumulative", &HEEP_cutflow41_dxyFirstPV_nCumulative, &b_HEEP_cutflow41_dxyFirstPV_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow41_dxyFirstPV_value", &HEEP_cutflow41_dxyFirstPV_value, &b_HEEP_cutflow41_dxyFirstPV_value);
   fChain->SetBranchAddress("HEEP_cutflow41_ID", &HEEP_cutflow41_ID, &b_HEEP_cutflow41_ID);
   fChain->SetBranchAddress("HEEP_cutflow41_ID_n", &HEEP_cutflow41_ID_n, &b_HEEP_cutflow41_ID_n);
   fChain->SetBranchAddress("HEEP_cutflow41_isolEMHadDepth1", &HEEP_cutflow41_isolEMHadDepth1, &b_HEEP_cutflow41_isolEMHadDepth1);
   fChain->SetBranchAddress("HEEP_cutflow41_isolEMHadDepth1_n", &HEEP_cutflow41_isolEMHadDepth1_n, &b_HEEP_cutflow41_isolEMHadDepth1_n);
   fChain->SetBranchAddress("HEEP_cutflow41_isolEMHadDepth1_nCumulative", &HEEP_cutflow41_isolEMHadDepth1_nCumulative, &b_HEEP_cutflow41_isolEMHadDepth1_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow41_isolEMHadDepth1_value", &HEEP_cutflow41_isolEMHadDepth1_value, &b_HEEP_cutflow41_isolEMHadDepth1_value);
   fChain->SetBranchAddress("HEEP_cutflow41_IsolPtTrks", &HEEP_cutflow41_IsolPtTrks, &b_HEEP_cutflow41_IsolPtTrks);
   fChain->SetBranchAddress("HEEP_cutflow41_IsolPtTrks_n", &HEEP_cutflow41_IsolPtTrks_n, &b_HEEP_cutflow41_IsolPtTrks_n);
   fChain->SetBranchAddress("HEEP_cutflow41_IsolPtTrks_nCumulative", &HEEP_cutflow41_IsolPtTrks_nCumulative, &b_HEEP_cutflow41_IsolPtTrks_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow41_IsolPtTrks_value", &HEEP_cutflow41_IsolPtTrks_value, &b_HEEP_cutflow41_IsolPtTrks_value);
   fChain->SetBranchAddress("HEEP_cutflow41_isolation", &HEEP_cutflow41_isolation, &b_HEEP_cutflow41_isolation);
   fChain->SetBranchAddress("HEEP_cutflow41_isolation_n", &HEEP_cutflow41_isolation_n, &b_HEEP_cutflow41_isolation_n);
   fChain->SetBranchAddress("HEEP_cutflow41_total", &HEEP_cutflow41_total, &b_HEEP_cutflow41_total);
   fChain->SetBranchAddress("HEEP_cutflow41_total_n", &HEEP_cutflow41_total_n, &b_HEEP_cutflow41_total_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_Et", &HEEP_cutflow50_25ns_Et, &b_HEEP_cutflow50_25ns_Et);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_Et_n", &HEEP_cutflow50_25ns_Et_n, &b_HEEP_cutflow50_25ns_Et_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_Et_nCumulative", &HEEP_cutflow50_25ns_Et_nCumulative, &b_HEEP_cutflow50_25ns_Et_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_Et_value", &HEEP_cutflow50_25ns_Et_value, &b_HEEP_cutflow50_25ns_Et_value);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_eta", &HEEP_cutflow50_25ns_eta, &b_HEEP_cutflow50_25ns_eta);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_eta_n", &HEEP_cutflow50_25ns_eta_n, &b_HEEP_cutflow50_25ns_eta_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_eta_nCumulative", &HEEP_cutflow50_25ns_eta_nCumulative, &b_HEEP_cutflow50_25ns_eta_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_eta_value", &HEEP_cutflow50_25ns_eta_value, &b_HEEP_cutflow50_25ns_eta_value);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_acceptance", &HEEP_cutflow50_25ns_acceptance, &b_HEEP_cutflow50_25ns_acceptance);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_acceptance_n", &HEEP_cutflow50_25ns_acceptance_n, &b_HEEP_cutflow50_25ns_acceptance_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_EcalDriven", &HEEP_cutflow50_25ns_EcalDriven, &b_HEEP_cutflow50_25ns_EcalDriven);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_EcalDriven_n", &HEEP_cutflow50_25ns_EcalDriven_n, &b_HEEP_cutflow50_25ns_EcalDriven_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_EcalDriven_nCumulative", &HEEP_cutflow50_25ns_EcalDriven_nCumulative, &b_HEEP_cutflow50_25ns_EcalDriven_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_EcalDriven_value", &HEEP_cutflow50_25ns_EcalDriven_value, &b_HEEP_cutflow50_25ns_EcalDriven_value);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_dEtaIn", &HEEP_cutflow50_25ns_dEtaIn, &b_HEEP_cutflow50_25ns_dEtaIn);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_dEtaIn_n", &HEEP_cutflow50_25ns_dEtaIn_n, &b_HEEP_cutflow50_25ns_dEtaIn_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_dEtaIn_nCumulative", &HEEP_cutflow50_25ns_dEtaIn_nCumulative, &b_HEEP_cutflow50_25ns_dEtaIn_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_dEtaIn_value", &HEEP_cutflow50_25ns_dEtaIn_value, &b_HEEP_cutflow50_25ns_dEtaIn_value);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_dPhiIn", &HEEP_cutflow50_25ns_dPhiIn, &b_HEEP_cutflow50_25ns_dPhiIn);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_dPhiIn_n", &HEEP_cutflow50_25ns_dPhiIn_n, &b_HEEP_cutflow50_25ns_dPhiIn_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_dPhiIn_nCumulative", &HEEP_cutflow50_25ns_dPhiIn_nCumulative, &b_HEEP_cutflow50_25ns_dPhiIn_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_dPhiIn_value", &HEEP_cutflow50_25ns_dPhiIn_value, &b_HEEP_cutflow50_25ns_dPhiIn_value);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_HOverE", &HEEP_cutflow50_25ns_HOverE, &b_HEEP_cutflow50_25ns_HOverE);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_HOverE_n", &HEEP_cutflow50_25ns_HOverE_n, &b_HEEP_cutflow50_25ns_HOverE_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_HOverE_nCumulative", &HEEP_cutflow50_25ns_HOverE_nCumulative, &b_HEEP_cutflow50_25ns_HOverE_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_HOverE_value", &HEEP_cutflow50_25ns_HOverE_value, &b_HEEP_cutflow50_25ns_HOverE_value);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_SigmaIetaIeta", &HEEP_cutflow50_25ns_SigmaIetaIeta, &b_HEEP_cutflow50_25ns_SigmaIetaIeta);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_SigmaIetaIeta_n", &HEEP_cutflow50_25ns_SigmaIetaIeta_n, &b_HEEP_cutflow50_25ns_SigmaIetaIeta_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_SigmaIetaIeta_nCumulative", &HEEP_cutflow50_25ns_SigmaIetaIeta_nCumulative, &b_HEEP_cutflow50_25ns_SigmaIetaIeta_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_SigmaIetaIeta_value", &HEEP_cutflow50_25ns_SigmaIetaIeta_value, &b_HEEP_cutflow50_25ns_SigmaIetaIeta_value);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_E1x5OverE5x5", &HEEP_cutflow50_25ns_E1x5OverE5x5, &b_HEEP_cutflow50_25ns_E1x5OverE5x5);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_E1x5OverE5x5_n", &HEEP_cutflow50_25ns_E1x5OverE5x5_n, &b_HEEP_cutflow50_25ns_E1x5OverE5x5_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_E1x5OverE5x5_nCumulative", &HEEP_cutflow50_25ns_E1x5OverE5x5_nCumulative, &b_HEEP_cutflow50_25ns_E1x5OverE5x5_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_E1x5OverE5x5_value", &HEEP_cutflow50_25ns_E1x5OverE5x5_value, &b_HEEP_cutflow50_25ns_E1x5OverE5x5_value);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_E2x5OverE5x5", &HEEP_cutflow50_25ns_E2x5OverE5x5, &b_HEEP_cutflow50_25ns_E2x5OverE5x5);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_E2x5OverE5x5_n", &HEEP_cutflow50_25ns_E2x5OverE5x5_n, &b_HEEP_cutflow50_25ns_E2x5OverE5x5_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_E2x5OverE5x5_nCumulative", &HEEP_cutflow50_25ns_E2x5OverE5x5_nCumulative, &b_HEEP_cutflow50_25ns_E2x5OverE5x5_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_E2x5OverE5x5_value", &HEEP_cutflow50_25ns_E2x5OverE5x5_value, &b_HEEP_cutflow50_25ns_E2x5OverE5x5_value);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_missingHits", &HEEP_cutflow50_25ns_missingHits, &b_HEEP_cutflow50_25ns_missingHits);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_missingHits_n", &HEEP_cutflow50_25ns_missingHits_n, &b_HEEP_cutflow50_25ns_missingHits_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_missingHits_nCumulative", &HEEP_cutflow50_25ns_missingHits_nCumulative, &b_HEEP_cutflow50_25ns_missingHits_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_missingHits_value", &HEEP_cutflow50_25ns_missingHits_value, &b_HEEP_cutflow50_25ns_missingHits_value);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_dxyFirstPV", &HEEP_cutflow50_25ns_dxyFirstPV, &b_HEEP_cutflow50_25ns_dxyFirstPV);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_dxyFirstPV_n", &HEEP_cutflow50_25ns_dxyFirstPV_n, &b_HEEP_cutflow50_25ns_dxyFirstPV_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_dxyFirstPV_nCumulative", &HEEP_cutflow50_25ns_dxyFirstPV_nCumulative, &b_HEEP_cutflow50_25ns_dxyFirstPV_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_dxyFirstPV_value", &HEEP_cutflow50_25ns_dxyFirstPV_value, &b_HEEP_cutflow50_25ns_dxyFirstPV_value);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_ID", &HEEP_cutflow50_25ns_ID, &b_HEEP_cutflow50_25ns_ID);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_ID_n", &HEEP_cutflow50_25ns_ID_n, &b_HEEP_cutflow50_25ns_ID_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_isolEMHadDepth1", &HEEP_cutflow50_25ns_isolEMHadDepth1, &b_HEEP_cutflow50_25ns_isolEMHadDepth1);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_isolEMHadDepth1_n", &HEEP_cutflow50_25ns_isolEMHadDepth1_n, &b_HEEP_cutflow50_25ns_isolEMHadDepth1_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_isolEMHadDepth1_nCumulative", &HEEP_cutflow50_25ns_isolEMHadDepth1_nCumulative, &b_HEEP_cutflow50_25ns_isolEMHadDepth1_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_isolEMHadDepth1_value", &HEEP_cutflow50_25ns_isolEMHadDepth1_value, &b_HEEP_cutflow50_25ns_isolEMHadDepth1_value);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_IsolPtTrks", &HEEP_cutflow50_25ns_IsolPtTrks, &b_HEEP_cutflow50_25ns_IsolPtTrks);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_IsolPtTrks_n", &HEEP_cutflow50_25ns_IsolPtTrks_n, &b_HEEP_cutflow50_25ns_IsolPtTrks_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_IsolPtTrks_nCumulative", &HEEP_cutflow50_25ns_IsolPtTrks_nCumulative, &b_HEEP_cutflow50_25ns_IsolPtTrks_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_IsolPtTrks_value", &HEEP_cutflow50_25ns_IsolPtTrks_value, &b_HEEP_cutflow50_25ns_IsolPtTrks_value);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_isolation", &HEEP_cutflow50_25ns_isolation, &b_HEEP_cutflow50_25ns_isolation);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_isolation_n", &HEEP_cutflow50_25ns_isolation_n, &b_HEEP_cutflow50_25ns_isolation_n);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_total", &HEEP_cutflow50_25ns_total, &b_HEEP_cutflow50_25ns_total);
   fChain->SetBranchAddress("HEEP_cutflow50_25ns_total_n", &HEEP_cutflow50_25ns_total_n, &b_HEEP_cutflow50_25ns_total_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_Et", &HEEP_cutflow50_50ns_Et, &b_HEEP_cutflow50_50ns_Et);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_Et_n", &HEEP_cutflow50_50ns_Et_n, &b_HEEP_cutflow50_50ns_Et_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_Et_nCumulative", &HEEP_cutflow50_50ns_Et_nCumulative, &b_HEEP_cutflow50_50ns_Et_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_Et_value", &HEEP_cutflow50_50ns_Et_value, &b_HEEP_cutflow50_50ns_Et_value);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_eta", &HEEP_cutflow50_50ns_eta, &b_HEEP_cutflow50_50ns_eta);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_eta_n", &HEEP_cutflow50_50ns_eta_n, &b_HEEP_cutflow50_50ns_eta_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_eta_nCumulative", &HEEP_cutflow50_50ns_eta_nCumulative, &b_HEEP_cutflow50_50ns_eta_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_eta_value", &HEEP_cutflow50_50ns_eta_value, &b_HEEP_cutflow50_50ns_eta_value);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_acceptance", &HEEP_cutflow50_50ns_acceptance, &b_HEEP_cutflow50_50ns_acceptance);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_acceptance_n", &HEEP_cutflow50_50ns_acceptance_n, &b_HEEP_cutflow50_50ns_acceptance_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_EcalDriven", &HEEP_cutflow50_50ns_EcalDriven, &b_HEEP_cutflow50_50ns_EcalDriven);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_EcalDriven_n", &HEEP_cutflow50_50ns_EcalDriven_n, &b_HEEP_cutflow50_50ns_EcalDriven_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_EcalDriven_nCumulative", &HEEP_cutflow50_50ns_EcalDriven_nCumulative, &b_HEEP_cutflow50_50ns_EcalDriven_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_EcalDriven_value", &HEEP_cutflow50_50ns_EcalDriven_value, &b_HEEP_cutflow50_50ns_EcalDriven_value);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_dEtaIn", &HEEP_cutflow50_50ns_dEtaIn, &b_HEEP_cutflow50_50ns_dEtaIn);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_dEtaIn_n", &HEEP_cutflow50_50ns_dEtaIn_n, &b_HEEP_cutflow50_50ns_dEtaIn_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_dEtaIn_nCumulative", &HEEP_cutflow50_50ns_dEtaIn_nCumulative, &b_HEEP_cutflow50_50ns_dEtaIn_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_dEtaIn_value", &HEEP_cutflow50_50ns_dEtaIn_value, &b_HEEP_cutflow50_50ns_dEtaIn_value);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_dPhiIn", &HEEP_cutflow50_50ns_dPhiIn, &b_HEEP_cutflow50_50ns_dPhiIn);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_dPhiIn_n", &HEEP_cutflow50_50ns_dPhiIn_n, &b_HEEP_cutflow50_50ns_dPhiIn_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_dPhiIn_nCumulative", &HEEP_cutflow50_50ns_dPhiIn_nCumulative, &b_HEEP_cutflow50_50ns_dPhiIn_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_dPhiIn_value", &HEEP_cutflow50_50ns_dPhiIn_value, &b_HEEP_cutflow50_50ns_dPhiIn_value);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_HOverE", &HEEP_cutflow50_50ns_HOverE, &b_HEEP_cutflow50_50ns_HOverE);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_HOverE_n", &HEEP_cutflow50_50ns_HOverE_n, &b_HEEP_cutflow50_50ns_HOverE_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_HOverE_nCumulative", &HEEP_cutflow50_50ns_HOverE_nCumulative, &b_HEEP_cutflow50_50ns_HOverE_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_HOverE_value", &HEEP_cutflow50_50ns_HOverE_value, &b_HEEP_cutflow50_50ns_HOverE_value);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_SigmaIetaIeta", &HEEP_cutflow50_50ns_SigmaIetaIeta, &b_HEEP_cutflow50_50ns_SigmaIetaIeta);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_SigmaIetaIeta_n", &HEEP_cutflow50_50ns_SigmaIetaIeta_n, &b_HEEP_cutflow50_50ns_SigmaIetaIeta_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_SigmaIetaIeta_nCumulative", &HEEP_cutflow50_50ns_SigmaIetaIeta_nCumulative, &b_HEEP_cutflow50_50ns_SigmaIetaIeta_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_SigmaIetaIeta_value", &HEEP_cutflow50_50ns_SigmaIetaIeta_value, &b_HEEP_cutflow50_50ns_SigmaIetaIeta_value);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_E1x5OverE5x5", &HEEP_cutflow50_50ns_E1x5OverE5x5, &b_HEEP_cutflow50_50ns_E1x5OverE5x5);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_E1x5OverE5x5_n", &HEEP_cutflow50_50ns_E1x5OverE5x5_n, &b_HEEP_cutflow50_50ns_E1x5OverE5x5_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_E1x5OverE5x5_nCumulative", &HEEP_cutflow50_50ns_E1x5OverE5x5_nCumulative, &b_HEEP_cutflow50_50ns_E1x5OverE5x5_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_E1x5OverE5x5_value", &HEEP_cutflow50_50ns_E1x5OverE5x5_value, &b_HEEP_cutflow50_50ns_E1x5OverE5x5_value);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_E2x5OverE5x5", &HEEP_cutflow50_50ns_E2x5OverE5x5, &b_HEEP_cutflow50_50ns_E2x5OverE5x5);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_E2x5OverE5x5_n", &HEEP_cutflow50_50ns_E2x5OverE5x5_n, &b_HEEP_cutflow50_50ns_E2x5OverE5x5_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_E2x5OverE5x5_nCumulative", &HEEP_cutflow50_50ns_E2x5OverE5x5_nCumulative, &b_HEEP_cutflow50_50ns_E2x5OverE5x5_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_E2x5OverE5x5_value", &HEEP_cutflow50_50ns_E2x5OverE5x5_value, &b_HEEP_cutflow50_50ns_E2x5OverE5x5_value);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_missingHits", &HEEP_cutflow50_50ns_missingHits, &b_HEEP_cutflow50_50ns_missingHits);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_missingHits_n", &HEEP_cutflow50_50ns_missingHits_n, &b_HEEP_cutflow50_50ns_missingHits_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_missingHits_nCumulative", &HEEP_cutflow50_50ns_missingHits_nCumulative, &b_HEEP_cutflow50_50ns_missingHits_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_missingHits_value", &HEEP_cutflow50_50ns_missingHits_value, &b_HEEP_cutflow50_50ns_missingHits_value);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_dxyFirstPV", &HEEP_cutflow50_50ns_dxyFirstPV, &b_HEEP_cutflow50_50ns_dxyFirstPV);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_dxyFirstPV_n", &HEEP_cutflow50_50ns_dxyFirstPV_n, &b_HEEP_cutflow50_50ns_dxyFirstPV_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_dxyFirstPV_nCumulative", &HEEP_cutflow50_50ns_dxyFirstPV_nCumulative, &b_HEEP_cutflow50_50ns_dxyFirstPV_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_dxyFirstPV_value", &HEEP_cutflow50_50ns_dxyFirstPV_value, &b_HEEP_cutflow50_50ns_dxyFirstPV_value);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_ID", &HEEP_cutflow50_50ns_ID, &b_HEEP_cutflow50_50ns_ID);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_ID_n", &HEEP_cutflow50_50ns_ID_n, &b_HEEP_cutflow50_50ns_ID_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_isolEMHadDepth1", &HEEP_cutflow50_50ns_isolEMHadDepth1, &b_HEEP_cutflow50_50ns_isolEMHadDepth1);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_isolEMHadDepth1_n", &HEEP_cutflow50_50ns_isolEMHadDepth1_n, &b_HEEP_cutflow50_50ns_isolEMHadDepth1_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_isolEMHadDepth1_nCumulative", &HEEP_cutflow50_50ns_isolEMHadDepth1_nCumulative, &b_HEEP_cutflow50_50ns_isolEMHadDepth1_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_isolEMHadDepth1_value", &HEEP_cutflow50_50ns_isolEMHadDepth1_value, &b_HEEP_cutflow50_50ns_isolEMHadDepth1_value);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_IsolPtTrks", &HEEP_cutflow50_50ns_IsolPtTrks, &b_HEEP_cutflow50_50ns_IsolPtTrks);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_IsolPtTrks_n", &HEEP_cutflow50_50ns_IsolPtTrks_n, &b_HEEP_cutflow50_50ns_IsolPtTrks_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_IsolPtTrks_nCumulative", &HEEP_cutflow50_50ns_IsolPtTrks_nCumulative, &b_HEEP_cutflow50_50ns_IsolPtTrks_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_IsolPtTrks_value", &HEEP_cutflow50_50ns_IsolPtTrks_value, &b_HEEP_cutflow50_50ns_IsolPtTrks_value);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_isolation", &HEEP_cutflow50_50ns_isolation, &b_HEEP_cutflow50_50ns_isolation);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_isolation_n", &HEEP_cutflow50_50ns_isolation_n, &b_HEEP_cutflow50_50ns_isolation_n);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_total", &HEEP_cutflow50_50ns_total, &b_HEEP_cutflow50_50ns_total);
   fChain->SetBranchAddress("HEEP_cutflow50_50ns_total_n", &HEEP_cutflow50_50ns_total_n, &b_HEEP_cutflow50_50ns_total_n);
   fChain->SetBranchAddress("HEEP_cutflow50_Et", &HEEP_cutflow50_Et, &b_HEEP_cutflow50_Et);
   fChain->SetBranchAddress("HEEP_cutflow50_Et_n", &HEEP_cutflow50_Et_n, &b_HEEP_cutflow50_Et_n);
   fChain->SetBranchAddress("HEEP_cutflow50_Et_nCumulative", &HEEP_cutflow50_Et_nCumulative, &b_HEEP_cutflow50_Et_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_Et_value", &HEEP_cutflow50_Et_value, &b_HEEP_cutflow50_Et_value);
   fChain->SetBranchAddress("HEEP_cutflow50_eta", &HEEP_cutflow50_eta, &b_HEEP_cutflow50_eta);
   fChain->SetBranchAddress("HEEP_cutflow50_eta_n", &HEEP_cutflow50_eta_n, &b_HEEP_cutflow50_eta_n);
   fChain->SetBranchAddress("HEEP_cutflow50_eta_nCumulative", &HEEP_cutflow50_eta_nCumulative, &b_HEEP_cutflow50_eta_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_eta_value", &HEEP_cutflow50_eta_value, &b_HEEP_cutflow50_eta_value);
   fChain->SetBranchAddress("HEEP_cutflow50_acceptance", &HEEP_cutflow50_acceptance, &b_HEEP_cutflow50_acceptance);
   fChain->SetBranchAddress("HEEP_cutflow50_acceptance_n", &HEEP_cutflow50_acceptance_n, &b_HEEP_cutflow50_acceptance_n);
   fChain->SetBranchAddress("HEEP_cutflow50_EcalDriven", &HEEP_cutflow50_EcalDriven, &b_HEEP_cutflow50_EcalDriven);
   fChain->SetBranchAddress("HEEP_cutflow50_EcalDriven_n", &HEEP_cutflow50_EcalDriven_n, &b_HEEP_cutflow50_EcalDriven_n);
   fChain->SetBranchAddress("HEEP_cutflow50_EcalDriven_nCumulative", &HEEP_cutflow50_EcalDriven_nCumulative, &b_HEEP_cutflow50_EcalDriven_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_EcalDriven_value", &HEEP_cutflow50_EcalDriven_value, &b_HEEP_cutflow50_EcalDriven_value);
   fChain->SetBranchAddress("HEEP_cutflow50_dEtaIn", &HEEP_cutflow50_dEtaIn, &b_HEEP_cutflow50_dEtaIn);
   fChain->SetBranchAddress("HEEP_cutflow50_dEtaIn_n", &HEEP_cutflow50_dEtaIn_n, &b_HEEP_cutflow50_dEtaIn_n);
   fChain->SetBranchAddress("HEEP_cutflow50_dEtaIn_nCumulative", &HEEP_cutflow50_dEtaIn_nCumulative, &b_HEEP_cutflow50_dEtaIn_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_dEtaIn_value", &HEEP_cutflow50_dEtaIn_value, &b_HEEP_cutflow50_dEtaIn_value);
   fChain->SetBranchAddress("HEEP_cutflow50_dPhiIn", &HEEP_cutflow50_dPhiIn, &b_HEEP_cutflow50_dPhiIn);
   fChain->SetBranchAddress("HEEP_cutflow50_dPhiIn_n", &HEEP_cutflow50_dPhiIn_n, &b_HEEP_cutflow50_dPhiIn_n);
   fChain->SetBranchAddress("HEEP_cutflow50_dPhiIn_nCumulative", &HEEP_cutflow50_dPhiIn_nCumulative, &b_HEEP_cutflow50_dPhiIn_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_dPhiIn_value", &HEEP_cutflow50_dPhiIn_value, &b_HEEP_cutflow50_dPhiIn_value);
   fChain->SetBranchAddress("HEEP_cutflow50_HOverE", &HEEP_cutflow50_HOverE, &b_HEEP_cutflow50_HOverE);
   fChain->SetBranchAddress("HEEP_cutflow50_HOverE_n", &HEEP_cutflow50_HOverE_n, &b_HEEP_cutflow50_HOverE_n);
   fChain->SetBranchAddress("HEEP_cutflow50_HOverE_nCumulative", &HEEP_cutflow50_HOverE_nCumulative, &b_HEEP_cutflow50_HOverE_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_HOverE_value", &HEEP_cutflow50_HOverE_value, &b_HEEP_cutflow50_HOverE_value);
   fChain->SetBranchAddress("HEEP_cutflow50_SigmaIetaIeta", &HEEP_cutflow50_SigmaIetaIeta, &b_HEEP_cutflow50_SigmaIetaIeta);
   fChain->SetBranchAddress("HEEP_cutflow50_SigmaIetaIeta_n", &HEEP_cutflow50_SigmaIetaIeta_n, &b_HEEP_cutflow50_SigmaIetaIeta_n);
   fChain->SetBranchAddress("HEEP_cutflow50_SigmaIetaIeta_nCumulative", &HEEP_cutflow50_SigmaIetaIeta_nCumulative, &b_HEEP_cutflow50_SigmaIetaIeta_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_SigmaIetaIeta_value", &HEEP_cutflow50_SigmaIetaIeta_value, &b_HEEP_cutflow50_SigmaIetaIeta_value);
   fChain->SetBranchAddress("HEEP_cutflow50_E1x5OverE5x5", &HEEP_cutflow50_E1x5OverE5x5, &b_HEEP_cutflow50_E1x5OverE5x5);
   fChain->SetBranchAddress("HEEP_cutflow50_E1x5OverE5x5_n", &HEEP_cutflow50_E1x5OverE5x5_n, &b_HEEP_cutflow50_E1x5OverE5x5_n);
   fChain->SetBranchAddress("HEEP_cutflow50_E1x5OverE5x5_nCumulative", &HEEP_cutflow50_E1x5OverE5x5_nCumulative, &b_HEEP_cutflow50_E1x5OverE5x5_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_E1x5OverE5x5_value", &HEEP_cutflow50_E1x5OverE5x5_value, &b_HEEP_cutflow50_E1x5OverE5x5_value);
   fChain->SetBranchAddress("HEEP_cutflow50_E2x5OverE5x5", &HEEP_cutflow50_E2x5OverE5x5, &b_HEEP_cutflow50_E2x5OverE5x5);
   fChain->SetBranchAddress("HEEP_cutflow50_E2x5OverE5x5_n", &HEEP_cutflow50_E2x5OverE5x5_n, &b_HEEP_cutflow50_E2x5OverE5x5_n);
   fChain->SetBranchAddress("HEEP_cutflow50_E2x5OverE5x5_nCumulative", &HEEP_cutflow50_E2x5OverE5x5_nCumulative, &b_HEEP_cutflow50_E2x5OverE5x5_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_E2x5OverE5x5_value", &HEEP_cutflow50_E2x5OverE5x5_value, &b_HEEP_cutflow50_E2x5OverE5x5_value);
   fChain->SetBranchAddress("HEEP_cutflow50_missingHits", &HEEP_cutflow50_missingHits, &b_HEEP_cutflow50_missingHits);
   fChain->SetBranchAddress("HEEP_cutflow50_missingHits_n", &HEEP_cutflow50_missingHits_n, &b_HEEP_cutflow50_missingHits_n);
   fChain->SetBranchAddress("HEEP_cutflow50_missingHits_nCumulative", &HEEP_cutflow50_missingHits_nCumulative, &b_HEEP_cutflow50_missingHits_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_missingHits_value", &HEEP_cutflow50_missingHits_value, &b_HEEP_cutflow50_missingHits_value);
   fChain->SetBranchAddress("HEEP_cutflow50_dxyFirstPV", &HEEP_cutflow50_dxyFirstPV, &b_HEEP_cutflow50_dxyFirstPV);
   fChain->SetBranchAddress("HEEP_cutflow50_dxyFirstPV_n", &HEEP_cutflow50_dxyFirstPV_n, &b_HEEP_cutflow50_dxyFirstPV_n);
   fChain->SetBranchAddress("HEEP_cutflow50_dxyFirstPV_nCumulative", &HEEP_cutflow50_dxyFirstPV_nCumulative, &b_HEEP_cutflow50_dxyFirstPV_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_dxyFirstPV_value", &HEEP_cutflow50_dxyFirstPV_value, &b_HEEP_cutflow50_dxyFirstPV_value);
   fChain->SetBranchAddress("HEEP_cutflow50_ID", &HEEP_cutflow50_ID, &b_HEEP_cutflow50_ID);
   fChain->SetBranchAddress("HEEP_cutflow50_ID_n", &HEEP_cutflow50_ID_n, &b_HEEP_cutflow50_ID_n);
   fChain->SetBranchAddress("HEEP_cutflow50_isolEMHadDepth1", &HEEP_cutflow50_isolEMHadDepth1, &b_HEEP_cutflow50_isolEMHadDepth1);
   fChain->SetBranchAddress("HEEP_cutflow50_isolEMHadDepth1_n", &HEEP_cutflow50_isolEMHadDepth1_n, &b_HEEP_cutflow50_isolEMHadDepth1_n);
   fChain->SetBranchAddress("HEEP_cutflow50_isolEMHadDepth1_nCumulative", &HEEP_cutflow50_isolEMHadDepth1_nCumulative, &b_HEEP_cutflow50_isolEMHadDepth1_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_isolEMHadDepth1_value", &HEEP_cutflow50_isolEMHadDepth1_value, &b_HEEP_cutflow50_isolEMHadDepth1_value);
   fChain->SetBranchAddress("HEEP_cutflow50_IsolPtTrks", &HEEP_cutflow50_IsolPtTrks, &b_HEEP_cutflow50_IsolPtTrks);
   fChain->SetBranchAddress("HEEP_cutflow50_IsolPtTrks_n", &HEEP_cutflow50_IsolPtTrks_n, &b_HEEP_cutflow50_IsolPtTrks_n);
   fChain->SetBranchAddress("HEEP_cutflow50_IsolPtTrks_nCumulative", &HEEP_cutflow50_IsolPtTrks_nCumulative, &b_HEEP_cutflow50_IsolPtTrks_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow50_IsolPtTrks_value", &HEEP_cutflow50_IsolPtTrks_value, &b_HEEP_cutflow50_IsolPtTrks_value);
   fChain->SetBranchAddress("HEEP_cutflow50_isolation", &HEEP_cutflow50_isolation, &b_HEEP_cutflow50_isolation);
   fChain->SetBranchAddress("HEEP_cutflow50_isolation_n", &HEEP_cutflow50_isolation_n, &b_HEEP_cutflow50_isolation_n);
   fChain->SetBranchAddress("HEEP_cutflow50_total", &HEEP_cutflow50_total, &b_HEEP_cutflow50_total);
   fChain->SetBranchAddress("HEEP_cutflow50_total_n", &HEEP_cutflow50_total_n, &b_HEEP_cutflow50_total_n);
   fChain->SetBranchAddress("HEEP_cutflow51_Et", &HEEP_cutflow51_Et, &b_HEEP_cutflow51_Et);
   fChain->SetBranchAddress("HEEP_cutflow51_Et_n", &HEEP_cutflow51_Et_n, &b_HEEP_cutflow51_Et_n);
   fChain->SetBranchAddress("HEEP_cutflow51_Et_nCumulative", &HEEP_cutflow51_Et_nCumulative, &b_HEEP_cutflow51_Et_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow51_Et_value", &HEEP_cutflow51_Et_value, &b_HEEP_cutflow51_Et_value);
   fChain->SetBranchAddress("HEEP_cutflow51_eta", &HEEP_cutflow51_eta, &b_HEEP_cutflow51_eta);
   fChain->SetBranchAddress("HEEP_cutflow51_eta_n", &HEEP_cutflow51_eta_n, &b_HEEP_cutflow51_eta_n);
   fChain->SetBranchAddress("HEEP_cutflow51_eta_nCumulative", &HEEP_cutflow51_eta_nCumulative, &b_HEEP_cutflow51_eta_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow51_eta_value", &HEEP_cutflow51_eta_value, &b_HEEP_cutflow51_eta_value);
   fChain->SetBranchAddress("HEEP_cutflow51_acceptance", &HEEP_cutflow51_acceptance, &b_HEEP_cutflow51_acceptance);
   fChain->SetBranchAddress("HEEP_cutflow51_acceptance_n", &HEEP_cutflow51_acceptance_n, &b_HEEP_cutflow51_acceptance_n);
   fChain->SetBranchAddress("HEEP_cutflow51_EcalDriven", &HEEP_cutflow51_EcalDriven, &b_HEEP_cutflow51_EcalDriven);
   fChain->SetBranchAddress("HEEP_cutflow51_EcalDriven_n", &HEEP_cutflow51_EcalDriven_n, &b_HEEP_cutflow51_EcalDriven_n);
   fChain->SetBranchAddress("HEEP_cutflow51_EcalDriven_nCumulative", &HEEP_cutflow51_EcalDriven_nCumulative, &b_HEEP_cutflow51_EcalDriven_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow51_EcalDriven_value", &HEEP_cutflow51_EcalDriven_value, &b_HEEP_cutflow51_EcalDriven_value);
   fChain->SetBranchAddress("HEEP_cutflow51_dEtaIn", &HEEP_cutflow51_dEtaIn, &b_HEEP_cutflow51_dEtaIn);
   fChain->SetBranchAddress("HEEP_cutflow51_dEtaIn_n", &HEEP_cutflow51_dEtaIn_n, &b_HEEP_cutflow51_dEtaIn_n);
   fChain->SetBranchAddress("HEEP_cutflow51_dEtaIn_nCumulative", &HEEP_cutflow51_dEtaIn_nCumulative, &b_HEEP_cutflow51_dEtaIn_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow51_dEtaIn_value", &HEEP_cutflow51_dEtaIn_value, &b_HEEP_cutflow51_dEtaIn_value);
   fChain->SetBranchAddress("HEEP_cutflow51_dPhiIn", &HEEP_cutflow51_dPhiIn, &b_HEEP_cutflow51_dPhiIn);
   fChain->SetBranchAddress("HEEP_cutflow51_dPhiIn_n", &HEEP_cutflow51_dPhiIn_n, &b_HEEP_cutflow51_dPhiIn_n);
   fChain->SetBranchAddress("HEEP_cutflow51_dPhiIn_nCumulative", &HEEP_cutflow51_dPhiIn_nCumulative, &b_HEEP_cutflow51_dPhiIn_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow51_dPhiIn_value", &HEEP_cutflow51_dPhiIn_value, &b_HEEP_cutflow51_dPhiIn_value);
   fChain->SetBranchAddress("HEEP_cutflow51_HOverE", &HEEP_cutflow51_HOverE, &b_HEEP_cutflow51_HOverE);
   fChain->SetBranchAddress("HEEP_cutflow51_HOverE_n", &HEEP_cutflow51_HOverE_n, &b_HEEP_cutflow51_HOverE_n);
   fChain->SetBranchAddress("HEEP_cutflow51_HOverE_nCumulative", &HEEP_cutflow51_HOverE_nCumulative, &b_HEEP_cutflow51_HOverE_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow51_HOverE_value", &HEEP_cutflow51_HOverE_value, &b_HEEP_cutflow51_HOverE_value);
   fChain->SetBranchAddress("HEEP_cutflow51_SigmaIetaIeta", &HEEP_cutflow51_SigmaIetaIeta, &b_HEEP_cutflow51_SigmaIetaIeta);
   fChain->SetBranchAddress("HEEP_cutflow51_SigmaIetaIeta_n", &HEEP_cutflow51_SigmaIetaIeta_n, &b_HEEP_cutflow51_SigmaIetaIeta_n);
   fChain->SetBranchAddress("HEEP_cutflow51_SigmaIetaIeta_nCumulative", &HEEP_cutflow51_SigmaIetaIeta_nCumulative, &b_HEEP_cutflow51_SigmaIetaIeta_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow51_SigmaIetaIeta_value", &HEEP_cutflow51_SigmaIetaIeta_value, &b_HEEP_cutflow51_SigmaIetaIeta_value);
   fChain->SetBranchAddress("HEEP_cutflow51_E1x5OverE5x5", &HEEP_cutflow51_E1x5OverE5x5, &b_HEEP_cutflow51_E1x5OverE5x5);
   fChain->SetBranchAddress("HEEP_cutflow51_E1x5OverE5x5_n", &HEEP_cutflow51_E1x5OverE5x5_n, &b_HEEP_cutflow51_E1x5OverE5x5_n);
   fChain->SetBranchAddress("HEEP_cutflow51_E1x5OverE5x5_nCumulative", &HEEP_cutflow51_E1x5OverE5x5_nCumulative, &b_HEEP_cutflow51_E1x5OverE5x5_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow51_E1x5OverE5x5_value", &HEEP_cutflow51_E1x5OverE5x5_value, &b_HEEP_cutflow51_E1x5OverE5x5_value);
   fChain->SetBranchAddress("HEEP_cutflow51_E2x5OverE5x5", &HEEP_cutflow51_E2x5OverE5x5, &b_HEEP_cutflow51_E2x5OverE5x5);
   fChain->SetBranchAddress("HEEP_cutflow51_E2x5OverE5x5_n", &HEEP_cutflow51_E2x5OverE5x5_n, &b_HEEP_cutflow51_E2x5OverE5x5_n);
   fChain->SetBranchAddress("HEEP_cutflow51_E2x5OverE5x5_nCumulative", &HEEP_cutflow51_E2x5OverE5x5_nCumulative, &b_HEEP_cutflow51_E2x5OverE5x5_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow51_E2x5OverE5x5_value", &HEEP_cutflow51_E2x5OverE5x5_value, &b_HEEP_cutflow51_E2x5OverE5x5_value);
   fChain->SetBranchAddress("HEEP_cutflow51_missingHits", &HEEP_cutflow51_missingHits, &b_HEEP_cutflow51_missingHits);
   fChain->SetBranchAddress("HEEP_cutflow51_missingHits_n", &HEEP_cutflow51_missingHits_n, &b_HEEP_cutflow51_missingHits_n);
   fChain->SetBranchAddress("HEEP_cutflow51_missingHits_nCumulative", &HEEP_cutflow51_missingHits_nCumulative, &b_HEEP_cutflow51_missingHits_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow51_missingHits_value", &HEEP_cutflow51_missingHits_value, &b_HEEP_cutflow51_missingHits_value);
   fChain->SetBranchAddress("HEEP_cutflow51_dxyFirstPV", &HEEP_cutflow51_dxyFirstPV, &b_HEEP_cutflow51_dxyFirstPV);
   fChain->SetBranchAddress("HEEP_cutflow51_dxyFirstPV_n", &HEEP_cutflow51_dxyFirstPV_n, &b_HEEP_cutflow51_dxyFirstPV_n);
   fChain->SetBranchAddress("HEEP_cutflow51_dxyFirstPV_nCumulative", &HEEP_cutflow51_dxyFirstPV_nCumulative, &b_HEEP_cutflow51_dxyFirstPV_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow51_dxyFirstPV_value", &HEEP_cutflow51_dxyFirstPV_value, &b_HEEP_cutflow51_dxyFirstPV_value);
   fChain->SetBranchAddress("HEEP_cutflow51_ID", &HEEP_cutflow51_ID, &b_HEEP_cutflow51_ID);
   fChain->SetBranchAddress("HEEP_cutflow51_ID_n", &HEEP_cutflow51_ID_n, &b_HEEP_cutflow51_ID_n);
   fChain->SetBranchAddress("HEEP_cutflow51_isolEMHadDepth1", &HEEP_cutflow51_isolEMHadDepth1, &b_HEEP_cutflow51_isolEMHadDepth1);
   fChain->SetBranchAddress("HEEP_cutflow51_isolEMHadDepth1_n", &HEEP_cutflow51_isolEMHadDepth1_n, &b_HEEP_cutflow51_isolEMHadDepth1_n);
   fChain->SetBranchAddress("HEEP_cutflow51_isolEMHadDepth1_nCumulative", &HEEP_cutflow51_isolEMHadDepth1_nCumulative, &b_HEEP_cutflow51_isolEMHadDepth1_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow51_isolEMHadDepth1_value", &HEEP_cutflow51_isolEMHadDepth1_value, &b_HEEP_cutflow51_isolEMHadDepth1_value);
   fChain->SetBranchAddress("HEEP_cutflow51_IsolPtTrks", &HEEP_cutflow51_IsolPtTrks, &b_HEEP_cutflow51_IsolPtTrks);
   fChain->SetBranchAddress("HEEP_cutflow51_IsolPtTrks_n", &HEEP_cutflow51_IsolPtTrks_n, &b_HEEP_cutflow51_IsolPtTrks_n);
   fChain->SetBranchAddress("HEEP_cutflow51_IsolPtTrks_nCumulative", &HEEP_cutflow51_IsolPtTrks_nCumulative, &b_HEEP_cutflow51_IsolPtTrks_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow51_IsolPtTrks_value", &HEEP_cutflow51_IsolPtTrks_value, &b_HEEP_cutflow51_IsolPtTrks_value);
   fChain->SetBranchAddress("HEEP_cutflow51_isolation", &HEEP_cutflow51_isolation, &b_HEEP_cutflow51_isolation);
   fChain->SetBranchAddress("HEEP_cutflow51_isolation_n", &HEEP_cutflow51_isolation_n, &b_HEEP_cutflow51_isolation_n);
   fChain->SetBranchAddress("HEEP_cutflow51_total", &HEEP_cutflow51_total, &b_HEEP_cutflow51_total);
   fChain->SetBranchAddress("HEEP_cutflow51_total_n", &HEEP_cutflow51_total_n, &b_HEEP_cutflow51_total_n);
   fChain->SetBranchAddress("HEEP_cutflow60_Et", &HEEP_cutflow60_Et, &b_HEEP_cutflow60_Et);
   fChain->SetBranchAddress("HEEP_cutflow60_Et_n", &HEEP_cutflow60_Et_n, &b_HEEP_cutflow60_Et_n);
   fChain->SetBranchAddress("HEEP_cutflow60_Et_nCumulative", &HEEP_cutflow60_Et_nCumulative, &b_HEEP_cutflow60_Et_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow60_Et_value", &HEEP_cutflow60_Et_value, &b_HEEP_cutflow60_Et_value);
   fChain->SetBranchAddress("HEEP_cutflow60_eta", &HEEP_cutflow60_eta, &b_HEEP_cutflow60_eta);
   fChain->SetBranchAddress("HEEP_cutflow60_eta_n", &HEEP_cutflow60_eta_n, &b_HEEP_cutflow60_eta_n);
   fChain->SetBranchAddress("HEEP_cutflow60_eta_nCumulative", &HEEP_cutflow60_eta_nCumulative, &b_HEEP_cutflow60_eta_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow60_eta_value", &HEEP_cutflow60_eta_value, &b_HEEP_cutflow60_eta_value);
   fChain->SetBranchAddress("HEEP_cutflow60_acceptance", &HEEP_cutflow60_acceptance, &b_HEEP_cutflow60_acceptance);
   fChain->SetBranchAddress("HEEP_cutflow60_acceptance_n", &HEEP_cutflow60_acceptance_n, &b_HEEP_cutflow60_acceptance_n);
   fChain->SetBranchAddress("HEEP_cutflow60_EcalDriven", &HEEP_cutflow60_EcalDriven, &b_HEEP_cutflow60_EcalDriven);
   fChain->SetBranchAddress("HEEP_cutflow60_EcalDriven_n", &HEEP_cutflow60_EcalDriven_n, &b_HEEP_cutflow60_EcalDriven_n);
   fChain->SetBranchAddress("HEEP_cutflow60_EcalDriven_nCumulative", &HEEP_cutflow60_EcalDriven_nCumulative, &b_HEEP_cutflow60_EcalDriven_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow60_EcalDriven_value", &HEEP_cutflow60_EcalDriven_value, &b_HEEP_cutflow60_EcalDriven_value);
   fChain->SetBranchAddress("HEEP_cutflow60_dEtaIn", &HEEP_cutflow60_dEtaIn, &b_HEEP_cutflow60_dEtaIn);
   fChain->SetBranchAddress("HEEP_cutflow60_dEtaIn_n", &HEEP_cutflow60_dEtaIn_n, &b_HEEP_cutflow60_dEtaIn_n);
   fChain->SetBranchAddress("HEEP_cutflow60_dEtaIn_nCumulative", &HEEP_cutflow60_dEtaIn_nCumulative, &b_HEEP_cutflow60_dEtaIn_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow60_dEtaIn_value", &HEEP_cutflow60_dEtaIn_value, &b_HEEP_cutflow60_dEtaIn_value);
   fChain->SetBranchAddress("HEEP_cutflow60_dPhiIn", &HEEP_cutflow60_dPhiIn, &b_HEEP_cutflow60_dPhiIn);
   fChain->SetBranchAddress("HEEP_cutflow60_dPhiIn_n", &HEEP_cutflow60_dPhiIn_n, &b_HEEP_cutflow60_dPhiIn_n);
   fChain->SetBranchAddress("HEEP_cutflow60_dPhiIn_nCumulative", &HEEP_cutflow60_dPhiIn_nCumulative, &b_HEEP_cutflow60_dPhiIn_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow60_dPhiIn_value", &HEEP_cutflow60_dPhiIn_value, &b_HEEP_cutflow60_dPhiIn_value);
   fChain->SetBranchAddress("HEEP_cutflow60_HOverE", &HEEP_cutflow60_HOverE, &b_HEEP_cutflow60_HOverE);
   fChain->SetBranchAddress("HEEP_cutflow60_HOverE_n", &HEEP_cutflow60_HOverE_n, &b_HEEP_cutflow60_HOverE_n);
   fChain->SetBranchAddress("HEEP_cutflow60_HOverE_nCumulative", &HEEP_cutflow60_HOverE_nCumulative, &b_HEEP_cutflow60_HOverE_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow60_HOverE_value", &HEEP_cutflow60_HOverE_value, &b_HEEP_cutflow60_HOverE_value);
   fChain->SetBranchAddress("HEEP_cutflow60_SigmaIetaIeta", &HEEP_cutflow60_SigmaIetaIeta, &b_HEEP_cutflow60_SigmaIetaIeta);
   fChain->SetBranchAddress("HEEP_cutflow60_SigmaIetaIeta_n", &HEEP_cutflow60_SigmaIetaIeta_n, &b_HEEP_cutflow60_SigmaIetaIeta_n);
   fChain->SetBranchAddress("HEEP_cutflow60_SigmaIetaIeta_nCumulative", &HEEP_cutflow60_SigmaIetaIeta_nCumulative, &b_HEEP_cutflow60_SigmaIetaIeta_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow60_SigmaIetaIeta_value", &HEEP_cutflow60_SigmaIetaIeta_value, &b_HEEP_cutflow60_SigmaIetaIeta_value);
   fChain->SetBranchAddress("HEEP_cutflow60_E1x5OverE5x5", &HEEP_cutflow60_E1x5OverE5x5, &b_HEEP_cutflow60_E1x5OverE5x5);
   fChain->SetBranchAddress("HEEP_cutflow60_E1x5OverE5x5_n", &HEEP_cutflow60_E1x5OverE5x5_n, &b_HEEP_cutflow60_E1x5OverE5x5_n);
   fChain->SetBranchAddress("HEEP_cutflow60_E1x5OverE5x5_nCumulative", &HEEP_cutflow60_E1x5OverE5x5_nCumulative, &b_HEEP_cutflow60_E1x5OverE5x5_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow60_E1x5OverE5x5_value", &HEEP_cutflow60_E1x5OverE5x5_value, &b_HEEP_cutflow60_E1x5OverE5x5_value);
   fChain->SetBranchAddress("HEEP_cutflow60_E2x5OverE5x5", &HEEP_cutflow60_E2x5OverE5x5, &b_HEEP_cutflow60_E2x5OverE5x5);
   fChain->SetBranchAddress("HEEP_cutflow60_E2x5OverE5x5_n", &HEEP_cutflow60_E2x5OverE5x5_n, &b_HEEP_cutflow60_E2x5OverE5x5_n);
   fChain->SetBranchAddress("HEEP_cutflow60_E2x5OverE5x5_nCumulative", &HEEP_cutflow60_E2x5OverE5x5_nCumulative, &b_HEEP_cutflow60_E2x5OverE5x5_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow60_E2x5OverE5x5_value", &HEEP_cutflow60_E2x5OverE5x5_value, &b_HEEP_cutflow60_E2x5OverE5x5_value);
   fChain->SetBranchAddress("HEEP_cutflow60_missingHits", &HEEP_cutflow60_missingHits, &b_HEEP_cutflow60_missingHits);
   fChain->SetBranchAddress("HEEP_cutflow60_missingHits_n", &HEEP_cutflow60_missingHits_n, &b_HEEP_cutflow60_missingHits_n);
   fChain->SetBranchAddress("HEEP_cutflow60_missingHits_nCumulative", &HEEP_cutflow60_missingHits_nCumulative, &b_HEEP_cutflow60_missingHits_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow60_missingHits_value", &HEEP_cutflow60_missingHits_value, &b_HEEP_cutflow60_missingHits_value);
   fChain->SetBranchAddress("HEEP_cutflow60_dxyFirstPV", &HEEP_cutflow60_dxyFirstPV, &b_HEEP_cutflow60_dxyFirstPV);
   fChain->SetBranchAddress("HEEP_cutflow60_dxyFirstPV_n", &HEEP_cutflow60_dxyFirstPV_n, &b_HEEP_cutflow60_dxyFirstPV_n);
   fChain->SetBranchAddress("HEEP_cutflow60_dxyFirstPV_nCumulative", &HEEP_cutflow60_dxyFirstPV_nCumulative, &b_HEEP_cutflow60_dxyFirstPV_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow60_dxyFirstPV_value", &HEEP_cutflow60_dxyFirstPV_value, &b_HEEP_cutflow60_dxyFirstPV_value);
   fChain->SetBranchAddress("HEEP_cutflow60_ID", &HEEP_cutflow60_ID, &b_HEEP_cutflow60_ID);
   fChain->SetBranchAddress("HEEP_cutflow60_ID_n", &HEEP_cutflow60_ID_n, &b_HEEP_cutflow60_ID_n);
   fChain->SetBranchAddress("HEEP_cutflow60_isolEMHadDepth1", &HEEP_cutflow60_isolEMHadDepth1, &b_HEEP_cutflow60_isolEMHadDepth1);
   fChain->SetBranchAddress("HEEP_cutflow60_isolEMHadDepth1_n", &HEEP_cutflow60_isolEMHadDepth1_n, &b_HEEP_cutflow60_isolEMHadDepth1_n);
   fChain->SetBranchAddress("HEEP_cutflow60_isolEMHadDepth1_nCumulative", &HEEP_cutflow60_isolEMHadDepth1_nCumulative, &b_HEEP_cutflow60_isolEMHadDepth1_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow60_isolEMHadDepth1_value", &HEEP_cutflow60_isolEMHadDepth1_value, &b_HEEP_cutflow60_isolEMHadDepth1_value);
   fChain->SetBranchAddress("HEEP_cutflow60_IsolPtTrks", &HEEP_cutflow60_IsolPtTrks, &b_HEEP_cutflow60_IsolPtTrks);
   fChain->SetBranchAddress("HEEP_cutflow60_IsolPtTrks_n", &HEEP_cutflow60_IsolPtTrks_n, &b_HEEP_cutflow60_IsolPtTrks_n);
   fChain->SetBranchAddress("HEEP_cutflow60_IsolPtTrks_nCumulative", &HEEP_cutflow60_IsolPtTrks_nCumulative, &b_HEEP_cutflow60_IsolPtTrks_nCumulative);
   fChain->SetBranchAddress("HEEP_cutflow60_IsolPtTrks_value", &HEEP_cutflow60_IsolPtTrks_value, &b_HEEP_cutflow60_IsolPtTrks_value);
   fChain->SetBranchAddress("HEEP_cutflow60_isolation", &HEEP_cutflow60_isolation, &b_HEEP_cutflow60_isolation);
   fChain->SetBranchAddress("HEEP_cutflow60_isolation_n", &HEEP_cutflow60_isolation_n, &b_HEEP_cutflow60_isolation_n);
   fChain->SetBranchAddress("HEEP_cutflow60_total", &HEEP_cutflow60_total, &b_HEEP_cutflow60_total);
   fChain->SetBranchAddress("HEEP_cutflow60_total_n", &HEEP_cutflow60_total_n, &b_HEEP_cutflow60_total_n);
   fChain->SetBranchAddress("Zee_n", &Zee_n, &b_Zee_n);
   fChain->SetBranchAddress("Zee_mass", &Zee_mass, &b_Zee_mass);
   fChain->SetBranchAddress("Zee_mass_HEEP", &Zee_mass_HEEP, &b_Zee_mass_HEEP);
   fChain->SetBranchAddress("Zee_i1", &Zee_i1, &b_Zee_i1);
   fChain->SetBranchAddress("Zee_i2", &Zee_i2, &b_Zee_i2);
   fChain->SetBranchAddress("Zee_highestMass", &Zee_highestMass, &b_Zee_highestMass);
   fChain->SetBranchAddress("Zmm_n", &Zmm_n, &b_Zmm_n);
   fChain->SetBranchAddress("Zmm_mass", &Zmm_mass, &b_Zmm_mass);
   fChain->SetBranchAddress("Zmm_i1", &Zmm_i1, &b_Zmm_i1);
   fChain->SetBranchAddress("Zmm_i2", &Zmm_i2, &b_Zmm_i2);
   fChain->SetBranchAddress("Zmm_highestMass", &Zmm_highestMass, &b_Zmm_highestMass);
   fChain->SetBranchAddress("Zem_n", &Zem_n, &b_Zem_n);
   fChain->SetBranchAddress("Zem_mass", &Zem_mass, &b_Zem_mass);
   fChain->SetBranchAddress("Zem_mass_HEEP", &Zem_mass_HEEP, &b_Zem_mass_HEEP);
   fChain->SetBranchAddress("Zem_i1", &Zem_i1, &b_Zem_i1);
   fChain->SetBranchAddress("Zem_i2", &Zem_i2, &b_Zem_i2);
   fChain->SetBranchAddress("Zem_highestMass", &Zem_highestMass, &b_Zem_highestMass);
   fChain->SetBranchAddress("trk_n", &trk_n, &b_trk_n);
   fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
   fChain->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_accept", &trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_accept, &b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_prescale", &trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_prescale, &b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_prescale);
   fChain->SetBranchAddress("trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltL1sL1DoubleEG2210_eta", &trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltL1sL1DoubleEG2210_eta, &b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltL1sL1DoubleEG2210_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltL1sL1DoubleEG2210_phi", &trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltL1sL1DoubleEG2210_phi, &b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltL1sL1DoubleEG2210_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg1TrackIsoFilter_eta", &trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg1TrackIsoFilter_eta, &b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg1TrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg1TrackIsoFilter_phi", &trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg1TrackIsoFilter_phi, &b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg1TrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg2TrackIsoFilter_eta", &trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg2TrackIsoFilter_eta, &b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg2TrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg2TrackIsoFilter_phi", &trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg2TrackIsoFilter_phi, &b_trig_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v2_hltEle24Ele22WPLooseGsfleg2TrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_accept", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_accept, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_prescale", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_prescale, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_prescale);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLPixelMatchFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLPixelMatchFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLPixelMatchFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLPixelMatchFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLPixelMatchFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLPixelMatchFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLNewPixelMatchFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLNewPixelMatchFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLNewPixelMatchFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLNewPixelMatchFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLNewPixelMatchFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltEle33CaloIdLNewPixelMatchFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEG33EtUnseededFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEG33EtUnseededFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEG33EtUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEG33EtUnseededFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEG33EtUnseededFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEG33EtUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLNewPixelMatchUnseededFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLNewPixelMatchUnseededFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLNewPixelMatchUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLNewPixelMatchUnseededFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLNewPixelMatchUnseededFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v3_hltDiEle33CaloIdLNewPixelMatchUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_accept", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_accept, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_prescale", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_prescale, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_prescale);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltL1sL1SingleEG40OR35OR30OR25OR20ORL1DoubleEG2210ORL1SingleJet200_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLPixelMatchFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLPixelMatchFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLPixelMatchFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLPixelMatchFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLPixelMatchFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLPixelMatchFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEG33EtUnseededFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEG33EtUnseededFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEG33EtUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEG33EtUnseededFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEG33EtUnseededFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEG33EtUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLPixelMatchUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_eta", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_eta, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_phi", &trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_phi, &b_trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_accept", &trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_accept, &b_trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_accept);
   fChain->SetBranchAddress("trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_prescale", &trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_prescale, &b_trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_eta", &trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_eta, &b_trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_eta);
   fChain->SetBranchAddress("trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_phi", &trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_phi, &b_trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_phi);
   fChain->SetBranchAddress("trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltSingleEle22WPLooseGsfTrackIsoFilter_eta", &trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltSingleEle22WPLooseGsfTrackIsoFilter_eta, &b_trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltSingleEle22WPLooseGsfTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltSingleEle22WPLooseGsfTrackIsoFilter_phi", &trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltSingleEle22WPLooseGsfTrackIsoFilter_phi, &b_trig_HLT_Ele22_eta2p1_WPLoose_Gsf_v2_hltSingleEle22WPLooseGsfTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_accept", &trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_accept, &b_trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_accept);
   fChain->SetBranchAddress("trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_prescale", &trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_prescale, &b_trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_eta", &trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_eta, &b_trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_eta);
   fChain->SetBranchAddress("trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_phi", &trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_phi, &b_trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG20erORSingleEG25_phi);
   fChain->SetBranchAddress("trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltSingleEle22WPTightGsfTrackIsoFilter_eta", &trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltSingleEle22WPTightGsfTrackIsoFilter_eta, &b_trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltSingleEle22WPTightGsfTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltSingleEle22WPTightGsfTrackIsoFilter_phi", &trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltSingleEle22WPTightGsfTrackIsoFilter_phi, &b_trig_HLT_Ele22_eta2p1_WPTight_Gsf_v2_hltSingleEle22WPTightGsfTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_accept", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_accept, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_accept);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_prescale", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_prescale, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltL1sL1SingleEG25OrSingleIsoEG30er_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltL1sL1SingleEG25OrSingleIsoEG30er_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltL1sL1SingleEG25OrSingleIsoEG30er_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltL1sL1SingleEG25OrSingleIsoEG30er_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltL1sL1SingleEG25OrSingleIsoEG30er_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltL1sL1SingleEG25OrSingleIsoEG30er_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtFilter_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtFilter_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtFilter_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtFilter_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8ClusterShapeFilter_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8ClusterShapeFilter_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8ClusterShapeFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8ClusterShapeFilter_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8ClusterShapeFilter_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8ClusterShapeFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HEFilter_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HEFilter_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HEFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HEFilter_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HEFilter_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HEFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EcalIsoFilter_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EcalIsoFilter_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EcalIsoFilter_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EcalIsoFilter_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HcalIsoFilter_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HcalIsoFilter_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HcalIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HcalIsoFilter_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HcalIsoFilter_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8HcalIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchFilter_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchFilter_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchFilter_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchFilter_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8OneOEMinusOneOPFilter_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8OneOEMinusOneOPFilter_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8OneOEMinusOneOPFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8OneOEMinusOneOPFilter_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8OneOEMinusOneOPFilter_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8OneOEMinusOneOPFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DetaFilter_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DetaFilter_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DetaFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DetaFilter_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DetaFilter_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DetaFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DphiFilter_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DphiFilter_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DphiFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DphiFilter_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DphiFilter_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8DphiFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8TrackIsoFilter_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8TrackIsoFilter_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8TrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8TrackIsoFilter_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8TrackIsoFilter_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8TrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtUnseededFilter_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtUnseededFilter_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtUnseededFilter_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtUnseededFilter_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8EtUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchUnseededFilter_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchUnseededFilter_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchUnseededFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchUnseededFilter_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchUnseededFilter_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8PixelMatchUnseededFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8Mass55Filter_eta", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8Mass55Filter_eta, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8Mass55Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8Mass55Filter_phi", &trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8Mass55Filter_phi, &b_trig_HLT_Ele30WP60_Ele8_Mass55_v2_hltEle30WP60Ele8Mass55Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_v2_accept", &trig_HLT_Ele23_WPLoose_Gsf_v2_accept, &b_trig_HLT_Ele23_WPLoose_Gsf_v2_accept);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_v2_prescale", &trig_HLT_Ele23_WPLoose_Gsf_v2_prescale, &b_trig_HLT_Ele23_WPLoose_Gsf_v2_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_v2_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_eta", &trig_HLT_Ele23_WPLoose_Gsf_v2_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_eta, &b_trig_HLT_Ele23_WPLoose_Gsf_v2_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_v2_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_phi", &trig_HLT_Ele23_WPLoose_Gsf_v2_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_phi, &b_trig_HLT_Ele23_WPLoose_Gsf_v2_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_v2_hltEle23WPLooseGsfTrackIsoFilter_eta", &trig_HLT_Ele23_WPLoose_Gsf_v2_hltEle23WPLooseGsfTrackIsoFilter_eta, &b_trig_HLT_Ele23_WPLoose_Gsf_v2_hltEle23WPLooseGsfTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_v2_hltEle23WPLooseGsfTrackIsoFilter_phi", &trig_HLT_Ele23_WPLoose_Gsf_v2_hltEle23WPLooseGsfTrackIsoFilter_phi, &b_trig_HLT_Ele23_WPLoose_Gsf_v2_hltEle23WPLooseGsfTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_accept", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_accept, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_accept);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_prescale", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_prescale, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_eta", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_eta, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_phi", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_phi, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltL1sL1SingleEG20ORSingleIsoEG18erORSingleIsoEG20_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltEle23WPLooseGsfTrackIsoFilter_eta", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltEle23WPLooseGsfTrackIsoFilter_eta, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltEle23WPLooseGsfTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltEle23WPLooseGsfTrackIsoFilter_phi", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltEle23WPLooseGsfTrackIsoFilter_phi, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltEle23WPLooseGsfTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltJetFilterEle23WPLoose_eta", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltJetFilterEle23WPLoose_eta, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltJetFilterEle23WPLoose_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltJetFilterEle23WPLoose_phi", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltJetFilterEle23WPLoose_phi, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltJetFilterEle23WPLoose_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltHCand80NoEle23WPLoose_eta", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltHCand80NoEle23WPLoose_eta, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltHCand80NoEle23WPLoose_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltHCand80NoEle23WPLoose_phi", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltHCand80NoEle23WPLoose_phi, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltHCand80NoEle23WPLoose_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMhtFilter0_eta", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMhtFilter0_eta, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMhtFilter0_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMhtFilter0_phi", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMhtFilter0_phi, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMhtFilter0_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMetFilter0_eta", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMetFilter0_eta, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMetFilter0_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMetFilter0_phi", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMetFilter0_phi, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltPFMetFilter0_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand80NoEle23WPLooseMET_eta", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand80NoEle23WPLooseMET_eta, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand80NoEle23WPLooseMET_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand80NoEle23WPLooseMET_phi", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand80NoEle23WPLooseMET_phi, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand80NoEle23WPLooseMET_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand70NoEle23WPLooseMHTIDTight_eta", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand70NoEle23WPLooseMHTIDTight_eta, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand70NoEle23WPLooseMHTIDTight_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand70NoEle23WPLooseMHTIDTight_phi", &trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand70NoEle23WPLooseMHTIDTight_phi, &b_trig_HLT_Ele23_WPLoose_Gsf_WHbbBoost_v1_hltWCand70NoEle23WPLooseMHTIDTight_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_v1_accept", &trig_HLT_Ele27_WPLoose_Gsf_v1_accept, &b_trig_HLT_Ele27_WPLoose_Gsf_v1_accept);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_v1_prescale", &trig_HLT_Ele27_WPLoose_Gsf_v1_prescale, &b_trig_HLT_Ele27_WPLoose_Gsf_v1_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_v1_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_eta", &trig_HLT_Ele27_WPLoose_Gsf_v1_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_eta, &b_trig_HLT_Ele27_WPLoose_Gsf_v1_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_v1_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_phi", &trig_HLT_Ele27_WPLoose_Gsf_v1_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_phi, &b_trig_HLT_Ele27_WPLoose_Gsf_v1_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_v1_hltEle27noerWPLooseGsfTrackIsoFilter_eta", &trig_HLT_Ele27_WPLoose_Gsf_v1_hltEle27noerWPLooseGsfTrackIsoFilter_eta, &b_trig_HLT_Ele27_WPLoose_Gsf_v1_hltEle27noerWPLooseGsfTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_v1_hltEle27noerWPLooseGsfTrackIsoFilter_phi", &trig_HLT_Ele27_WPLoose_Gsf_v1_hltEle27noerWPLooseGsfTrackIsoFilter_phi, &b_trig_HLT_Ele27_WPLoose_Gsf_v1_hltEle27noerWPLooseGsfTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_accept", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_accept, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_accept);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_prescale", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_prescale, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_eta", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_eta, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_phi", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_phi, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltL1sL1SingleIsoEG20OrSingleEG20OrSingleEG25_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltEle27noerWPLooseGsfTrackIsoFilter_eta", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltEle27noerWPLooseGsfTrackIsoFilter_eta, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltEle27noerWPLooseGsfTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltEle27noerWPLooseGsfTrackIsoFilter_phi", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltEle27noerWPLooseGsfTrackIsoFilter_phi, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltEle27noerWPLooseGsfTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltJetFilterEle27WPLoose_eta", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltJetFilterEle27WPLoose_eta, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltJetFilterEle27WPLoose_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltJetFilterEle27WPLoose_phi", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltJetFilterEle27WPLoose_phi, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltJetFilterEle27WPLoose_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltHCand80NoEle27WPLoose_eta", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltHCand80NoEle27WPLoose_eta, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltHCand80NoEle27WPLoose_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltHCand80NoEle27WPLoose_phi", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltHCand80NoEle27WPLoose_phi, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltHCand80NoEle27WPLoose_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMhtFilter0_eta", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMhtFilter0_eta, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMhtFilter0_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMhtFilter0_phi", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMhtFilter0_phi, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMhtFilter0_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMetFilter0_eta", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMetFilter0_eta, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMetFilter0_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMetFilter0_phi", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMetFilter0_phi, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltPFMetFilter0_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand80NoEle27WPLooseMET_eta", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand80NoEle27WPLooseMET_eta, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand80NoEle27WPLooseMET_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand80NoEle27WPLooseMET_phi", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand80NoEle27WPLooseMET_phi, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand80NoEle27WPLooseMET_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand70NoEle27WPLooseMHTIDTight_eta", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand70NoEle27WPLooseMHTIDTight_eta, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand70NoEle27WPLooseMHTIDTight_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand70NoEle27WPLooseMHTIDTight_phi", &trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand70NoEle27WPLooseMHTIDTight_phi, &b_trig_HLT_Ele27_WPLoose_Gsf_WHbbBoost_v2_hltWCand70NoEle27WPLooseMHTIDTight_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_accept", &trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_accept, &b_trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_accept);
   fChain->SetBranchAddress("trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_prescale", &trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_prescale, &b_trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_eta", &trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_eta, &b_trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_phi", &trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_phi, &b_trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltEle27WPLooseGsfTrackIsoFilter_eta", &trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltEle27WPLooseGsfTrackIsoFilter_eta, &b_trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltEle27WPLooseGsfTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltEle27WPLooseGsfTrackIsoFilter_phi", &trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltEle27WPLooseGsfTrackIsoFilter_phi, &b_trig_HLT_Ele27_eta2p1_WPLoose_Gsf_v2_hltEle27WPLooseGsfTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_accept", &trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_accept, &b_trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_accept);
   fChain->SetBranchAddress("trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_prescale", &trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_prescale, &b_trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_eta", &trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_eta, &b_trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_phi", &trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_phi, &b_trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG22erOrSingleEG25_phi);
   fChain->SetBranchAddress("trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltEle27WPTightGsfTrackIsoFilter_eta", &trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltEle27WPTightGsfTrackIsoFilter_eta, &b_trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltEle27WPTightGsfTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltEle27WPTightGsfTrackIsoFilter_phi", &trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltEle27WPTightGsfTrackIsoFilter_phi, &b_trig_HLT_Ele27_eta2p1_WPTight_Gsf_v2_hltEle27WPTightGsfTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_accept", &trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_accept, &b_trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_accept);
   fChain->SetBranchAddress("trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_prescale", &trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_prescale, &b_trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG30erOrSingleEG40_eta", &trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG30erOrSingleEG40_eta, &b_trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG30erOrSingleEG40_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG30erOrSingleEG40_phi", &trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG30erOrSingleEG40_phi, &b_trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltL1sL1SingleIsoEG30erOrSingleEG40_phi);
   fChain->SetBranchAddress("trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltEle32WPTightGsfTrackIsoFilter_eta", &trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltEle32WPTightGsfTrackIsoFilter_eta, &b_trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltEle32WPTightGsfTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltEle32WPTightGsfTrackIsoFilter_phi", &trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltEle32WPTightGsfTrackIsoFilter_phi, &b_trig_HLT_Ele32_eta2p1_WPTight_Gsf_v2_hltEle32WPTightGsfTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_accept", &trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_accept, &b_trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_accept);
   fChain->SetBranchAddress("trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_prescale", &trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_prescale, &b_trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_eta", &trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_eta, &b_trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_eta);
   fChain->SetBranchAddress("trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_phi", &trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_phi, &b_trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_phi);
   fChain->SetBranchAddress("trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter_eta", &trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter_eta, &b_trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter_phi", &trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter_phi, &b_trig_HLT_Ele105_CaloIdVT_GsfTrkIdT_v3_hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_accept", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_accept, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_accept);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_prescale", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_prescale, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG2210_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG2210_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG2210_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG2210_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG2210_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG2210_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_accept", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_accept, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_accept);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_prescale", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_prescale, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG1510_eta", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG1510_eta, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG1510_eta);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG1510_phi", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG1510_phi, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltL1sL1DoubleEG1510_phi);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter_eta", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter_eta, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter_phi", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter_phi, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG2210_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG2210_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG2210_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG2210_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG2210_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG2210_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG1510_eta", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG1510_eta, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG1510_eta);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG1510_phi", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG1510_phi, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1DoubleEG1510_phi);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi", &trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi, &b_trig_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_accept", &trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_accept, &b_trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_accept);
   fChain->SetBranchAddress("trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_prescale", &trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_prescale, &b_trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG20_eta", &trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG20_eta, &b_trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG20_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG20_phi", &trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG20_phi, &b_trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG20_phi);
   fChain->SetBranchAddress("trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter_eta", &trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter_eta, &b_trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter_phi", &trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter_phi, &b_trig_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v3_hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_accept", &trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_accept, &b_trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_accept);
   fChain->SetBranchAddress("trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_prescale", &trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_prescale, &b_trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltL1sL1SingleEG15_eta", &trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltL1sL1SingleEG15_eta, &b_trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltL1sL1SingleEG15_eta);
   fChain->SetBranchAddress("trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltL1sL1SingleEG15_phi", &trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltL1sL1SingleEG15_phi, &b_trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltL1sL1SingleEG15_phi);
   fChain->SetBranchAddress("trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter_eta", &trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter_eta, &b_trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter_phi", &trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter_phi, &b_trig_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v2_hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept", &trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept, &b_trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_accept);
   fChain->SetBranchAddress("trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale", &trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale, &b_trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG15ORSingleEG10_eta", &trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG15ORSingleEG10_eta, &b_trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG15ORSingleEG10_eta);
   fChain->SetBranchAddress("trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG15ORSingleEG10_phi", &trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG15ORSingleEG10_phi, &b_trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltL1sL1SingleEG15ORSingleEG10_phi);
   fChain->SetBranchAddress("trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter_eta", &trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter_eta, &b_trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter_phi", &trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter_phi, &b_trig_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v3_hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter_phi);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_accept", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_accept, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_accept);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_prescale", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_prescale, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_prescale);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_eta", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_eta, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_eta);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_phi", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_phi, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltL1sL1SingleEG40ORL1SingleEG35ORL1SingleJet200_phi);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_eta", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_eta, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_eta);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_phi", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_phi, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_v2_hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter_phi);
   fChain->SetBranchAddress("trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_accept", &trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_accept, &b_trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_accept);
   fChain->SetBranchAddress("trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_prescale", &trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_prescale, &b_trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_prescale);
   fChain->SetBranchAddress("trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltL1sL1SingleIsoEG22erOrSingleIsoEG20OrSingleEG25OrSingleEG20_eta", &trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltL1sL1SingleIsoEG22erOrSingleIsoEG20OrSingleEG25OrSingleEG20_eta, &b_trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltL1sL1SingleIsoEG22erOrSingleIsoEG20OrSingleEG25OrSingleEG20_eta);
   fChain->SetBranchAddress("trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltL1sL1SingleIsoEG22erOrSingleIsoEG20OrSingleEG25OrSingleEG20_phi", &trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltL1sL1SingleIsoEG22erOrSingleIsoEG20OrSingleEG25OrSingleEG20_phi, &b_trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltL1sL1SingleIsoEG22erOrSingleIsoEG20OrSingleEG25OrSingleEG20_phi);
   fChain->SetBranchAddress("trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltEle23WPLooseGsfTrackIsoFilterForEleStream_eta", &trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltEle23WPLooseGsfTrackIsoFilterForEleStream_eta, &b_trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltEle23WPLooseGsfTrackIsoFilterForEleStream_eta);
   fChain->SetBranchAddress("trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltEle23WPLooseGsfTrackIsoFilterForEleStream_phi", &trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltEle23WPLooseGsfTrackIsoFilterForEleStream_phi, &b_trig_AlCa_Ele23_WPVeryLoose_Gsf_v1_hltEle23WPLooseGsfTrackIsoFilterForEleStream_phi);
   Notify();
}

Bool_t IIHEAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void IIHEAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t IIHEAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef IIHEAnalysis_cxx
