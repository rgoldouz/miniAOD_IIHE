#include "UserCode/IIHETree/interface/IIHEModuleData.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleData::IIHEModuleData(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig){
  electronCollectionLabel_                  = iConfig.getParameter<edm::InputTag>("electronsBeforeGSFixCollection") ;
  electronCollectionToken_                  = iC.consumes<View<pat::Electron> > (electronCollectionLabel_);
  ETThreshold_ = iConfig.getUntrackedParameter<double>("electrons_ETThreshold", 0.0 ) ;

  METCollectionLabel_                       = iConfig.getParameter<edm::InputTag>("METsMuEGCleanCollection") ; 
  METCollectionToken_                       = iC.consumes<View<pat::MET> > (METCollectionLabel_);

  particleFlowEGammaGSFixedCollectionLabel_ = iConfig.getParameter<edm::InputTag>("particleFlowEGammaGSFixedCollection") ;
  particleFlowEGammaGSFixedCollectionToken_ = iC.consumes<bool> (particleFlowEGammaGSFixedCollectionLabel_);

  pfcandidateCollectionLabel_               = iConfig.getParameter<edm::InputTag>("discardedMuonCollection") ;
  pfcandidateCollectionToken_               = iC.consumes<View<pat::PackedCandidate> > (pfcandidateCollectionLabel_);

  ecalMultiAndGSGlobalRecHitEBLabel_        = iConfig.getParameter<edm::InputTag>("ecalMultiAndGSGlobalRecHitEBCollection") ;
  ecalMultiAndGSGlobalRecHitEBToken_        = iC.consumes<edm::EDCollection<DetId>>(ecalMultiAndGSGlobalRecHitEBLabel_);

  triggerResultsLabel_                      = iConfig.getParameter<edm::InputTag>("TriggerResults") ;
  triggerResultsToken_                      = iC.consumes<edm::TriggerResults>(triggerResultsLabel_);

}
IIHEModuleData::~IIHEModuleData(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleData::beginJob(){
  setBranchType(kVectorFloat) ;
  addBranch("gsf_bGSfix_deltaPhiSuperClusterTrackAtVtx") ;
  addBranch("gsf_bGSfix_full5x5_sigmaIetaIeta");
  addBranch("gsf_bGSfix_superClusterEnergy") ;
  addBranch("gsf_bGSfix_hadronicOverEm") ;
  addBranch("gsf_bGSfix_full5x5_e5x5") ;
  addBranch("gsf_bGSfix_full5x5_e1x5") ;
  addBranch("gsf_bGSfix_full5x5_e2x5Max") ;
  addBranch("gsf_bGSfix_full5x5_sigmaIetaIeta");
  addBranch("gsf_bGSfix_dr03EcalRecHitSumEt") ;
  addBranch("gsf_bGSfix_dr03HcalDepth1TowerSumEt") ;
  addBranch("gsf_bGSfix_deltaEtaSeedClusterTrackAtVtx") ;
  addBranch("gsf_bGSfix_caloEnergy") ;
  addBranch("gsf_bGSfix_theta") ;
  addBranch("gsf_bGSfix_sc_eta") ;
  addBranch("gsf_bGSfix_phi") ;
  addBranch("gsf_bGSfix_sc_phi") ;

  setBranchType(kVectorBool) ;
  addBranch("gsf_bGSfix_ecaldrivenSeed"   ) ;

  setBranchType(kVectorInt) ;
  addBranch("gsf_bGSfix_nLostInnerHits"   ) ;

  setBranchType(kFloat) ;
  addBranch("MET_pfMetMuEGClean_et"   ) ;
  addBranch("MET_pfMetMuEGClean_phi"  ) ;

  addBranch("ev_particleFlowEGammaGSFixed", kBool) ;
  addBranch("ev_ecalMultiAndGSGlobalRecHitEB", kBool) ;
  addBranch("ev_duplicateMuons", kBool) ;
}

// ------------ method called to for each event  ------------
void IIHEModuleData::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){



  edm::Handle<edm::View<pat::Electron>>     electronCollection_ ;
  iEvent.getByToken( electronCollectionToken_ , electronCollection_) ;

  edm::Handle<edm::View<pat::MET> > METCollection_;
  iEvent.getByToken(METCollectionToken_, METCollection_);

  edm::Handle<bool> particleFlowEGammaGSFixedCollection_ ;
  iEvent.getByToken(particleFlowEGammaGSFixedCollectionToken_, particleFlowEGammaGSFixedCollection_) ;

  edm::Handle<edm::EDCollection<DetId>> ecalMultiAndGSGlobalRecHitEB_;
  iEvent.getByToken(ecalMultiAndGSGlobalRecHitEBToken_, ecalMultiAndGSGlobalRecHitEB_) ;

  edm::Handle<TriggerResults> triggerResultsCollection_ ;
  iEvent.getByToken(triggerResultsToken_, triggerResultsCollection_);

  store("ev_ecalMultiAndGSGlobalRecHitEB", ecalMultiAndGSGlobalRecHitEB_.isValid()) ;

  bool particleFlowEGammaGSFixed = *particleFlowEGammaGSFixedCollection_ ;
  store("ev_particleFlowEGammaGSFixed", particleFlowEGammaGSFixed) ;

  if (METCollection_.isValid()) {
    const pat::MET *pfMET = 0;
    pfMET = &(METCollection_->front());
    store("MET_pfMetMuEGClean_et"   , pfMET->et()    ) ;
    store("MET_pfMetMuEGClean_phi"  , pfMET->phi()   ) ;
  }

  for( unsigned int i = 0 ; i < electronCollection_->size() ; i++ ) {
    Ptr<pat::Electron> gsfiter = electronCollection_->ptrAt( i );

    float ET = gsfiter->caloEnergy()*sin(2.*atan(exp(-1.*gsfiter->superCluster()->eta()))) ;
    if(ET<ETThreshold_ && gsfiter->pt()<ETThreshold_) continue ;

    store("gsf_bGSfix_deltaPhiSuperClusterTrackAtVtx", gsfiter->deltaPhiSuperClusterTrackAtVtx()) ;
    store("gsf_bGSfix_full5x5_sigmaIetaIeta"         , gsfiter->full5x5_sigmaIetaIeta()) ;
    store("gsf_bGSfix_nLostInnerHits"                , gsfiter->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)) ;
    store("gsf_bGSfix_hadronicOverEm"                , gsfiter->hadronicOverEm()                ) ;
    store("gsf_bGSfix_superClusterEnergy"            , gsfiter->superCluster()->energy()        ) ;
    store("gsf_bGSfix_full5x5_e5x5"                  , gsfiter->full5x5_e5x5()) ;
    store("gsf_bGSfix_full5x5_e1x5"                  , gsfiter->full5x5_e1x5()) ;
    store("gsf_bGSfix_full5x5_e2x5Max"               , gsfiter->full5x5_e2x5Max()) ;
    store("gsf_bGSfix_full5x5_sigmaIetaIeta"         , gsfiter->full5x5_sigmaIetaIeta()) ;
    store("gsf_bGSfix_dr03EcalRecHitSumEt"           , gsfiter->dr03EcalRecHitSumEt()           ) ;
    store("gsf_bGSfix_dr03HcalDepth1TowerSumEt"      , gsfiter->dr03HcalDepth1TowerSumEt()      ) ;
    store("gsf_bGSfix_ecaldrivenSeed"                , gsfiter->ecalDrivenSeed()                ) ;
    store("gsf_bGSfix_deltaEtaSeedClusterTrackAtVtx" , gsfiter->deltaEtaSeedClusterTrackAtVtx() ) ;
    store("gsf_bGSfix_caloEnergy"                    , gsfiter->caloEnergy()                    ) ;
    store("gsf_bGSfix_theta"                         , gsfiter->theta()                         ) ;
    store("gsf_bGSfix_sc_eta"                        , gsfiter->superCluster()->eta()                    ) ;
    store("gsf_bGSfix_phi"                           , gsfiter->phi()                           ) ;
    store("gsf_bGSfix_sc_phi"        , gsfiter->superCluster()->phi()                    ) ;
  }


}
void IIHEModuleData::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleData::beginEvent(){}
void IIHEModuleData::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleData::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleData);
