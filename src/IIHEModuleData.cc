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

  METCollectionLabel_                       = iConfig.getParameter<edm::InputTag>("METsMuEGCleanCollection") ; 
  METCollectionToken_                       = iC.consumes<View<pat::MET> > (METCollectionLabel_);

  particleFlowEGammaGSFixedCollectionLabel_ = iConfig.getParameter<edm::InputTag>("particleFlowEGammaGSFixedCollection") ;
  particleFlowEGammaGSFixedCollectionToken_ = iC.consumes<bool> (particleFlowEGammaGSFixedCollectionLabel_);

  pfcandidateCollectionLabel_               = iConfig.getParameter<edm::InputTag>("discardedMuonCollection") ;
  pfcandidateCollectionToken_               = iC.consumes<View<pat::PackedCandidate> > (pfcandidateCollectionLabel_);

  ecalMultiAndGSGlobalRecHitEBLabel_        = iConfig.getParameter<edm::InputTag>("ecalMultiAndGSGlobalRecHitEBCollection") ;
  ecalMultiAndGSGlobalRecHitEBToken_        = iC.consumes<View<edm::EDCollection<DetId>>>(ecalMultiAndGSGlobalRecHitEBLabel_);

}
IIHEModuleData::~IIHEModuleData(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleData::beginJob(){
}

// ------------ method called to for each event  ------------
void IIHEModuleData::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

}
void IIHEModuleData::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleData::beginEvent(){}
void IIHEModuleData::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleData::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleData);
