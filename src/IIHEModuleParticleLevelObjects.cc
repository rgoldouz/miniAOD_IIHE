#include "UserCode/IIHETree/interface/IIHEModuleParticleLevelObjects.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleParticleLevelObjects::IIHEModuleParticleLevelObjects(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig){
  particleLevelJetsToken_ = iC.consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("particleLevelJetsCollection"));
  particleLevelBJetsToken_ = iC.consumes<std::vector<reco::GenJet>> (iConfig.getParameter<edm::InputTag>("particleLevelBJetsCollection"));
  particleLevelak1DressedMuonsToken_ = iC.consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("particleLevelak1DressedMuonsCollection"));
  particleLevelak1DressedElectronsToken_ = iC.consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("particleLevelak1DressedElectronsCollection"));
  particleLevelNeutrinosToken_ = iC.consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("particleLevelNeutrinosCollection"));
}
IIHEModuleParticleLevelObjects::~IIHEModuleParticleLevelObjects(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleParticleLevelObjects::beginJob(){
  setBranchType(kVectorFloat) ;
  addBranch("pl_jet_pt"   ) ;
  addBranch("pl_jet_eta"   ) ;
  addBranch("pl_jet_phi"   ) ;
  addBranch("pl_bjet_pt"   ) ;
  addBranch("pl_bjet_eta"   ) ;
  addBranch("pl_bjet_phi"   ) ;
  addBranch("pl_ele_pt"   ) ;
  addBranch("pl_ele_eta"   ) ;
  addBranch("pl_ele_phi"   ) ;
  addBranch("pl_mu_pt"   ) ;
  addBranch("pl_mu_eta"   ) ;
  addBranch("pl_mu_phi"   ) ;
  addBranch("pl_nu_pt"   ) ;
  addBranch("pl_nu_eta"   ) ;
  addBranch("pl_nu_phi"   ) ;
}

// ------------ method called to for each event  ------------
void IIHEModuleParticleLevelObjects::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle<std::vector<reco::GenJet>> particleLevelJetsHandle_;
  iEvent.getByToken(particleLevelJetsToken_, particleLevelJetsHandle_);

  edm::Handle<std::vector<reco::GenJet>> particleLevelBJetsHandle_;
  iEvent.getByToken(particleLevelBJetsToken_, particleLevelBJetsHandle_);

  edm::Handle<std::vector<reco::GenJet>> particleLevelak1DressedMuonsHandle_;
  iEvent.getByToken(particleLevelak1DressedMuonsToken_, particleLevelak1DressedMuonsHandle_);

  edm::Handle<std::vector<reco::GenJet>> particleLevelak1DressedElectronsHandle_;
  iEvent.getByToken(particleLevelak1DressedElectronsToken_, particleLevelak1DressedElectronsHandle_);

  edm::Handle<std::vector<reco::GenParticle>> particleLevelNeutrinosHandle_;
  iEvent.getByToken(particleLevelNeutrinosToken_, particleLevelNeutrinosHandle_);

  for (size_t j = 0; j < particleLevelJetsHandle_->size();++j){
    store("pl_jet_pt"     , particleLevelJetsHandle_->at(j).pt()) ;
    store("pl_jet_eta"    , particleLevelJetsHandle_->at(j).eta()) ;
    store("pl_jet_phi"    , particleLevelJetsHandle_->at(j).phi()) ;
  }

  for (size_t j = 0; j < particleLevelBJetsHandle_->size();++j){
    store("pl_bjet_pt"     , particleLevelBJetsHandle_->at(j).pt()) ;
    store("pl_bjet_eta"    , particleLevelBJetsHandle_->at(j).eta()) ;
    store("pl_bjet_phi"    , particleLevelBJetsHandle_->at(j).phi()) ;
  }

  for (size_t j = 0; j < particleLevelak1DressedMuonsHandle_->size();++j){
    store("pl_mu_pt"     , particleLevelak1DressedMuonsHandle_->at(j).pt()) ;
    store("pl_mu_eta"    , particleLevelak1DressedMuonsHandle_->at(j).eta()) ;
    store("pl_mu_phi"    , particleLevelak1DressedMuonsHandle_->at(j).phi()) ;
  }


  for (size_t j = 0; j < particleLevelak1DressedElectronsHandle_->size();++j){
    store("pl_ele_pt"     , particleLevelak1DressedElectronsHandle_->at(j).pt()) ;
    store("pl_ele_eta"    , particleLevelak1DressedElectronsHandle_->at(j).eta()) ;
    store("pl_ele_phi"    , particleLevelak1DressedElectronsHandle_->at(j).phi()) ;
  }

  for (size_t j = 0; j < particleLevelNeutrinosHandle_->size();++j){
    store("pl_nu_pt"     , particleLevelNeutrinosHandle_->at(j).pt()) ;
    store("pl_nu_eta"    , particleLevelNeutrinosHandle_->at(j).eta()) ;
    store("pl_nu_phi"    , particleLevelNeutrinosHandle_->at(j).phi()) ;
  }

}
void IIHEModuleParticleLevelObjects::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleParticleLevelObjects::beginEvent(){}
void IIHEModuleParticleLevelObjects::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleParticleLevelObjects::endJob(){
}

DEFINE_FWK_MODULE(IIHEModuleParticleLevelObjects);
