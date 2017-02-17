#include "UserCode/IIHETree/interface/IIHEModuleTau.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleTau::IIHEModuleTau(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig){
}
IIHEModuleTau::~IIHEModuleTau(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleTau::beginJob(){
}

// ------------ method called to for each event  ------------
void IIHEModuleTau::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
}
void IIHEModuleTau::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleTau::beginEvent(){}
void IIHEModuleTau::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleTau::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleTau);
