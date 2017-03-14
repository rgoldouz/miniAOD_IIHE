#include "UserCode/IIHETree/interface/IIHEModuleMET.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h" 
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

//////////////////////////////////////////////////////////////////////////////////////////
//                             IIHEMETVariable classes                            
//////////////////////////////////////////////////////////////////////////////////////////
IIHEMETVariableBase::IIHEMETVariableBase(std::string prefix, std::string name, int type){
  name_       = name ;
  branchName_ = prefix + "_" + name_ ;
  branchType_ = type ;
}
bool IIHEMETVariableBase::addBranch(IIHEAnalysis* analysis){
  return analysis->addBranch(branchName_, branchType_) ;
}

IIHEMETVariableInt::IIHEMETVariableInt(std::string prefix, std::string name):
IIHEMETVariableBase(prefix, name, kVectorInt){
  reset() ;
}
void IIHEMETVariableInt::store(IIHEAnalysis* analysis){
  analysis->store(BranchName(), value_) ;
}

IIHEMETVariableFloat::IIHEMETVariableFloat(std::string prefix, std::string name):
IIHEMETVariableBase(prefix, name, kVectorFloat){
  reset() ;
}
void IIHEMETVariableFloat::store(IIHEAnalysis* analysis){
  analysis->store(BranchName(), value_ ) ;
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                  IIHEMET class                                 
//////////////////////////////////////////////////////////////////////////////////////////
IIHEMETWrapper::IIHEMETWrapper(std::string prefix){
  prefix_ = prefix ;

  et_               = new IIHEMETVariableInt  (prefix_, "et"                 ) ;
  phi_              = new IIHEMETVariableInt  (prefix_, "phi"                ) ;
  significance_     = new IIHEMETVariableInt  (prefix_, "significance"       ) ;

  variables_.push_back((IIHEMETVariableBase*) et_                   ) ;
  variables_.push_back((IIHEMETVariableBase*) phi_                  ) ;
  variables_.push_back((IIHEMETVariableBase*) significance_         ) ;
}

void IIHEMETWrapper::addBranches(IIHEAnalysis* analysis){
  for(unsigned int i=0 ; i<variables_.size() ; ++i){
    variables_.at(i)->addBranch(analysis) ;
  }
}
void IIHEMETWrapper::reset(){
  for(unsigned int i=0 ; i<variables_.size() ; ++i){
    variables_.at(i)->reset() ;
  }
}

void IIHEMETWrapper::fill(pat::MET MET){
  et_        ->fill(MET.et()                ) ;
  phi_       ->fill(MET.phi()               ) ;
  significance_   ->fill(MET.metSignificance()   ) ;
}
void IIHEMETWrapper::store(IIHEAnalysis* analysis){
  for(unsigned int i=0 ; i<variables_.size() ; ++i){
    variables_.at(i)->store(analysis) ;
  }
}


//////////////////////////////////////////////////////////////////////////////////////////
//                                  Main IIHEMuonModule                                 
//////////////////////////////////////////////////////////////////////////////////////////


IIHEModuleMET::IIHEModuleMET(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig),
  metWrapper_(new IIHEMETWrapper("MET_T1"))
{
  pfMETToken_ =  iC.consumes<View<pat::MET> > (iConfig.getParameter<edm::InputTag>("METCollection"));
}
IIHEModuleMET::~IIHEModuleMET(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleMET::beginJob(){
  setBranchType(kFloat) ;
  IIHEAnalysis* analysis = parent_ ;
  addBranch("MET_gen"   ) ;
  metWrapper_->addBranches(analysis) ;
}

// ------------ method called to for each event  ------------
void IIHEModuleMET::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::Handle<edm::View<pat::MET> > pfMETHandle;
  iEvent.getByToken(pfMETToken_, pfMETHandle);
  IIHEAnalysis* analysis = parent_ ;
  metWrapper_->reset() ;
  if (pfMETHandle.isValid()) {
    const pat::MET *pfMET = 0;
    pfMET = &(pfMETHandle->front());
    store("MET_gen"   , pfMET->genMET()     ) ;
    metWrapper_->fill(pfMETHandle->front()) ;
    metWrapper_ ->store(analysis) ;
  }
}
void IIHEModuleMET::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleMET::beginEvent(){}
void IIHEModuleMET::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleMET::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleMET);
