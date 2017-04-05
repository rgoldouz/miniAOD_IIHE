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

  et_               = new IIHEMETVariableFloat  (prefix_, "et"                 ) ;
  phi_              = new IIHEMETVariableFloat  (prefix_, "phi"                ) ;
  significance_     = new IIHEMETVariableFloat  (prefix_, "significance"       ) ;

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
//                                  Main IIHEMETModule                                 
//////////////////////////////////////////////////////////////////////////////////////////


IIHEModuleMET::IIHEModuleMET(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig),
  metnominalWrapper_(new IIHEMETWrapper("MET_nominal")),
  metWrapper_(new IIHEMETWrapper("MET")),
  metT1Wrapper_(new IIHEMETWrapper("MET_T1")),
  metT1JetEnDownWrapper_(new IIHEMETWrapper("MET_T1JetEnDown")),
  metT1JetEnUpWrapper_(new IIHEMETWrapper("MET_T1JetEnUp")),
  metT1SmearJetEnDownWrapper_(new IIHEMETWrapper("MET_T1SmearJetEnDown")),
  metT1SmearJetEnUpWrapper_(new IIHEMETWrapper("MET_T1SmearJetEnUp")),
  metT1SmearJetResDownWrapper_(new IIHEMETWrapper("MET_T1SmearJetResDown")),
  metT1SmearJetResUpWrapper_(new IIHEMETWrapper("MET_T1SmearJetResUp")),
  metTxyWrapper_(new IIHEMETWrapper("MET_Txy")),
  metFinalWrapper_(new IIHEMETWrapper("MET_FinalCollection"))
{
  pfMETToken_                               =  iC.consumes<View<pat::MET> > (iConfig.getParameter<edm::InputTag>("METCollection"));
  patPFMetCollectionToken_                  =  iC.consumes<View<pat::MET> > (iConfig.getParameter<edm::InputTag>("patPFMetCollection"));
  patPFMetT1CollectionToken_                =  iC.consumes<View<pat::MET> > (iConfig.getParameter<edm::InputTag>("patPFMetT1Collection"));
  patPFMetT1JetEnDownCollectionToken_       =  iC.consumes<View<pat::MET> > (iConfig.getParameter<edm::InputTag>("patPFMetT1JetEnDownCollection"));
  patPFMetT1JetEnUpCollectionToken_         =  iC.consumes<View<pat::MET> > (iConfig.getParameter<edm::InputTag>("patPFMetT1JetEnUpCollection"));
  patPFMetT1SmearJetEnDownCollectionToken_  =  iC.consumes<View<pat::MET> > (iConfig.getParameter<edm::InputTag>("patPFMetT1SmearJetEnDownCollection"));
  patPFMetT1SmearJetEnUpCollectionToken_    =  iC.consumes<View<pat::MET> > (iConfig.getParameter<edm::InputTag>("patPFMetT1SmearJetEnUpCollection"));
  patPFMetT1SmearJetResDownCollectionToken_ =  iC.consumes<View<pat::MET> > (iConfig.getParameter<edm::InputTag>("patPFMetT1SmearJetResDownCollection"));
  patPFMetT1SmearJetResUpCollectionToken_   =  iC.consumes<View<pat::MET> > (iConfig.getParameter<edm::InputTag>("patPFMetT1SmearJetResUpCollection"));
  patPFMetTxyToken_                         =  iC.consumes<View<pat::MET> > (iConfig.getParameter<edm::InputTag>("patPFMetTxyCollection"));
  patPFMetFinalCollectionToken_             =  iC.consumes<View<pat::MET> > (iConfig.getParameter<edm::InputTag>("patPFMetFinalCollection"));
  isMC_ = iConfig.getUntrackedParameter<bool>("isMC") ;
}
IIHEModuleMET::~IIHEModuleMET(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleMET::beginJob(){
  setBranchType(kFloat) ;
  IIHEAnalysis* analysis = parent_ ;
  addBranch("MET_gen"   ) ;
  metnominalWrapper_->addBranches(analysis) ;
  metWrapper_->addBranches(analysis) ;
  metT1Wrapper_->addBranches(analysis) ;
  metT1JetEnDownWrapper_->addBranches(analysis) ;
  metT1JetEnUpWrapper_->addBranches(analysis) ;
  metT1SmearJetEnDownWrapper_->addBranches(analysis) ;
  metT1SmearJetEnUpWrapper_->addBranches(analysis) ;
  metT1SmearJetResDownWrapper_->addBranches(analysis) ;
  metT1SmearJetResUpWrapper_->addBranches(analysis) ;
  metTxyWrapper_->addBranches(analysis) ;
  metFinalWrapper_->addBranches(analysis) ;
}

// ------------ method called to for each event  ------------
void IIHEModuleMET::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::Handle<edm::View<pat::MET> > pfMETHandle_;
  iEvent.getByToken(pfMETToken_, pfMETHandle_);

  edm::Handle<edm::View<pat::MET> > patPFMetCollectionHandle_;
  iEvent.getByToken(patPFMetCollectionToken_, patPFMetCollectionHandle_);

  edm::Handle<edm::View<pat::MET> > patPFMetT1CollectionHandle_;
  iEvent.getByToken(patPFMetT1CollectionToken_, patPFMetT1CollectionHandle_);

  edm::Handle<edm::View<pat::MET> > patPFMetT1JetEnDownCollectionHandle_;
  iEvent.getByToken(patPFMetT1JetEnDownCollectionToken_, patPFMetT1JetEnDownCollectionHandle_);

  edm::Handle<edm::View<pat::MET> > patPFMetT1JetEnUpCollectionHandle_;
  iEvent.getByToken(patPFMetT1JetEnUpCollectionToken_, patPFMetT1JetEnUpCollectionHandle_);

  edm::Handle<edm::View<pat::MET> > patPFMetT1SmearJetEnDownCollectionHandle_;
  iEvent.getByToken(patPFMetT1SmearJetEnDownCollectionToken_, patPFMetT1SmearJetEnDownCollectionHandle_);

  edm::Handle<edm::View<pat::MET> > patPFMetT1SmearJetEnUpCollectionHandle_;
  iEvent.getByToken(patPFMetT1SmearJetEnUpCollectionToken_, patPFMetT1SmearJetEnUpCollectionHandle_);

  edm::Handle<edm::View<pat::MET> > patPFMetT1SmearJetResUpCollectionHandle_;
  iEvent.getByToken(patPFMetT1SmearJetResUpCollectionToken_, patPFMetT1SmearJetResUpCollectionHandle_);

  edm::Handle<edm::View<pat::MET> > patPFMetT1SmearJetResDownCollectionHandle_;
  iEvent.getByToken(patPFMetT1SmearJetResDownCollectionToken_, patPFMetT1SmearJetResDownCollectionHandle_);

  edm::Handle<edm::View<pat::MET> > patPFMetTxyHandle_;
  iEvent.getByToken(patPFMetTxyToken_, patPFMetTxyHandle_);

  edm::Handle<edm::View<pat::MET> > patPFMetFinalCollectionHandle_;
  iEvent.getByToken(patPFMetFinalCollectionToken_, patPFMetFinalCollectionHandle_);


  metnominalWrapper_->reset() ;
  metWrapper_->reset() ;
  metT1Wrapper_->reset() ;
  metT1JetEnDownWrapper_->reset() ;
  metT1JetEnUpWrapper_->reset() ;
  metT1SmearJetEnDownWrapper_->reset() ;
  metT1SmearJetEnUpWrapper_->reset() ;
  metT1SmearJetResDownWrapper_->reset() ;
  metT1SmearJetResUpWrapper_->reset() ;
  metTxyWrapper_->reset() ;
  metFinalWrapper_->reset() ;

  IIHEAnalysis* analysis = parent_ ;  
  Ptr<pat::MET> pfMET = pfMETHandle_->ptrAt( 0 );
  store("MET_gen"   , pfMET->genMET()     ) ;

  metnominalWrapper_->fill(pfMETHandle_->front()) ;
  metnominalWrapper_->store(analysis) ;

  metWrapper_->fill(patPFMetCollectionHandle_->front()) ;
  metWrapper_ ->store(analysis) ;

  metT1Wrapper_->fill(patPFMetT1CollectionHandle_->front()) ;
  metT1Wrapper_->store(analysis) ;

  metT1JetEnDownWrapper_->fill(patPFMetT1JetEnDownCollectionHandle_->front()) ;
  metT1JetEnDownWrapper_->store(analysis) ;

  metT1JetEnUpWrapper_->fill(patPFMetT1JetEnUpCollectionHandle_->front()) ;
  metT1JetEnUpWrapper_->store(analysis) ;

  metT1SmearJetEnDownWrapper_->fill(patPFMetT1SmearJetEnDownCollectionHandle_->front()) ;
  metT1SmearJetEnDownWrapper_->store(analysis) ;

  metT1SmearJetEnUpWrapper_->fill(patPFMetT1SmearJetEnUpCollectionHandle_->front()) ;
  metT1SmearJetEnUpWrapper_->store(analysis) ;

  metT1SmearJetResDownWrapper_->fill(patPFMetT1SmearJetResDownCollectionHandle_->front()) ;
  metT1SmearJetResDownWrapper_->store(analysis) ;

  metT1SmearJetResUpWrapper_->fill(patPFMetT1SmearJetResUpCollectionHandle_->front()) ;
  metT1SmearJetResUpWrapper_->store(analysis) ;

  metTxyWrapper_->fill(patPFMetTxyHandle_->front()) ;
  metTxyWrapper_->store(analysis) ;

  metFinalWrapper_->fill(patPFMetFinalCollectionHandle_->front()) ;
  metFinalWrapper_->store(analysis) ;
}
void IIHEModuleMET::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleMET::beginEvent(){}
void IIHEModuleMET::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleMET::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleMET);
