#include "UserCode/IIHETree/interface/IIHEModuleJet.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

//////////////////////////////////////////////////////////////////////////////////////////
////                             IIHEJetVariable classes
//////////////////////////////////////////////////////////////////////////////////////////
IIHEJetVariableBase::IIHEJetVariableBase(std::string prefix, std::string name, int type){
  name_       = name ;
  branchName_ = prefix + "_" + name_ ;
  branchType_ = type ;
}
bool IIHEJetVariableBase::addBranch(IIHEAnalysis* analysis){
  return analysis->addBranch(branchName_, branchType_) ;
}

IIHEJetVariableInt::IIHEJetVariableInt(std::string prefix, std::string name):
IIHEJetVariableBase(prefix, name, kVectorInt){
  reset() ;
}
void IIHEJetVariableInt::store(IIHEAnalysis* analysis){
  analysis->store(BranchName(), value_) ;
}

IIHEJetVariableFloat::IIHEJetVariableFloat(std::string prefix, std::string name):
IIHEJetVariableBase(prefix, name, kVectorFloat){
  reset() ;
}
void IIHEJetVariableFloat::store(IIHEAnalysis* analysis){
  analysis->store(BranchName(), value_ ) ;
}

//////////////////////////////////////////////////////////////////////////////////////////
////                                  IIHEJet class
//////////////////////////////////////////////////////////////////////////////////////////
IIHEJetWrapper::IIHEJetWrapper(std::string prefix){
  prefix_ = prefix ;

  pt_               = new IIHEJetVariableFloat  (prefix_, "pt"                 ) ;
  energy_              = new IIHEJetVariableFloat  (prefix_, "energy"                ) ;

  variables_.push_back((IIHEJetVariableBase*) pt_                   ) ;
  variables_.push_back((IIHEJetVariableBase*) energy_                  ) ;
}

void IIHEJetWrapper::addBranches(IIHEAnalysis* analysis){
  for(unsigned int i=0 ; i<variables_.size() ; ++i){
    variables_.at(i)->addBranch(analysis) ;
  }
}
void IIHEJetWrapper::reset(){
  for(unsigned int i=0 ; i<variables_.size() ; ++i){
    variables_.at(i)->reset() ;
  }
}

void IIHEJetWrapper::fill(pat::Jet Jet){
  pt_        ->fill(Jet.et()                ) ;
  energy_       ->fill(Jet.energy()               ) ;
}
void IIHEJetWrapper::store(IIHEAnalysis* analysis){
  for(unsigned int i=0 ; i<variables_.size() ; ++i){
    variables_.at(i)->store(analysis) ;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////
////                                  Main IIHEJetModule
////////////////////////////////////////////////////////////////////////////////////////////

IIHEModuleJet::IIHEModuleJet(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): 
  IIHEModule(iConfig),
  jetEnUpWrapper_(new IIHEJetWrapper("jet_EnUp")),
  jetEnDownWrapper_(new IIHEJetWrapper("jet_EnDown")),
  jetSmearedWrapper_(new IIHEJetWrapper("jet_Smeared")),
  jetSmearedJetResUpWrapper_(new IIHEJetWrapper("jet_SmearedJetResUp")),
  jetSmearedJetResDownWrapper_(new IIHEJetWrapper("jet_SmearedJetResDown")){
  pfJetLabel_                  =  iConfig.getParameter<edm::InputTag>("JetCollection");
  pfJetToken_                  =  iC.consumes<View<pat::Jet> > (pfJetLabel_);
  pfJetLabelSmeared_           =  iConfig.getParameter<edm::InputTag>("JetCollectionSmeared");
  pfJetTokenSmeared_           =  iC.consumes<View<pat::Jet> > (pfJetLabelSmeared_);
  pfJetLabelEnUp_              =  iConfig.getParameter<edm::InputTag>("JetCollectionEnUp");
  pfJetTokenEnUp_              =  iC.consumes<View<pat::Jet> > (pfJetLabelEnUp_);
  pfJetLabelEnDown_            =  iConfig.getParameter<edm::InputTag>("JetCollectionEnDown");
  pfJetTokenEnDown_            =  iC.consumes<View<pat::Jet> > (pfJetLabelEnDown_);
  pfJetLabelSmearedJetResUp_   =  iConfig.getParameter<edm::InputTag>("JetCollectionSmearedJetResUp");
  pfJetTokenSmearedJetResUp_   =  iC.consumes<View<pat::Jet> > (pfJetLabelSmearedJetResUp_);
  pfJetLabelSmearedJetResDown_ =  iConfig.getParameter<edm::InputTag>("JetCollectionSmearedJetResDown");
  pfJetTokenSmearedJetResDown_ =  iC.consumes<View<pat::Jet> > (pfJetLabelSmearedJetResDown_);

  ETThreshold_ = iConfig.getUntrackedParameter<double>("jetPtThreshold") ;
  isMC_ = iConfig.getUntrackedParameter<bool>("isMC") ;
}
IIHEModuleJet::~IIHEModuleJet(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleJet::beginJob(){
  addBranch("jet_n", kUInt);

  setBranchType(kVectorFloat);
  addBranch("jet_px");
  addBranch("jet_py");
  addBranch("jet_pz");
  addBranch("jet_pt");
  addBranch("jet_eta");
  addBranch("jet_theta");
  addBranch("jet_phi");
  addBranch("jet_energy");
  addBranch("jet_mass");
  addBranch("jet_chargedEmEnergyFraction");
  addBranch("jet_neutralHadronEnergyFraction");
  addBranch("jet_neutralEmEnergyFraction");
  addBranch("jet_chargedHadronEnergyFraction");
  addBranch("jet_muonEnergyFraction");
  setBranchType(kVectorInt);
  addBranch("jet_chargedMultiplicity");
  addBranch("jet_neutralMultiplicity");
  addBranch("jet_partonFlavour");
  addBranch("jet_hadronFlavour");

  setBranchType(kVectorFloat);
  addBranch("jet_CSV");
  addBranch("jet_CSVv2");
  addBranch("jet_CvsL");
  addBranch("jet_CvsB");
  setBranchType(kVectorBool);
  addBranch("jet_isJetIDLoose");
  addBranch("jet_isJetIDTight");
  addBranch("jet_isJetIDTightLepVeto");

  IIHEAnalysis* analysis = parent_ ;
  if (isMC_){
    jetEnUpWrapper_->addBranches(analysis) ;
    jetEnDownWrapper_->addBranches(analysis) ;
    jetSmearedWrapper_->addBranches(analysis) ;
    jetSmearedJetResUpWrapper_->addBranches(analysis) ;
    jetSmearedJetResDownWrapper_->addBranches(analysis) ;
  }
}

// ------------ method called to for each event  ------------
void IIHEModuleJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::Handle<edm::View<pat::Jet> > pfJetHandle_;
  iEvent.getByToken(pfJetToken_, pfJetHandle_);

  edm::Handle<edm::View<pat::Jet> > pfJetHandleSmeared_;
  iEvent.getByToken(pfJetTokenSmeared_, pfJetHandleSmeared_);

  edm::Handle<edm::View<pat::Jet> > pfJetHandleEnUp_;
  iEvent.getByToken(pfJetTokenEnUp_, pfJetHandleEnUp_);

  edm::Handle<edm::View<pat::Jet> > pfJetHandleEnDown_;
  iEvent.getByToken(pfJetTokenEnDown_, pfJetHandleEnDown_);

  edm::Handle<edm::View<pat::Jet> > pfJetHandleSmearedJetResUp_;
  iEvent.getByToken(pfJetTokenSmearedJetResUp_, pfJetHandleSmearedJetResUp_);

  edm::Handle<edm::View<pat::Jet> > pfJetHandleSmearedJetResDown_;
  iEvent.getByToken(pfJetTokenSmearedJetResDown_, pfJetHandleSmearedJetResDown_);

  jetEnUpWrapper_->reset();
  jetEnDownWrapper_->reset();
  jetSmearedWrapper_->reset();
  jetSmearedJetResUpWrapper_->reset();
  jetSmearedJetResDownWrapper_->reset();

  IIHEAnalysis* analysis = parent_ ;

  store("jet_n", (unsigned int) pfJetHandle_ -> size() );
  for ( unsigned int i = 0; i <pfJetHandle_->size(); ++i) {
    Ptr<pat::Jet> pfjet = pfJetHandle_->ptrAt( i );

//cout<<i<<"  "<<pfjet->pt()<<" Smeared   "<<pfjetSmeared->pt()<<"  JES UP  "<<(pfjetEnUp->pt()- pfjet->pt())/pfjet->pt()<<"  JES DOWN  "<<(pfjetEnDown->pt() - pfjet->pt())/pfjet->pt()<<endl;
    if(pfjet->pt() < ETThreshold_) continue ;

    store("jet_px"    , pfjet->px()) ;
    store("jet_py"    , pfjet->py()) ;
    store("jet_pz"    , pfjet->pz()) ;
    store("jet_pt"    , pfjet->pt()) ;
    store("jet_eta"   , pfjet->eta()) ;
    store("jet_theta" , pfjet->theta()) ;
    store("jet_phi"   , pfjet->phi()) ;
    store("jet_energy", pfjet->energy()) ;
    store("jet_mass"  , pfjet->mass()) ;
    store("jet_neutralHadronEnergyFraction"         ,pfjet->neutralHadronEnergyFraction());
    store("jet_neutralEmEnergyFraction"             ,pfjet->neutralEmEnergyFraction());
    store("jet_chargedHadronEnergyFraction"         ,pfjet->chargedHadronEnergyFraction());
    store("jet_muonEnergyFraction"                  ,pfjet->muonEnergyFraction());
    store("jet_chargedEmEnergyFraction"             ,pfjet->chargedEmEnergyFraction());
    store("jet_chargedMultiplicity"                 ,pfjet->chargedMultiplicity());
    store("jet_neutralMultiplicity"                 ,pfjet->neutralMultiplicity());
    store("jet_partonFlavour"                       ,pfjet->partonFlavour());
    store("jet_hadronFlavour"                       ,pfjet->hadronFlavour());

    store("jet_CSV",pfjet->bDiscriminator("combinedSecondaryVertexBJetTags"));
    store("jet_CSVv2",pfjet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    store("jet_CvsL",pfjet->bDiscriminator("pfCombinedCvsLJetTags"));
    store("jet_CvsB",pfjet->bDiscriminator("pfCombinedCvsBJetTags"));

    float eta = pfjet->eta();
    float NHF = pfjet->neutralHadronEnergyFraction();
    float NEMF = pfjet->neutralEmEnergyFraction();
    float CHF = pfjet->chargedHadronEnergyFraction();
    float MUF = pfjet->muonEnergyFraction();
    float CEMF = pfjet->chargedEmEnergyFraction();
    int NumConst = pfjet->chargedMultiplicity()+pfjet->neutralMultiplicity();
    int NumNeutralParticles =pfjet->neutralMultiplicity();
    int CHM = pfjet->chargedMultiplicity(); 

    bool isJetIDLoose;
    bool isJetIDTight;
    bool isJetIDTightLepVeto;

    if( fabs(eta) <= 2.7){
      isJetIDLoose = ((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4));
      isJetIDTight = ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4));
      isJetIDTightLepVeto = ((NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(eta)>2.4));
      }
    else if (abs(eta)>2.7 && abs(eta)<=3.0 ){
      isJetIDLoose = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 );
      isJetIDTight = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 );
      isJetIDTightLepVeto = false;
      }
    else{
      isJetIDLoose = (NEMF<0.90 && NumNeutralParticles>10 );
      isJetIDTight = (NEMF<0.90 && NumNeutralParticles>10 ) ;
      isJetIDTightLepVeto = false;
    }
    store("jet_isJetIDLoose"                 ,isJetIDLoose);
    store("jet_isJetIDTight"                 ,isJetIDTight);
    store("jet_isJetIDTightLepVeto"          ,isJetIDTightLepVeto);

    if (isMC_){
      Ptr<pat::Jet> pfjetEnUp = pfJetHandleEnUp_->ptrAt( i );
      Ptr<pat::Jet> pfjetEnDown = pfJetHandleEnDown_->ptrAt( i );
      Ptr<pat::Jet> pfjetSmeared = pfJetHandleSmeared_->ptrAt( i );
      Ptr<pat::Jet> pfjetSmearedJetResUp = pfJetHandleSmearedJetResUp_->ptrAt( i );
      Ptr<pat::Jet> pfjetSmearedJetResDown = pfJetHandleSmearedJetResDown_->ptrAt( i );
      jetEnUpWrapper_->fill(pfJetHandleEnUp_->at(i)) ;
      jetEnUpWrapper_ ->store(analysis) ;
      jetEnDownWrapper_->fill(pfJetHandleEnDown_->at(i)) ;
      jetEnDownWrapper_ ->store(analysis) ;
      jetSmearedWrapper_->fill(pfJetHandleSmeared_->at(i));
      jetSmearedWrapper_->store(analysis) ;
      jetSmearedJetResUpWrapper_->fill(pfJetHandleSmearedJetResUp_->at(i));
      jetSmearedJetResUpWrapper_->store(analysis) ;
      jetSmearedJetResDownWrapper_->fill(pfJetHandleSmearedJetResDown_->at(i));
      jetSmearedJetResDownWrapper_->store(analysis) ;
    }

  }


}
void IIHEModuleJet::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleJet::beginEvent(){}
void IIHEModuleJet::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleJet::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleJet);
