#include "UserCode/IIHETree/interface/IIHEModuleJet.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleJet::IIHEModuleJet(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC): IIHEModule(iConfig){
  pfJetLabel_ =  iConfig.getParameter<edm::InputTag>("JetCollection");
  pfJetToken_ =  iC.consumes<View<pat::Jet> > (pfJetLabel_);
  ETThreshold_ = iConfig.getUntrackedParameter<double>("jetPtThreshold") ;
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
  setBranchType(kVectorFloat);
  addBranch("jet_CSV");
  addBranch("jet_CSVv2");
  addBranch("jet_CvsL");
  addBranch("jet_CvsB");
  setBranchType(kVectorBool);
  addBranch("jet_isJetIDLoose");
  addBranch("jet_isJetIDTight");
  addBranch("jet_isJetIDTightLepVeto");
}

// ------------ method called to for each event  ------------
void IIHEModuleJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::Handle<edm::View<pat::Jet> > pfJetHandle_;
  iEvent.getByToken(pfJetToken_, pfJetHandle_);

  store("jet_n", (unsigned int) pfJetHandle_ -> size() );
  for ( unsigned int i = 0; i <pfJetHandle_->size(); ++i) {
    Ptr<pat::Jet> pfjet = pfJetHandle_->ptrAt( i );
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
      isJetIDLoose = (NEMF<0.90 && NumNeutralParticles>2 );
      isJetIDTight = (NEMF<0.90 && NumNeutralParticles>2 );
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
  }
}
void IIHEModuleJet::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleJet::beginEvent(){}
void IIHEModuleJet::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleJet::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleJet);
