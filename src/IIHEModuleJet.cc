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
  addBranch("jet_neutralHadronEnergyFraction");
  addBranch("jet_neutralEmEnergyFraction");
  addBranch("jet_chargedHadronEnergyFraction");
  addBranch("jet_muonEnergyFraction");
  setBranchType(kVectorInt);
  addBranch("jet_chargedMultiplicity");
  addBranch("jet_neutralMultiplicity");
  setBranchType(kVectorFloat);
  addBranch("jet_TCHE");
  addBranch("jet_TCHP");
  addBranch("jet_JTP");
  addBranch("jet_JBTP ");
  addBranch("jet_SSV");
  addBranch("jet_CSV ");
  addBranch("jet_MSV");
  addBranch("jet_IPM");
  addBranch("jet_SET");
  addBranch("jet_SMT");
  addBranch("jet_SMNoIPT");  



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

    store("jet_TCHE",pfjet->bDiscriminator("trackCountingHighEffBJetTags"));
    store("jet_TCHP",pfjet->bDiscriminator("trackCountingHighPurBJetTags"));
    store("jet_JTP",pfjet->bDiscriminator("jetProbabilityBJetTags"));
    store("jet_JBTP",pfjet->bDiscriminator("jetBProbabilityBJetTags"));
    store("jet_SSV",pfjet->bDiscriminator("simpleSecondaryVertexBJetTags"));
    store("jet_CSV",pfjet->bDiscriminator("combinedSecondaryVertexBJetTags"));
    store("jet_MSV",pfjet->bDiscriminator("combinedSecondaryVertexMVABJetTags"));
    store("jet_IPM",pfjet->bDiscriminator("impactParameterMVABJetTags"));
    store("jet_SET",pfjet->bDiscriminator("softElectronBJetTags"));
    store("jet_SMT",pfjet->bDiscriminator("softMuonBJetTags"));
    store("jet_SMNoIPT",pfjet->bDiscriminator("softMuonNoIPBJetTags"));


  }
}
void IIHEModuleJet::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleJet::beginEvent(){}
void IIHEModuleJet::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleJet::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleJet);
