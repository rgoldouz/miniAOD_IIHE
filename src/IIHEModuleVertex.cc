#include "UserCode/IIHETree/interface/IIHEModuleVertex.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/LorentzVector.h"
#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleVertex::IIHEModuleVertex(const edm::ParameterSet& iConfig , edm::ConsumesCollector && iC): IIHEModule(iConfig){
  pfToken_ = iC.consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"));
  electronToken_ = iC.consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronCollection"));
}
IIHEModuleVertex::~IIHEModuleVertex(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleVertex::beginJob(){
  addBranch("pv_n", kUInt) ;
  setBranchType(kVectorFloat) ;
  addBranch("pv_x") ;
  addBranch("pv_y") ;
  addBranch("pv_z") ;
  addBranch("pv_isValid", kVectorBool) ;
  addBranch("pv_normalizedChi2", kVectorFloat) ;
  addBranch("pv_ndof", kVectorFloat) ;

  setBranchType(kVectorInt) ;
  addBranch("pf_pdgId") ;
  addBranch("pf_fromPV") ;
  setBranchType(kVectorFloat) ;
  addBranch("pf_px") ;
  addBranch("pf_py") ;
  addBranch("pf_pz") ;
  addBranch("pf_mass") ;
  addBranch("pf_dz") ;

}

// ------------ method called to for each event  ------------
void IIHEModuleVertex::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  // Get the beamspot from the Event:
  // The beamspot is passed to the IIHEAnalysis class, so we call it from parent_
  // Don't forget to declare IIHEModuleVertex as a friend of IIHEAnalysis!

  // Retrieve primary vertex collection
  const reco::VertexCollection* pvcoll = parent_->getPrimaryVertices() ;

  store("pv_n", (unsigned int) pvcoll->size()) ;

 for (vector<reco::Vertex>::const_iterator pvIt = pvcoll->begin(); pvIt != pvcoll->end(); ++pvIt) {
    store("pv_x"             , pvIt->x()                ) ;
    store("pv_y"             , pvIt->y()                ) ;   
    store("pv_z"             , pvIt->z()                ) ;  
    store("pv_isValid"       , pvIt->isValid()          ) ;
    store("pv_ndof"          , pvIt->ndof()        ) ;
    store("pv_normalizedChi2", pvIt->normalizedChi2()   ) ;
  }

  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);

  std::vector<math::XYZTLorentzVector> elp4s;
  for (const pat::Electron &el : *electrons) {
    elp4s.push_back(el.p4());
  }

  bool isGoodToSave = false;
  bool escape = false ;
  for(unsigned int i1=0 ; i1<elp4s.size() ; ++i1){
    if (escape) break;
    for(unsigned int i2=i1+1 ; i2<elp4s.size() ; ++i2){
      math::XYZTLorentzVector Zeep4 = elp4s .at(i1) + elp4s .at(i2) ;
      float mZee = Zeep4.M() ;
      if (mZee>500) {
        isGoodToSave = true;
        escape = true;
        break;
      }
    }
  }

  
    
  //for more info please look at the https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#Packed_ParticleFlow_Candidates
  if (isGoodToSave){
    edm::Handle<pat::PackedCandidateCollection> pfs;
    iEvent.getByToken(pfToken_, pfs);
    for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
      const pat::PackedCandidate &pf = (*pfs)[i];
      if (pf.fromPV() <= 2) continue;
      store("pf_pdgId"          ,pf.pdgId());
      store("pf_fromPV"         ,pf.fromPV());
      store("pf_px"             ,pf.px());
      store("pf_py"             ,pf.py());
      store("pf_pz"             ,pf.pz());
      store("pf_mass"           ,pf.mass());
      store("pf_dz"             ,pf.dz());
    }
  }
}

void IIHEModuleVertex::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleVertex::beginEvent(){}
void IIHEModuleVertex::endEvent(){}

// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleVertex::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleVertex);
