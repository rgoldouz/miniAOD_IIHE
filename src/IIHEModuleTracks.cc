#include "UserCode/IIHETree/interface/IIHEModuleTracks.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleTracks::IIHEModuleTracks(const edm::ParameterSet& iConfig): IIHEModule(iConfig){}
IIHEModuleTracks::~IIHEModuleTracks(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleTracks::beginJob(){
  addBranch("trk_n", kUInt) ;
  setBranchType(kVectorFloat) ;
  addBranch("trk_pt") ;
  addBranch("trk_eta") ;
  addBranch("trk_phi") ;
}

// ------------ method called to for each event  ------------
void IIHEModuleTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
/*
// Retrieve primary vertex collection
   const reco::Vertex* pvcoll = parent_->getPrimaryVertices();
// Retrieve electon collection
   reco::GsfElectronCollection electrons = parent_->getElectronCollection() ;


  float mEl = 0.000511;
  std::vector<TLorentzVector> elp4s ;
//boolian to find if there is a pair of electron with mass > 500 GeV
  bool isgoodtosave = false;  
  int ele1index = 0;
  int ele2index = 0;
  float dr=999;
  reco::Vertex::const_iterator PrimaryVertex = pvcoll->begin();


  for(reco::GsfElectronCollection::const_iterator gsfiter=electrons.begin() ; gsfiter!=electrons.end() ; ++gsfiter){
  float px = gsfiter->px() ;
  float py = gsfiter->py() ;
  float pz = gsfiter->pz() ;
  float E = sqrt(mEl*mEl+px*px+py*py+pz*pz) ;
  elp4s.push_back(TLorentzVector(px, py, pz, E)) ;
  }
  for(unsigned int i1=0 ; i1<elp4s.size() ; ++i1){
    for(unsigned int i2=i1+1 ; i2<elp4s.size() ; ++i2){
      TLorentzVector Zeep4 = elp4s .at(i1) + elp4s .at(i2) ;
      float mZee = Zeep4.M() ;
      if (mZee>500) isgoodtosave=true;
    }
  }
  

  if (isgoodtosave && pvcoll->size()>0){
    for(reco::Vertex::const_iterator pvIt = pvcoll->begin(); pvIt!=pvcoll->end(); ++pvIt){
    float dxy1 = electrons[ele1index].gsfTrack()->dxy(pvIt->position());
    float dxy2 = electrons[ele2index].gsfTrack()->dxy(pvIt->position());
    float dz1 =  electrons[ele1index].gsfTrack()->dz(pvIt->position());
    float dz2 =  electrons[ele2index].gsfTrack()->dz(pvIt->position());
    if (sqrt(pow(dxy1,2)+pow(dz1,2)) + sqrt(pow(dxy2,2)+pow(dz2,2)) < dr) {
      dr = sqrt(pow(dxy1,2)+pow(dz1,2)) + sqrt(pow(dxy2,2)+pow(dz2,2));
      PrimaryVertex = pvIt;
      }
    }
    store("trk_n"     , (unsigned int) PrimaryVertex->nTracks()) ;
    for(reco::Vertex::trackRef_iterator it=PrimaryVertex->tracks_begin();  it!=PrimaryVertex->tracks_end(); ++it) {
    store("trk_pt" , (*it)->pt());
    store("trk_eta", (*it)->eta());
    store("trk_phi", (*it)->phi());
    }
  }
*/
}

void IIHEModuleTracks::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleTracks::beginEvent(){}
void IIHEModuleTracks::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleTracks::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleTracks);
