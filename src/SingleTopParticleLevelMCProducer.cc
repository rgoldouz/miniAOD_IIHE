
// Top MC particle-level object producer
// This produces sets of particle-level objects (leptons, jets, neutrinos) with different definitions
// There are different possibilities for defining particle-level objects:
// 1) Leptons: one should use only leptons stemming from resonances (W,Z,a*) and not hadrons
// 1.1) One has the option to "dress" leptons with photons, i.e. defining collections of leptons and photons to use ak1 jet clustering onto
// 2) Jets: should be made out of the non-leptons and non-photons for dressing.
// 3) b-Jets: one possibility is to define b-jets clustering the "rescaled b-hadrons" instead of their decay products. For this definition, a b-jet will be a jet that has a rescaled b-hadrons amongst its component
// 3.1) Another possibility is to use the "rescaled b-hadron descendants", i.e. not modify the jet clustering but just going back to the b-decay products for identification.


#define DEBUG    0 // 0=false
//#define DEBUG    1 // 0=false

#include "PhysicsTools/PatAlgos/plugins/PATJetProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackProbabilityTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackCountingTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "PhysicsTools/JetMCUtils/interface/CandMCTag.h"
#include "TVector3.h"
////////Part to be made official
#include "UserCode/IIHETree/interface/SingleTopParticleLevelMCProducer.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include <vector>
#include <memory>
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"

class isFJBHadron : public fastjet::PseudoJet::UserInfoBase{
public:
  isFJBHadron(const bool & isBh):isB_(isBh){;}
  const bool isB() const {return isB_;}
protected:
  bool isB_;
};

//using namespace pat;
typedef boost::shared_ptr<fastjet::ClusterSequence>  ClusterSequencePtr;
typedef boost::shared_ptr<fastjet::JetDefinition>    JetDefPtr;

SingleTopParticleLevelMCProducer::SingleTopParticleLevelMCProducer(const edm::ParameterSet& iConfig) 
{
  // initialize the configurables
  genParticlesSrc_ = consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>( "genParticlesSource" ));
  genJetsSrc_ = consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>( "genJetsSource" ));
  
  if (iConfig.exists("ignoreParticleIDs"))
    setIgnoredParticles(iConfig.getParameter<std::vector<unsigned int> >("ignoreParticleIDs"));
    setExcludeFromResonancePids(iConfig.getParameter<std::vector<unsigned int> >("excludeFromResonancePids"));
  
  doDressedCollections_          = iConfig.getUntrackedParameter<bool>("doDressedCollections",true); 
  doBaseCollections_             = iConfig.getUntrackedParameter<bool>("doBaseCollections",true); 
  doBDescendentJetSelection_     = iConfig.getUntrackedParameter<bool>("doBDescendentJetSelection",true);
  doRescaledBHadronJetSelection_ = iConfig.getUntrackedParameter<bool>("doRescaledBHadronJetSelection",false);

  useJetsNoNu_     = iConfig.getUntrackedParameter<bool>("useJetsNoNu",false); 
  jetDefinition_   = iConfig.getUntrackedParameter<std::string>("jetDefinition","ak4");
  minJetPt_        = iConfig.getUntrackedParameter<double>("minJetPt",3.0);
  rParam_          = iConfig.getUntrackedParameter<double>("rParam",0.4);
  rescaling_       = iConfig.getUntrackedParameter<double>("rescaling",1e-18);
  
  if(doBaseCollections_){
  //Part 1: particle-level defined "bare" objects, no dressing:
  //Leptons from resonances ("bare" leptons):
    produces<std::vector<reco::GenParticle> >("particleLevelLeptons").setBranchAlias("particleLevelLeptons");

  //Jets using the plain resonance definition:
    produces <reco::GenParticleRefVector> ("genParticlesForJetsNoRescaledBHadronsNoDressLeptons").setBranchAlias("genParticlesForJetsNoRescaledBHadronsNoDressLeptons");
    produces <reco::GenParticleRefVector> ("genParticlesForJetsRescaledBHadronsNoDressLeptons").setBranchAlias("genParticlesForJetsRescaledBHadronsNoDressLeptons");
  }
  
  if(doDressedCollections_){
    //Part 2: particles for dressing:
    //Leptons and photons for dressing:
    produces <reco::GenParticleRefVector> ("genLeptonsForDressing").setBranchAlias("genEGammaForDressing");
    produces <reco::GenParticleRefVector> ("genEGammaForDressing").setBranchAlias("genEGammaForDressing");
    produces <reco::GenParticleRefVector> ("genMuGammaForDressing").setBranchAlias("genMuGammaForDressing");
 
    // using the lepton dressing definition
    //Jets: Adding Rescaled B-Hadrons for particle-level b-jet definition:
    produces <reco::GenParticleRefVector> ("genParticlesForJetsNoRescaledBHadrons").setBranchAlias("genParticlesForJetsNoRescaledBHadrons");
    
    //Jets without rescaled b-hadrons:   
    produces <reco::GenParticleRefVector> ("genParticlesForJetsRescaledBHadrons").setBranchAlias("genParticlesForJetsRescaledBHadrons");
  }


  //Part 3: selecting jets using the b-descendent definition
  if(doBDescendentJetSelection_ || doRescaledBHadronJetSelection_){
    //Selection of genJets, using the "b-descendent" definition:
    produces<std::vector<reco::GenJet > >("particleLevelJets").setBranchAlias("particleLevelJets");
    produces<std::vector<reco::GenJet > >("particleLevelBJets").setBranchAlias("particleLevelBJets");
    //Note: if there is a rescaled b-hadron, it will automatically be a b-descendent.
  }
  
  if(doBaseCollections_){
    //Neutrinos from resonances and MET = SUM of nu-4 momenta
    produces<std::vector<reco::GenParticle> >("particleLevelNeutrinos").setBranchAlias("particleLevelNeutrinos");
    //    produces<std::vector<reco::NamedCompositeCandidate> >("particleLevelMET").setBranchAlias("particleLevelMET");
  }
  
  size_t nignored = ignoreParticleIDs.size();
  for (size_t ig = 0;  ig < nignored;++ig) {
//    std::cout << "ignored# "<< ig <<":"<< ignoreParticleIDs.at(ig)<<std::endl;
  }
  
  size_t nexcluded = excludeFromResonancePids.size();
  for (size_t ig = 0;  ig < nexcluded;++ig) {
//    std::cout << "ignored# "<< ig <<":"<< excludeFromResonancePids.at(ig)<<std::endl;
  }
}

void SingleTopParticleLevelMCProducer::produce(edm::Event & iEvent, const edm::EventSetup & iEventSetup){
  using namespace std;
  using namespace edm;
  using namespace reco;

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesSrc_, genParticles) ;

  edm::Handle<std::vector<reco::GenJet> > genJets;
  iEvent.getByToken(genJetsSrc_, genJets);

  
  std::auto_ptr< std::vector< reco::GenParticle > > particleLevelLeptons( new std::vector<reco::GenParticle> );
  std::auto_ptr<reco::GenParticleRefVector> genLeptonsForDressing (new reco::GenParticleRefVector);
  std::auto_ptr<reco::GenParticleRefVector> genEGammaForDressing (new reco::GenParticleRefVector);
  std::auto_ptr<reco::GenParticleRefVector> genMuGammaForDressing (new reco::GenParticleRefVector);
  std::auto_ptr< std::vector< reco::GenJet > > particleLevelJets( new std::vector<reco::GenJet> );
  std::auto_ptr< std::vector< reco::GenJet > > particleLevelBJets( new std::vector<reco::GenJet> );
  std::auto_ptr<reco::GenParticleRefVector> genParticlesForJetsRescaledBHadronsNoDressLeptons(new reco::GenParticleRefVector);
  std::auto_ptr<reco::GenParticleRefVector> genParticlesForJetsNoRescaledBHadronsNoDressLeptons(new reco::GenParticleRefVector);
  std::auto_ptr<reco::GenParticleRefVector> genParticlesForJetsRescaledBHadrons(new reco::GenParticleRefVector);
  std::auto_ptr<reco::GenParticleRefVector> genParticlesForJetsNoRescaledBHadrons(new reco::GenParticleRefVector);
  std::auto_ptr< std::vector< reco::GenParticle > > particleLevelNeutrinos( new std::vector<reco::GenParticle> );

  std::map<const reco::GenParticle*,size_t> particlePtrIdxMap;
  std::map<reco::GenParticle,size_t> bParticlePtrIdxMap;

  ParticleVector particles,bHadParticles;

  for (reco::GenParticleCollection::const_iterator iter=genParticles->begin();iter!=genParticles->end();++iter){
    particles.push_back(&*iter);
    particlePtrIdxMap[&*iter] = (iter - genParticles->begin());
  }
  std::sort(particles.begin(), particles.end());
  unsigned int size = particles.size();
  ParticleBitmap invalidFromResonance(size,false);
  ParticleBitmap selectedForNeutrinos(size, false), invalidForNeutrinos(size, false);
  ParticleBitmap selectedForBJets(size, false), invalidForBJets(size, false);
  ParticleBitmap selectedForJets(size, false), invalidForJets(size, false);
  ParticleBitmap selectedForLeptons(size, false), invalidForLeptons(size, false);
  ParticleBitmap selectedForMuGamma(size, false), invalidForMuGamma(size, false);
  ParticleBitmap selectedForEGamma(size, false), invalidForEGamma(size, false);
    
  for(unsigned int i = 0; i < size; i++) {
    const reco::GenParticle *particle = particles[i];
    if (isBHadronRescaled(particle,particles))bHadParticles.push_back(particle);
  }
  
  for(unsigned int i = 0; i < size; i++) {
    const reco::GenParticle *particle = particles[i];
    int pdgid = particle->pdgId();

    bool selected = false; //bool invalid = false;
    
    selectedForMuGamma[i]= false; selectedForEGamma[i]= false; selectedForJets[i]= false; selectedForLeptons[i]= false; selectedForBJets[i]=false;
    invalidForMuGamma[i]= false; invalidForEGamma[i]= false; invalidForJets[i]= false; invalidForLeptons[i]= false; invalidForBJets[i]=false;
    
    //Basic selection of final state particles: not always good, e.g. rescaled b-hadrons have daughters.
    if (particle->status() == 1 ) selected = true;
    bool isRescaledBHadron        = isBHadronRescaled(particle,particles);
    bool isRescaledBHadronProduct = isBHadronRescaledDecayProduct(particle,particles,bHadParticles);
    ResonanceState rstate         = isFromResonance(particle,particles,invalidFromResonance);
    bool fromHadron               = isFromHadron(particle,particles);
    bool isnotfromresonance       = (!rstate);
    
    if(selected && !isnotfromresonance)invalidFromResonance[i]=true;
    
    //If comes from a rescaled b hadron: remove it from lepton list
    if(isRescaledBHadronProduct && false) std::cout << " " <<std::endl;
    
    //if is a prompt lepton
    if(selected && !fromHadron) {//Note: resonances do not include taus at this stage!
//            std::cout<< "selected particle id "<< pdgid<< std::endl;
      //      if(!isnotfromresonance ){
      if(fabs(pdgid)==13 ) selectedForMuGamma[i]= true;
      if(fabs(pdgid)==11 ) selectedForEGamma[i]= true;
	//      }
      if(fabs(pdgid)==22 ) {selectedForMuGamma[i]= true;selectedForEGamma[i]= true;}
      selectedForLeptons[i]= selectedForMuGamma[i] || selectedForEGamma[i];
      
      if(fabs(pdgid)== 12 || fabs(pdgid)==14 || fabs(pdgid)==16){ selectedForNeutrinos[i]=true;
	//	std::cout << " neutrino found"<<std::endl;
      }
    }
    else{
      invalidForMuGamma[i]= true; invalidForEGamma[i]= true;
      invalidForLeptons[i]= true; invalidForNeutrinos[i]= true;
    }

	

    
    //if it's a prompt lepton it should not be in a jet
    if(selectedForLeptons[i]){ 
      invalidForJets[i] = true; 
      invalidForBJets[i] =true;
    }

   
    //   std::cout << pdgid<<" isignored? "<< isIgnored(pdgid)<< " is from resonance? "<< !isnotfromresonance <<std::endl;
    selected = selected && !(isIgnored(pdgid) && particle->pt()>=0);    

    //Jet selection: remove resonances
    if(fromHadron){// we need to  choose the ones that are not from hadrons. If you ask "! fro resonance" you'll get the leptons from taus
      if(selected) selectedForJets[i] = true;
      else invalidForJets[i] = true; 
    
      if (useJetsNoNu_){
	if(selected &&  (fabs(pdgid)== 12 || fabs(pdgid)==14 || fabs(pdgid)==16)){
	  selectedForJets[i] = false;
	  invalidForJets[i] = true;
	  selectedForNeutrinos[i]=false;
	  invalidForNeutrinos[i]=true;
	}
      }
      //bHadrons: adding rescaled hadron definition
      //if((selected && !isRescaledBHadronProduct))selectedForBJets[i] = true;
      if(isRescaledBHadron) selectedForBJets[i] = true;
      else invalidForBJets[i]= true;
      
    }
    else {
      invalidForJets[i] = true; 
      invalidForBJets[i] = true; 
    }
  }
  
  for(size_t idx = 0; idx < size; ++idx){   
    const reco::GenParticle *particle = particles[idx];
    if (!selectedForJets[idx] || invalidForJets[idx]){
      continue;
    }
    edm::Ref<reco::GenParticleCollection> particleRef(genParticles, particlePtrIdxMap[particle]);
    genParticlesForJetsNoRescaledBHadrons->push_back(particleRef);
    
  }
	
  for(size_t idx = 0; idx < size; ++idx){   
    const reco::GenParticle *particle = particles[idx];
    if (!selectedForBJets[idx] || invalidForBJets[idx]){
      continue;
    }
    edm::Ref<reco::GenParticleCollection> particleRef(genParticles,particlePtrIdxMap[particle]);
    genParticlesForJetsRescaledBHadrons->push_back(particleRef);
  }
  //Add rescaled b-hadrons:
  //  for(size_t idx = 0; idx < sizeb; ++idx){   
  //    reco::GenParticle particle = bHadParticlesToAddC[idx];
  //    edm::Ref<reco::GenParticleCollection> dynamic_cast< edm::Ref<reco::GenParticleCollection> > particle.ref();
  //    // particleRef(bHadParticlesToAddC,bParticlePtrIdxMap[particle]);
  //    genParticlesForJetsRescaledBHadrons->push_back(particleRef);
  //}
  
  //Lepton Collection Loop
  for(size_t idx = 0; idx < size; ++idx){   
    const reco::GenParticle *particle = particles[idx];
//    if (particle->pdgId()== 22 && particle->pt()>30 && particle->status()==3) std::cout<< "GEN photon   "<<particle->pt()<<std::endl;
//    if (particle->pdgId()== 22 && particle->pt()>30 && particle->status()==1) {
//std::cout<< "reco photon   "<<particle->pt()<<std::endl;
//std::cout<< isFromHadron(particle,particles)<<std::endl;
//}
    if (!selectedForLeptons[idx] || invalidForLeptons[idx]){
      continue;
    }
    edm::Ref<reco::GenParticleCollection> particleRef(genParticles,particlePtrIdxMap[particle]);
    genLeptonsForDressing->push_back(particleRef);
    particleLevelLeptons->push_back(*particle);
//std::cout<< particle->pdgId()<<std::endl;
  
}

  //Neutrino Collection Loop
  for(size_t idx = 0; idx < size; ++idx){   
    const reco::GenParticle *particle = particles[idx];
    if (!selectedForNeutrinos[idx] || invalidForNeutrinos[idx]){
      continue;
    }
    particleLevelNeutrinos->push_back(*particle);
  }
  //EGamma Collection Loop
  for(size_t idx = 0; idx < size; ++idx){   
    const reco::GenParticle *particle = particles[idx];
    if (!selectedForEGamma[idx] || invalidForEGamma[idx]){
      continue;
    }
    edm::Ref<reco::GenParticleCollection> particleRef(genParticles,particlePtrIdxMap[particle]);
    genEGammaForDressing->push_back(particleRef);
  }
  
  //MuGamma Collection Loop
  for(size_t idx = 0; idx < size; ++idx){   
    const reco::GenParticle *particle = particles[idx];
    if (!selectedForMuGamma[idx] || invalidForMuGamma[idx]){
      continue;
    }
    edm::Ref<reco::GenParticleCollection> particleRef(genParticles,particlePtrIdxMap[particle]);
    genMuGammaForDressing->push_back(particleRef);
  }
  
  for (reco::GenParticleCollection::const_iterator t = genParticles->begin (); t != genParticles->end (); ++t){
    ;    
  }
  if(doBDescendentJetSelection_){
    
    //    std::cout << "#genJets collection is "<<  genJets->size()<<std::endl;
    for (size_t j = 0; j < genJets->size();++j){
//      std::vector<const reco::Candidate*> constituents = genJets->at(j).getJetConstituentsQuick();
      //      std::cout << "ak5GenJet #"<<j+1<<" #constituents is "<< constituents.size()<<std::endl;
      bool isBJet = false;
      for(size_t c = 0; c < genJets->at(j).numberOfDaughters();++c){
	//	std::cout << "constituent #"<< c+1 << " pdgId "<< constituents.at(c)->pdgId() <<std::endl;
	
	const reco::GenParticle *particle = dynamic_cast <const reco::GenParticle*>(genJets->at(j).daughter(c));  
	bool isBHadronDP=isBHadronRescaledDecayProduct(particle,particles,bHadParticles);
	//     if(isBHadron(constituents.at(c)->pdgId()))isBJet = true;
	if(isBHadronDP)isBJet = true;
	
	//	std::cout << "is b hadron? "<< isBHadron(constituents.at(c)->pdgId())<< std::endl;
	//	std::cout << "has rescaled b-hadron product? "<< isBHadronDP <<std::endl;
      }
      if (isBJet)particleLevelBJets->push_back(genJets->at(j));
      particleLevelJets->push_back(genJets->at(j));
      //      std::cout << " is it a b jet? "<< isBJet<<std::endl;
    }
  }
  if( doRescaledBHadronJetSelection_){
    std::vector<fastjet::PseudoJet> fjInputs;
    for (size_t j = 0; j < genJets->size();++j){
      if (fabs(genJets->at(j).eta()) >5.5) continue;
      for(size_t c = 0; c < genJets->at(j).numberOfDaughters();++c){
        if(genJets->at(j).daughter(c)->pt() > 0){
          fjInputs.push_back(fastjet::PseudoJet(genJets->at(j).daughter(c)->px(),genJets->at(j).daughter(c)->py(),genJets->at(j).daughter(c)->pz(),genJets->at(j).daughter(c)->energy()));
  	}
      }
    }
    size_t sizeb = bHadParticles.size();
    for (size_t b = 0; b < sizeb;++b){
      const reco::GenParticle *particle = dynamic_cast <const reco::GenParticle*>(bHadParticles.at(b));
      if(particle->pt()!=0  && fabs(particle->eta())<5.5){
	fastjet::PseudoJet pj = (fastjet::PseudoJet(particle->px(),particle->py(),particle->pz(),particle->energy()))*rescaling_;
	const bool isBh= true;
	pj.set_user_info(new isFJBHadron(isBh));
	fjInputs.push_back(pj);
	  
      }
    }  
    JetDefPtr  fjJetDefinition_= JetDefPtr( new fastjet::JetDefinition(fastjet::antikt_algorithm, rParam_) );
    ClusterSequencePtr fjClusterSeq_ = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs, *fjJetDefinition_ ) );
    std::vector<fastjet::PseudoJet> inclusiveJets = fastjet::sorted_by_pt( fjClusterSeq_->inclusive_jets(minJetPt_) );
    
    std::vector<int> reclusteredIndices;
    //    std::cout << " gen jet size "<< genJets->size() << " recl size "<< inclusiveJets.size() <<std::endl; 
//    if( inclusiveJets.size() < genJets->size() )std::cout<< " less recl. jets than original!" <<std::endl;
    matchReclusteredJets(genJets,inclusiveJets,reclusteredIndices);
   size_t jInd = 0;
    for( size_t j = 0;j< genJets->size();++j){
      if(fabs(genJets->at(j).eta()) >5.5) continue;
      bool isBJet = false;
      if( reclusteredIndices.at(jInd) < 0 ){
	//	std::cout<< " jet "<<j<<" has no recl jet "<<std::endl;
	continue;
	  }
      else{
	//	std::cout << "studying recl jets for j ="<< j<<std::endl;
	std::vector<fastjet::PseudoJet> constituents = fastjet::sorted_by_pt( inclusiveJets.at(reclusteredIndices.at(jInd)).constituents() );
        ++jInd;
	for(std::vector<fastjet::PseudoJet>::const_iterator it = constituents.begin(); it != constituents.end(); ++it){
	  if( !it->has_user_info() ) continue; // skip if not a "ghost"
	  else{
	    if(it->user_info<isFJBHadron>().isB()){
	      //	      std::cout << "fastjet #" << j << "is b jet! "<<std::endl;
	      isBJet=true;
	    }
	  }
	}
      }
      if (isBJet)particleLevelBJets->push_back(genJets->at(j));
      particleLevelJets->push_back(genJets->at(j));
    }
  }
  
  //cout<<"AT THE END"<<endl;
  // iEvent.put(particleLevelTops,"particleLevelTops");
 if(doBaseCollections_){
   iEvent.put(particleLevelLeptons,"particleLevelLeptons");
   iEvent.put(genParticlesForJetsNoRescaledBHadronsNoDressLeptons,"genParticlesForJetsNoRescaledBHadronsNoDressLeptons");
   iEvent.put(genParticlesForJetsRescaledBHadronsNoDressLeptons,"genParticlesForJetsRescaledBHadronsNoDressLeptons");
 }
 
 
 if(doDressedCollections_){
  iEvent.put(genLeptonsForDressing,"genLeptonsForDressing");
  iEvent.put(genEGammaForDressing,"genEGammaForDressing");
  iEvent.put(genMuGammaForDressing,"genMuGammaForDressing");
  iEvent.put(genParticlesForJetsRescaledBHadrons,"genParticlesForJetsRescaledBHadrons");
  iEvent.put(genParticlesForJetsNoRescaledBHadrons,"genParticlesForJetsNoRescaledBHadrons");
 }

 if(doBaseCollections_){
   iEvent.put(particleLevelNeutrinos,"particleLevelNeutrinos");
   //   iEvent.put(particleLevelMET, "particleLevelMET");
 }  
 if(doBDescendentJetSelection_ || doRescaledBHadronJetSelection_){
   iEvent.put(particleLevelBJets,"particleLevelBJets");
   iEvent.put(particleLevelJets,"particleLevelJets");
 }
 
}


SingleTopParticleLevelMCProducer::ResonanceState SingleTopParticleLevelMCProducer::isFromResonance( const reco::GenParticle * particle, const ParticleVector &allParticles, ParticleBitmap &invalid){

  unsigned int partIdx(const SingleTopParticleLevelMCProducer::ParticleVector &p,
			      const reco::GenParticle *particle);

  int id = particle->pdgId();
  unsigned int idx = partIdx(allParticles, particle);
  //  std::cout << " id is "<< id <<std::endl;
  
  if (isResonance(id) && (particle->status() == 23 || particle->status() == 22) ){  
    //    std::cout << " id "<< id <<" is direct resonance "<<std::endl; 
    return kDirect;
  }
  if (invalid[idx]) return kIndirect;
  
  if(!isIgnored(id) && isParton(id)) {
    //    std::cout << " isignored or isparton kNo "<<std::endl;
    return kNo;
  }  
  unsigned int nMo=particle->numberOfMothers();
  if (!nMo)
    {//  std::cout << " has no mom kNo "<<std::endl;

      return kNo;
    }
  for(unsigned int i=0;i<nMo;++i){
    ResonanceState result = isFromResonance(dynamic_cast<const reco::GenParticle*>(particle->mother(i)),allParticles,invalid);
    switch(result){
    case kNo: 
      //      std::cout << " mother " <<i <<"is kNo "<<std::endl;
      break;
    case kDirect:
      if (dynamic_cast<const reco::GenParticle*>(particle->mother(i))->pdgId()==id || isResonance(id))
	{
	  //  std::cout << " mother "<< i << " id "<< dynamic_cast<const reco::GenParticle*>(particle->mother(i))->pdgId() <<" is resonance " <<std::endl;
	  return kDirect;
	}
      if(!isExcludedFromResonance(id))break;
    case kIndirect: 
      return kIndirect;
    }
     
  }
  //  std::cout << "reached end kNo "<<std::endl;
  return kNo;
}

bool SingleTopParticleLevelMCProducer::isParton(int pdgId)
{
  pdgId = (pdgId > 0 ? pdgId : -pdgId) % 10000;
  return (pdgId > 0 && pdgId < 6) ||
    pdgId == 9 || pdgId == 21;
  // tops are not considered "regular" partons
  // but taus eventually are (since they may hadronize later)
   
}
bool SingleTopParticleLevelMCProducer::isIgnored(int pdgId)
{
  
  pdgId = pdgId > 0 ? pdgId : -pdgId;
  std::vector<unsigned int>::const_iterator pos =
    std::lower_bound(ignoreParticleIDs.begin(),
		     ignoreParticleIDs.end(),
		     (unsigned int)pdgId);
  return pos != ignoreParticleIDs.end() && *pos == (unsigned int)pdgId;
  
}

bool SingleTopParticleLevelMCProducer::isExcludedFromResonance(int pdgId)
{
  pdgId = pdgId > 0 ? pdgId : -pdgId;
  std::vector<unsigned int>::const_iterator pos =
    std::lower_bound(excludeFromResonancePids.begin(),
		     excludeFromResonancePids.end(),
		     (unsigned int)pdgId);
  return pos != excludeFromResonancePids.end() && *pos == (unsigned int)pdgId;
  
}

unsigned int partIdx(const SingleTopParticleLevelMCProducer::ParticleVector &p,
			    const reco::GenParticle *particle){
  SingleTopParticleLevelMCProducer::ParticleVector::const_iterator pos =
    std::lower_bound(p.begin(), p.end(), particle);
  if (pos == p.end() || *pos != particle)
    throw cms::Exception("CorruptedData")
      << "reco::GenEvent corrupted: Unlisted particles"
      " in decay tree." << std::endl;
  
  return pos - p.begin();
}

bool SingleTopParticleLevelMCProducer::isResonance(int pdgId)
{
  // gauge bosons and tops
  pdgId = (pdgId > 0 ? pdgId : -pdgId) % 10000;
  return (pdgId > 21 && pdgId <= 42) || pdgId == 6 || pdgId == 7 || pdgId == 8 ;  //BUG! was 21. 22=gamma..
}

bool SingleTopParticleLevelMCProducer::isPromptLepton( const reco::GenParticle * particle, const ParticleVector &allParticles){
  return true;
}

bool SingleTopParticleLevelMCProducer::isBHadronRescaled( const reco::GenParticle * particle, const ParticleVector &allParticles){

  if(isBHadron(particle->pdgId())){
    size_t ndaughters = particle->numberOfDaughters();
    bool hasBHadronDaughters = false;
    //    std::cout<< " inside b rescaled func, particle id "<< particle->pdgId() << " pt "<< particle->pt()<< " daughters "<< ndaughters<<std::endl;
    for (size_t d = 0; d < ndaughters;++d ){
      //#if DEBUG
      //std::cout<< " daughterd "<< d << " id "<< particle->daughter(d)->pdgId()<<std::endl; 
      //#endif
    if(isBHadron(particle->daughter(d)->pdgId())) hasBHadronDaughters = true;
    }
    if(!hasBHadronDaughters && particle->pt()>5.0){
      //#if DEBUG
      //      std::cout << " is rescaled b "<<std::endl;
      //#endif
      return true;
    }
  }
  return false;
}

bool SingleTopParticleLevelMCProducer::isHadron( int pdgId ){
  pdgId = (pdgId > 0 ? pdgId : -pdgId) % 10000;
  return (pdgId > 100 && pdgId < 900) || (pdgId > 1000 && pdgId < 9000);
}

bool SingleTopParticleLevelMCProducer::isFromHadron( const reco::GenParticle * particle, const ParticleVector &allParticles){
  bool result = false;
  int pdgid = particle->pdgId();
  size_t nmothers = particle->numberOfMothers();
  result = isHadron(pdgid)|| isParton(pdgid);
  //  std::cout<< "id "<< pdgid << " nmothers "<< nmothers << " result "<<result<<std::endl;
//if(abs(pdgid) ==13)std::cout << result <<std::endl;
 
  if(result) return true;
  for (size_t m = 0; m < nmothers;++m ){
    int motherid =  dynamic_cast<const reco::GenParticle*>(particle->mother(m))->pdgId();
    //    std::cout<< " mother m " << m << " id "<< motherid<<std::endl;
    //    if(motherid ==24)std::cout << " mother status " << dynamic_cast<const reco::GenParticle*>(particle->mother(m))->status()<< " mother from had? " << isFromHadron(dynamic_cast<const reco::GenParticle*>(particle->mother(m)), allParticles)<<std::endl;

//if (abs(pdgid) ==22&& particle->pt()>30) std::cout <<" particle->status()     "<< particle->status() <<"    particle->pdgId()    " <<particle->pdgId() <<"     motherid   "<<motherid<<std::endl;

    if( isHadron (motherid)|| isParton(motherid)) { result = true; break;}
    
    if( ((motherid == pdgid) || (abs(motherid)==15)) && dynamic_cast<const reco::GenParticle*>(particle->mother(m))->status()!=23){//Hadrons in taus are also considered
//if(abs(pdgid) ==22)std::cout << "GHABL   "<<result <<std::endl;   

      result = result || isFromHadron(dynamic_cast<const reco::GenParticle*>(particle->mother(m)), allParticles);
//if(abs(pdgid) ==22)std::cout << "BAAD   "<<result <<std::endl;
    }
  }
//if(abs(pdgid) ==13)std::cout << result <<std::endl;
  return result;
}


bool SingleTopParticleLevelMCProducer::isBHadronRescaledDecayProduct( const reco::GenParticle * particle, const ParticleVector &allParticles, const ParticleVector &rescaledBHadrons){
  bool result = false;
  size_t nrescaled = rescaledBHadrons.size();
  size_t nmothers = particle->numberOfMothers();
  //  std::cout<< " inside b rescaled product func, particle id "<< particle->pdgId() << " pt "<< particle->pt()<< " nmothers " << nmothers << std::endl; 
  for(size_t r=0;r<nrescaled;++r){
    if( rescaledBHadrons.at(r)==particle) return true;
  } 
  
  for (size_t m = 0; m < nmothers;++m ){
    // std::cout<< " mother m " << m << std::endl;
    result = result || isBHadronRescaledDecayProduct(dynamic_cast <const reco::GenParticle*>(particle->mother(m)), allParticles, rescaledBHadrons);
    if(result) return result;
  }
  
  return result;
}

bool SingleTopParticleLevelMCProducer::isBHadron(int pdgid){
  //Bottom mesons
  int code1;
  int code2;
  bool tmpHasBottom = false;
  code1 = (int)( ( abs(pdgid ) / 100)%10 );
  code2 = (int)( ( abs(pdgid ) /1000)%10 );
  if ( code1 == 5 || code2 == 5) tmpHasBottom = true;
  return tmpHasBottom;

  if (fabs(pdgid)==511) return true;
  if (fabs(pdgid)==521) return true;
  if (fabs(pdgid)==10511) return true;
  if (fabs(pdgid)==10521) return true;


  if (fabs(pdgid)==513) return true;
  if (fabs(pdgid)==523) return true;
  if (fabs(pdgid)==20513) return true;
  if (fabs(pdgid)==20523) return true;


  if (fabs(pdgid)==515) return true;
  if (fabs(pdgid)==525) return true;

  if (fabs(pdgid)==531) return true;
  if (fabs(pdgid)==10531) return true;

  if (fabs(pdgid)==533) return true;
  if (fabs(pdgid)==10533) return true;
  if (fabs(pdgid)==20533) return true;

  if (fabs(pdgid)==535) return true;

  if (fabs(pdgid)==541) return true;
  if (fabs(pdgid)==10541) return true;

  if (fabs(pdgid)==543) return true;
  if (fabs(pdgid)==10543) return true;
  if (fabs(pdgid)==20543) return true;

  if (fabs(pdgid)==545) return true;

  //BBbar mesons

  if (fabs(pdgid)==551) return true;
  if (fabs(pdgid)==10551) return true;
  if (fabs(pdgid)==100551) return true;
  if (fabs(pdgid)==110551) return true;
  if (fabs(pdgid)==200551) return true;
  if (fabs(pdgid)==210551) return true;


  if (fabs(pdgid)==553) return true;
  if (fabs(pdgid)==10553) return true;

  if (fabs(pdgid)==553) return true;
  if (fabs(pdgid)==10553) return true;
  if (fabs(pdgid)==20553) return true;
  if (fabs(pdgid)==30553) return true;

  if (fabs(pdgid)==100553) return true;
  if (fabs(pdgid)==110553) return true;
  if (fabs(pdgid)==120553) return true;
  if (fabs(pdgid)==130553) return true;

  if (fabs(pdgid)==200553) return true;
  if (fabs(pdgid)==210553) return true;
  if (fabs(pdgid)==220553) return true;

  if (fabs(pdgid)==300553) return true;

  if (fabs(pdgid)==9000553) return true;
  if (fabs(pdgid)==9010553) return true;

  if (fabs(pdgid)==555) return true;
  if (fabs(pdgid)==10555) return true;
  if (fabs(pdgid)==20555) return true;
  if (fabs(pdgid)==100555) return true;
  if (fabs(pdgid)==110555) return true;
  if (fabs(pdgid)==120555) return true;
  if (fabs(pdgid)==200555) return true;

  if (fabs(pdgid)==557) return true;
  if (fabs(pdgid)==10557) return true;

  //Bottom baryons
  
  if (fabs(pdgid)==5122) return true;
  if (fabs(pdgid)==5112) return true;
  if (fabs(pdgid)==5212) return true;
  if (fabs(pdgid)==5222) return true;
  if (fabs(pdgid)==5114) return true;
  if (fabs(pdgid)==5214) return true;
  if (fabs(pdgid)==5224) return true;
  if (fabs(pdgid)==5132) return true;
  if (fabs(pdgid)==5232) return true;
  

  if (fabs(pdgid)==5312) return true;
  if (fabs(pdgid)==5322) return true;
  if (fabs(pdgid)==5314) return true;
  if (fabs(pdgid)==5324) return true;
  if (fabs(pdgid)==5332) return true;
  if (fabs(pdgid)==5334) return true;
  
  if (fabs(pdgid)==5142) return true;
  if (fabs(pdgid)==5242) return true;
  if (fabs(pdgid)==5412) return true;
  if (fabs(pdgid)==5422) return true;
  if (fabs(pdgid)==5414) return true;
  if (fabs(pdgid)==5424) return true;
  if (fabs(pdgid)==5342) return true;
  if (fabs(pdgid)==5432) return true;
  if (fabs(pdgid)==5434) return true;
  if (fabs(pdgid)==5442) return true;
  if (fabs(pdgid)==5444) return true;

  if (fabs(pdgid)==5512) return true;
  if (fabs(pdgid)==5522) return true;
  if (fabs(pdgid)==5514) return true;
  if (fabs(pdgid)==5524) return true;
  if (fabs(pdgid)==5532) return true;
  if (fabs(pdgid)==5534) return true;
  if (fabs(pdgid)==5542) return true;
  if (fabs(pdgid)==5544) return true;
  if (fabs(pdgid)==5554) return true;
  
  return false;
;}

void SingleTopParticleLevelMCProducer::setIgnoredParticles(const std::vector<unsigned int> &particleIDs){
  ignoreParticleIDs = particleIDs;
  std::sort(ignoreParticleIDs.begin(), ignoreParticleIDs.end());
   }
 
void SingleTopParticleLevelMCProducer::setExcludeFromResonancePids(const std::vector<unsigned int> &particleIDs)
 {
      excludeFromResonancePids = particleIDs;
      std::sort( excludeFromResonancePids.begin(), excludeFromResonancePids.end());
       }

bool SingleTopParticleLevelMCProducer::isRadiation(int status){
  return false;
;}

bool SingleTopParticleLevelMCProducer::isFinalStateParticle(int pdgid){
  return false;
;}

void SingleTopParticleLevelMCProducer::matchReclusteredJets(const edm::Handle<std::vector<reco::GenJet> >& jets, const std::vector<fastjet::PseudoJet>& reclusteredJets, std::vector<int>& matchedIndices) {
  std::vector<bool> matchedLocks(reclusteredJets.size(),false);
  for(size_t j=0; j<jets->size(); ++j){
    //    std::cout << "jet j "<< j << " pt "<< jets->at(j).pt()<<std::endl;
    double matchedDR2 = 1e9;
    int matchedIdx = -1;
        if(fabs(jets->at(j).eta()) >5.5) continue;
    for(size_t rj=0; rj<reclusteredJets.size(); ++rj)
      {
	//	std::cout << "recl jet rj "<< rj << " pt "<< reclusteredJets.at(rj).pt()<<std::endl;
	if( matchedLocks.at(rj) ) continue; // skip jets that have already been matched
	
	double tempDR2 = reco::deltaR2( jets->at(j).eta(), jets->at(j).phi(), reclusteredJets.at(rj).eta(), reclusteredJets.at(rj).phi() );
	if( tempDR2 < matchedDR2 )
	  {
	    matchedDR2 = tempDR2;
	    matchedIdx = rj;
	    //	    std::cout << " they match!" <<std::endl;
	  }
      }
    
    if( matchedIdx>=0 )
      {
	if ( matchedDR2 > rParam_*rParam_ ) {
	  edm::LogError("JetMatchingFailed") << "Matched reclustered jet " << matchedIdx << " and original jet " << j <<" are separated by dR=" << sqrt(matchedDR2) << " which is greater than the jet size R=" << rParam_ << ".\n"
					     << "This is not expected so please check that the jet algorithm and jet size match those used for the original jet collection.";
std::cout<<"original jet "<<jets->at(j).rapidity()<<"  "<<jets->at(j).eta()<<std::endl;
std::cout<<"reclustered jet "<<reclusteredJets.at(matchedIdx).rapidity()<<"  "<<reclusteredJets.at(matchedIdx).eta()<<std::endl;
	}
	else
	  matchedLocks.at(matchedIdx) = true;
      }
    else
      edm::LogError("JetMatchingFailed") << "Matching reclustered to original jets failed. Please check that the jet algorithm and jetsize match those used for the original jet collection.";
    
    matchedIndices.push_back(matchedIdx);
  }
}


SingleTopParticleLevelMCProducer::~SingleTopParticleLevelMCProducer(){;}

DEFINE_FWK_MODULE( SingleTopParticleLevelMCProducer );

//  LocalWords:  MCtopsQuarkBar
