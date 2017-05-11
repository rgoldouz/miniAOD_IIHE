#ifndef IIHE_Particle_Level_Producer_h
#define IIHE_Particle_Level_Producer_h


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include <FWCore/Framework/interface/Run.h>
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"
#include "PhysicsTools/PatAlgos/interface/PATUserDataHelper.h"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"

//class JetFlavourIdentifier;

  class IIHEParticleLevelMCProducer : public edm::EDProducer {

    public:
    typedef std::vector<bool>     ParticleBitmap;
    typedef std::vector<const reco::GenParticle*> ParticleVector;
        
    explicit IIHEParticleLevelMCProducer(const edm::ParameterSet & iConfig);
    ~IIHEParticleLevelMCProducer();
    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    enum ResonanceState {
           kNo = 0,
           kDirect =1,
	   kIndirect = 2
    };

    bool isFinalStateParticle(int status);
    bool isRadiation(int status);

    bool isPromptLepton( const reco::GenParticle * particle, const ParticleVector &allParticles);
    bool isParton(int pdgid);

    bool isIgnored(int pdgid);
    bool isExcludedFromResonance(int pdgid);

    bool isHadron(int pdgid);
    bool isFromHadron(const reco::GenParticle * particle, const ParticleVector &allParticles);
   
    bool isResonance(int pdgid);
    ResonanceState isFromResonance( const reco::GenParticle * particle, const ParticleVector &allParticles, ParticleBitmap &map);
    
    bool isBHadron(int pdgid);
    bool isBHadronRescaled( const reco::GenParticle * particle, const ParticleVector &allParticles);
    bool isBHadronRescaledDecayProduct( const reco::GenParticle * particle, const ParticleVector &allParticles, const ParticleVector &rescaledBHadrons);
    //const reco::GenParticle * rescaleParticle(const reco::GenParticle*);
    //       static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

    void setExcludeFromResonancePids(const std::vector<unsigned int> &particleIDs);
    void setIgnoredParticles(const std::vector<unsigned int> &particleIDs);
    void matchReclusteredJets(const edm::Handle<std::vector<reco::GenJet> >& jets, const std::vector<fastjet::PseudoJet>& reclusteredJets, std::vector<int>& matchedIndices) ;

  private:

    //InputTags
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesSrc_;
    edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsSrc_; 
    std::string hadronizer_, channel_,jetDefinition_;
    std::vector<unsigned int> ignoreParticleIDs;
    std::vector<unsigned int> excludeFromResonancePids;    
    double minJetPt_,rParam_,rescaling_;
    bool doDressedCollections_,doRescaledBHadronJetSelection_,doBaseCollections_,doBDescendentJetSelection_, useJetsNoNu_;

  };
#endif
