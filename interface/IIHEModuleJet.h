#ifndef UserCode_IIHETree_IIHEModuleJet_h
#define UserCode_IIHETree_IIHEModuleJet_h

#include "UserCode/IIHETree/interface/IIHEModule.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

class IIHEJetVariableBase{
public:
  IIHEJetVariableBase(std::string, std::string, int) ;
  ~IIHEJetVariableBase(){} ;

  const std::string       Name(){ return       name_ ; }
  const std::string BranchName(){ return branchName_ ; }
  const int         branchType(){ return branchType_ ; }

  virtual void reset(){} ;
  virtual void store(IIHEAnalysis*){} ;
  bool addBranch(IIHEAnalysis*) ;

private:
  std::string       name_ ;
  std::string branchName_ ;
  int         branchType_ ;
};

class IIHEJetVariableInt: IIHEJetVariableBase{
public:
  IIHEJetVariableInt(std::string, std::string) ;
  ~IIHEJetVariableInt(){} ;
  void reset(){ value_ = -999 ; }
  void fill(int value){ value_ = value ; }
  void store(IIHEAnalysis*) ;
private:
  int value_ ;
};

class IIHEJetVariableFloat: IIHEJetVariableBase{
public:
  IIHEJetVariableFloat(std::string, std::string) ;
  ~IIHEJetVariableFloat(){} ;
  void reset(){ value_ = -999.0 ; }
  void fill(float value){ value_ = value ; }
  void store(IIHEAnalysis*) ;
private:
  float value_ ;
};

class IIHEJetWrapper{
public:
  explicit IIHEJetWrapper(std::string);
  ~IIHEJetWrapper(){};

  void addBranches(IIHEAnalysis*) ;
  void reset() ;
  void fill(pat::Jet) ;
  void store(IIHEAnalysis*) ;

  private:
  int type_ ;
  std::string prefix_ ;
  IIHEJetVariableFloat*   pt_        ;
  IIHEJetVariableFloat*   energy_        ;
  std::vector<IIHEJetVariableBase*> variables_ ;
};
// class decleration

class IIHEModuleJet : public IIHEModule {
public:
  explicit IIHEModuleJet(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC);
  explicit IIHEModuleJet(const edm::ParameterSet& iConfig): IIHEModule(iConfig){};
  ~IIHEModuleJet();
  
  void   pubBeginJob(){   beginJob() ; } ;
  void pubBeginEvent(){ beginEvent() ; } ;
  void   pubEndEvent(){   endEvent() ; } ;
  virtual void pubAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ analyze(iEvent, iSetup) ; } ;
  
  virtual void beginEvent() ;
  virtual void endEvent() ;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);

private:
  edm::InputTag pfJetLabel_;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetToken_;
  edm::InputTag pfJetLabelSmeared_;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetTokenSmeared_;
  edm::InputTag pfJetLabelEnUp_;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetTokenEnUp_;
  edm::InputTag pfJetLabelEnDown_;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetTokenEnDown_;
  edm::InputTag pfJetLabelSmearedJetResUp_;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetTokenSmearedJetResUp_;
  edm::InputTag pfJetLabelSmearedJetResDown_;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetTokenSmearedJetResDown_;

  IIHEJetWrapper*  jetEnUpWrapper_;
  IIHEJetWrapper*  jetEnDownWrapper_;
  IIHEJetWrapper*  jetSmearedWrapper_;
  IIHEJetWrapper*  jetSmearedJetResUpWrapper_;
  IIHEJetWrapper*  jetSmearedJetResDownWrapper_;


  float ETThreshold_ ;
  bool isMC_;
};
#endif
