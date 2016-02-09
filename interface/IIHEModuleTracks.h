#ifndef UserCode_IIHETree_IIHEModuleTracks_h
#define UserCode_IIHETree_IIHEModuleTracks_h

#include "UserCode/IIHETree/interface/IIHEModule.h"

// class decleration
class IIHEModuleTracks : public IIHEModule {
public:
  explicit IIHEModuleTracks(const edm::ParameterSet& iConfig);
  ~IIHEModuleTracks();
  
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
};
#endif
