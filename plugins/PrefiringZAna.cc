// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/transform.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "HLTrigger/HLTcore/interface/HLTConfigData.h"

#include "TTree.h"
#include "TMath.h"

#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"
#include "TRegexp.h"

// == reco::LorentzVector
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> LorentzVector;

namespace {
  struct EventStruct {
    Long64_t run;
    Long64_t lumi;
    Long64_t event;
    int bunchCrossing;
    
    LorentzVector tag_electron;
    std::vector<int> L1EG_bx;
    std::vector<LorentzVector> L1EG_p4;
    std::vector<int> L1EG_iso;
  };

}

//class PrefiringZAna : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
class PrefiringZAna : public edm::EDAnalyzer  {

  public:
    explicit PrefiringZAna(const edm::ParameterSet&);
    ~PrefiringZAna();
  
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
  virtual void beginJob() ;//override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);//override;
  virtual void endJob() ;//override;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<BXVector<l1t::EGamma>> l1egToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken_;
  edm::EDGetTokenT<edm::View<reco::GsfElectron> > gsfElectronsCollection_;

  std::string processName_;
  bool loadTriggersFromHLT_;
  int verbose_;
  std::vector<std::string>   triggerNames_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventTag_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  const edm::EDGetTokenT<trigger::TriggerEvent> triggerEventToken_;

  edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
   int evtCounter;
  std::map<std::string, bool> triggerInMenu;

  HLTConfigProvider hltConfig_;

  unsigned int triggerIndex_;
  unsigned int moduleIndex_;
  string moduleLabel_;
  vector<string> moduleLabels_;

  TTree * tree_;
  EventStruct event_;
  std::vector<double> id;
  std::vector<double> pt;
  std::vector<double> eta;
  std::vector<double> phi;
  std::vector<double> mass;
};

PrefiringZAna::PrefiringZAna(const edm::ParameterSet& iConfig):
  l1egToken_(consumes<BXVector<l1t::EGamma>>(iConfig.getParameter<edm::InputTag>("l1egSrc"))),
  processName_(iConfig.getParameter<std::string>("processName")),
  loadTriggersFromHLT_(iConfig.getUntrackedParameter<bool>("loadTriggersFromHLT",false)),
  triggerNames_(iConfig.getParameter<std::vector<std::string> >("triggerNames")),
  triggerResultsTag_(iConfig.getParameter<edm::InputTag>("triggerResults")),
  triggerEventTag_(iConfig.getParameter<edm::InputTag>("triggerEvent")),
  triggerResultsToken_(consumes<edm::TriggerResults>(triggerResultsTag_)),
  triggerEventToken_(consumes<trigger::TriggerEvent>(triggerEventTag_))
{

  gsfElectronsCollection_ = consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("gsfElectronLabel"));
  edm::Service<TFileService> fs;
  verbose_ = 0;
  tree_ = fs->make<TTree>("tree","Event Summary");
  tree_->Branch("run", &event_.run);
  tree_->Branch("lumi", &event_.lumi);
  tree_->Branch("event", &event_.event);
  tree_->Branch("bunchCrossing", &event_.bunchCrossing);
  tree_->Branch("tag_electron", &event_.tag_electron);
  tree_->Branch("L1EG_bx", &event_.L1EG_bx);
  tree_->Branch("L1EG_p4", &event_.L1EG_p4);
  tree_->Branch("L1EG_iso", &event_.L1EG_iso);
}


PrefiringZAna::~PrefiringZAna()
{
}

void
PrefiringZAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  event_.run = iEvent.run();
  event_.lumi = iEvent.luminosityBlock();
  event_.event = iEvent.id().event();
  event_.bunchCrossing = iEvent.bunchCrossing();

  Handle<BXVector<l1t::EGamma>> l1egHandle;
  iEvent.getByToken(l1egToken_, l1egHandle);

  iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
  iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);

  Handle<edm::View<reco::GsfElectron> > gsfElectronsHandle;
  iEvent.getByToken(gsfElectronsCollection_, gsfElectronsHandle);

  for(unsigned int itrig=0; itrig<triggerNames_.size(); itrig++){

    std::map<std::string,bool>::iterator inMenu = triggerInMenu.find(triggerNames_[itrig]);
    if (inMenu==triggerInMenu.end()){ continue; }
    triggerIndex_ = hltConfig_.triggerIndex(triggerNames_[itrig]);
    const unsigned int mIndex = triggerResultsHandle_->index(triggerIndex_);
    
    // Results from TriggerEvent product - Attention: must look only for
    // modules actually run in this path for this event!
    for (unsigned int j=0; j<=mIndex; ++j) {
      // check whether the module is packed up in TriggerEvent product
      string trigFilterIndex = hltConfig_.moduleLabels(triggerIndex_).at(j); //this is simple to put into a loop to get all triggers...
        const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(trigFilterIndex,"",processName_)));
	if (filterIndex<triggerEventHandle_->sizeFilters()) {
	  const trigger::Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
	  const trigger::Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
	  const unsigned int nI(VIDS.size());
	  const unsigned int nK(KEYS.size());
	  assert(nI==nK);
	  const unsigned int n(max(nI,nK));
	  
	  const trigger::TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());

	  for (unsigned int i=0; i!=n; ++i) {
	    const trigger::TriggerObject& TO(TOC[KEYS[i]]);
	    //This check prevents grabbing the L1 trigger object (VIDS < 0), and finds the max trigger pt within all trigger collections
	    if(VIDS[i]>0){ // && pt<TO.pt()){
	      if(verbose_){
		std::cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] 
			  << TO.id() << " " << TO.pt() << " " << TO.et() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass() 
			  << std::endl;
	      }
	      id.push_back(TO.id());
	      pt.push_back(TO.pt());
	      eta.push_back(TO.eta());
	      phi.push_back(TO.phi());
	      mass.push_back(TO.mass());
	    }
	  }

	}
    }
  }

 for (edm::View<reco::GsfElectron>::const_iterator ele = gsfElectronsHandle->begin(); ele != gsfElectronsHandle->end(); ++ele) {
   double dphi;
   for (unsigned loop=0; loop < id.size(); loop++) {
     dphi = ele->phi()-phi.at(loop);
   if (fabs(dphi)>TMath::Pi())
     { dphi = dphi < 0? (TMath::TwoPi()) + dphi : dphi - TMath::TwoPi(); }
   double deltaR = sqrt(pow((ele->eta()-eta.at(loop)),2) + pow(dphi,2));
      if (ele->pt()>20 && deltaR<0.1){
     //std::cout << loop<<" match " << ele->pt()<< " "<< pt.at(loop)<< std::endl;
       event_.tag_electron = ele->p4();
   }
   }
 }


  event_.L1EG_bx.clear();
  event_.L1EG_p4.clear();
  event_.L1EG_iso.clear();

  auto readBx = [&] (const BXVector<l1t::EGamma>& egVect, int bx) {
    for (auto itL1=l1egHandle->begin(bx); itL1!=l1egHandle->end(bx); ++itL1) {
      event_.L1EG_bx.push_back(bx);
      event_.L1EG_p4.push_back(itL1->p4());
      event_.L1EG_iso.push_back(itL1->hwIso());
    }
  };

  readBx(*l1egHandle, -2);
  readBx(*l1egHandle, -1);
  readBx(*l1egHandle,  0);
  readBx(*l1egHandle,  1);
  readBx(*l1egHandle,  2);

  tree_->Fill();
}
// ------------ method called once each job just before starting event loop  ------------
void
PrefiringZAna::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
PrefiringZAna::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
PrefiringZAna::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{


  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      std::vector<std::string> activeHLTPathsInThisEvent = hltConfig_.triggerNames();
     
      triggerInMenu.clear(); 

      if(loadTriggersFromHLT_){
	      triggerNames_ = hltConfig_.triggerNames(); 
	      //if(!evtCounter) nt_.reserve(triggerNames_.size());
	      for(unsigned int isize=0; isize<triggerNames_.size(); isize++){
		      std::size_t versionLoc = triggerNames_.at(isize).find_last_of("_v");
		      std::string genericTrg = triggerNames_.at(isize).substr(0,versionLoc+1);
		      if(versionLoc < triggerNames_.at(isize).size()) triggerNames_.at(isize) = genericTrg;
 		      //if(!evtCounter) nt_[isize] = fs->make<TTree>(triggerNames_.at(isize).c_str(),Form("trigger %d",isize));
	      }
	      evtCounter++;
      }
      for(unsigned int itrig=0; itrig<triggerNames_.size(); itrig++){
	      for (std::vector<std::string>::const_iterator iHLT = activeHLTPathsInThisEvent.begin(); iHLT != activeHLTPathsInThisEvent.end(); ++iHLT){
		//matching with regexp filter name. More than 1 matching filter is allowed so trig versioning is transparent to analyzer
		if (TString(*iHLT).Contains(TRegexp(TString(triggerNames_[itrig])))){
		  triggerInMenu[*iHLT] = true;
		  triggerNames_[itrig] = TString(*iHLT);
		}
	      }
      }
      for(unsigned int itrig=0; itrig<triggerNames_.size(); itrig++){
	      std::map<std::string,bool>::iterator inMenu = triggerInMenu.find(triggerNames_[itrig]);
	      if (inMenu==triggerInMenu.end()) {
		      cout << "<HLT Object Analyzer> Warning! Trigger " << triggerNames_[itrig] << " not found in HLTMenu. Skipping..." << endl;
	      }
      }
      if(verbose_){
	hltConfig_.dump("ProcessName");
	hltConfig_.dump("GlobalTag");
	hltConfig_.dump("TableName");
	hltConfig_.dump("Streams");
	hltConfig_.dump("Datasets");
	hltConfig_.dump("PrescaleTable");
	hltConfig_.dump("ProcessPSet");
      }
    }
  } else {
    cout << "HLTObjectAnalyzer::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }

}
// ------------ method called when ending the processing of a run  ------------
void
PrefiringZAna::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
PrefiringZAna::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
PrefiringZAna::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

void
PrefiringZAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(PrefiringZAna);
