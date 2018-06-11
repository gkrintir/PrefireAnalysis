// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/transform.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "TTree.h"

#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"
// == reco::LorentzVector
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> LorentzVector;

namespace {
  struct EventStruct {
    Long64_t run;
    Long64_t lumi;
    Long64_t event;
    int bunchCrossing;
    
    std::vector<LorentzVector> jet_p4;
    std::vector<float> jet_neutralEmFrac;
    std::vector<float> jet_neutralHadFrac;

    std::vector<int> L1EG_bx;
    std::vector<LorentzVector> L1EG_p4;
    std::vector<int> L1EG_iso;

    std::vector<int> L1GtBx;
  };
  struct JRA{

    static const int MAXJETS = 1000;
    static const int MAXTRACKS = 5000;
    static const int MAXHLTBITS = 5000;
    static const int MAXBFRAG = 500;

    int nref;

    float rawpt[MAXJETS];
    float jtpt[MAXJETS];
    float jteta[MAXJETS];
    float jtphi[MAXJETS];
    float jty[MAXJETS];
    float jtpu[MAXJETS];
    float jtm[MAXJETS];
    float jtarea[MAXJETS];

    float trackMax[MAXJETS];
    float trackSum[MAXJETS];
    int trackN[MAXJETS];

    float chargedMax[MAXJETS];
    float chargedSum[MAXJETS];
    int chargedN[MAXJETS];

    float photonMax[MAXJETS];
    float photonSum[MAXJETS];
    int photonN[MAXJETS];

    float trackHardSum[MAXJETS];
    float chargedHardSum[MAXJETS];
    float photonHardSum[MAXJETS];

    int trackHardN[MAXJETS];
    int chargedHardN[MAXJETS];
    int photonHardN[MAXJETS];

    float neutralMax[MAXJETS];
    float neutralSum[MAXJETS];
    int neutralN[MAXJETS];

    float eMax[MAXJETS];
    float eSum[MAXJETS];
    int eN[MAXJETS];

    float muMax[MAXJETS];
    float muSum[MAXJETS];
    int muN[MAXJETS];


  };
}

class PrefiringJetAna : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit PrefiringJetAna(const edm::ParameterSet&);
    ~PrefiringJetAna();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    // virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    // virtual void endJob() override;

    PFJetIDSelectionFunctor looseJetIdSelector_{PFJetIDSelectionFunctor::WINTER16, PFJetIDSelectionFunctor::LOOSE};
  //pat::strbitset hasLooseId_;

   edm::EDGetTokenT<reco::JetView> jetToken_;
   edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidateLabel_;
   edm::EDGetTokenT<reco::TrackCollection> trackTag_;
    edm::EDGetTokenT<BXVector<l1t::EGamma>> l1egToken_;
    edm::EDGetTokenT<BXVector<GlobalAlgBlk>> l1GtToken_;

    TTree * tree_;
    EventStruct event_;
    JRA jets_;
};

PrefiringJetAna::PrefiringJetAna(const edm::ParameterSet& iConfig):
  jetToken_(consumes<reco::JetView>(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  //tagJetCut_(iConfig.getParameter<std::string>("tagJetCut")),
  l1egToken_(consumes<BXVector<l1t::EGamma>>(iConfig.getParameter<edm::InputTag>("l1egSrc"))),
  l1GtToken_(consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("l1GtSrc")))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  pfCandidateLabel_ = consumes<reco::PFCandidateCollection> (iConfig.getUntrackedParameter<edm::InputTag>("pfCandidateLabel",edm::InputTag("particleFlow")));
  trackTag_ = consumes<reco::TrackCollection> (iConfig.getParameter<edm::InputTag>("trackTag"));
  //hasLooseId_ = looseJetIdSelector_.getBitTemplate();

  tree_ = fs->make<TTree>("tree","Event Summary");
  tree_->Branch("run", &event_.run);
  tree_->Branch("lumi", &event_.lumi);
  tree_->Branch("event", &event_.event);
  tree_->Branch("bunchCrossing", &event_.bunchCrossing);
  tree_->Branch("jet_p4", &event_.jet_p4);
  //tree_->Branch("jet_neutralEmFrac", &event_.jet_neutralEmFrac);
  //tree_->Branch("jet_neutralHadFrac", &event_.jet_neutralHadFrac);
  tree_->Branch("nref",&jets_.nref,"nref/I");
  tree_->Branch("rawpt",jets_.rawpt,"rawpt[nref]/F");
  tree_->Branch("jteta",jets_.jteta,"jteta[nref]/F");
  tree_->Branch("jty",jets_.jty,"jty[nref]/F");
  tree_->Branch("jtphi",jets_.jtphi,"jtphi[nref]/F");
  tree_->Branch("jtpu",jets_.jtpu,"jtpu[nref]/F");
  tree_->Branch("jtm",jets_.jtm,"jtm[nref]/F");
  tree_->Branch("jtarea",jets_.jtarea,"jtarea[nref]/F");
  tree_->Branch("trackMax", jets_.trackMax,"trackMax[nref]/F");
  tree_->Branch("trackSum", jets_.trackSum,"trackSum[nref]/F");
  tree_->Branch("trackN", jets_.trackN,"trackN[nref]/I");
  tree_->Branch("trackHardSum", jets_.trackHardSum,"trackHardSum[nref]/F");
  tree_->Branch("trackHardN", jets_.trackHardN,"trackHardN[nref]/I");
  
  tree_->Branch("chargedMax", jets_.chargedMax,"chargedMax[nref]/F");
  tree_->Branch("chargedSum", jets_.chargedSum,"chargedSum[nref]/F");
  tree_->Branch("chargedN", jets_.chargedN,"chargedN[nref]/I");
  tree_->Branch("chargedHardSum", jets_.chargedHardSum,"chargedHardSum[nref]/F");
  tree_->Branch("chargedHardN", jets_.chargedHardN,"chargedHardN[nref]/I");
  
  tree_->Branch("photonMax", jets_.photonMax,"photonMax[nref]/F");
  tree_->Branch("photonSum", jets_.photonSum,"photonSum[nref]/F");
  tree_->Branch("photonN", jets_.photonN,"photonN[nref]/I");
  tree_->Branch("photonHardSum", jets_.photonHardSum,"photonHardSum[nref]/F");
  tree_->Branch("photonHardN", jets_.photonHardN,"photonHardN[nref]/I");
  
  tree_->Branch("neutralMax", jets_.neutralMax,"neutralMax[nref]/F");
  tree_->Branch("neutralSum", jets_.neutralSum,"neutralSum[nref]/F");
  tree_->Branch("neutralN", jets_.neutralN,"neutralN[nref]/I");
  
  tree_->Branch("eMax", jets_.eMax,"eMax[nref]/F");
  tree_->Branch("eSum", jets_.eSum,"eSum[nref]/F");
  tree_->Branch("eN", jets_.eN,"eN[nref]/I");
  
  tree_->Branch("muMax", jets_.muMax,"muMax[nref]/F");
  tree_->Branch("muSum", jets_.muSum,"muSum[nref]/F");
  tree_->Branch("muN", jets_.muN,"muN[nref]/I");
    
  tree_->Branch("L1EG_bx", &event_.L1EG_bx);
  tree_->Branch("L1EG_p4", &event_.L1EG_p4);
  tree_->Branch("L1EG_iso", &event_.L1EG_iso);
  tree_->Branch("L1GtBx", &event_.L1GtBx);
}


PrefiringJetAna::~PrefiringJetAna()
{
}

void
PrefiringJetAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  event_.run = iEvent.run();
  event_.lumi = iEvent.luminosityBlock();
  event_.event = iEvent.id().event();
  event_.bunchCrossing = iEvent.bunchCrossing();

  Handle<reco::JetView> jetHandle;
  iEvent.getByToken(jetToken_, jetHandle);

  Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByToken(pfCandidateLabel_,pfCandidates);

  Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackTag_,tracks);

  Handle<BXVector<l1t::EGamma>> l1egHandle;
  iEvent.getByToken(l1egToken_, l1egHandle);

  Handle<BXVector<GlobalAlgBlk>> l1GtHandle;
  iEvent.getByToken(l1GtToken_, l1GtHandle);


  event_.jet_p4.clear();
  event_.jet_neutralEmFrac.clear();
  event_.jet_neutralHadFrac.clear();

  jets_.nref = 0;

  for(unsigned int j = 0; j < jetHandle->size(); ++j){
   const reco::Jet& jet = (*jetHandle)[j];
    jets_.muMax[jets_.nref] = 0;
    jets_.muSum[jets_.nref] = 0;
    jets_.muN[jets_.nref] = 0;
    
    jets_.eMax[jets_.nref] = 0;
    jets_.eSum[jets_.nref] = 0;
    jets_.eN[jets_.nref] = 0;
    
    jets_.neutralMax[jets_.nref] = 0;
    jets_.neutralSum[jets_.nref] = 0;
    jets_.neutralN[jets_.nref] = 0;
    
    jets_.photonMax[jets_.nref] = 0;
    jets_.photonSum[jets_.nref] = 0;
    jets_.photonN[jets_.nref] = 0;
    jets_.photonHardSum[jets_.nref] = 0;
    jets_.photonHardN[jets_.nref] = 0;
    
    jets_.chargedMax[jets_.nref] = 0;
    jets_.chargedSum[jets_.nref] = 0;
    jets_.chargedN[jets_.nref] = 0;
    jets_.chargedHardSum[jets_.nref] = 0;
    jets_.chargedHardN[jets_.nref] = 0;
    
    jets_.trackMax[jets_.nref] = 0;
    jets_.trackSum[jets_.nref] = 0;
    jets_.trackN[jets_.nref] = 0;
    jets_.trackHardSum[jets_.nref] = 0;
    jets_.trackHardN[jets_.nref] = 0;
    
    
    jets_.rawpt[jets_.nref]=jet.pt();
    jets_.jteta[jets_.nref] = jet.eta();
    jets_.jtphi[jets_.nref] = jet.phi();
    jets_.jty[jets_.nref] = jet.eta();
    jets_.jtpu[jets_.nref] = jet.pileup();
    jets_.jtm[jets_.nref] = jet.mass();
    jets_.jtarea[jets_.nref] = jet.jetArea();
    event_.jet_p4.push_back( jet.p4() );
    //event_.jet_neutralEmFrac.push_back( jetNew.neutralEmEnergyFraction() );
    //event_.jet_neutralHadFrac.push_back( jetNew.neutralHadronEnergyFraction() );


    for(unsigned int icand = 0; icand < tracks->size(); ++icand){
      const reco::Track& track = (*tracks)[icand];
      
	double dr = deltaR(jet,track);
	if(dr < 0.4){
	  double ptcand = track.pt();
	  jets_.trackSum[jets_.nref] += ptcand;
	  jets_.trackN[jets_.nref] += 1;
	  if(ptcand > 5.0){
	    jets_.trackHardSum[jets_.nref] += ptcand;
	    jets_.trackHardN[jets_.nref] += 1;
	    
	  }
	  if(ptcand > jets_.trackMax[jets_.nref]) jets_.trackMax[jets_.nref] = ptcand;
	}

	for(unsigned int icand = 0; icand < pfCandidates->size(); ++icand){
	  const reco::PFCandidate& track = (*pfCandidates)[icand];
        double dr = deltaR(jet,track);
        if(dr < 0.4){
	  double ptcand = track.pt();
	  int pfid = track.particleId();

	  switch(pfid){
	  case 1:
	    jets_.chargedSum[jets_.nref] += ptcand;
	    jets_.chargedN[jets_.nref] += 1;
	    if(ptcand > 5.0){
	      jets_.chargedHardSum[jets_.nref] += ptcand;
	      jets_.chargedHardN[jets_.nref] += 1;
	    }
	    if(ptcand > jets_.chargedMax[jets_.nref]) jets_.chargedMax[jets_.nref] = ptcand;
	    break;

	  case 2:
	    jets_.eSum[jets_.nref] += ptcand;
	    jets_.eN[jets_.nref] += 1;
	    if(ptcand > jets_.eMax[jets_.nref]) jets_.eMax[jets_.nref] = ptcand;
	    break;

	  case 3:
	    jets_.muSum[jets_.nref] += ptcand;
	    jets_.muN[jets_.nref] += 1;
	    if(ptcand > jets_.muMax[jets_.nref]) jets_.muMax[jets_.nref] = ptcand;
	    break;

	  case 4:
	    jets_.photonSum[jets_.nref] += ptcand;
	    jets_.photonN[jets_.nref] += 1;
	    if(ptcand > 5.0){
	      jets_.photonHardSum[jets_.nref] += ptcand;
	      jets_.photonHardN[jets_.nref] += 1;
	    }
	    if(ptcand > jets_.photonMax[jets_.nref]) jets_.photonMax[jets_.nref] = ptcand;
	    break;

	  case 5:
	    jets_.neutralSum[jets_.nref] += ptcand;
	    jets_.neutralN[jets_.nref] += 1;
	    if(ptcand > jets_.neutralMax[jets_.nref]) jets_.neutralMax[jets_.nref] = ptcand;
	    break;

	  default:
	    break;

	  }
	}
	}}
    jets_.nref++;
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

  event_.L1GtBx.clear();
  event_.L1GtBx.push_back(l1GtHandle->begin(-2)->getFinalOR());
  event_.L1GtBx.push_back(l1GtHandle->begin(-1)->getFinalOR());
  event_.L1GtBx.push_back(l1GtHandle->begin( 0)->getFinalOR());
  event_.L1GtBx.push_back(l1GtHandle->begin( 1)->getFinalOR());
  event_.L1GtBx.push_back(l1GtHandle->begin( 2)->getFinalOR());

  tree_->Fill();
}


void
PrefiringJetAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(PrefiringJetAna);
