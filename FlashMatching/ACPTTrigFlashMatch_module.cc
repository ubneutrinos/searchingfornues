////////////////////////////////////////////////////////////////////////
// Class:       ACPTTrigFlashMatch
// Plugin Type: analyzer (art v3_00_00)
// File:        ACPTTrigFlashMatch_module.cc
//
// Generated at Sun Jan 13 18:33:43 2019 by David Caratelli using cetskelgen
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Utilities/make_tool.h"
#include "FlashMatchingToolBase_tool.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "TTree.h"

// art TOOLS
#include "art/Utilities/ToolMacros.h"


class ACPTTrigFlashMatch;


class ACPTTrigFlashMatch : public art::EDAnalyzer {
public:
  explicit ACPTTrigFlashMatch(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ACPTTrigFlashMatch(ACPTTrigFlashMatch const&) = delete;
  ACPTTrigFlashMatch(ACPTTrigFlashMatch&&) = delete;
  ACPTTrigFlashMatch& operator=(ACPTTrigFlashMatch const&) = delete;
  ACPTTrigFlashMatch& operator=(ACPTTrigFlashMatch&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  TTree* _tree;

  int _run, _sub, _evt;

  float _score;
  int _acpt;
  TTree* _evt_tree;
  float _best_score;
  float _acpt_score;
  std::vector<float> _score_v;
  int _score_acpt_idx;
  int _best_acpt;
  float _trk_start_x, _trk_start_y, _trk_start_z;
  float _trk_end_x, _trk_end_y, _trk_end_z;
  float _trk_start_x_acpt, _trk_start_y_acpt, _trk_start_z_acpt;
  float _trk_end_x_acpt, _trk_end_y_acpt, _trk_end_z_acpt;
  std::vector<float> _peSpectrum, _peHypothesis;
  std::vector<float> _peSpectrum_acpt, _peHypothesis_acpt;

  std::string fPFPproducer, fTrackproducer, fSpacePointproducer, fT0producer;
  
  bool fOnlyTagged;
  
  std::unique_ptr<flashmatch::FlashMatchingToolBase> _flashmatchTool;             ///< The slice id tool

};


ACPTTrigFlashMatch::ACPTTrigFlashMatch(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fPFPproducer         = p.get<std::string>("PFPproducer");
  fTrackproducer       = p.get<std::string>("Trackproducer");
  fSpacePointproducer  = p.get<std::string>("SpacePointproducer");
  fT0producer          = p.get<std::string>("T0producer");
  fOnlyTagged          = p.get<bool       >("OnlyTagged");

  const fhicl::ParameterSet& flashconfig = p.get<fhicl::ParameterSet>("SliceTool");  
  
  _flashmatchTool = art::make_tool<flashmatch::FlashMatchingToolBase>(flashconfig);

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("flashmatch","flashmatching tree");
  _tree->Branch("_run",&_run,"run/I");
  _tree->Branch("_sub",&_sub,"sub/I");
  _tree->Branch("_evt",&_evt,"evt/I");
  _tree->Branch("_score",&_score,"score/F");
  _tree->Branch("_acpt",&_acpt,"acpt/I");
  _tree->Branch("_trk_start_x",&_trk_start_x,"trk_start_x/F");
  _tree->Branch("_trk_start_y",&_trk_start_y,"trk_start_y/F");
  _tree->Branch("_trk_start_z",&_trk_start_z,"trk_start_z/F");
  _tree->Branch("_trk_end_x",&_trk_end_x,"trk_end_x/F");
  _tree->Branch("_trk_end_y",&_trk_end_y,"trk_end_y/F");
  _tree->Branch("_trk_end_z",&_trk_end_z,"trk_end_z/F");
  _tree->Branch("peSpectrum"  ,"std::vector<float>",&_peSpectrum);
  _tree->Branch("peHypothesis","std::vector<float>",&_peHypothesis);

  _evt_tree = tfs->make<TTree>("evtflashmatch","evt flashmatching tree");
  _evt_tree->Branch("_best_score",&_best_score,"best_score/F");
  _evt_tree->Branch("_acpt_score",&_acpt_score,"acpt_score/F");
  _evt_tree->Branch("_best_acpt",&_best_acpt,"best_acpt/I");
  _evt_tree->Branch("_score_acpt_idx",&_score_acpt_idx,"score_acpt_idx/I");
  _evt_tree->Branch("_score_v","std::vector<float>",&_score_v);
  _evt_tree->Branch("_trk_start_x_acpt",&_trk_start_x_acpt,"trk_start_x_acpt/F");
  _evt_tree->Branch("_trk_start_y_acpt",&_trk_start_y_acpt,"trk_start_y_acpt/F");
  _evt_tree->Branch("_trk_start_z_acpt",&_trk_start_z_acpt,"trk_start_z_acpt/F");
  _evt_tree->Branch("_trk_end_x_acpt",&_trk_end_x_acpt,"trk_end_x_acpt/F");
  _evt_tree->Branch("_trk_end_y_acpt",&_trk_end_y_acpt,"trk_end_y_acpt/F");
  _evt_tree->Branch("_trk_end_z_acpt",&_trk_end_z_acpt,"trk_end_z_acpt/F");
  //_evt_tree->Branch("peSpectrum_acpt"  ,"std::vector<float>",&_peSpectrum_acpt);
  //_evt_tree->Branch("peHypothesis_acpt","std::vector<float>",&_peHypothesis_acpt);

}

void ACPTTrigFlashMatch::analyze(art::Event const& e)
{

  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();

  // grab pfp in event
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);

  // grab tracks in event
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track> >(fTrackproducer);
  
  // grab T0 association (to find ACPT tagged tracks)
  art::FindManyP<anab::T0> trk_t0_assn_v(trk_h, e, fT0producer);

  // grab tracks associated to pfparticles
  art::FindManyP<recob::Track> pfp_trk_assn_v(pfp_h, e, fPFPproducer);
  
  // grab spacepoints associated with PFParticles
  art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fPFPproducer);
  
  auto const& spacepoint_h = e.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointproducer);
  
  art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, e, fSpacePointproducer);

  _best_score = 10000.;
  _acpt_score = -1; // this way we know if a MuCS track was tagged in the event

  size_t nbad = 0;

  // have we tagged an ACPT track in this event?
  bool acpttagged = false;

  _score_v.clear();

  // loop through PFParticles
  for (unsigned int p=0; p < pfp_h->size(); p++){

    //std::cout << "3D NEW TRACK" << std::endl;
    _peSpectrum.clear();
    _peHypothesis.clear();
    
    auto const& pfp = pfp_h->at(p);
    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p); 
   auto pfp_key = pfp_ptr.key();

    // only primary PFParticles
    if (pfp.IsPrimary() == false) continue;
    
    // grab the associated track
    const std::vector< art::Ptr<recob::Track> >& trk_ptr_v = pfp_trk_assn_v.at(pfp_key);

    // are there assocaited tracks? if no skip
    if (trk_ptr_v.size() != 1) continue;

    // grab the track key
    auto trk_key = trk_ptr_v.at(0).key();

    const art::Ptr<recob::Track> trk_ptr(trk_h, trk_key);

    _trk_start_x = trk_ptr->Vertex().X();
    _trk_start_y = trk_ptr->Vertex().Y();
    _trk_start_z = trk_ptr->Vertex().Z();
    _trk_end_x = trk_ptr->End().X();
    _trk_end_y = trk_ptr->End().Y();
    _trk_end_z = trk_ptr->End().Z();


    if (trk_ptr->Length() < 20.) { nbad += 1; continue; }

    if ( (_trk_start_x < -10) || (_trk_start_x > 270) ) { nbad += 1; continue; }

    if ( (_trk_end_x < -10)   || (_trk_end_x > 270)   ) { nbad += 1; continue; }

    // associations
    const std::vector< art::Ptr<recob::SpacePoint> >& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(pfp_key);
    const std::vector< art::Ptr<anab::T0> >& t0_ptr_v = trk_t0_assn_v.at(trk_key);
    
    _acpt = t0_ptr_v.size();

    //std::cout << "\t\t tagged?" << std::endl;

    if (fOnlyTagged && (_acpt != 1) ) continue;
    //if (!fOnlyTagged && (_acpt != 0) ) continue;

    //std::cout << "\t\t tagged!!!" << std::endl;

    
    std::vector< art::Ptr<recob::Hit> > hit_ptr_v;
    
    for (size_t sp=0; sp < spacepoint_ptr_v.size(); sp++) {
      auto const& spkey = spacepoint_ptr_v.at(sp).key();
      //std::cout << "3D XYZ " << spacepoint_ptr_v.at(sp)->XYZ()[0] << " " << spacepoint_ptr_v.at(sp)->XYZ()[1] << " " << spacepoint_ptr_v.at(sp)->XYZ()[2] << std::endl;
      const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = spacepoint_hit_assn_v.at( spkey );
      for (size_t h=0; h < this_hit_ptr_v.size(); h++) {
	hit_ptr_v.push_back( this_hit_ptr_v.at( h ) );
      }// for all hits associated to this spacepoint
    }// fpr all spacepoints
    
    
    _score = _flashmatchTool->ClassifyTrack(e, spacepoint_ptr_v, hit_ptr_v, _peSpectrum, _peHypothesis);

    _score_v.push_back(_score);


    if ( (_score < _best_score) && (_score > 0) ) { _best_score = _score; _best_acpt = _acpt; }

    if (_acpt) { 
      acpttagged = true;
      _score_acpt_idx = _score_v.size()-1;
      _acpt_score = _score; 
      _trk_end_x_acpt = _trk_end_x;
      _trk_start_x_acpt = _trk_start_x;
      _trk_end_y_acpt = _trk_end_y;
      _trk_start_y_acpt = _trk_start_y;
      _trk_end_z_acpt = _trk_end_z;
      _trk_start_z_acpt = _trk_start_z;
    }
    
    _tree->Fill();
    
  }// for all tracks

  if (acpttagged == true)
    _evt_tree->Fill();

  return;
}

void ACPTTrigFlashMatch::beginJob()
{
  // Implementation of optional member function here.
}

void ACPTTrigFlashMatch::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ACPTTrigFlashMatch)
