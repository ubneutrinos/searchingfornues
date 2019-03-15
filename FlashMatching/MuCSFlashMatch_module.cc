////////////////////////////////////////////////////////////////////////
// Class:       MuCSFlashMatch
// Plugin Type: analyzer (art v3_00_00)
// File:        MuCSFlashMatch_module.cc
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


class MuCSFlashMatch;


class MuCSFlashMatch : public art::EDAnalyzer {
public:
  explicit MuCSFlashMatch(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MuCSFlashMatch(MuCSFlashMatch const&) = delete;
  MuCSFlashMatch(MuCSFlashMatch&&) = delete;
  MuCSFlashMatch& operator=(MuCSFlashMatch const&) = delete;
  MuCSFlashMatch& operator=(MuCSFlashMatch&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  TTree* _tree;
  float _score;
  int _mucs;
  int _through;

  TTree* _evt_tree;

  int _run, _sub, _evt;

  float _best_score;
  float _mucs_score;
  std::vector<float> _score_v;
  int _score_mucs_idx;
  int _best_mucs;
  float _trk_start_x, _trk_start_y, _trk_start_z;
  float _trk_end_x, _trk_end_y, _trk_end_z;
  float _trk_start_x_mucs, _trk_start_y_mucs, _trk_start_z_mucs;
  float _trk_end_x_mucs, _trk_end_y_mucs, _trk_end_z_mucs;
  std::vector<float> _peSpectrum, _peHypothesis;
  std::vector<float> _peSpectrum_mucs, _peHypothesis_mucs;

  std::string fPFPproducer, fTrackproducer, fSpacePointproducer, fCTagproducer;
  
  bool fOnlyTagged;
  
  std::unique_ptr<flashmatch::FlashMatchingToolBase> _flashmatchTool;             ///< The slice id tool

  bool ThroughGoing(const float& xs, const float& ys, const float& zs, const float& xe, const float& ye, const float& ze);

};


MuCSFlashMatch::MuCSFlashMatch(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fPFPproducer         = p.get<std::string>("PFPproducer");
  fTrackproducer       = p.get<std::string>("Trackproducer");
  fSpacePointproducer  = p.get<std::string>("SpacePointproducer");
  fCTagproducer        = p.get<std::string>("CTagproducer");
  fOnlyTagged          = p.get<bool       >("OnlyTagged");

  const fhicl::ParameterSet& flashconfig = p.get<fhicl::ParameterSet>("SliceTool");  
  
  _flashmatchTool = art::make_tool<flashmatch::FlashMatchingToolBase>(flashconfig);

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("flashmatch","flashmatching tree");
  _tree->Branch("_run",&_run,"run/I");
  _tree->Branch("_sub",&_sub,"sub/I");
  _tree->Branch("_evt",&_evt,"evt/I");
  _tree->Branch("_score",&_score,"score/F");
  _tree->Branch("_mucs",&_mucs,"mucs/I");
  _tree->Branch("_through",&_through,"through/I");
  _tree->Branch("_trk_start_x",&_trk_start_x,"trk_start_x/F");
  _tree->Branch("_trk_start_y",&_trk_start_y,"trk_start_y/F");
  _tree->Branch("_trk_start_z",&_trk_start_z,"trk_start_z/F");
  _tree->Branch("_trk_end_x",&_trk_end_x,"trk_end_x/F");
  _tree->Branch("_trk_end_y",&_trk_end_y,"trk_end_y/F");
  _tree->Branch("_trk_end_z",&_trk_end_z,"trk_end_z/F");
  _tree->Branch("peSpectrum"  ,"std::vector<float>",&_peSpectrum);
  _tree->Branch("peHypothesis","std::vector<float>",&_peHypothesis);

  _evt_tree = tfs->make<TTree>("evtflashmatch","evt flashmatching tree");
  _evt_tree->Branch("_run",&_run,"run/I");
  _evt_tree->Branch("_sub",&_sub,"sub/I");
  _evt_tree->Branch("_evt",&_evt,"evt/I");
  _evt_tree->Branch("_best_score",&_best_score,"best_score/F");
  _evt_tree->Branch("_mucs_score",&_mucs_score,"mucs_score/F");
  _evt_tree->Branch("_best_mucs",&_best_mucs,"best_mucs/I");
  _evt_tree->Branch("_score_mucs_idx",&_score_mucs_idx,"score_mucs_idx/I");
  _evt_tree->Branch("_score_v","std::vector<float>",&_score_v);
  _evt_tree->Branch("_trk_start_x_mucs",&_trk_start_x_mucs,"trk_start_x_mucs/F");
  _evt_tree->Branch("_trk_start_y_mucs",&_trk_start_y_mucs,"trk_start_y_mucs/F");
  _evt_tree->Branch("_trk_start_z_mucs",&_trk_start_z_mucs,"trk_start_z_mucs/F");
  _evt_tree->Branch("_trk_end_x_mucs",&_trk_end_x_mucs,"trk_end_x_mucs/F");
  _evt_tree->Branch("_trk_end_y_mucs",&_trk_end_y_mucs,"trk_end_y_mucs/F");
  _evt_tree->Branch("_trk_end_z_mucs",&_trk_end_z_mucs,"trk_end_z_mucs/F");
  _evt_tree->Branch("peSpectrum_mucs"  ,"std::vector<float>",&_peSpectrum_mucs);
  _evt_tree->Branch("peHypothesis_mucs","std::vector<float>",&_peHypothesis_mucs);

}

void MuCSFlashMatch::analyze(art::Event const& e)
{

  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();

  // grab pfp in event
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);

  // grab tracks in event
  auto const& trk_h = e.getValidHandle<std::vector<recob::Track> >(fTrackproducer);
  
  // grab cosmic score association (to find MuCS tagged tracks)
  art::FindManyP<anab::CosmicTag> trk_ctag_assn_v(trk_h, e, fCTagproducer);

  // grab tracks associated to pfparticles
  art::FindManyP<recob::Track> pfp_trk_assn_v(pfp_h, e, fPFPproducer);
  
  // grab spacepoints associated with PFParticles
  art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fPFPproducer);
  
  auto const& spacepoint_h = e.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointproducer);
  
  art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, e, fSpacePointproducer);

  _best_score = 10000.;
  _mucs_score = -1; // this way we know if a MuCS track was tagged in the event

  size_t nbad = 0;

  // have we tagged an MUCS track in this event?
  bool mucstagged = false;
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
    const std::vector< art::Ptr<anab::CosmicTag> >& ctag_ptr_v = trk_ctag_assn_v.at(trk_key);
    
    _mucs = ctag_ptr_v.size();
    _through = 0;

    std::cout << " this track is mucs-associated!" << std::endl;

    if (_mucs) { // check that it is through-oging)
      if (ThroughGoing( _trk_start_x, _trk_start_y, _trk_start_z, _trk_end_x, _trk_end_y, _trk_end_z ) == false) {
	std::cout << "MUCS tagged track is not through-going" << std::endl;
	_through = 1;
	continue;
      }
    }// if a MuCS tagged track

    if (fOnlyTagged && (_mucs != 1) ) continue;
    //if (!fOnlyTagged && (_mucs != 0) ) continue;

    
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
    
    std::cout << "new score " << _score << std::endl;

    //if ( (_peSpectrum.size() < 232) || (_peHypothesis.size() < 32) ) continue;

    /*
    float PEhypo = 0;
    float PEreco = 0;
    for (size_t pmt=0; pmt < 32; pmt++){
      auto H = _peHypothesis[pmt];
      auto O = _peHypothesis[pmt+200];
      if ( (O < 10) && (H < 15.) ) 
	continue;
      PEhypo += H;
      PEreco += O;
    }
    float deltaPE = fabs(PEhypo-PEreco)/(PEhypo+PEreco);
    if (deltaPE > 0.5)
      continue;
    */
      

    _score_v.push_back(_score);

    if ( (_score < _best_score) && (_score > 0) ) { _best_score = _score; _best_mucs = _mucs; }
    if (_mucs) { 
      mucstagged = true;
      _score_mucs_idx = _score_v.size()-1;
      _peHypothesis_mucs = _peHypothesis;
      _peSpectrum_mucs = _peSpectrum;
      _mucs_score = _score; 
      _trk_end_x_mucs = _trk_end_x;
      _trk_start_x_mucs = _trk_start_x;
      _trk_end_y_mucs = _trk_end_y;
      _trk_start_y_mucs = _trk_start_y;
      _trk_end_z_mucs = _trk_end_z;
      _trk_start_z_mucs = _trk_start_z;
    }
    
    _tree->Fill();
    
  }// for all tracks

  if (mucstagged == true)
    _evt_tree->Fill();

  return;
}

void MuCSFlashMatch::beginJob()
{
  // Implementation of optional member function here.
}

void MuCSFlashMatch::endJob()
{
  // Implementation of optional member function here.
}

bool MuCSFlashMatch::ThroughGoing(const float& xs, const float& ys, const float& zs, const float& xe, const float& ye, const float& ze) {

  if (ys < 90) return false;
  if ( (xe < 250.) && (ye > -90) ) return false;
  return true;
}

DEFINE_ART_MODULE(MuCSFlashMatch)
