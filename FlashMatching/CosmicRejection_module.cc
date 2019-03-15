////////////////////////////////////////////////////////////////////////
// Class:       CosmicRejection
// Plugin Type: analyzer (art v2_11_03)
// File:        CosmicRejection_module.cc
//
// Generated at Sun Oct 28 14:38:04 2018 by Wouter Van de pontseele using cetskelgen
// from cetlib version v3_03_01.
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

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardataobj/AnalysisBase/T0.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "TTree.h"

// art TOOLS
#include "art/Utilities/ToolMacros.h"

class CosmicRejection;


class CosmicRejection : public art::EDProducer {
public:
  explicit CosmicRejection(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicRejection(CosmicRejection const &) = delete;
  CosmicRejection(CosmicRejection &&) = delete;
  CosmicRejection & operator = (CosmicRejection const &) = delete;
  CosmicRejection & operator = (CosmicRejection &&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  std::string fPFPproducer, fSpacePointproducer;

  std::unique_ptr<flashmatch::FlashMatchingToolBase> _flashmatchTool;             ///< The slice id tool

  void AddDaughters(const art::Ptr<recob::PFParticle>& pfp,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v);

  bool ObviousCosmic(const std::vector<art::Ptr<recob::Track> >& pfp_track_assn_v);

  bool TrackInTime(const art::Ptr<recob::Track>& pfp_track_assn_v);

  // Declare member data here.
  TTree* _tree;
  float _score;
  float _best_score;
  float _xe, _ye, _ze, _xs, _ys, _zs;
  int _obvious;
  std::vector<float> _peSpectrum, _peHypothesis;

  // PFP map
  std::map<unsigned int, unsigned int> _pfpmap;

};

DEFINE_ART_MODULE(CosmicRejection)

CosmicRejection::CosmicRejection(fhicl::ParameterSet const & p)
  :
  EDProducer(p)  // ,
 // More initializers here.
{

  produces< std::vector<anab::T0> >();
  produces< art::Assns<recob::PFParticle, anab::T0> >();

  fPFPproducer  = p.get<std::string>("PFPproducer");
  fSpacePointproducer  = p.get<std::string>("SpacePointproducer");

  const fhicl::ParameterSet& flashconfig = p.get<fhicl::ParameterSet>("SliceTool");  
  
  _flashmatchTool = art::make_tool<flashmatch::FlashMatchingToolBase>(flashconfig);
  //_flashmatchTool->configure(flashconfig);

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("flashmatch","flashmatching tree");
  _tree->Branch("_score",&_score,"score/F");
  _tree->Branch("_xe",&_xe,"xe/F");
  _tree->Branch("_ye",&_ye,"ye/F");
  _tree->Branch("_ze",&_ze,"ze/F");
  _tree->Branch("_xs",&_xs,"xs/F");
  _tree->Branch("_ys",&_ys,"ys/F");
  _tree->Branch("_zs",&_zs,"zs/F");
  _tree->Branch("_best_score",&_best_score,"best_score/F");
  _tree->Branch("_obvious",&_obvious,"obvious/I");
  _tree->Branch("peSpectrum"  ,"std::vector<float>",&_peSpectrum);
  _tree->Branch("peHypothesis","std::vector<float>",&_peHypothesis);

}

void CosmicRejection::produce(art::Event& e)
{

  std::unique_ptr< std::vector< anab::T0 > >                   cosmicTagTrackVector ( new std::vector<anab::T0> );
  std::unique_ptr< art::Assns<recob::PFParticle, anab::T0 > >  assnOutCosmicTagTrack( new art::Assns<recob::PFParticle, anab::T0> );

  // grab PFParticles in event
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);
  
  // grab spacepoints associated with PFParticles
  art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fPFPproducer);
  // grab tracks associated with PFParticles
  art::FindManyP<recob::Track> pfp_track_assn_v(pfp_h, e, fPFPproducer);
  
  // grab associated metadata
  art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfp_h, e, fPFPproducer);    
  
  auto const& spacepoint_h = e.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointproducer);
  
  art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, e, fSpacePointproducer);
  
  /*
  // ADDITION FROM PETRILLO
  e.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointproducer);
  
  // grab the hits associated to the PFParticles
  auto pfp_hit_assn_v = lar::FindManyInChainP<recob::Hit, recob::SpacePoint>::find(pfp_h, e, fPFPproducer);
  */
  
  std::cout << "There are " << pfp_h->size() << " pfparticles in the event " << std::endl;

  _best_score = 1e5;
  _obvious = 0;
  
  // fill map: pfparticle Self() -> index/key
  _pfpmap.clear();
  for (unsigned int p=0; p < pfp_h->size(); p++)
    _pfpmap[pfp_h->at(p).Self()] = p;

  for (unsigned int p=0; p < pfp_h->size(); p++){
    
    auto const& pfp = pfp_h->at(p);
    
    // start from primary PFParticles
    if (pfp.IsPrimary() == false) continue;
    
    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
    
    // now build vectors of PFParticles, space-points, and hits for this slice
    std::vector<recob::PFParticle> pfp_v;
    std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
    std::vector<std::vector<art::Ptr<recob::SpacePoint>>> spacepoint_v_v;
    std::vector<std::vector<art::Ptr<recob::Hit>>> hit_v_v;
    
    std::cout << "creating PFP hierarchy." << std::endl;
    AddDaughters(pfp_ptr, pfp_h, pfp_ptr_v);
    std::cout << "There are " << pfp_ptr_v.size() << " PFParticles in this hierarchy " << std::endl << std::endl;
    
    // go through these pfparticles and fill info needed for matching
    for (size_t i=0; i < pfp_ptr_v.size(); i++) {
      
      auto key = pfp_ptr_v.at(i).key();
      recob::PFParticle pfp = *pfp_ptr_v.at(i);
      
      pfp_v.push_back(pfp);
      
      //auto const& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
      const std::vector< art::Ptr<recob::SpacePoint> >& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
      
      std::vector< art::Ptr<recob::Hit> > hit_ptr_v;
      
      for (size_t sp=0; sp < spacepoint_ptr_v.size(); sp++) {
	auto const& spkey = spacepoint_ptr_v.at(sp).key();
	const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = spacepoint_hit_assn_v.at( spkey );
	for (size_t h=0; h < this_hit_ptr_v.size(); h++) {
	  hit_ptr_v.push_back( this_hit_ptr_v.at( h ) );
	}// for all hits associated to this spacepoint
      }// fpr all spacepoints
      
      spacepoint_v_v.push_back( spacepoint_ptr_v );
      hit_v_v.push_back( hit_ptr_v );
      
    }// for all pfp pointers

    // ignore events where the primary is not a muon [we are looking for obvious cosmics!]
    if (pfp_track_assn_v.at(p).size() != 1) continue;

    auto const trk = pfp_track_assn_v.at(p).at(0);

    // tracks too short should not be considered...
    if (trk->Length() < 10) continue;

    // do not attempt flash-matching with tracks not in time with the beam drift window
    if ( TrackInTime(trk) == false ) continue;

    // ready to call flash-matching
    _score = _flashmatchTool->ClassifySlice(e, pfp_ptr_v, spacepoint_v_v, hit_v_v, _peSpectrum, _peHypothesis);
    if (_score <= 0) continue;

    if (_score < _best_score) { 
      _best_score = _score; 
      // is this hierarchy an obvious cosmic?
      _obvious = 0;
      if ( ObviousCosmic(pfp_track_assn_v.at(p)) ) _obvious = 1;
    }// if the best score yet
    
  }// for all PFParticles

  _tree->Fill();

  e.put(std::move(cosmicTagTrackVector) );
  e.put(std::move(assnOutCosmicTagTrack));

  return;
}

bool CosmicRejection::TrackInTime(const art::Ptr<recob::Track>& trk) {

  auto vtx = trk->Vertex();
  auto end = trk->End();

  if (vtx.X() < -10 || vtx.X() > 260.) return false;
  if (end.X() < -10 || end.X() > 260.) return false;

  return true;
}// is the track in time with the beam?

bool CosmicRejection::ObviousCosmic(const std::vector<art::Ptr<recob::Track> >& pfp_track_assn_v) {

  // if no track is associated to the primary PFP, this is not an obvious cosmic
  if (pfp_track_assn_v.size() == 0) return false;

  if (pfp_track_assn_v.size() != 1) {
    std::cout << "ERROR. More then one track associated to this PFP!!" << std::endl;
    return false;
  }

  auto const trk = pfp_track_assn_v.at(0);

  auto vtx = trk->Vertex();
  auto end = trk->End();

  _xe = end.X();
  _ye = end.Y();
  _ze = end.Z();

  _xs = vtx.X();
  _ys = vtx.Y();
  _zs = vtx.Z();

  if ( (vtx.X() < 10 || vtx.X() > 245 || vtx.Y() < -100. || vtx.Y() > 100. || vtx.Z() < 10 || vtx.Z() > 1000.) and
       (end.X() < 10 || end.X() > 245 || end.Y() < -100. || end.Y() > 100. || end.Z() < 10 || end.Z() > 1000.) )
    return true;

  return false;
}

void CosmicRejection::AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v) {
  
  auto daughters = pfp_ptr->Daughters();
  
  pfp_v.push_back(pfp_ptr);
  
  std::cout << "\t PFP w/ PdgCode " << pfp_ptr->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;
  
  for(auto const& daughterid : daughters) {

    if (_pfpmap.find(daughterid) == _pfpmap.end()) {
      std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
      continue;
    }
    
    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, _pfpmap.at(daughterid) );
    
    AddDaughters(pfp_ptr, pfp_h, pfp_v);
    
  }// for all daughters
  
  return;
}
