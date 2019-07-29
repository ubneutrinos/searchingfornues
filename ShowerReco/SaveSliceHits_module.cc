////////////////////////////////////////////////////////////////////////
// Class:       SaveSliceHits
// Plugin Type: producer (art v3_01_02)
// File:        SaveSliceHits_module.cc
//
// Generated at Fri Jun 21 13:31:12 2019 by David Caratelli using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

class SaveSliceHits;


class SaveSliceHits : public art::EDProducer {
public:
  explicit SaveSliceHits(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SaveSliceHits(SaveSliceHits const&) = delete;
  SaveSliceHits(SaveSliceHits&&) = delete;
  SaveSliceHits& operator=(SaveSliceHits const&) = delete;
  SaveSliceHits& operator=(SaveSliceHits&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;

  art::InputTag fHitproducer, fClusterproducer, fPfpproducer, fSliceproducer;
  float fMinHitCharge; // ADC, Hit->Integral()

  void addDaughter(const recob::PFParticle pfp,
		   art::ValidHandle<std::vector<recob::PFParticle> > pfp_h,
		   art::FindManyP<recob::Cluster> pfp_clus_assn_v,
		   art::FindManyP<recob::Hit> clus_hit_assn_v,
		   std::vector<unsigned int> &PfpHitIdx_v);
  
};


SaveSliceHits::SaveSliceHits(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{

  fHitproducer     = p.get< art::InputTag > ("Hitproducer"    );
  fClusterproducer = p.get< art::InputTag > ("Clusterproducer");
  fPfpproducer     = p.get< art::InputTag > ("Pfpproducer"    );
  fSliceproducer   = p.get< art::InputTag > ("Sliceproducer"  );
  fMinHitCharge    = p.get< float         >("MinHitCharge", 0.); // ADC for Hit.Integral()

  produces<std::vector<recob::Hit> >();
}

void SaveSliceHits::produce(art::Event& e)
{

  std::cout << "****************** NEW EVNT ****************************8" << std::endl;


  // produce output Hits
  std::unique_ptr< std::vector<recob::Hit> > Hit_v(new std::vector<recob::Hit> );  
  
  // get event hits
  auto const& hit_h = e.getValidHandle<std::vector<recob::Hit> >(fHitproducer);
  // get a handle to the pfparticles
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPfpproducer);
  // grab slice associated to slices
  art::FindManyP<recob::Slice> pfp_slice_assn_v(pfp_h, e, fPfpproducer);
  // grab slices themselves
  auto const& slice_h = e.getValidHandle<std::vector<recob::Slice> >(fSliceproducer);
  // grab hits associated to slices
  art::FindManyP<recob::Hit> slice_hit_assn_v(slice_h, e, fSliceproducer);
  // grab clusters associated with PFParticles
  art::FindManyP<recob::Cluster> pfp_clus_assn_v(pfp_h, e, fClusterproducer);
  // grab clusters themselves
  auto const& cluster_h = e.getValidHandle<std::vector<recob::Cluster> >(fClusterproducer);
  // get hits associated to clusters
  art::FindManyP<recob::Hit> clus_hit_assn_v(cluster_h, e, fClusterproducer);
  
  // vector of hit indices in the slice for hits belonging to a PFParticle
  std::vector<unsigned int> PfpHitIdx_v;
  // vector of hit indices for the slice
  std::vector<unsigned int> SliceHitIdx_v;
  

  _pfpmap.clear();

  size_t neutrinos = 0;

  
  // buid pfp self -> index map
  std::map<size_t, size_t> PFP_self_key_map;
  for (unsigned int p=0; p < pfp_h->size(); p++) {
    auto pfpself = pfp_h->at(p).Self();
    _pfpmap[pfpself] = p;
  }

  size_t slicekey; // key for neutrino slice
  
  // loop through PFParticles
  for (size_t p=0; p < pfp_h->size(); p++) {
    
    auto pfp = pfp_h->at(p);
    
    //  find neutrino candidate
    if (pfp.IsPrimary() == false) continue; // only want primaries
    
    auto PDG = fabs(pfp.PdgCode());
    if ( (PDG != 12) && (PDG != 14) ) continue; // only want the neutrino

    neutrinos += 1;
    
    // grab slice associated to neutrino PFParticle
    auto pfp_slice_ass = pfp_slice_assn_v.at( p );
    for (size_t si=0; si < pfp_slice_ass.size(); si++) {
      // grab slice index
      slicekey = pfp_slice_ass[si].key();
    }// for all slices associated to PFP
    
    addDaughter(pfp, pfp_h, pfp_clus_assn_v, clus_hit_assn_v, PfpHitIdx_v);
    
  }// for all PFParticles
  
  if (neutrinos != 1) {
    e.put(std::move(Hit_v));  
    return;
  }

  // grab slice -> hit ass vector
  auto slice_hit_ass = slice_hit_assn_v.at(slicekey);
  // loop through slice hits and add to display
  for (size_t slicehitidx = 0; slicehitidx < slice_hit_ass.size(); slicehitidx++) {
    
    //auto hit = *(slice_hit_ass.at(slicehitidx));
    SliceHitIdx_v.push_back( slice_hit_ass.at(slicehitidx).key() );
    
  }// for all slice-hits

  
  // loop through hits and save those that are not associated to a PFParticle
  for (unsigned int hitidx : SliceHitIdx_v) {
    bool matchedtopfp = false;
    for (auto const& pfphitidx : PfpHitIdx_v) {
      if (pfphitidx == hitidx) {
	matchedtopfp = true;
	break;
      }// if matched to PFP hit index
    }// for all PFP hit indices
    if (matchedtopfp == false) {
      if (hit_h->at(hitidx).Integral() > fMinHitCharge)
	Hit_v->emplace_back(hit_h->at(hitidx));
    }// if not matched
  }// for all slice hits

  std::cout << "adding " << Hit_v->size() << " hits to event" << std::endl;
  std::cout << std::endl;

  e.put(std::move(Hit_v));  
  
  return;
}

  void SaveSliceHits::addDaughter(const recob::PFParticle pfp,
				  art::ValidHandle<std::vector<recob::PFParticle> > pfp_h,
				  art::FindManyP<recob::Cluster> pfp_clus_assn_v,
				  art::FindManyP<recob::Hit> clus_hit_assn_v,
				  std::vector<unsigned int> &PfpHitIdx_v) {
  
  auto daughters = pfp.Daughters();
  
  for(auto const& daughterid : daughters) {
    
    if (_pfpmap.find(daughterid) == _pfpmap.end()) {
      std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
      continue;
    }
    
    const auto daughter = pfp_h->at( _pfpmap[daughterid] );

    size_t pfphitnum = 0;
    
    // find cluster associated to this PFP
    auto const& ass_clus_v = pfp_clus_assn_v.at( _pfpmap[daughterid] );
    // for each cluster, get associated hits
    for (size_t c=0; c < ass_clus_v.size(); c++) {
      auto const& clus_key = ass_clus_v[c].key();
      auto const& ass_hits = clus_hit_assn_v.at(clus_key);
      pfphitnum += ass_hits.size();
      for (size_t h=0; h < ass_hits.size(); h++) {
	//std::cout << " \t pfp hit has key " << ass_hits[h].key() << std::endl;
	PfpHitIdx_v.push_back( ass_hits[h].key() );
      }// for all hits associated to cluster
    }// for all clusters associated to PFP
    
    std::cout << "PFP has " << pfphitnum << " hits associated" << std::endl;

    addDaughter(daughter, pfp_h, pfp_clus_assn_v, clus_hit_assn_v, PfpHitIdx_v);
    
  }// for all daughters

  return;
}


void SaveSliceHits::beginJob()
{
  // Implementation of optional member function here.
}

void SaveSliceHits::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(SaveSliceHits)
