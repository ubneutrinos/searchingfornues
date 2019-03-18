#ifndef BACTRACKINGFUNCS_H
#define BACTRACKINGFUNCS_H

// services for detector properties
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "lardata/Utilities/FindManyInChainP.h"

namespace searchingfornues {

// shift coordinates for truth particles according to SCE offsets + time offsets
void ApplyDetectorOffsets(const float _vtx_t, const float _vtx_x, const float _vtx_y, const float _vtx_z,
                          float& _xtimeoffset, float& _xsceoffset, float& _ysceoffset, float& _zsceoffset) {

  auto const& detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
  auto const& detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  double g4Ticks = detClocks->TPCG4Time2Tick(_vtx_t) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
  _xtimeoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);

  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto offset = SCE->GetPosOffsets(geo::Point_t(_vtx_x,_vtx_y,_vtx_z));
  _xsceoffset = offset.X();
  _ysceoffset = offset.Y();
  _zsceoffset = offset.Z();

}

// BackTrack a single hit collection (i.e. from a PFParticle)
art::Ptr<simb::MCParticle> getAssocMCParticle(art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& hittruth,
					      const std::vector<art::Ptr<recob::Hit> >& hits,
					      float& purity, float& completeness) {

  // store total charge from hits
  float pfpcharge = 0; // total hit charge from clusters
  float maxcharge = 0; // charge backtracked to best match
  
  //credit: Wes Ketchum
  std::unordered_map<int,double> trkide;
  std::unordered_map<int,float> trkq;
  double maxe=-1, tote=0;
  art::Ptr<simb::MCParticle> maxp_me; //pointer for the particle match we will calculate
  //simb::MCParticle* maxp_me; //pointer for the particle match we will calculate
  for (auto h : hits) {
    pfpcharge += h->Integral();
    //const std::vector<const simb::MCParticle*> particle_vec = hittruth.at(h.key());
    std::vector<art::Ptr<simb::MCParticle> > particle_vec = hittruth.at(h.key());
    //auto particle_vec = hittruth.at(h.key());
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec = hittruth.data(h.key());;
    //loop over particles
    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
      trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id
      trkq  [ particle_vec[i_p]->TrackId() ] += h->Integral()  * match_vec[i_p]->ideFraction; //store hit integral associated to this hit
      tote += match_vec[i_p]->energy; //calculate total energy deposited
      if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
	maxe = trkide[ particle_vec[i_p]->TrackId() ];
	maxp_me = particle_vec[i_p];
	maxcharge = trkq[ particle_vec[i_p]->TrackId() ];
      }
    }//end loop over particles per hit
  }

  purity       = maxcharge / pfpcharge;
  completeness = 0;

  return maxp_me;
}

}

#endif
