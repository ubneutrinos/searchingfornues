#ifndef BACTRACKINGFUNCS_H
#define BACTRACKINGFUNCS_H

// services for detector properties
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"

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


 // struct storing the info from MCTrack and MCShower to perform backtracking
 struct BtPart {
  public:
   BtPart(const int pdg_, const float px_, const float py_, const float pz_, const float e_, const std::vector<unsigned int>& tids_)
   : pdg(pdg_), px(px_), py(py_), pz(pz_), e(e_), tids(tids_) {}
   BtPart(const int pdg_, const float px_, const float py_, const float pz_, const float e_, const unsigned int tid_)
   : pdg(pdg_), px(px_), py(py_), pz(pz_), e(e_) { tids.push_back(tid_); }
   int pdg;
   float px,py,pz,e;
   std::vector<unsigned int> tids;
   int nhits=0;
 };

 std::vector<BtPart> initBacktrackingParticleVec(const std::vector<sim::MCShower>& inputMCShower,
						 const std::vector<sim::MCTrack>& inputMCTrack,
						 const std::vector<recob::Hit>& inputHits,
						 const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& assocMCPart) {
    //
    // Notes for reference: in this way we have the full history of electrons and photons
    // however, we do not attach delta-rays nor michel electrons to muon (which is probably good),
    // but also we do not attach the scatter products to pi+ nor neutrons.
    // Note that the only non-primary particles stored are photons (also dalitz electrons!) from pi0 (which are not stored).
    std::vector<BtPart> btparts_v;
    for (auto mcs : inputMCShower) {
      if (mcs.Process()=="primary" || (mcs.MotherPdgCode()==111 && mcs.Process()=="Decay" && mcs.MotherProcess()=="primary")) {
	btparts_v.push_back(BtPart(mcs.PdgCode(), mcs.Start().Momentum().Px()*0.001, mcs.Start().Momentum().Py()*0.001,
				   mcs.Start().Momentum().Pz()*0.001, mcs.Start().Momentum().E()*0.001, mcs.DaughterTrackID()));
      }
    }
    for (auto mct : inputMCTrack) {
      if (mct.Process()=="primary") {
	btparts_v.push_back(BtPart(mct.PdgCode(), mct.Start().Momentum().Px()*0.001, mct.Start().Momentum().Py()*0.001,
				   mct.Start().Momentum().Pz()*0.001, mct.Start().Momentum().E()*0.001, mct.TrackID()));
      }
    }
    // Now let's fill the nhits member using all input hits
    for (unsigned int ih=0; ih<inputHits.size(); ih++) {
      auto assmcp = assocMCPart->at(ih);
      auto assmdt = assocMCPart->data(ih);
      for (unsigned int ia=0; ia<assmcp.size(); ++ia) {
	auto mcp = assmcp[ia];
	auto amd = assmdt[ia];
	if (amd->isMaxIDE!=1) continue;
	for (auto& btp : btparts_v) {
	  if (std::find(btp.tids.begin(),btp.tids.end(),mcp->TrackId())!=btp.tids.end()) {
	    btp.nhits++;
	  }
	}
      }
    }
    return btparts_v;
 }

 // BackTrack a single hit collection (i.e. from a PFParticle)
 // backtrack based on #hits, not on energy... could be done but it's a bit more complicated
  int getAssocBtPart(const std::vector<art::Ptr<recob::Hit> >& hits,
		     const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& assocMCPart,
		     const std::vector<BtPart>& btpartsv,
		     float& purity, float& completeness) {
   //
   std::vector<unsigned int> bthitsv(btpartsv.size(),0);
   //
   for (unsigned int ih=0; ih<hits.size(); ih++) {
     art::Ptr<recob::Hit> hitp = hits[ih];
     auto assmcp = assocMCPart->at(hitp.key());
     auto assmdt = assocMCPart->data(hitp.key());
     for (unsigned int ia=0; ia<assmcp.size(); ++ia) {
       auto mcp = assmcp[ia];
       auto amd = assmdt[ia];
       if (amd->isMaxIDE!=1) continue;
       for (unsigned int ib=0; ib<btpartsv.size(); ++ib) {
	 auto& btp = btpartsv[ib];
	 if (std::find(btp.tids.begin(),btp.tids.end(),mcp->TrackId())!=btp.tids.end()) {
	   bthitsv[ib]++;
	 }
       }
     }
   }
   purity = 0.;
   completeness = 0.;
   unsigned int maxel = (std::max_element(bthitsv.begin(),bthitsv.end()) - bthitsv.begin());
   if (maxel == bthitsv.size()) return -1;
   if (bthitsv[maxel]==0) return -1;
   //
   purity       = float(bthitsv[maxel])/float(hits.size());
   completeness = float(bthitsv[maxel])/float(btpartsv[maxel].nhits);
   //
   return maxel;
 }


}

#endif
