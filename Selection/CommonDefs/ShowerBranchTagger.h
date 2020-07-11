#ifndef SHOWERBRANCHTAGGER_H
#define SHOWERBRANCHTAGGER_H

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TMatrixDSymEigen.h" 

namespace searchingfornues
{

  /**
   * @brief go from 3D coordinates to 2D coordinates in order to compare 3D reco with hits
   * @input pt3d -> 3d point to be projected
   * @input pl -> which plane are we on?
   * @input wire2cm -> wire 2 cm conversion
   * @input time2cm -> time 2 cm conversion
   * @output wirecm -> pt3d wire coordinate in cm
   * @output timecm -> pt3d time coordinate in cm
   */
  void Project3Dto2D(const TVector3& pt3d, const int& pl,
		     const float& wire2cm, const float& time2cm,
		     float& wirecm, float& timecm) {

    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    
    wirecm = geom->WireCoordinate(pt3d[1],pt3d[2],geo::PlaneID(0,0,pl)) * wire2cm;
    timecm = pt3d[0];

    return;
  }

  /**
   * @brief get hit wire/time in cm
   * @input recob::hit 
   * @output hitwire -> hit wire coordinate in cm
   * @output hittime -> hit time coordinate in cm
   */
  void GetHitWireTime(const art::Ptr<recob::Hit> &hit, 
		      const float& wire2cm, const float& time2cm,		      
		      float& hitwire, float& hittime) {

    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();

    hitwire = hit->WireID().Wire * wire2cm;
    hittime = (hit->PeakTime() - detp->TriggerOffset())  * time2cm;

    return;
  }

  /**
   * @brief given a 3D pt and a 2D hit get their distance in 2D on the plane
   * @input pt3d -> 3d point to be projected
   * @input hit -> hit
   * @input wire2cm -> wire 2 cm conversion
   * @input time2cm -> time 2 cm conversion
   * @return 2d distance [cm]
   */
  float HitPtDistance(const TVector3& pt3d, const art::Ptr<recob::Hit> &hit,
		      const float& wire2cm, const float& time2cm) {

    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    // what plane are we on?
    auto pl = hit->WireID().Plane;

    float ptwire, pttime;
    Project3Dto2D(pt3d,pl,wire2cm,time2cm,ptwire,pttime);
    
    float hitwire = hit->WireID().Wire * wire2cm;
    float hittime = (hit->PeakTime() - detp->TriggerOffset())  * time2cm;
    
    float distance = sqrt( (ptwire-hitwire)*(ptwire-hitwire) + (pttime-hittime)*(pttime-hittime) );

    return distance;
  }
  
  /**
   * @brief given a 2D cluster and the vertex, how well aligned are they?
   * @input pt3d -> 3d point to be projected
   * @input hits -> hit vector for cluster
   * @input wire2cm -> wire 2 cm conversion
   * @input time2cm -> time 2 cm conversion
   * @return dot product
   */
  float ClusterVtxAlignment(const TVector3& pt3d, const std::vector<art::Ptr<recob::Hit>> &hits,
		      const float& wire2cm, const float& time2cm) {

    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    if (hits.size() == 0) return 0;

    // what plane are we on?
    auto pl = hits[0]->WireID().Plane;

    // find closest hit
    size_t closesthitidx = 0;
    float dmin = 1e6;
    for (size_t i=0; i < hits.size(); i++) {
      auto vtxdistance = HitPtDistance(pt3d,hits[i],wire2cm,time2cm);
      if (vtxdistance < dmin) { dmin = vtxdistance; closesthitidx = i; }
    }// for all hits

    // now calculate two, 2D vectors:
    // (1) pt3d -> closest hit
    // (2) charge-weighted 2D vector from cluster points to closest point

    // (1)
    float pt2dW, pt2dT; // 2D coordinates of 3D point on plane
    Project3Dto2D(pt3d,pl,wire2cm,time2cm,pt2dW,pt2dT);
    
    float Shitwire = (hits[closesthitidx])->WireID().Wire * wire2cm;
    float Shittime = ((hits[closesthitidx])->PeakTime() - detp->TriggerOffset())  * time2cm;

    float vtx2start_w = (Shitwire-pt2dW);
    float vtx2start_t = (Shittime-pt2dT);

    // (2)
    
    float charge = 0;
    float dwire = 0;
    float dtime = 0;
    for (size_t i=0; i < hits.size(); i++) {
      if (i == closesthitidx) continue;
      float hitwire = (hits[i])->WireID().Wire * wire2cm;
      float hittime = ((hits[i])->PeakTime() - detp->TriggerOffset())  * time2cm;
      
      dwire += (hitwire-Shitwire) * (hits[i])->Integral();
      dtime += (hittime-Shittime) * (hits[i])->Integral();
      charge += (hits[i])->Integral();

    }// for all hits
    dwire /= charge;
    dtime /= charge;

    if (charge==0 || ( (dwire*dwire) + (dtime*dtime) )==0 || ( (vtx2start_w*vtx2start_w) + (vtx2start_t*vtx2start_t) )==0) return std::numeric_limits<float>::max();

    // calculate dot product
    float dot = (dwire * vtx2start_w) + (dtime * vtx2start_t);
    dot /= sqrt( (dwire*dwire) + (dtime*dtime) );
    dot /= sqrt( (vtx2start_w*vtx2start_w) + (vtx2start_t*vtx2start_t) );
    
    return dot;
  }


  /**
   * @brief 2D direction of charge-weighted cluster hits w.r.t. vertex 
   * @input pt3d -> 3d point for vertex
   * @input hits -> hit vector for cluster
   * @input wire2cm -> wire 2 cm conversion
   * @input time2cm -> time 2 cm conversion
   * @return direction [angle in radians w.r.t. vertical]
   */
  float ClusterVtxDirection(const TVector3& pt3d, const std::vector<art::Ptr<recob::Hit>> &hits,
		      const float& wire2cm, const float& time2cm) {

    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    if (hits.size() == 0) return 0;
    
    // what plane are we on?
    auto pl = hits[0]->WireID().Plane;
    
    // (1) project vertex to 2D
    float pt2dW, pt2dT; // 2D coordinates of 3D point on plane
    Project3Dto2D(pt3d,pl,wire2cm,time2cm,pt2dW,pt2dT);
    
    // (2) charge-weighted direction vector from projected vtx
    
    float charge = 0;
    float dwire = 0;
    float dtime = 0;
    for (size_t i=0; i < hits.size(); i++) {
      float hitwire = (hits[i])->WireID().Wire * wire2cm;
      float hittime = ((hits[i])->PeakTime() - detp->TriggerOffset())  * time2cm;
      
      dwire += (hitwire-pt2dW) * (hits[i])->Integral();
      dtime += (hittime-pt2dT) * (hits[i])->Integral();
      charge += (hits[i])->Integral();

    }// for all hits
    dwire /= charge;
    dtime /= charge;

    // calculate the direction w.r.t. vertical
    float angle = (180./3.1415) * atan(fabs(dtime)/fabs(dwire));
    if ( (dtime >= 0) && (dwire >= 0) ) 
      angle = 90. - angle;
    if ( (dtime < 0) && (dwire >= 0) ) 
      angle += 90;
    if ( (dtime < 0) && (dwire < 0) ) 
      angle = 180. + (90. - angle);
    if ( (dtime >= 0) && (dwire < 0) ) 
      angle = 270. + angle;
    
    return angle;
  }


  /**
   * @brief get hits from a given handle (gaushit) for a different collection of hits
   * @input hits -> hits associated to cluster, for example
   * @input gaushit_h -> hit handle from gaushit (for example)
   * @return gaushit hits identical to the cluster hits
   */
  std::vector<art::Ptr<recob::Hit>> getGaussHits(const std::vector<art::Ptr<recob::Hit>> &hits,
						 const art::ValidHandle<std::vector<recob::Hit> > gaushit_h) {
    
    std::vector<art::Ptr<recob::Hit> > gaushit_v;
    
    for (size_t h1=0; h1 < hits.size(); h1++) {
      
      auto protonhit = hits.at(h1);
      
      for (size_t h2=0; h2 < gaushit_h->size(); h2++) {
	
	auto gaushit = gaushit_h->at(h2);
	
	// if idntical, add to output vector
	if ( (fabs(protonhit->PeakTime() - gaushit.PeakTime()) < 0.001) &&
	     (fabs(protonhit->WireID().Wire - gaushit.WireID().Wire) == 0) )
	  
	  gaushit_v.push_back( art::Ptr<recob::Hit>(gaushit_h,h2) );
	
      }// for gaushits
    }// for proton cluster's hits
    
    return gaushit_v;
  }
  
  /**
   * @brief get the 2D dot product between a cluster and the shower
   * @input gammaWire -> average charge-weighted wire of cluster
   * @input gammaTime -> average charge-weighted tick of cluster
   * @input pl        -> plane on which cluster was found
   * @input showerVtx -> 3D vertex of shower
   * @input showerDir -> 3D direction of shower
   * @input wire2cm   -> conversion from wire to cm
   * @output dot -> dot product between cluster and shower direction
   * @output d2d -> 2D distance between average cluster position and shower vertex
   */
  void GammaDot(const float& gammaWire, const float& gammaTime,const int& pl,
		const TVector3& showerVtx, const TVector3& showerDir,
		const float wire2cm,
		float &dot, float& d2d) {
    
    float Vtxwire, Vtxtime;
    Project3Dto2D(showerVtx,pl,wire2cm,0.,Vtxwire,Vtxtime);

    float Dirwire, Dirtime;
    Project3Dto2D(showerVtx,pl,wire2cm,0.,Dirwire,Dirtime);
    
    TVector3 showerDir2D(Dirwire,Dirtime,0.);
    TVector3 gammaDir2D(gammaWire-Vtxwire,gammaTime-Vtxtime,0.);
    
    d2d = sqrt( ((gammaWire - Vtxwire) * (gammaWire - Vtxwire)) +
		((gammaTime - Vtxtime) * (gammaTime - Vtxtime)) );
    
    dot = showerDir2D.Dot(gammaDir2D);
    dot /= showerDir2D.Mag();
    dot /= gammaDir2D.Mag();
    
    return;
  }

  /**
   * @brief get Principal Component Analysis components for 2D hits
   * @input hits -> cluster hits
   * @input wire2cm -> wire to cm conversion
   * @input time2cm -> tick to cm conversion
   * @output eigenVal -> set of eigenvalues for cluster
   * @output eigenVec -> eigenvectors for cluster
   */
  void  PCA(const std::vector<art::Ptr<recob::Hit>> &hits, 
	    const float wire2cm, const float time2cm,
	    TVectorD& eigenVal, TMatrixD& eigenVec) {
    
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    // x -> wire
    // y -> time
    
    // find average w,t
    float wavg = 0;
    float tavg = 0;
    for (size_t h=0; h < hits.size(); h++) {
      
      auto hit = hits.at(h);
      auto wire = hit->WireID().Wire * wire2cm;
      auto time = ( hit->PeakTime() - detp->TriggerOffset() )   * time2cm;
      
      wavg += wire;
      tavg += time;
      
    }// for all his
    
    wavg /= hits.size();
    tavg /= hits.size();
    
    float norm = 1. / hits.size();
    
    // build matrix
  TMatrixDSym HitMatrix(2);

  for (size_t h=0; h < hits.size(); h++) {

    auto hit = hits.at(h);
    auto wire = hit->WireID().Wire * wire2cm;
    auto time = ( hit->PeakTime() - detp->TriggerOffset() )   * time2cm;
    
    double x = wire - wavg;
    double y = time - tavg;
    
    HitMatrix(0,0) += x*x*norm;
    HitMatrix(0,1) += x*y*norm;
    HitMatrix(1,0) += x*y*norm;
    HitMatrix(1,1) += y*y*norm;
    
  }// for all hits
  
  const TMatrixDSymEigen me(HitMatrix);
  eigenVal = me.GetEigenValues();
  eigenVec = me.GetEigenVectors();
  
  std::cout << "Matrix contents : [ " << HitMatrix(0,0) << ", " << HitMatrix(1,0) << " ], [ "<< HitMatrix(0,1) << ", " << HitMatrix(1,1) << " ] " << std::endl;

  for (int i=0; i<2; ++i)
    std::cout << "\t eigenvalue " << eigenVal(i) << " has eigenvector [" << eigenVec(0, i) << "," << eigenVec(1, i) << " ]" << std::endl;

  return;
  }



  /**
   * @brief Merge 2D clusters in the event with shower, if appropriate
   * @input event -> full event record
   * @input fClusterproducer -> proucer for 2D clusters associated to slice un-clustered hits
   * @input fHitproducer -> producer for default 2D hits in the event (gaushit)
   * @input shr_pfp_idx -> index of main shower candidate PFP
   * @input slice_pfp_v -> all PFParticles in slice
   */
  void Merge2DClusters(art::Event const& e,
		       const art::InputTag fClusterproducer,
		       const art::InputTag fHitproducer,
		       const size_t shr_pfp_idx,
		       const std::vector<searchingfornues::ProxyPfpElem_t>& slice_pfp_v) {

    // get detector specific properties
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    float _wire2cm = geom->WirePitch(0,0,0);
    float _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

    // get shower candidate
    if (slice_pfp_v.size() <= shr_pfp_idx) return ;
    auto ass_shr_v = slice_pfp_v[shr_pfp_idx].get<recob::Shower>();
    if (ass_shr_v.size() != 1) return;
    std::cout << "\t shower @ idx " << shr_pfp_idx << std::endl;
    auto shr = ass_shr_v[0];
    // get shower 3D direction and starting point
    auto shrVtx = shr->ShowerStart(); // xyz coordinate [TVector3]
    auto shrDir = shr->Direction(); // unit vector [TVector3]
    
    // grab clusters themselves
    auto const& cluster_h = e.getValidHandle<std::vector<recob::Cluster> >(fClusterproducer);
    // get hits associated to clusters
    art::FindManyP<recob::Hit> clus_hit_assn_v(cluster_h, e, fClusterproducer);
    // grab hits themselves
    auto const& hit_h = e.getValidHandle<std::vector<recob::Hit> >(fHitproducer);
    
    // loop through clusters
    for (size_t c=0; c < cluster_h->size(); c++) {
      
      //auto clus = cluster_h->at(c);
      
      // get associated hits
      auto clus_hit_v = clus_hit_assn_v.at( c );
      
      if (clus_hit_v.size() < 10) continue;
      
      // is the proton isolated?
      //if (IsProtonIsolated(clus_hit_v,hit_h) == false) continue;
      
      // create vector of gaushits corresponding to new proton hits
      auto gaushit_hit_v = searchingfornues::getGaussHits(clus_hit_v, hit_h);
      
      //auto clus_tmin = clus.StartTick() * _time2cm;
      //auto clus_tmax = clus.EndTick()   * _time2cm;
      
      auto plane = (gaushit_hit_v.at(0))->WireID().Plane;
      //auto nhit  = gaushit_hit_v.size();
      
      float shrdot = 1e6;
      float shrdist = 1e6;
      float charge = 0;
      
      
      float gammaWire = 0;
      float gammaTime = 0;
      
      for (size_t hi=0; hi < gaushit_hit_v.size(); hi++) {
	auto hit = gaushit_hit_v.at(hi);
	gammaWire += hit->WireID().Wire * _wire2cm * hit->Integral();
	gammaTime += (hit->PeakTime() - detp->TriggerOffset())  * _time2cm * hit->Integral();
	charge += hit->Integral();
      }
      gammaWire /= charge;
      gammaTime /= charge;
      
      std::cout << "Gamma [wire,time] -> [ " << gammaWire << ", " << gammaTime << " ]"  << std::endl;
      
      searchingfornues::GammaDot(gammaWire, gammaTime, plane, shrVtx, shrDir, _wire2cm, shrdot, shrdist);
      
      TVectorD eigenVal(2);
      TMatrixD eigenVec(2,2);
      PCA(clus_hit_v, _wire2cm, _time2cm, eigenVal, eigenVec);
      
      //float eigenratio = eigenVal(0) / eigenVal(1);
      
      // get dot-product between principal eigenvector and shower direction
      //float eigendot = searchingfornues::EigenDot(_plane, ShrDir, eigenVec(0,0), eigenVec(0,1));
      
    }
    
    return;
  }
  
  
  /**
   * @brief find PFParticles associated to the shower turnk
   * @input shr_pfp_idx -> index of main shower candidate PFP
   * @input slice_pfp_v -> all PFParticles in slice
   * @output vtxcandidate -> new candidate vertex based on track-merging
   * @output trktrunk_pfp_idx -> which PFP was merged as trunk, if any?
   * @return bool for whether anything was merged
   */
  bool FindShowerTrunk(const size_t shr_pfp_idx,
		       const std::vector<searchingfornues::ProxyPfpElem_t>& slice_pfp_v,
		       TVector3& vtxcandidate,
		       size_t& trktrunk_pfp_idx,
		       float& bestmgdot,
		       float& bestmgdist) {

    //    float fShrVtxTrkDistMax = 0; // ToDo -> add as input
    //float fShrTrkDotMin     = 0; // ToDo -> add as input

    // do any of the slices have a track-like segment aligned with the shower?
    
    if (slice_pfp_v.size() <= shr_pfp_idx) return false;
    
    auto ass_shr_v = slice_pfp_v[shr_pfp_idx].get<recob::Shower>();
    
    if (ass_shr_v.size() != 1) return false;
    
    std::cout << "\t shower @ idx " << shr_pfp_idx << std::endl;
    
    auto shr = ass_shr_v[0];
    // get shower 3D direction and starting point
    auto shrVtx = shr->ShowerStart(); // xyz coordinate [TVector3]
    auto shrDir = shr->Direction(); // unit vector [TVector3]
    
    // keep track of possible matching track
    // track with best alignment in dot product is the candidate
    float bestdot = 0.;
    float bestdist = 1e6;
    
    for (size_t p=0; p < slice_pfp_v.size(); p++) {
      
      // skip the pfparticle associated to the shower itself
      if (p == shr_pfp_idx) continue;
      
      // shower trunk needs to be track-like to be merged by this algorithm
      // filter out non-track like PFPs
      
      auto ass_trk_v = slice_pfp_v[p].get<recob::Track>();
      if (ass_trk_v.size() != 1) continue;
      
      auto trk = ass_trk_v[0];
      
      // track start direction / start point / end point
      TVector3 trkDir(trk->StartDirection().X(), trk->StartDirection().Y(), trk->StartDirection().Z());
      TVector3 trkVtx(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
      TVector3 trkEnd(trk->End().X(),   trk->End().Y(),   trk->End().Z());
      
      // compatibility?
      double TrkShrDot = trkDir.Dot(shrDir);
      // distance between shower start and track start/end
      // whicever is smallest
      double TrkShrDist = (shrVtx - trkVtx).Mag();
      if ( ((shrVtx - trkEnd).Mag() ) < TrkShrDist )
	TrkShrDist = (shrVtx - trkEnd).Mag();
      
      std::cout << "Compare new track..." << std::endl;
      std::cout << "TrkShr Dist is " << TrkShrDist << std::endl;
      std::cout << "TrkShr Dot  is " << TrkShrDot << std::endl;
      std::cout << std::endl;
      
      //      if ( (TrkShrDist < fShrVtxTrkDistMax) && (fabs(TrkShrDot) > fShrTrkDotMin) ) {
	
	// if this is the most aligned track
	if ( fabs(TrkShrDot) > fabs(bestdot) ) {
	  trktrunk_pfp_idx = p;
	  bestdot  = TrkShrDot;
	  bestdist = TrkShrDist;
	  // vertex candidate is track start/end depending on sign
	  // of dot-product.
	  if (bestdot > 0) { vtxcandidate = trkVtx; }
	  if (bestdot < 0) { vtxcandidate = trkEnd; }
	  // also save index of track to be mergd
	}
	
	// }// if agrees within user-defined specs
      
    }// for all PFParticles
    
    bestmgdot=bestdot;
    bestmgdist = bestdist;
    // if bestdot != 0 -> means we found a compatible match
    if (bestdot == 0)
      return false;
    
    return true;
  }// end FindShowerTrunk
  
  
  /**
   * @brief find PFParticles associated to the shower turnk
   * @input shr_pfp_idx -> index of main shower candidate PFP
   * @input slice_pfp_v -> all PFParticles in slice
   * @output trktrunk_pfp_idx_v -> which PFPs were merged as shower-branches, if any?
   * @return bool for whether anything was merged
   */
  bool FindShowerBranches(const size_t shr_pfp_idx,
			  const std::vector<searchingfornues::ProxyPfpElem_t>& slice_pfp_v,
			  std::vector<size_t>& branch_pfp_idx_v) {

    float fMinBranchConeDot = 0; // ToDo -> add as input
    float fMinBranchVtxDist = 0; // ToDo -> add as input

    // do any slices lie in the cone of an upstream shower?
    
    if (slice_pfp_v.size() <= shr_pfp_idx) return false;
    
    auto ass_shr_v = slice_pfp_v[shr_pfp_idx].get<recob::Shower>();
    
    if (ass_shr_v.size() != 1) return false;
    
    std::cout << "\t shower @ idx " << shr_pfp_idx << std::endl;
    
    auto shr = ass_shr_v[0];
    // get shower 3D direction and starting point
    auto shrVtx = shr->ShowerStart(); // xyz coordinate [TVector3]
    auto shrDir = shr->Direction(); // unit vector [TVector3]
    
    for (size_t p=0; p < slice_pfp_v.size(); p++) {
      
      // skip the pfparticle associated to the shower itself
      if (p == shr_pfp_idx) continue;
      
      // find start point and direction of downstream PFParticles
      // regardless of whether shower/track like
      TVector3 BranchVtx;
      TVector3 BranchDir;
      
      if (slice_pfp_v[p].get<recob::Track>().size() == 1) {
	
	auto trk = slice_pfp_v[p].get<recob::Track>()[0];
	
	BranchDir = TVector3(trk->StartDirection().X(), trk->StartDirection().Y(), trk->StartDirection().Z());
	BranchVtx = TVector3(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
	
      }// if associated to a track
      else if (slice_pfp_v[p].get<recob::Shower>().size() == 1) {
	
	auto shrB = slice_pfp_v[p].get<recob::Shower>()[0];
	
	BranchVtx = shrB->ShowerStart();
	BranchDir = shrB->Direction();
	
      }// if associated to a shower
      
      // is the branch candidate vertex compatible with the shower cone?
      
      // (1) angle between shower direction and shr start -> branch start
      TVector3 Branch2ShrDir = (BranchVtx-shrVtx).Unit();
      double branchDot = Branch2ShrDir.Dot( shrDir );
      
      std::cout << " Branch finding dot  : " << branchDot << std::endl;
      
      
      if (branchDot < fMinBranchConeDot) continue;
      
      // (2) make sure branch starts enough downstream
      double BranchDist = (BranchVtx-shrVtx).Mag();
      
      std::cout << " Branch finding dist : " << BranchDist << std::endl;
      
      if ( BranchDist < fMinBranchVtxDist ) continue;
      
      // made it this far -> add candidate branch
      branch_pfp_idx_v.push_back( p );
      
    }// for all PFParticles
    
    
    // if bestdot != 0 -> means we found a compatible match
    if (branch_pfp_idx_v.size() == 0)
      return false;
    
    return true;
  }// end FindShowerTrunk
  
} // namespace searchingfornues

#endif
