////////////////////////////////////////////////////////////////////////
// Class:       SecondShowerPurity
// Plugin Type: analyzer (art v3_01_02)
// File:        SecondShowerPurity_module.cc
//
// Generated at Sun Jun 23 07:31:37 2019 by David Caratelli using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "../Selection/CommonDefs/Typedefs.h"
#include "../Selection/CommonDefs/TrackShowerScoreFuncs.h"
#include "../Selection/CommonDefs/BacktrackingFuncs.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"
#include "TMatrixDSymEigen.h"

class SecondShowerPurity;


class SecondShowerPurity : public art::EDAnalyzer {
public:
  explicit SecondShowerPurity(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SecondShowerPurity(SecondShowerPurity const&) = delete;
  SecondShowerPurity(SecondShowerPurity&&) = delete;
  SecondShowerPurity& operator=(SecondShowerPurity const&) = delete;
  SecondShowerPurity& operator=(SecondShowerPurity&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  art::InputTag fClusterproducer, fHTproducer, fHitproducer, fPFPproducer, fSHRproducer, fMCSproducer;

  /**
   * @brief function to builf a map linking PFParticle index to Self() attribute
   *
   * @input handle to event pfparticle record
   */
  void BuildPFPMap(const searchingfornues::ProxyPfpColl_t& pfp_pxy_col);

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;  
  
  /**
   * @brief build PFParticle hierarchy (i.e. slice) from parent [recursive function]
   *
   * @input pfp_pxy : parent pfparticle proxy for which to add daughters
   * @input pfp_pxy_col : evnt PFP proxy collection
   * @input slice_v : passed by reference, slice containing all PFParticles in hierarchy
   *
   */
  void AddDaughters(const searchingfornues::ProxyPfpElem_t& pfp_pxy,
		    const searchingfornues::ProxyPfpColl_t& pfp_pxy_col,
		    std::vector<searchingfornues::ProxyPfpElem_t>& slice_v);
  
  float _wire2cm, _time2cm;

  // TTree
  TTree* _tree;
  int _evt, _sub, _run;
  int _nhit;
  int _plane;
  float _charge;
  float _purity, _shrpurity;
  float _shrdot;
  float _shrdist;
  float _shrenergy1;
  float _shrenergy2;
  float _elecenergy;
  int _btpdg, _btshrpdg; // PDG of backtracked particle
  float _btmom, _btshrmom;
  int _nshr;
  float _clus_tmin, _clus_tmax;
  float _twoshrdot;
  float _eigenratio, _eigendot;
  float _gammaemin;

  // event TTree

  void ResetTTree();

  std::vector<art::Ptr<recob::Hit>> getGaussHits(const std::vector<art::Ptr<recob::Hit>> &hits,
						 const art::ValidHandle<std::vector<recob::Hit> > gaushit_h);

  bool IsProtonIsolated(const std::vector<art::Ptr<recob::Hit>> &hits,
			const art::ValidHandle<std::vector<recob::Hit> > gaushit_h);

  void GammaDot(const float& gammaWire, const float& gammaTime,const int& pl,
		const TVector3& showerVtx, const TVector3& showerDir,
		float &dot, float& d2d);

  void  PCA(const std::vector<art::Ptr<recob::Hit>> &hits, TVectorD& eigenVal, TMatrixD& eigenVec);

  float EigenDot(const int& pl,const TVector3& ShowerDir,const float& gammaWireDir,const float& gammaTimeDir);
  
};


SecondShowerPurity::SecondShowerPurity(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

  fClusterproducer = p.get< art::InputTag > ("Clusterproducer"  );
  fHTproducer      = p.get< art::InputTag > ("HTproducer","gaushitTruthMatch");
  fHitproducer     = p.get< art::InputTag > ("Hitproducer");
  fMCSproducer     = p.get< art::InputTag > ("MCSproducer","mcreco");
  fPFPproducer     = p.get< art::InputTag > ("PFPproducer");
  fSHRproducer     = p.get< art::InputTag > ("SHRproducer");

  art::ServiceHandle<art::TFileService> tfs;

  _tree = tfs->make<TTree>("_tree", "Second Shower Tagging");
  _tree->Branch("evt",&_evt,"evt/I");
  _tree->Branch("sub",&_sub,"sub/I");
  _tree->Branch("run",&_run,"run/I");
  _tree->Branch("nhit",&_nhit,"nhit/I");
  _tree->Branch("nshr",&_nshr,"nshr/I");
  _tree->Branch("btpdg",&_btpdg,"btpdg/I");
  _tree->Branch("btmom",&_btmom,"btmom/F");
  _tree->Branch("btshrpdg",&_btshrpdg,"btshrpdg/I");
  _tree->Branch("btshrmom",&_btshrmom,"btshrmom/F");
  _tree->Branch("plane",&_plane,"plane/I");
  _tree->Branch("charge",&_charge,"charge/F");
  _tree->Branch("purity",&_purity,"purity/F");
  _tree->Branch("shrpurity",&_shrpurity,"shrpurity/F");
  _tree->Branch("shrdot",&_shrdot,"shrdot/F");
  _tree->Branch("twoshrdot",&_twoshrdot,"twoshrdot/F");
  _tree->Branch("shrdist",&_shrdist,"shrdist/F");
  _tree->Branch("shrenergy1",&_shrenergy1,"shrenergy1/F");
  _tree->Branch("shrenergy2",&_shrenergy2,"shrenergy2/F");
  _tree->Branch("elecenergy",&_elecenergy,"elecenergy/F");
  _tree->Branch("clus_tmin",&_clus_tmin,"clus_tmin/F");
  _tree->Branch("clus_tmax",&_clus_tmax,"clus_tmax/F");
  _tree->Branch("eigenratio",&_eigenratio,"eigenratio/F");
  _tree->Branch("eigendot",&_eigendot,"eigendot/F");
  _tree->Branch("gammaemin",&_gammaemin,"gammaemin/F");

  /*
  _evt_tree = tfs->make<TTree>("_event_tree", "Second Shower Tagging Event Tree");
  _evt_tree->Branch("shrenergy",&_shrenergy,"shrenergy/F");
  _evt_tree->Branch("elecenergy",&_elecenergy,"elecenergy/F");
  */

  // get detector specific properties
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _wire2cm = geom->WirePitch(0,0,0);
  _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
  
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void SecondShowerPurity::ResetTTree() {

  _nhit = 0;
  _nshr = 0;
  _btpdg  = 0;
  _btmom = 0;
  _plane = -1;
  _charge = 0;
  _purity = -1;
  _shrdot = -1;
  _shrdist = -1;
  _shrenergy1 = 0;
  _shrenergy2 = 0;
  _elecenergy = 0;
  _btshrpdg = 0;
  _btshrmom = 0;
  _shrpurity = -1;
  _twoshrdot = 0;
  _eigenratio = _eigendot = 0;

}
  
void SecondShowerPurity::analyze(art::Event const& e)
{

  ResetTTree();

  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();

  // grab clusters themselves
  auto const& cluster_h = e.getValidHandle<std::vector<recob::Cluster> >(fClusterproducer);
  // get hits associated to clusters
  art::FindManyP<recob::Hit> clus_hit_assn_v(cluster_h, e, fClusterproducer);
  // grab hits themselves
  auto const& hit_h = e.getValidHandle<std::vector<recob::Hit> >(fHitproducer);

  // load hits <-> truth associations
  auto const& assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>(hit_h, e, fHTproducer));
  
  // load MCTruth
  //auto const& mcparticle_h = e.getValidHandle<std::vector<simb::MCParticle> >(fMCPproducer);
  auto const& mcshower_h   = e.getValidHandle<std::vector<sim::MCShower> >(fMCSproducer);
  

  // grab PFParticles in event
  searchingfornues::ProxyPfpColl_t const& pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle> >(e,fPFPproducer,
													    proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
													    proxy::withAssociated<recob::Cluster>(fPFPproducer),
													    proxy::withAssociated<recob::Slice>(fPFPproducer),
													    proxy::withAssociated<recob::Track>(fPFPproducer),
													    proxy::withAssociated<recob::Vertex>(fPFPproducer),
													    proxy::withAssociated<recob::PCAxis>(fPFPproducer),
													    proxy::withAssociated<recob::Shower>(fSHRproducer));
  
  // grab cluster proxy for PFParticles
  searchingfornues::ProxyClusColl_t const& clus_proxy = proxy::getCollection<std::vector<recob::Cluster> >(e, fPFPproducer, proxy::withAssociated<recob::Hit>(fPFPproducer));
  
  // build PFParticle map  for this event
  BuildPFPMap(pfp_proxy);


  // vector of BtPart for backtracking
  std::vector<searchingfornues::BtPart> Gamma_v;
  
  // loop through MCParticles and identify gammas
  float _gammaemin = 1e6;
  for ( unsigned int im=0; im< mcshower_h->size(); ++im ) {

    const auto& mcs = mcshower_h->at(im);
    
    if ( ( (mcs.PdgCode() == 22) && (mcs.MotherPdgCode() == 111) && (mcs.Process() == "Decay") && (mcs.MotherProcess() == "primary") ) || (fabs(mcs.PdgCode()) == 11) ) 
      
      Gamma_v.push_back( searchingfornues::BtPart(mcs.PdgCode(),
						  mcs.Start().Momentum().Px() * 0.001,
						  mcs.Start().Momentum().Py() * 0.001,
						  mcs.Start().Momentum().Pz() * 0.001,
						  mcs.Start().Momentum().E() * 0.001, 
						  mcs.DaughterTrackID()) );

    if (mcs.PdgCode() == 11) { _elecenergy = mcs.Start().Momentum().E(); }

    if (mcs.PdgCode() == 22 && (mcs.MotherPdgCode() == 111 && mcs.Process() == "Decay" && mcs.MotherProcess() == "primary") )   { 
      if (mcs.DetProfile().E() < _gammaemin) { _gammaemin = mcs.DetProfile().E(); }
    }// if pi0 gamma
    
  }// for all MCParticles


  // collect PFParticle hierarchy originating from this neutrino candidate
  std::vector<searchingfornues::ProxyPfpElem_t> slice_pfp_v;

  // grab candidate electron direction
  for (const searchingfornues::ProxyPfpElem_t& pfp_pxy : pfp_proxy) {
    //  find neutrino candidate
    if (pfp_pxy->IsPrimary() == false) continue;
    
    auto PDG = fabs(pfp_pxy->PdgCode());
    
    if (PDG != 12) continue;

    AddDaughters(pfp_pxy, pfp_proxy, slice_pfp_v);
    
  }// for all proxys
  
  TVector3 ShowerVtx1, ShowerDir1;
  TVector3 ShowerVtx2, ShowerDir2;
  _shrenergy1 = 0;
  _shrenergy2 = 0;
  
  // go through slice and find shower-like particles
  for (size_t p=0; p < slice_pfp_v.size(); p++) {
    
    auto pfp_pxy = slice_pfp_v[p];
    
    auto ass_shr_v = pfp_pxy.get<recob::Shower>();
    
    if (ass_shr_v.size() != 1) continue;
    
    auto trkshrscore = searchingfornues::GetTrackShowerScore(pfp_pxy);
    
    if (trkshrscore > 0.5) continue;
    
    auto shr = ass_shr_v[0];
    
    if (shr->Energy()[2] > _shrenergy1) {
    _shrenergy1 = shr->Energy()[2];
    ShowerVtx1 = shr->ShowerStart();
    ShowerDir1 = shr->Direction();  
    }    

    else if (shr->Energy()[2] > _shrenergy2) {
    _shrenergy2 = shr->Energy()[2];
    ShowerVtx2 = shr->ShowerStart();
    ShowerDir2 = shr->Direction();  
    }    

    // backtrack this shower candidate
    // get hits associated to this PFParticle through the clusters
    std::vector<art::Ptr<recob::Hit> > pfp_hit_v;
    auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();
    if (clus_pxy_v.size() != 0) {
      for (auto ass_clus : clus_pxy_v) {
	// get cluster proxy
	const auto& clus = clus_proxy[ass_clus.key()];
	auto clus_hit_v = clus.get<recob::Hit>();
	for (const auto& hit : clus_hit_v)
	  pfp_hit_v.push_back(hit);
      }// for all clusters associated to PFP
    }
    float completeness = 0;
    _btshrpdg = 0;
    _btshrmom = 0;
    _shrpurity = 0;
    auto btpartidx = getAssocBtPart(pfp_hit_v,assocMCPart,Gamma_v,_shrpurity,completeness);
    if (btpartidx > 0) {
      _btshrpdg = Gamma_v.at(btpartidx).pdg;
      _btshrmom = Gamma_v.at(btpartidx).e;
    }
    
  }// for all proxys in slice

  _twoshrdot = ShowerDir1.Dot(ShowerDir2);
  _twoshrdot /= ShowerDir1.Mag();
  _twoshrdot /= ShowerDir2.Mag();

  // if no shower nue candidate -> skip event
  if (_shrenergy1 == 0)
    return;

  size_t nclusters = 0;
  
  // loop through clusters
  for (size_t c=0; c < cluster_h->size(); c++) {
    auto clus = cluster_h->at(c);
    
    // get associated hits
    auto clus_hit_v = clus_hit_assn_v.at( c );

    if (clus_hit_v.size() < 10) continue;

    nclusters += 1;
    
    float completeness = 0;

    // is the proton isolated?
    //if (IsProtonIsolated(clus_hit_v,hit_h) == false) continue;
    
    // create vector of gaushits corresponding to new proton hits
    auto gaushit_hit_v = getGaussHits(clus_hit_v, hit_h);
    
    _clus_tmin = clus.StartTick() * _time2cm;
    _clus_tmax = clus.EndTick()   * _time2cm;
    
    _plane = (gaushit_hit_v.at(0))->WireID().Plane;
    _nhit  = gaushit_hit_v.size();
    _charge = 0.;
    _purity = 0.;
    
    _shrdot = 1e6;
    _shrdist = 1e6;
    
    float gammaWire = 0;
    float gammaTime = 0;

    for (size_t hi=0; hi < gaushit_hit_v.size(); hi++) {
      auto hit = gaushit_hit_v.at(hi);
      gammaWire += hit->WireID().Wire * _wire2cm * hit->Integral();
      gammaTime += (hit->PeakTime() - detp->TriggerOffset())  * _time2cm * hit->Integral();
      _charge += hit->Integral();
    }
    gammaWire /= _charge;
    gammaTime /= _charge;

    std::cout << "Gamma [wire,time] -> [ " << gammaWire << ", " << gammaTime << " ]"  << std::endl;
    
    GammaDot(gammaWire, gammaTime, _plane, ShowerVtx1, ShowerDir1, _shrdot, _shrdist);

    TVectorD eigenVal(2);
    TMatrixD eigenVec(2,2);
    PCA(clus_hit_v, eigenVal, eigenVec);

    _eigenratio = eigenVal(0) / eigenVal(1);

    // get dot-product between principal eigenvector and shower direction
    _eigendot = EigenDot(_plane, ShowerDir1, eigenVec(0,0), eigenVec(0,1));
    
    auto btpartidx = getAssocBtPart(gaushit_hit_v,assocMCPart,Gamma_v,_purity,completeness);
    if (btpartidx < 0) {
      _btpdg = 0;
      _btmom = 0;
    }
    else {
      _btpdg = Gamma_v.at(btpartidx).pdg; 
      _btmom = Gamma_v.at(btpartidx).e;
    }
    
    _tree->Fill();
    
    std::cout << "GAMMA Found " << gaushit_hit_v.size() << " matching original " << clus_hit_v.size() << " proton hits" << std::endl;
    std::cout << "GAMMA This cluster has " << gaushit_hit_v.size() << " hits and a purity/completenss of [ "
	      << _purity << ", " << completeness << " ]" << std::endl;
    std::cout << "GAMMA" << std::endl;
    
  }// for all clusters

  if (nclusters == 0) { _tree->Fill(); }
  
  return;
}
  
  
void SecondShowerPurity::GammaDot(const float& gammaWire, const float& gammaTime,const int& pl,
				  const TVector3& showerVtx, const TVector3& showerDir,
				  float &dot, float& d2d) {
  
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  
  auto Vtxwire = geom->WireCoordinate(showerVtx[1],showerVtx[2],geo::PlaneID(0,0,pl)) * _wire2cm;
  auto Vtxtime = showerVtx[0];
  
  auto Dirwire = geom->WireCoordinate(showerDir[1],showerDir[2],geo::PlaneID(0,0,pl)) * _wire2cm;
  auto Dirtime = showerDir[0];

  std::cout << "Shower Dir [x,y,z] -> [ " << showerDir[0] << ", " << showerDir[1] << ", " << showerDir[2] << " ]"  << std::endl;
  std::cout << "Shower Vtx [x,y,z] -> [ " << showerVtx[0] << ", " << showerVtx[1] << ", " << showerVtx[2] << " ]"  << std::endl;

  std::cout << "Shower Dir [wire,time] -> [ " << Dirwire << ", " << Dirtime << " ]"  << std::endl;
  std::cout << "Shower Vtx [wire,time] -> [ " << Vtxwire << ", " << Vtxtime << " ]"  << std::endl;


  
  TVector3 showerDir2D(Dirwire,Dirtime,0.);
  TVector3 gammaDir2D(gammaWire-Vtxwire,gammaTime-Vtxtime,0.);
  
  d2d = sqrt( ((gammaWire - Vtxwire) * (gammaWire - Vtxwire)) +
	      ((gammaTime - Vtxtime) * (gammaTime - Vtxtime)) );
  
  dot = showerDir2D.Dot(gammaDir2D);
  dot /= showerDir2D.Mag();
  dot /= gammaDir2D.Mag();
  
  std::cout << "dot product : " << dot << std::endl;
  std::cout << std::endl;

  return;
}

float SecondShowerPurity::EigenDot(const int& pl,const TVector3& ShowerDir,const float& gammaWireDir,const float& gammaTimeDir) { 

  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  
  auto Dirwire = geom->WireCoordinate(ShowerDir[1],ShowerDir[2],geo::PlaneID(0,0,pl)) * _wire2cm;
  auto Dirtime = ShowerDir[0];

  float dot = Dirwire * gammaWireDir + Dirtime * gammaTimeDir;
  dot /= sqrt( Dirwire * Dirwire + Dirtime * Dirtime);
  dot /= sqrt( gammaWireDir * gammaWireDir + gammaTimeDir * gammaTimeDir);

  return dot;
}

std::vector<art::Ptr<recob::Hit>> SecondShowerPurity::getGaussHits(const std::vector<art::Ptr<recob::Hit>> &hits,
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

void  SecondShowerPurity::PCA(const std::vector<art::Ptr<recob::Hit>> &hits, TVectorD& eigenVal, TMatrixD& eigenVec) {

  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  // x -> wire
  // y -> time

  // find average w,t
  float wavg = 0;
  float tavg = 0;
  for (size_t h=0; h < hits.size(); h++) {
    
    auto hit = hits.at(h);
    auto wire = hit->WireID().Wire * _wire2cm;
    auto time = ( hit->PeakTime() - detp->TriggerOffset() )   * _time2cm;

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
    auto wire = hit->WireID().Wire * _wire2cm;
    auto time = ( hit->PeakTime() - detp->TriggerOffset() )   * _time2cm;

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

bool SecondShowerPurity::IsProtonIsolated(const std::vector<art::Ptr<recob::Hit>> &hits,
				       const art::ValidHandle<std::vector<recob::Hit> > gaushit_h) {

  float TimeMin = 1e6;
  float TimeMax = 0;
  float WireMin = 1e4;
  float WireMax = 0;

  std::vector< std::pair<float,float> > protonHitCoordinates_v;

  for (size_t h1=0; h1 < hits.size(); h1++) {
    
    auto protonhit = hits.at(h1);
    auto wire = protonhit->WireID().Wire * _wire2cm;
    auto time = protonhit->PeakTime() * _time2cm;

    protonHitCoordinates_v.push_back( std::make_pair(wire,time) );

    if (time < TimeMin) { TimeMin = time; }
    if (time > TimeMax) { TimeMax = time; }
    if (wire < WireMin) { WireMin = wire; }
    if (wire > WireMax) { WireMax = wire; }

  }// end proton loop

  size_t nclose = 0;

  // loop through gauss hits

  for (size_t h2=0; h2 < gaushit_h->size(); h2++) {
    
      auto gaushit = gaushit_h->at(h2);

      auto gTime = gaushit.PeakTime() * _time2cm;
      auto gWire = gaushit.WireID().Wire * _wire2cm;

      if (gTime < (TimeMin - 5) ) continue;
      if (gTime > (TimeMax + 5) ) continue;
      if (gWire < (WireMin - 5) ) continue;
      if (gWire > (WireMax + 5) ) continue;

      float d2dmin = 1e4;

      // made it this far, the hit is in a 5 cm box. compute 2D distance
      for (size_t p=0; p < protonHitCoordinates_v.size(); p++) {

	auto coord = protonHitCoordinates_v.at(p);
	float d2d = sqrt( (coord.first - gWire)*(coord.first - gWire) + (coord.second - gTime)*(coord.second - gTime) );

	if (d2d < d2dmin) { d2dmin = d2d; }
	
      }// for all proton hits

      // if hits are overlapping -> skip
      if (d2dmin < 1e-3) continue;
      
      if (d2dmin < 2.0) { nclose += 1; }
      
  }// for all gaushits

  if (nclose < 3) return true;

  return false;
}

 void SecondShowerPurity::BuildPFPMap(const searchingfornues::ProxyPfpColl_t& pfp_pxy_col) {
   
   _pfpmap.clear();
   
   unsigned int p=0;
   for (const auto& pfp_pxy : pfp_pxy_col) {
     _pfpmap[pfp_pxy->Self()] = p;
     p++;
   }
   
   return;
 }// BuildPFPMap
 
 void SecondShowerPurity::AddDaughters(const searchingfornues::ProxyPfpElem_t& pfp_pxy,
				    const searchingfornues::ProxyPfpColl_t& pfp_pxy_col,
				    std::vector<searchingfornues::ProxyPfpElem_t>& slice_v) {
   
   auto daughters = pfp_pxy->Daughters();
   
   slice_v.push_back(pfp_pxy);
   
   for(auto const& daughterid : daughters) {
     
     if (_pfpmap.find(daughterid) == _pfpmap.end()) {
       std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
       continue;
     }
     
     // const art::Ptr<recob::PFParticle> pfp_pxy(pfp_pxy_col, _pfpmap.at(daughterid) );
     auto pfp_pxy2 = pfp_pxy_col.begin();
     for (size_t j=0; j<_pfpmap.at(daughterid); ++j) ++pfp_pxy2;
     // const T& pfp_pxy2 = (pfp_pxy_col.begin()+_pfpmap.at(daughterid));
     
     AddDaughters(*pfp_pxy2, pfp_pxy_col, slice_v);
     
   }// for all daughters
   
   return;
 }// AddDaughters
 

void SecondShowerPurity::beginJob()
{
  // Implementation of optional member function here.
}

void SecondShowerPurity::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(SecondShowerPurity)
