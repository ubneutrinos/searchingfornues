////////////////////////////////////////////////////////////////////////
// Class:       ProtonTruthStudies
// Plugin Type: analyzer (art v3_01_02)
// File:        ProtonTruthStudies_module.cc
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

class ProtonTruthStudies;


class ProtonTruthStudies : public art::EDAnalyzer {
public:
  explicit ProtonTruthStudies(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ProtonTruthStudies(ProtonTruthStudies const&) = delete;
  ProtonTruthStudies(ProtonTruthStudies&&) = delete;
  ProtonTruthStudies& operator=(ProtonTruthStudies const&) = delete;
  ProtonTruthStudies& operator=(ProtonTruthStudies&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  art::InputTag fClusterproducer, fHTproducer, fHitproducer, fMCSproducer, fPFPproducer, fSHRproducer;

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
  float _energy1, _energy2; // proton KE [MeV]
  float _pfppur1, _clspur1, _pfppur2, _clspur2;
  int _hitpfp1, _hitpfp2, _hitcls1, _hitcls2;

  void ResetTTree();

  std::vector<art::Ptr<recob::Hit>> getGaussHits(const std::vector<art::Ptr<recob::Hit>> &hits,
						 const art::ValidHandle<std::vector<recob::Hit> > gaushit_h);

  int IsProtonIsolated(const std::vector<art::Ptr<recob::Hit>> &hits,
		       const art::ValidHandle<std::vector<recob::Hit> > gaushit_h);

  void ProtonDot(const float& protonWire, const float& protonTime,const int& pl,
		 const TVector3& showerVtx, const TVector3& showerDir,
		 float &dot, float& d2d);

  float IoU(const int& idx, const float& time_min, const float& time_max, const int& pl,
	    const std::map< size_t, std::pair<float,float> >& protonIdx2TimeSpanMap,
	    const std::map< size_t, int >& protonIdx2PlaneMap,
	    const std::map< size_t, float >& protonIdx2ChargeMap,
	    size_t& matchedidx);
    
  
};

void ProtonTruthStudies::ResetTTree() {

  _evt = 0;
  _sub = 0;
  _run = 0;
  _energy1 = _energy2 = 0;
  _clspur1 = _clspur2 = 0;
  _hitpfp1 = _hitpfp2 = 0;
  _pfppur1 = _pfppur2 = 0;

}


ProtonTruthStudies::ProtonTruthStudies(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

  fClusterproducer = p.get< art::InputTag > ("Clusterproducer"  );
  fHTproducer      = p.get< art::InputTag > ("HTproducer","gaushitTruthMatch");
  fHitproducer     = p.get< art::InputTag > ("Hitproducer");
  fMCSproducer     = p.get< art::InputTag > ("MCSproducer");
  fPFPproducer     = p.get< art::InputTag > ("PFPproducer");
  fSHRproducer     = p.get< art::InputTag > ("SHRproducer");

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree", "Proton Tagging");
  _tree->Branch("evt",&_evt,"evt/I");
  _tree->Branch("sub",&_sub,"sub/I");
  _tree->Branch("run",&_run,"run/I");
  _tree->Branch("energy1",&_energy1,"energy1/F");
  _tree->Branch("pfppur1",&_pfppur1,"pfppur1/F");
  _tree->Branch("clspur1",&_clspur1,"clspur1/F");
  _tree->Branch("hitpfp1",&_hitpfp1,"hitpfp1/I");
  _tree->Branch("hitcls1",&_hitcls1,"hitcls1/I");
  _tree->Branch("energy2",&_energy2,"energy2/F");
  _tree->Branch("pfppur2",&_pfppur2,"pfppur2/F");
  _tree->Branch("clspur2",&_clspur2,"clspur2/F");
  _tree->Branch("hitpfp2",&_hitpfp2,"hitpfp2/I");
  _tree->Branch("hitcls2",&_hitcls2,"hitcls2/I");

  // get detector specific properties
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _wire2cm = geom->WirePitch(0,0,0);
  _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
  
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}
  
void ProtonTruthStudies::analyze(art::Event const& e)
{

  ResetTTree();

  //auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  
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
  
  // load MCShower
  auto const& mcshower_h   = e.getValidHandle<std::vector<sim::MCShower> >(fMCSproducer);

  // grab PFParticles in event
  searchingfornues::ProxyPfpColl_t const& pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle> >(e,fPFPproducer,
													    proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
													    proxy::withAssociated<recob::Cluster>(fPFPproducer),
													    proxy::withAssociated<recob::Slice>(fPFPproducer),
													    proxy::withAssociated<recob::Track>(fPFPproducer),
													    proxy::withAssociated<recob::Vertex>(fPFPproducer),
													    proxy::withAssociated<recob::PCAxis>(fPFPproducer),
													    proxy::withAssociated<recob::Shower>(fSHRproducer),
                                                                                                            proxy::withAssociated<recob::SpacePoint>(fSHRproducer));
  
  searchingfornues::ProxyClusColl_t const& clus_proxy = proxy::getCollection<std::vector<recob::Cluster> >(e, fPFPproducer, proxy::withAssociated<recob::Hit>(fPFPproducer));
  
  // build PFParticle map  for this event
  BuildPFPMap(pfp_proxy);

  // loop through MCParticles and identify protons
  // fill a map for each MCParticle proton to the backtracked PFParticle or Cluster [value indicates <KE, purity> ]
  std::map<unsigned int, std::pair<float, float> > Gamma2PFPMap;
  std::map<unsigned int, int > Gamma2PFPMapNHit;
  std::map<unsigned int, std::pair<float, float> > Gamma2CLSMap;
  std::map<unsigned int, int > Gamma2CLSMapNHit;
  // map linking MCParticle TrackId to index
  std::map<unsigned int, unsigned int> MCParticleTID2IDXMap;

  int gamma1IDX = -1;
  int gamma2IDX = -1;
  float gamma1E = 0;
  float gamma2E = 0;

  // vector of BtPart for backtracking
  std::vector<searchingfornues::BtPart> Gamma_v;
  
  // loop through MCParticles and identify gammas
  for ( unsigned int im=0; im< mcshower_h->size(); ++im ) {

    const auto& mcs = mcshower_h->at(im);
    
    if ( ( (mcs.PdgCode() == 22) && (mcs.MotherPdgCode() == 111) && (mcs.Process() == "Decay") && (mcs.MotherProcess() == "primary") ) || (fabs(mcs.PdgCode()) == 11) ) {
      
      Gamma_v.push_back( searchingfornues::BtPart(mcs.PdgCode(),
						  mcs.Start().Momentum().Px(),
						  mcs.Start().Momentum().Py(),
						  mcs.Start().Momentum().Pz(),
						  mcs.DetProfile().E(), 
						  mcs.DaughterTrackID()) );

      float gammaE = mcs.DetProfile().E();

      std::cout << "gamma with energy " << gammaE << " has idx " << im << std::endl;

      std::cout << "IDX " << im << " assigned energy " << gammaE << std::endl;
      
      Gamma2PFPMap[im] = std::make_pair(gammaE,0.);
      Gamma2CLSMap[im] = std::make_pair(gammaE,0.);

      Gamma2PFPMapNHit[im] = 0;
      Gamma2CLSMapNHit[im] = 0;

      if (gammaE > gamma2E) {
	std::cout << "\t update gamma2" << std::endl;
	gamma2E = gammaE;
	gamma2IDX = im;
      }

      if (gammaE > gamma1E) {
	std::cout << "\t update gamma1 & gamma2" << std::endl;
	gamma2E = gamma1E;
	gamma2IDX = gamma1IDX;
	gamma1E = gammaE;
	gamma1IDX = im;
      }

      MCParticleTID2IDXMap[mcs.DaughterTrackID().at(0)] = im;
      
      std::cout << "MCS Energy : " << gammaE << ", E1 : " << gamma1E << " & idx " << gamma1IDX << ". E2 : " << gamma2E << " & idx " << gamma2IDX << std::endl;
      
    }// if pi0 gamma
    
  }// for all MCParticles

  std::cout << "Found " << Gamma2PFPMap.size() << " truth gammas" << std::endl;

  // found neutrinos?
  bool recoslice = false;

  // loop through reconstructed PFParticles
  for (const searchingfornues::ProxyPfpElem_t& pfp_pxy : pfp_proxy) {

    auto PDG = fabs(pfp_pxy->PdgCode());

    if ( (PDG == 12) || (PDG == 14) ) {
      recoslice = true;
      continue;
    }
    
    // get hits associated to this PFParticle through the clusters
    std::vector<art::Ptr<recob::Hit> > hit_v;
    auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();
    if (clus_pxy_v.size() != 0) {
      for (auto ass_clus : clus_pxy_v) {
	// get cluster proxy
	const auto& clus = clus_proxy[ass_clus.key()];
	auto clus_hit_v = clus.get<recob::Hit>();
	for (const auto& hit : clus_hit_v)
	  hit_v.push_back(hit);
      }// for all clusters associated to PFP
      
    }// if there are associated clusters
    

    // are there associated hits?
    if (hit_v.size() == 0) continue;
    

    float purity = 0.;
    float completeness = 0.;

    auto btpart = getAssocBtPart(hit_v,assocMCPart,Gamma_v,purity,completeness);

    purity *= hit_v.size();

    if (btpart >= 0) {
      
      auto backtrackedPart = Gamma_v.at(btpart);
      
      if (backtrackedPart.tids.size() > 0) {
	
	auto tid = backtrackedPart.tids.at(0);

	// is purity for this proton improved? if so update backtracked pfparticle
	if ( Gamma2PFPMap[ MCParticleTID2IDXMap[ tid ] ].second < purity ) {
	  
	  Gamma2PFPMap[ MCParticleTID2IDXMap[ tid ] ] = std::make_pair( backtrackedPart.e , purity );

	  Gamma2PFPMapNHit[ MCParticleTID2IDXMap[ tid ] ] = hit_v.size();

	  //std::cout << "Matched proton with PFP : TrackID " << tid << " and purity " << purity << " and completeness " << completeness << std::endl;

	}
	
      }// if there is at least one TrackID associated to this MCParticle
      
    }// if we backtracked to something
   
  }// for all PFPs 

  if (recoslice == false) return;

  // loop through tagged clusters in slice
  for (size_t c=0; c < cluster_h->size(); c++) {
    
    auto clus = cluster_h->at(c);
    
    // get associated hits
    auto clus_hit_v = clus_hit_assn_v.at( c );
    
    // create vector of gaushits corresponding to new proton hits
    auto gaushit_hit_v = getGaussHits(clus_hit_v, hit_h);

    float purity = 0;
    float completeness = 0;

    auto btpart = getAssocBtPart(gaushit_hit_v,assocMCPart,Gamma_v,purity,completeness);

    purity *= gaushit_hit_v.size();

    if (btpart >= 0) {
      
      auto backtrackedPart = Gamma_v.at(btpart);
      
      if (backtrackedPart.tids.size() > 0) {
	
	auto tid = backtrackedPart.tids.at(0);

	// is purity for this proton improved? if so update backtracked cluster
	if ( Gamma2CLSMap[ MCParticleTID2IDXMap[ tid ] ].second < purity ) {
	  
	  Gamma2CLSMap[ MCParticleTID2IDXMap[ tid ] ] = std::make_pair( backtrackedPart.e , purity );

	  Gamma2CLSMapNHit[ MCParticleTID2IDXMap[ tid ] ] = gaushit_hit_v.size();
	  
	  //std::cout << "Matched proton with CLUSTER : TrackID " << tid << " and purity " << purity << " and completeness " << completeness << std::endl;

	}
	
      }// if there is at least one TrackID associated to this MCParticle
      
    }// if we backtracked to something

  }// loop through all clusters in slice

  /*
  std::cout << "gamma1IDX is " << gamma1IDX << " and gamma2IDX is " << gamma2IDX << std::endl;

  std::cout << "e pfp map element 0 : " << Gamma2PFPMap[0].first << std::endl;
  std::cout << "e cls map element 0 : " << Gamma2CLSMap[0].first << std::endl;
  std::cout << "e pfp map element 1 : " << Gamma2PFPMap[1].first << std::endl;
  std::cout << "e cls map element 1 : " << Gamma2CLSMap[1].first << std::endl;
  */

  if (gamma1IDX >= 0) {
    
    // get PFP idx , pur pair
    auto PFPmatch1 = Gamma2PFPMap[gamma1IDX];
    // get CLS idx, pur pair
    auto CLSmatch1 = Gamma2CLSMap[gamma1IDX];
    
    std::cout << "Gamma with energy " << PFPmatch1.first
	      << " has PFP purity of " << PFPmatch1.second
	      << " and CLS purity of " << CLSmatch1.second << std::endl;
    
    _energy1 = PFPmatch1.first;
    _hitpfp1 = Gamma2PFPMapNHit[gamma1IDX];
    _hitcls1 = Gamma2CLSMapNHit[gamma1IDX];
    _pfppur1 = _clspur1 = 0;
    if (_hitpfp1 > 0)
      _pfppur1 = PFPmatch1.second / _hitpfp1;
    if (_hitcls1 > 0)
      _clspur1 = CLSmatch1.second / _hitcls1;

  }

  if (gamma2IDX >= 0) {
    
    // get PFP idx , pur pair
    auto PFPmatch2 = Gamma2PFPMap[gamma2IDX];
    // get CLS idx, pur pair
    auto CLSmatch2 = Gamma2CLSMap[gamma2IDX];
    
    std::cout << "Gamma with energy " << PFPmatch2.first
	      << " has PFP purity of " << PFPmatch2.second
	      << " and CLS purity of " << CLSmatch2.second << std::endl;
    
    _energy2 = PFPmatch2.first;
    _hitpfp2 = Gamma2PFPMapNHit[gamma2IDX];
    _hitcls2 = Gamma2CLSMapNHit[gamma2IDX];
    _pfppur2 = _clspur2 = 0;
    if (_hitpfp2 > 0)
      _pfppur2 = PFPmatch2.second / _hitpfp2;
    if (_hitcls2 > 0)
      _clspur2 = CLSmatch2.second / _hitcls2;

  }

  
  _tree->Fill();

return;
}

void ProtonTruthStudies::ProtonDot(const float& protonWire, const float& protonTime,const int& pl,
				   const TVector3& showerVtx, const TVector3& showerDir,
				   float &dot, float& d2d) {
  
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  
  std::cout << "3D shower dir : [ " << showerDir[0] << ", " << showerDir[1] << ", " << showerDir[2] << " ]" << std::endl;
  std::cout << "3D shower vtx : [ " << showerVtx[0] << ", " << showerVtx[1] << ", " << showerVtx[2] << " ]" << std::endl;
  
  auto Vtxwire = geom->WireCoordinate(showerVtx[1],showerVtx[2],geo::PlaneID(0,0,pl)) * _wire2cm;
  auto Vtxtime = showerVtx[0];
  
  auto Dirwire = geom->WireCoordinate(showerDir[1],showerDir[2],geo::PlaneID(0,0,pl)) * _wire2cm;
  auto Dirtime = showerDir[0];
  
  TVector3 showerDir2D(Dirwire,Dirtime,0.);
  TVector3 protonDir2D(protonWire-Vtxwire,protonTime-Vtxtime,0.);

  std::cout << "2D shower dir : [ " << Dirwire << ", " << Dirtime << " ]" << std::endl;
  std::cout << "2D shower vtx : [ " << Vtxwire << ", " << Vtxtime << " ]" << std::endl;

  std::cout << "2D proton pos : [ " << protonWire << ", " << protonTime << " ]" << std::endl;
  
  std::cout << std::endl;
  
  d2d = sqrt( ((protonWire - Vtxwire) * (protonWire - Vtxwire)) +
	      ((protonTime - Vtxtime) * (protonTime - Vtxtime)) );
  
  dot  = showerDir2D.Dot(protonDir2D);
  dot /= showerDir2D.Mag();
  dot /= protonDir2D.Mag();
  
  return;
}

float ProtonTruthStudies::IoU(const int& idx, const float& time_min, const float& time_max, const int& pl,
			   const std::map< size_t, std::pair<float,float> >& protonIdx2TimeSpanMap,
			   const std::map< size_t, int >& protonIdx2PlaneMap,
			   const std::map< size_t, float >& protonIdx2ChargeMap,
			   size_t& matchedidx) {

  float ioumin = 1e4;
  float dtmin = 1e4;
  matchedidx = 999;

  for (auto const& clus : protonIdx2PlaneMap) {
    
    auto c2idx = clus.first;

    if (clus.second == pl) continue; // cannot match to cluster on the same plane!
    
    float c1time = (time_min+time_max)/2.;
    float c2timespanmin = protonIdx2TimeSpanMap.at(c2idx).first;
    float c2timespanmax = protonIdx2TimeSpanMap.at(c2idx).second;
    float c2time = (c2timespanmin + c2timespanmax) / 2.;
    
    float iou = c1time - c2time;
    float dt = fabs(iou);

    if (dt < dtmin) { dtmin = dt; matchedidx = c2idx; ioumin = iou; }

  }// for all clusters

  return ioumin;
}// IoU

std::vector<art::Ptr<recob::Hit>> ProtonTruthStudies::getGaussHits(const std::vector<art::Ptr<recob::Hit>> &hits,
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

int ProtonTruthStudies::IsProtonIsolated(const std::vector<art::Ptr<recob::Hit>> &hits,
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

  return nclose;
}

 void ProtonTruthStudies::BuildPFPMap(const searchingfornues::ProxyPfpColl_t& pfp_pxy_col) {
   
   _pfpmap.clear();
   
   unsigned int p=0;
   for (const auto& pfp_pxy : pfp_pxy_col) {
     _pfpmap[pfp_pxy->Self()] = p;
     p++;
   }
   
   return;
 }// BuildPFPMap
 
 void ProtonTruthStudies::AddDaughters(const searchingfornues::ProxyPfpElem_t& pfp_pxy,
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
 

void ProtonTruthStudies::beginJob()
{
  // Implementation of optional member function here.
}

void ProtonTruthStudies::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ProtonTruthStudies)
