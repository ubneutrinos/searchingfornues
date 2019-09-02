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
  art::InputTag fClusterproducer, fHTproducer, fHitproducer, fMCPproducer, fPFPproducer, fSHRproducer;

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
  float _energy; // proton KE [MeV]
  float _pfppur, _clspur;
  int _ismaxproton; // is this the highest energy proton in the event?

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

}


ProtonTruthStudies::ProtonTruthStudies(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

  fClusterproducer = p.get< art::InputTag > ("Clusterproducer"  );
  fHTproducer      = p.get< art::InputTag > ("HTproducer","gaushitTruthMatch");
  fHitproducer     = p.get< art::InputTag > ("Hitproducer");
  fMCPproducer     = p.get< art::InputTag > ("MCPproducer");
  fPFPproducer     = p.get< art::InputTag > ("PFPproducer");
  fSHRproducer     = p.get< art::InputTag > ("SHRproducer");

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree", "Proton Tagging");
  _tree->Branch("evt",&_evt,"evt/I");
  _tree->Branch("sub",&_sub,"sub/I");
  _tree->Branch("run",&_run,"run/I");
  _tree->Branch("energy",&_energy,"energy/F");
  _tree->Branch("pfppur",&_pfppur,"pfppur/F");
  _tree->Branch("clspur",&_clspur,"clspur/F");
  _tree->Branch("ismaxproton",&_ismaxproton,"ismaxproton/I");

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
  
  // load MCTruth
  auto const& mcparticle_h = e.getValidHandle<std::vector<simb::MCParticle> >(fMCPproducer);
  

  // grab PFParticles in event
  searchingfornues::ProxyPfpColl_t const& pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle> >(e,fPFPproducer,
													    proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
													    proxy::withAssociated<recob::Cluster>(fPFPproducer),
													    proxy::withAssociated<recob::Slice>(fPFPproducer),
													    proxy::withAssociated<recob::Track>(fPFPproducer),
													    proxy::withAssociated<recob::Vertex>(fPFPproducer),
													    proxy::withAssociated<recob::PCAxis>(fPFPproducer),
													    proxy::withAssociated<recob::Shower>(fSHRproducer));
  
  searchingfornues::ProxyClusColl_t const& clus_proxy = proxy::getCollection<std::vector<recob::Cluster> >(e, fPFPproducer, proxy::withAssociated<recob::Hit>(fPFPproducer));
  
  // build PFParticle map  for this event
  BuildPFPMap(pfp_proxy);

  // loop through MCParticles and identify protons
  // fill a map for each MCParticle proton to the backtracked PFParticle or Cluster [value indicates <KE, purity> ]
  std::map<unsigned int, std::pair<float, float> > Proton2PFPMap;
  std::map<unsigned int, std::pair<float, float> > Proton2CLSMap;
  // map linking MCParticle TrackId to index
  std::map<unsigned int, unsigned int> MCParticleTID2IDXMap;

  // vector of BtPart for backtracking
  std::vector<searchingfornues::BtPart> Proton_v;

  size_t maxprotonIDX = 0;
  float maxprotonKE = 0.;

  for ( unsigned int im=0; im< mcparticle_h->size(); ++im ) {
    
    const auto& mcp = mcparticle_h->at(im);

    if ( (mcp.PdgCode() == 2212) && (mcp.StatusCode() == 1) && (mcp.Process() == "primary") ) {
      


      Proton_v.push_back( searchingfornues::BtPart(mcp.PdgCode(),
						   mcp.Momentum(0).Px() ,
						   mcp.Momentum(0).Py() ,
						   mcp.Momentum(0).Pz() ,
						   mcp.Momentum(0).E() ,  // in MeV, Kinetic Energy
						   mcp.TrackId()) );

      float ke = (mcp.Momentum(0).E() - 0.938277) * 1000.;
      
      Proton2PFPMap[im] = std::make_pair(ke,0.);
      Proton2CLSMap[im] = std::make_pair(ke,0.);

      if (ke > maxprotonKE) {
	maxprotonIDX = im;
	maxprotonKE  = ke;
      }
      
      MCParticleTID2IDXMap[mcp.TrackId()] = im;
      
      std::cout << "True proton w/ energy : " << mcp.Momentum(0).E() 
		<< " \t TrackID : " << mcp.TrackId() 
		<< std::endl;
      
    }// if a proton
    
  }// for all MCParticles

  std::cout << "Found " << Proton2PFPMap.size() << " truth protons" << std::endl;

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
    auto btpart = getAssocBtPart(hit_v,assocMCPart,Proton_v,purity,completeness);

    if (btpart >= 0) {
      
      auto backtrackedPart = Proton_v.at(btpart);
      
      if (backtrackedPart.tids.size() > 0) {
	
	auto tid = backtrackedPart.tids.at(0);

	// is purity for this proton improved? if so update backtracked pfparticle
	if ( Proton2PFPMap[ MCParticleTID2IDXMap[ tid ] ].second < purity ) {
	  
	  Proton2PFPMap[ MCParticleTID2IDXMap[ tid ] ] = std::make_pair( 1000.*(backtrackedPart.e-0.938277) , purity );

	  std::cout << "Matched proton with PFP : TrackID " << tid << " and purity " << purity << std::endl;

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

    auto btpart = getAssocBtPart(gaushit_hit_v,assocMCPart,Proton_v,purity,completeness);

    if (btpart >= 0) {
      
      auto backtrackedPart = Proton_v.at(btpart);
      
      if (backtrackedPart.tids.size() > 0) {
	
	auto tid = backtrackedPart.tids.at(0);

	// is purity for this proton improved? if so update backtracked cluster
	if ( Proton2CLSMap[ MCParticleTID2IDXMap[ tid ] ].second < purity ) {
	  
	  Proton2CLSMap[ MCParticleTID2IDXMap[ tid ] ] = std::make_pair( 1000.*(backtrackedPart.e-0.938277) , purity );
	  
	  std::cout << "Matched proton with CLUSTER : TrackID " << tid << " and purity " << purity << std::endl;

	}
	
      }// if there is at least one TrackID associated to this MCParticle
      
    }// if we backtracked to something

  }// loop through all clusters in slice

  // loop through all protons
  for (auto const& mapelem : Proton2PFPMap) {

    // get PFP idx , pur pair
    auto PFPmatch = mapelem.second;
    // get CLS idx, pur pair
    auto CLSmatch = Proton2CLSMap[mapelem.first];

    std::cout << "Proton with energy " << PFPmatch.first
	      << " has PFP purity of " << PFPmatch.second
	      << " and CLS purity of " << CLSmatch.second << std::endl;

    _energy = PFPmatch.first;
    _pfppur = PFPmatch.second;
    _clspur = CLSmatch.second;

    _ismaxproton = 0;
    if (mapelem.first == maxprotonIDX)
      _ismaxproton = 1;
    
    _tree->Fill();

  }
  
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
