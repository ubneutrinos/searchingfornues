////////////////////////////////////////////////////////////////////////
// Class:       ProtonHitPurity
// Plugin Type: analyzer (art v3_01_02)
// File:        ProtonHitPurity_module.cc
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

class ProtonHitPurity;


class ProtonHitPurity : public art::EDAnalyzer {
public:
  explicit ProtonHitPurity(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ProtonHitPurity(ProtonHitPurity const&) = delete;
  ProtonHitPurity(ProtonHitPurity&&) = delete;
  ProtonHitPurity& operator=(ProtonHitPurity const&) = delete;
  ProtonHitPurity& operator=(ProtonHitPurity&&) = delete;

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
  int _nhit;
  int _plane;
  float _charge;
  float _purity;
  float _shrdot;
  float _shrdist;
  float _shrenergy;
  int _nshr;

  std::vector<art::Ptr<recob::Hit>> getGaussHits(const std::vector<art::Ptr<recob::Hit>> &hits,
						 const art::ValidHandle<std::vector<recob::Hit> > gaushit_h);

  bool IsProtonIsolated(const std::vector<art::Ptr<recob::Hit>> &hits,
			const art::ValidHandle<std::vector<recob::Hit> > gaushit_h);

  void ProtonDot(const float& protonWire, const float& protonTime,const int& pl,
		 const TVector3& showerVtx, const TVector3& showerDir,
		 float &dot, float& d2d);
  
};


ProtonHitPurity::ProtonHitPurity(fhicl::ParameterSet const& p)
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
  _tree->Branch("nhit",&_nhit,"nhit/I");
  _tree->Branch("nshr",&_nshr,"nshr/I");
  _tree->Branch("plane",&_plane,"plane/I");
  _tree->Branch("charge",&_charge,"charge/F");
  _tree->Branch("purity",&_purity,"purity/F");
  _tree->Branch("shrdot",&_shrdot,"shrdot/F");
  _tree->Branch("shrdist",&_shrdist,"shrdist/F");
  _tree->Branch("shrenergy",&_shrenergy,"shrenergy/F");

  // get detector specific properties
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _wire2cm = geom->WirePitch(0,0,0);
  _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
  
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}
  
void ProtonHitPurity::analyze(art::Event const& e)
{

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
  
  // build PFParticle map  for this event
  BuildPFPMap(pfp_proxy);

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
  
  TVector3 ShowerVtx, ShowerDir;
  _shrenergy = 0;
  
  // go through slice and find shower-like particles
  for (size_t p=0; p < slice_pfp_v.size(); p++) {
    
    auto pfp_pxy = slice_pfp_v[p];
    
    auto ass_shr_v = pfp_pxy.get<recob::Shower>();
    
    if (ass_shr_v.size() != 1) continue;
    
    auto trkshrscore = searchingfornues::GetTrackShowerScore(pfp_pxy);
    
    if (trkshrscore > 0.5) continue;
    
    auto shr = ass_shr_v[0];
    
    if (shr->Energy()[2] < _shrenergy) continue;
    
    _shrenergy = shr->Energy()[2];
    
    ShowerVtx = shr->ShowerStart();
    ShowerDir = shr->Direction();
    
  }// for all proxys in slice
  
  // vector of BtPart for backtracking
  std::vector<searchingfornues::BtPart> Proton_v;
  
  // loop through MCParticles and identify protons
  for ( unsigned int im=0; im< mcparticle_h->size(); ++im ) {
    const auto& mcp = mcparticle_h->at(im);
    if (mcp.PdgCode() == 2212) 
      
      Proton_v.push_back( searchingfornues::BtPart(mcp.PdgCode(),
						   mcp.Momentum(0).Px() * 0.001,
						   mcp.Momentum(0).Py() * 0.001,
						   mcp.Momentum(0).Pz() * 0.001,
						   mcp.Momentum(0).E() * 0.001, 
						   mcp.TrackId()) );
    
  }// for all MCParticles
  
  std::cout << "PROTON there are " << Proton_v.size() << " protons in the event" << std::endl;
  
  // loop through clusters
  for (size_t c=0; c < cluster_h->size(); c++) {
    
    auto clus = cluster_h->at(c);
    
    // get associated hits
    auto clus_hit_v = clus_hit_assn_v.at( c );
    
    float completeness = 0;
    
    // is the proton isolated?
    if (IsProtonIsolated(clus_hit_v,hit_h) == false) continue;
    
    // create vector of gaushits corresponding to new proton hits
    auto gaushit_hit_v = getGaussHits(clus_hit_v, hit_h);
    
    //if ((gaushit_hit_v.at(0))->WireID().Plane != 2) continue;
    
    _plane = (gaushit_hit_v.at(0))->WireID().Plane;
    _nhit  = gaushit_hit_v.size();
    _charge = 0.;
    _purity = 0.;
    
    _shrdot = 1e6;
    _shrdist = 1e6;
    
    float protonWire = (gaushit_hit_v.at(0))->WireID().Wire * _wire2cm;
    float protonTime = (gaushit_hit_v.at(0))->PeakTime()     * _time2cm;
    
    if (_shrenergy > 0)
      ProtonDot(protonWire, protonTime, _plane, ShowerVtx, ShowerDir, _shrdot, _shrdist);
    
    for (size_t ph=0; ph < gaushit_hit_v.size(); ph++) 
      _charge += (gaushit_hit_v.at(ph))->Integral();
    
    getAssocBtPart(gaushit_hit_v,assocMCPart,Proton_v,_purity,completeness);
    
    _tree->Fill();
    
    std::cout << "PROTON Found " << gaushit_hit_v.size() << " matching original " << clus_hit_v.size() << " proton hits" << std::endl;
    std::cout << "PROTON This cluster has " << gaushit_hit_v.size() << " hits and a purity/completenss of [ "
	      << _purity << ", " << completeness << " ]" << std::endl;
    std::cout << "PROTON" << std::endl;
    
  }// for all MCparticles
  
  return;
}
  
  
void ProtonHitPurity::ProtonDot(const float& protonWire, const float& protonTime,const int& pl,
				const TVector3& showerVtx, const TVector3& showerDir,
				float &dot, float& d2d) {
  
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  
  auto Vtxwire = geom->WireCoordinate(showerVtx[1],showerVtx[2],geo::PlaneID(0,0,pl)) * _wire2cm;
  auto Vtxtime = showerVtx[0];
  
  auto Dirwire = geom->WireCoordinate(showerDir[1],showerDir[2],geo::PlaneID(0,0,pl)) * _wire2cm;
  auto Dirtime = showerDir[0];
  
  TVector3 showerDir2D(Dirwire,Dirtime,0.);
  TVector3 protonDir2D(Vtxwire-protonWire,Vtxtime-protonTime,0.);
  
  d2d = sqrt( ((protonWire - Vtxwire) * (protonWire - Vtxwire)) +
	      ((protonTime - Vtxtime) * (protonTime - Vtxtime)) );
  
  dot = showerDir2D.Dot(protonDir2D);
  dot /= showerDir2D.Mag();
  dot /= protonDir2D.Mag();
  
  return;
}

std::vector<art::Ptr<recob::Hit>> ProtonHitPurity::getGaussHits(const std::vector<art::Ptr<recob::Hit>> &hits,
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

bool ProtonHitPurity::IsProtonIsolated(const std::vector<art::Ptr<recob::Hit>> &hits,
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

 void ProtonHitPurity::BuildPFPMap(const searchingfornues::ProxyPfpColl_t& pfp_pxy_col) {
   
   _pfpmap.clear();
   
   unsigned int p=0;
   for (const auto& pfp_pxy : pfp_pxy_col) {
     _pfpmap[pfp_pxy->Self()] = p;
     p++;
   }
   
   return;
 }// BuildPFPMap
 
 void ProtonHitPurity::AddDaughters(const searchingfornues::ProxyPfpElem_t& pfp_pxy,
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
 

void ProtonHitPurity::beginJob()
{
  // Implementation of optional member function here.
}

void ProtonHitPurity::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ProtonHitPurity)
