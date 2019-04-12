#ifndef SELECTION_SELECTIONEXAMPLE_CXX
#define SELECTION_SELECTIONEXAMPLE_CXX

#include <iostream>
#include "SelectionToolBase.h"

#include "nusimdata/SimulationBase/MCTruth.h"

namespace selection
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       Pi0Selection
    // File:        Pi0Selection.cc
    //
    //              A basic selection example
    //
    // Configuration parameters:
    //
    // TBD
    //
    // Created by David Caratelli (davidc@fnal.gov) on 01/30/2019
    //
    ////////////////////////////////////////////////////////////////////////
    
  class Pi0Selection : public SelectionToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    Pi0Selection(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~Pi0Selection(){};
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset);
    
    /**
     * @brief Selection function
     */
    bool selectEvent(art::Event const& e,
		     const std::vector<ProxyPfpElem_t>& pfp_pxy_v);

    /**
cioe'     * @brief set branches for TTree
     */
    void setBranches(TTree* _tree);

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree* _tree);
    
  private:

    std::pair<double,double> VtxCompatibility(const TVector3& nuvtx, const TVector3& shrvtx, const TVector3& shrdir);

    // get shower score
    float GetTrackShowerScore(const ProxyPfpElem_t& pfp_pxy);

    /**
     * @brief calculate PFP energy based on hits associated to clusters
     */
    template <typename T> float PFPEnergy(const T& ass_clus_v);

    /**
     * @brief read truth info associated to pi0
     */
    void ReadTruth(art::Event const& e);

    void Reset();

    // TTree variables
    
    // truth variables
    int _ispi0;
    TVector3 _mcgamma0_mom, _mcgamma1_mom;
    float _mcgamma0_e, _mcgamma0_px, _mcgamma0_py, _mcgamma0_pz;
    float _mcgamma1_e, _mcgamma1_px, _mcgamma1_py, _mcgamma1_pz;
    float _mcrcdot0, _mcrcdot1; // dot product between MC gamma and RC gamma. Each MC gamma matched to best RC gamma
    float _mcrce0, _mcrce1; // energy of reconstructed photon associated to the MC one

    // reco variables
    int _nshower;
    int _ntrack;
    int _ngamma;
    float _radlen1, _radlen2;
    float _dot1, _dot2;
    float _energy1, _energy2;
    float _dedx1, _dedx2;
    float _shrscore1, _shrscore2;
    float _gammadot;
    float _mass;
    float _rc_vtx_x, _rc_vtx_y, _rc_vtx_z; // reco neutrino vertex
    

    // module-specific settings
    bool _onlyshower;   // should we use only showers to reconstruct pi0s?
    float _dmin;        // what is the minimum distance of the trk/shr vertex to the neutrino vertex?
    float _dotmin;      // maximum dot product between shower direction and vtx->start vector
    float _trkshrscore; // score on which to cut for track/shower classification
    
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  Pi0Selection::Pi0Selection(const fhicl::ParameterSet& pset)
  {
    configure(pset);
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void Pi0Selection::configure(fhicl::ParameterSet const & pset)
  {
    _onlyshower  = pset.get< bool >  ("onlyshower" );
    _dotmin      = pset.get< float > ("dotmin"     );
    _dmin        = pset.get< float > ("dmin"       );
    _trkshrscore = pset.get< float > ("trkshrscore");
}
  
  //----------------------------------------------------------------------------
  /// selectEvent
  ///
  /// Arguments:
  ///
  /// art::Event
  /// slice track pointer vector
  /// slice shower pointer vector
  ///
  bool Pi0Selection::selectEvent(art::Event const& e,
				 const std::vector<ProxyPfpElem_t>& pfp_pxy_v)
  {

    TVector3 nuvtx;
    Double_t xyz[3] = {};

    Reset();
    
    if (!fData)
      ReadTruth(e);

    TVector3 gammadir1, gammadir2;

    // vector of pfp indices for showers
    std::vector<float> shr_energy_v;
    // vector of shower energies
    std::vector<short> pfp_idx_v;

    short pfp_ctr = 0;

    // loop over all PFPs and get metadata, find the neutrino and save its vertex
    for (const auto& pfp_pxy : pfp_pxy_v) {

      pfp_ctr += 1;

      auto PDG = fabs(pfp_pxy->PdgCode());

      if ( (PDG == 12) || (PDG == 14) ) {

	// grab vertex
	auto vtx = pfp_pxy.get<recob::Vertex>();
	if (vtx.size() != 1) {
	  std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
	  return false;
	}
	
	// save vertex to array
	vtx.at(0)->XYZ(xyz);
	nuvtx = TVector3(xyz[0],xyz[1],xyz[2]);

	_rc_vtx_x = nuvtx.X();
	_rc_vtx_y = nuvtx.Y();
	_rc_vtx_z = nuvtx.Z();

      }// if neutrino PFP

      else { // if not neutrino PFP

	// grab shower/track score
	auto trkshrscore = GetTrackShowerScore(pfp_pxy);

	auto nshr = pfp_pxy.get<recob::Shower>().size();
	auto ntrk = pfp_pxy.get<recob::Track>().size();

	std::cout << "DAVIDC : there are " << nshr << " showers and " << ntrk << " tracks associated to this PFP w/ trackScore " << trkshrscore << std::endl;

	// 1 -> track-like
	if (trkshrscore > _trkshrscore)  continue;

	if (nshr != 1) continue;

	auto const& shr = pfp_pxy.get<recob::Shower>().at(0);

	auto energy = shr->Energy()[2];
	
	auto vtxcompat = VtxCompatibility(nuvtx, shr->ShowerStart(), shr->Direction());

	// if blank result, continue
	if ( (vtxcompat.first == -1) && (vtxcompat.second == -1) ) continue;
	// if too close to vertex or too mis-algined, continue
	if ( (vtxcompat.second < _dmin) || (vtxcompat.first < _dotmin) ) continue;

	shr_energy_v.push_back( energy );
	pfp_idx_v.push_back( pfp_ctr - 1 );

      }// if not the neutrino PFP

    }// for all PFP

    // if we did not find any shower, return
    if (shr_energy_v.size() == 0) return false;

    // if a single shower, backtrack and quit
    if (shr_energy_v.size() == 1) {

      auto shr = pfp_pxy_v.at(pfp_idx_v[0]).get<recob::Shower>().at(0);
      
      float dot = shr->Direction().Dot(_mcgamma0_mom);
      if (dot > _mcrcdot0) {
	_mcrcdot0 = dot;
	_mcrce0   = shr->Energy()[2];
      }// if new most aligned shower
      dot = shr->Direction().Dot(_mcgamma1_mom);
      if (dot > _mcrcdot1) {
	_mcrcdot1 = dot;
	_mcrce1   = shr->Energy()[2];
      }// if new most aligned shower

      return false;
    }// if a single shower

    // if two or more, sort by energy
    double e1 = 0; // energy of highest energy reco shower
    short  i1 = 0; // pfp index for highest energy reco shower
    double e2 = 0; // energy of 2nd highest energy reco shower
    short  i2 = 0; // pfp index for 2nd highest energy reco shower

    // find highest energy shower
    for (size_t n0=0; n0 < shr_energy_v.size(); n0++) {
      if (shr_energy_v[n0] > e1) {
	e1 = shr_energy_v[n0];
	i1 = pfp_idx_v[n0];
      } 
    }
    // find second highest energy shower
    for (size_t n0=0; n0 < shr_energy_v.size(); n0++) {
      if (pfp_idx_v[n0] == i1) continue;
      if (shr_energy_v[n0] > e2) {
	e2 = shr_energy_v[n0];
	i2 = pfp_idx_v[n0];
      }
    }

    auto shr1 = pfp_pxy_v.at(i1).get<recob::Shower>().at(0);
    auto shr2 = pfp_pxy_v.at(i2).get<recob::Shower>().at(0);

    _shrscore1 = GetTrackShowerScore(pfp_pxy_v.at(i1));
    _shrscore2 = GetTrackShowerScore(pfp_pxy_v.at(i2));

    auto vtxcompat1 = VtxCompatibility(nuvtx, shr1->ShowerStart(), shr1->Direction());

    _radlen1  = vtxcompat1.second;
    _dot1     = vtxcompat1.first;
    _energy1  = shr1->Energy()[2];
    _dedx1    = shr1->dEdx()[2];
    
    auto vtxcompat2 = VtxCompatibility(nuvtx, shr2->ShowerStart(), shr2->Direction());

    _radlen2 = vtxcompat2.second;
    _dot2    = vtxcompat2.first;
    _energy2 = shr2->Energy()[2];
    _dedx2   = shr2->dEdx()[2];
    
    _gammadot = shr1->Direction().Dot(shr2->Direction());
    _mass = sqrt( 2 * _energy1 * _energy2 * (1 - _gammadot ) );

    // backtracking for shower 1
    float dot1 = shr1->Direction().Dot(_mcgamma0_mom);
    if (dot1 > _mcrcdot0) {
      _mcrcdot0 = dot1;
      _mcrce0   = shr1->Energy()[2];
    }// if new most aligned shower
    dot1 = shr1->Direction().Dot(_mcgamma1_mom);
    if (dot1 > _mcrcdot1) {
      _mcrcdot1 = dot1;
      _mcrce1   = shr1->Energy()[2];
    }// if new most aligned shower

    // backtracking for shower 1
    float dot2 = shr2->Direction().Dot(_mcgamma0_mom);
    if (dot2 > _mcrcdot0) {
      _mcrcdot0 = dot2;
      _mcrce0   = shr2->Energy()[2];
    }// if new most aligned shower
    dot2 = shr2->Direction().Dot(_mcgamma1_mom);
    if (dot2 > _mcrcdot1) {
      _mcrcdot1 = dot2;
      _mcrce1   = shr2->Energy()[2];
    }// if new most aligned shower
    
    return true;
  }
  
  
  void Pi0Selection::setBranches(TTree* _tree) {

    // truth
    _tree->Branch("_mcgamma0_e" ,&_mcgamma0_e ,"mcgamma0_e/F" );
    _tree->Branch("_mcgamma0_px",&_mcgamma0_px,"mcgamma0_px/F");
    _tree->Branch("_mcgamma0_py",&_mcgamma0_py,"mcgamma0_py/F");
    _tree->Branch("_mcgamma0_pz",&_mcgamma0_pz,"mcgamma0_pz/F");
    _tree->Branch("_mcrcdot0",&_mcrcdot0,"mcrcdot0/F");
    _tree->Branch("_mcrce0",&_mcrce0,"mcrce0/F");
    _tree->Branch("_mcgamma1_e" ,&_mcgamma1_e ,"mcgamma1_e/F" );
    _tree->Branch("_mcgamma1_px",&_mcgamma1_px,"mcgamma1_px/F");
    _tree->Branch("_mcgamma1_py",&_mcgamma1_py,"mcgamma1_py/F");
    _tree->Branch("_mcgamma1_pz",&_mcgamma1_pz,"mcgamma1_pz/F");
    _tree->Branch("_mcrcdot1",&_mcrcdot1,"mcrcdot1/F");
    _tree->Branch("_mcrce1",&_mcrce1,"mcrce1/F");

    // reco
    _tree->Branch("_nshower",&_nshower,"nshower/I");
    _tree->Branch("_ntrack" ,&_ntrack ,"ntrack/I" );
    _tree->Branch("_ngamma",&_ngamma,"ngamma/I");
    _tree->Branch("_radlen1",&_radlen1,"radlen1/F");
    _tree->Branch("_radlen2",&_radlen2,"radlen2/F");
    _tree->Branch("_dot1",&_dot1,"dot1/F");
    _tree->Branch("_dot2",&_dot2,"dot2/F");
    _tree->Branch("_energy1",&_energy1,"energy1/F");
    _tree->Branch("_energy2",&_energy2,"energy2/F");
    _tree->Branch("_dedx1",&_dedx1,"dedx1/F");
    _tree->Branch("_dedx2",&_dedx2,"dedx2/F");
    _tree->Branch("_shrscore1",&_shrscore1,"shrscore1/F");
    _tree->Branch("_shrscore2",&_shrscore2,"shrscore2/F");
    _tree->Branch("_gammadot",&_gammadot,"gammadot/F");
    _tree->Branch("_mass",&_mass,"mass/F");
    _tree->Branch("_rc_vtx_x",&_rc_vtx_x,"rc_vtx_x/F");
    _tree->Branch("_rc_vtx_y",&_rc_vtx_y,"rc_vtx_y/F");
    _tree->Branch("_rc_vtx_z",&_rc_vtx_z,"rc_vtx_z/F");

    return;
  }

  void Pi0Selection::resetTTree(TTree* _tree) {

    _nshower = std::numeric_limits<int>::min();
    _ntrack  = std::numeric_limits<int>::min();

    return;
  }

  std::pair<double,double> Pi0Selection::VtxCompatibility(const TVector3& nuvtx, const TVector3& shrvtx, const TVector3& shrdir) {

    /*
      auto shr_v = pfp_pxy.get<recob::Shower>();
      auto trk_v = pfp_pxy.get<recob::Track>(;)
      
      if ( (shr_v.size() + trk_v.size()) != 1) {
      std::cout << "\t there are " << shr_v.size() << " showers associated to this PFP" << std::endl;
      std::cout << "\t there are " << trk_v.size() << " tracks  associated to this PFP" << std::endl;
      std::cout << "ERROR. PFP associated with != (shr+trk)." << std::endl;
      return std::make_pair(-1,-1);
      }
      
      if (shr_v.size() == 1) {
      
      auto shr = shr_v.at(0);
    */
    
    // grab shower start point and direction
    //auto shrvtx = shr->ShowerStart();
    //auto shrdir = shr->Direction();
    //gammadir = shrdir;
    
    // assess compatibility
    auto nuvtx2shrvtx = (shrvtx - nuvtx).Unit();
    auto shrdirnormed = shrdir.Unit();
    
    double dot  = nuvtx2shrvtx.Dot(shrdirnormed);
    double dist = (nuvtx-shrvtx).Mag();
    
    return std::make_pair(dot,dist);
    
    /*
      }// associated with a shower
      
      if (trk_v.size() == 1) {
      
      auto trk = trk_v.at(0);
      
      auto trkvtx = trk->Vertex();
      TVector3 trkvtx3(trkvtx.x(), trkvtx.y(), trkvtx.z());
      auto trkdir = trk->VertexDirection();
      TVector3 trkdir3(trkdir.x(), trkdir.y(), trkdir.z());
      gammadir = trkdir3;
      
      auto nuvtx2trkvtx = (trkvtx3 - nuvtx).Unit();
      auto trkdirnormed = trkdir3.Unit();
      
      double dot  = nuvtx2trkvtx.Dot(trkdirnormed);
      double dist = (nuvtx-trkvtx3).Mag();
      
      return std::make_pair(dot,dist);
      
      }// associated with a track
      
      return std::make_pair(-1,-1);
    */
    
  }// end of vertex compatibility

  template <typename T> float Pi0Selection::PFPEnergy(const T& ass_clus_v) {

    float energy0 = 0;
    float energy1 = 0;
    float energy2 = 0;
    
    for (auto ass_clus : ass_clus_v) {

      //auto clus = clus_v[ass_clus.key()];
      
      std::cout << "cluster integral : " << ass_clus->Integral() << std::endl;

      if (ass_clus->View() == 0) { energy0 = ass_clus->Integral(); }
      if (ass_clus->View() == 1) { energy1 = ass_clus->Integral(); }
      if (ass_clus->View() == 2) { energy2 = ass_clus->Integral(); }

    }// for all clusters

    if (energy2 != 0) return energy2;
    return (energy1+energy0)/2.;
  }// calculte PFP energy based on associated hit charge

  void Pi0Selection::ReadTruth(art::Event const& e) {

    auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");
    auto const& mcp_h = e.getValidHandle<std::vector<simb::MCParticle> >("largeant");
    
    _ispi0 = 0;
    
    auto mct = mct_h->at(0);
    size_t npart = mct.NParticles();
    
    int pi0trkid = -1;
    double pi0_vtx_x, pi0_vtx_y, pi0_vtx_z;
    
    for (size_t i=0; i < npart; i++){
      auto const& part = mct.GetParticle(i);
      if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
	pi0_vtx_x = part.Trajectory().X(0);
	pi0_vtx_y = part.Trajectory().Y(0);
	pi0_vtx_z = part.Trajectory().Z(0);
	pi0trkid = part.TrackId();
	std::cout << "DAVIDC found a pi0! " << std::endl;
	break;
      }
    }
    
    if (pi0trkid < 0) return;
    
    // how many photons from the pi0 did we find?
    int ngamma = 0;
    
    // go through MCParticles and find photons
    for (size_t i=0; i < mcp_h->size(); i++) {
      
      auto mcp = mcp_h->at(i);
      
      if (mcp.PdgCode() != 22) continue;
      
      double x = mcp.Trajectory().X(0);
      double y = mcp.Trajectory().Y(0);
      double z = mcp.Trajectory().Z(0);
      double d = sqrt( ( (pi0_vtx_x - x) * (pi0_vtx_x - x) ) +
		       ( (pi0_vtx_y - y) * (pi0_vtx_y - y) ) +
		       ( (pi0_vtx_z - z) * (pi0_vtx_z - z) ) );
      
      if ( d < 0.01 ) {
	
	std::cout << "DAVIDC ngamma is " << ngamma << std::endl;
	
	if (ngamma == 0) {
	  std::cout << "DAVIDC gamma0 " << std::endl;
	  _mcgamma0_e = mcp.E(0);
	  _mcgamma0_mom = mcp.Momentum(0).Vect().Unit();
	  _mcgamma0_px = _mcgamma0_mom.X();
	  _mcgamma0_py = _mcgamma0_mom.Y();
	  _mcgamma0_pz = _mcgamma0_mom.Z();
	}
	
	if (ngamma == 1) {
	  std::cout << "DAVIDC gamma1 " << std::endl;
	  _mcgamma1_e = mcp.E(0);
	  _mcgamma1_mom = mcp.Momentum(0).Vect().Unit();
	  _mcgamma1_px = _mcgamma1_mom.X();
	  _mcgamma1_py = _mcgamma1_mom.Y();
	  _mcgamma1_pz = _mcgamma1_mom.Z();
	}
	
	ngamma += 1; // update photon count
	
	std::cout << "DAVIDC found " << mcp.PdgCode() << " with pi0 parent " << std::endl;
	std::cout << "DAVIDC found " << mcp.PdgCode() << " with pi0 parent " << std::endl;
	
      }// if a pi0 induced gamma
      
    }// for all mcparticles
    
    /*
      if ( (foundShowers == true) && fPDG == 111) {
      
    size_t idx_1 = 0;
    size_t idx_2 = 0;
    size_t n_found = 0;
    for (size_t i=0; i < mcs_h->size(); i++){
      auto const& mcs = mcs_h->at(i);
      // distance from vertex                                                                
      double x = mcs.Start().X();
      double y = mcs.Start().Y();
      double z = mcs.Start().Z();
      double d = sqrt( ( (_mc_vtx_x - x) * (_mc_vtx_x - x) ) +
		       ( (_mc_vtx_y - y) * (_mc_vtx_y - y) ) +
		       ( (_mc_vtx_z - z) * (_mc_vtx_z - z) ) );
      if ( d < 0.01 ){
	if (n_found == 0){
	  idx_1 = i;
	  n_found += 1;
	}
	else if (n_found == 1){
	  idx_2 = i;
	  n_found += 1;
	}
	else
	  n_found += 1;
      }// if mother is a Pi0   
    }// for all MCShowers                                                               
    
    
    size_t idxLARGE = idx_1;
    size_t idxSMALL = idx_2;
    
    if (mcs_h->at(idx_1).Start().E() < mcs_h->at(idx_2).Start().E() )
      { idxLARGE = idx_2; idxSMALL = idx_1; }
    
    auto const& mcshr1 = mcs_h->at(idxLARGE);
    auto const& mcshr2 = mcs_h->at(idxSMALL);
    
    _pi0_e  = mcshr1.MotherEnd().E();
    
    _mc_shr1_e  = mcshr1.Start().E();
    _mc_shr1_edep  = mcshr1.DetProfile().E();
    _mc_shr1_x  = mcshr1.DetProfile().X();
    _mc_shr1_y  = mcshr1.DetProfile().Y();
    _mc_shr1_z  = mcshr1.DetProfile().Z();

    double mom1 = mcshr1.Start().Momentum().Vect().Mag();    
    _mc_shr1_px = mcshr1.Start().Px() / mom1;
    _mc_shr1_py = mcshr1.Start().Py() / mom1;
    _mc_shr1_pz = mcshr1.Start().Pz() / mom1;
    
    _mc_shr2_e  = mcshr2.Start().E();
    _mc_shr2_edep  = mcshr2.DetProfile().E();
    _mc_shr2_x  = mcshr2.DetProfile().X();
    _mc_shr2_y  = mcshr2.DetProfile().Y();
    _mc_shr2_z  = mcshr2.DetProfile().Z();

    double mom2 = mcshr2.Start().Momentum().Vect().Mag();    
    _mc_shr2_px = mcshr2.Start().Px() / mom2;
    _mc_shr2_py = mcshr2.Start().Py() / mom2;
    _mc_shr2_pz = mcshr2.Start().Pz() / mom2;
  */

  return;
  }

  float Pi0Selection::GetTrackShowerScore(const ProxyPfpElem_t& pfp_pxy) {

    const auto& pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();

    if (pfParticleMetadataList.size() == 0) 
      return 1;

    for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j) {
      
      const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
      auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
      if (!pfParticlePropertiesMap.empty()) {
	for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	  if (it->first == "TrackScore")
	    return it->second;
	}// for map elements
      }// if pfp metadata map not empty
    }// for list

    return 1;
  }

  void Pi0Selection::Reset() {

    _ispi0 = 0;
    _mcgamma0_e = 0;
    _mcgamma0_px = 0;
    _mcgamma0_py = 0;
    _mcgamma0_pz = 0;
    _mcgamma1_e = 0;
    _mcgamma1_px = 0;
    _mcgamma1_py = 0;
    _mcgamma1_pz = 0;
    _mcrce0 = 0;
    _mcrce1 = 0;
    _mcrcdot0 = -1;
    _mcrcdot1 = -1;

    _dot1 = -1;
    _dot2 = -1;
    _radlen1 = -1;
    _radlen2 = -1;
    _energy1 = -1;
    _energy2 = -1;
    _dedx1 = -1;
    _dedx2 = -1;
    _ngamma = 0;
    _ntrack = 0;
    _nshower = 0;
    _shrscore1 = -1;
    _shrscore2 = -1;
    _gammadot = -1;
    _mass = -1;
    
    return;
  }// end of reset function

  
  DEFINE_ART_CLASS_TOOL(Pi0Selection)
} // namespace selection

#endif
