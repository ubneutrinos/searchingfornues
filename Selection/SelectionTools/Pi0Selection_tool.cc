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

    // obtain track fit dedx in first 4 cm from track calo
    void TrackFitdEdx(const searchingfornues::ProxyCaloElem_t& trk,
		      float& dedxU, float& dedxV ,float& dedxY);

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
    float _vx1, _vy1, _vz1; // momentum of first shower
    float _px1, _py1, _pz1; // direction of first shower
    float _vx2, _vy2, _vz2; // momentum of second shower
    float _px2, _py2, _pz2; // direction of second shower
    float _radlen1, _radlen2;
    float _dot1, _dot2;
    float _energy1_Y, _energy2_Y;
    float _dedx1_Y, _dedx2_Y;
    float _dedx1_fit_Y, _dedx2_fit_Y;
    float _energy1_V, _energy2_V;
    float _dedx1_V, _dedx2_V;
    float _dedx1_fit_V, _dedx2_fit_V;
    float _energy1_U, _energy2_U;
    float _dedx1_U, _dedx2_U;
    float _dedx1_fit_U, _dedx2_fit_U;
    float _shrscore1, _shrscore2;
    float _gammadot;
    float _mass_Y, _mass_V, _mass_U;
    float _rc_vtx_x, _rc_vtx_y, _rc_vtx_z; // reco neutrino vertex
    

    // module-specific settings
    bool _onlyshower;   // should we use only showers to reconstruct pi0s?
    float _dmin;        // what is the minimum distance of the trk/shr vertex to the neutrino vertex?
    float _dotmin;      // maximum dot product between shower direction and vtx->start vector
    float _trkshrscore; // score on which to cut for track/shower classification

    art::InputTag fTRKproducer;
    art::InputTag fCALproducer;
    
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

    fTRKproducer = pset.get< art::InputTag > ("TRKproducer", "");
    fCALproducer = pset.get< art::InputTag > ("CALproducer", "");

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

    searchingfornues::ProxyCaloColl_t const* tkcalo_proxy = NULL;
    if (fTRKproducer!="") {
      tkcalo_proxy = new searchingfornues::ProxyCaloColl_t( proxy::getCollection<std::vector<recob::Track> >(e,fTRKproducer,proxy::withAssociated<anab::Calorimetry>(fCALproducer)) );
    }
    
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
	//auto ntrk = pfp_pxy.get<recob::Track>().size();

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
    
    if (tkcalo_proxy!=NULL) {
      for (const searchingfornues::ProxyCaloElem_t& tk : *tkcalo_proxy) {
	// find track with ID matching the pfp index (this convention apparently works only for shower fits...)
	if (tk->ID()==int(pfp_pxy_v[i1].index())) {
	  TrackFitdEdx(tk, _dedx1_fit_U, _dedx1_fit_V, _dedx1_fit_Y);
	}// if track matches shower index -> this is the track-fitted to the shower
      }// for all track fits to showers
    }// if track-fits to showers exist
    
    _radlen1    = vtxcompat1.second;
    _dot1       = vtxcompat1.first;
    _energy1_Y  = shr1->Energy()[2];
    _dedx1_Y    = shr1->dEdx()[2];
    _energy1_V  = shr1->Energy()[1];
    _dedx1_V    = shr1->dEdx()[1];
    _energy1_U  = shr1->Energy()[0];
    _dedx1_U    = shr1->dEdx()[0];
    _vx1        = shr1->ShowerStart()[0];
    _vy1        = shr1->ShowerStart()[1];
    _vz1        = shr1->ShowerStart()[2];
    _px1        = shr1->Direction()[0];
    _py1        = shr1->Direction()[1];
    _pz1        = shr1->Direction()[2];
    
    auto vtxcompat2 = VtxCompatibility(nuvtx, shr2->ShowerStart(), shr2->Direction());

    if (tkcalo_proxy!=NULL) {
      for (const searchingfornues::ProxyCaloElem_t& tk : *tkcalo_proxy) {
	// find track with ID matching the pfp index (this convention apparently works only for shower fits...)
	if (tk->ID()==int(pfp_pxy_v[i2].index())) {
	  TrackFitdEdx(tk, _dedx2_fit_U, _dedx2_fit_V, _dedx2_fit_Y);
	}// if track matches shower index -> this is the track-fitted to the shower
      }// for all track fits to showers
    }// if track-fits to showers exist

    _radlen2   = vtxcompat2.second;
    _dot2      = vtxcompat2.first;
    _energy2_Y = shr2->Energy()[2];
    _dedx2_Y   = shr2->dEdx()[2];
    _energy2_V = shr2->Energy()[1];
    _dedx2_V   = shr2->dEdx()[1];
    _energy2_U = shr2->Energy()[0];
    _dedx2_U   = shr2->dEdx()[0];
    _vx2        = shr2->ShowerStart()[0];
    _vy2        = shr2->ShowerStart()[1];
    _vz2        = shr2->ShowerStart()[2];
    _px2        = shr2->Direction()[0];
    _py2        = shr2->Direction()[1];
    _pz2        = shr2->Direction()[2];
    
    _gammadot = shr1->Direction().Dot(shr2->Direction());
    _mass_Y = sqrt( 2 * _energy1_Y * _energy2_Y * (1 - _gammadot ) );
    _mass_V = sqrt( 2 * _energy1_V * _energy2_V * (1 - _gammadot ) );
    _mass_U = sqrt( 2 * _energy1_U * _energy2_U * (1 - _gammadot ) );

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
    _tree->Branch("mcgamma0_e" ,&_mcgamma0_e ,"mcgamma0_e/F" );
    _tree->Branch("mcgamma0_px",&_mcgamma0_px,"mcgamma0_px/F");
    _tree->Branch("mcgamma0_py",&_mcgamma0_py,"mcgamma0_py/F");
    _tree->Branch("mcgamma0_pz",&_mcgamma0_pz,"mcgamma0_pz/F");
    _tree->Branch("mcrcdot0",&_mcrcdot0,"mcrcdot0/F");
    _tree->Branch("mcrce0",&_mcrce0,"mcrce0/F");
    _tree->Branch("mcgamma1_e" ,&_mcgamma1_e ,"mcgamma1_e/F" );
    _tree->Branch("mcgamma1_px",&_mcgamma1_px,"mcgamma1_px/F");
    _tree->Branch("mcgamma1_py",&_mcgamma1_py,"mcgamma1_py/F");
    _tree->Branch("mcgamma1_pz",&_mcgamma1_pz,"mcgamma1_pz/F");
    _tree->Branch("mcrcdot1",&_mcrcdot1,"mcrcdot1/F");
    _tree->Branch("mcrce1",&_mcrce1,"mcrce1/F");

    // reco
    _tree->Branch("nshower",&_nshower,"nshower/I");
    _tree->Branch("ntrack" ,&_ntrack ,"ntrack/I" );
    _tree->Branch("ngamma",&_ngamma,"ngamma/I");
    _tree->Branch("radlen1",&_radlen1,"radlen1/F");
    _tree->Branch("radlen2",&_radlen2,"radlen2/F");
    _tree->Branch("dot1",&_dot1,"dot1/F");
    _tree->Branch("dot2",&_dot2,"dot2/F");
    _tree->Branch("energy1_Y",&_energy1_Y,"energy1_Y/F");
    _tree->Branch("energy2_Y",&_energy2_Y,"energy2_Y/F");
    _tree->Branch("dedx1_Y",&_dedx1_Y,"dedx1_Y/F");
    _tree->Branch("dedx2_Y",&_dedx2_Y,"dedx2_Y/F");
    _tree->Branch("dedx1_fit_Y",&_dedx1_fit_Y,"dedx1_fit_Y/F");
    _tree->Branch("dedx2_fit_Y",&_dedx2_fit_Y,"dedx2_fit_Y/F");
    _tree->Branch("energy1_V",&_energy1_V,"energy1_V/F");
    _tree->Branch("energy2_V",&_energy2_V,"energy2_V/F");
    _tree->Branch("dedx1_V",&_dedx1_V,"dedx1_V/F");
    _tree->Branch("dedx2_V",&_dedx2_V,"dedx2_V/F");
    _tree->Branch("dedx1_fit_V",&_dedx1_fit_V,"dedx1_fit_V/F");
    _tree->Branch("dedx2_fit_V",&_dedx2_fit_V,"dedx2_fit_V/F");
    _tree->Branch("energy1_U",&_energy1_U,"energy1_U/F");
    _tree->Branch("energy2_U",&_energy2_U,"energy2_U/F");
    _tree->Branch("dedx1_U",&_dedx1_U,"dedx1_U/F");
    _tree->Branch("dedx2_U",&_dedx2_U,"dedx2_U/F");
    _tree->Branch("dedx1_fit_U",&_dedx1_fit_U,"dedx1_fit_U/F");
    _tree->Branch("dedx2_fit_U",&_dedx2_fit_U,"dedx2_fit_U/F");
    _tree->Branch("shrscore1",&_shrscore1,"shrscore1/F");
    _tree->Branch("shrscore2",&_shrscore2,"shrscore2/F");
    _tree->Branch("gammadot",&_gammadot,"gammadot/F");
    _tree->Branch("mass_Y",&_mass_Y,"mass_Y/F");
    _tree->Branch("mass_V",&_mass_V,"mass_V/F");
    _tree->Branch("mass_U",&_mass_U,"mass_U/F");
    _tree->Branch("rc_vtx_x",&_rc_vtx_x,"rc_vtx_x/F");
    _tree->Branch("rc_vtx_y",&_rc_vtx_y,"rc_vtx_y/F");
    _tree->Branch("rc_vtx_z",&_rc_vtx_z,"rc_vtx_z/F");
    _tree->Branch("vx1",&_vx1,"vx1/F");
    _tree->Branch("vy1",&_vy1,"vy1/F");
    _tree->Branch("vz1",&_vz1,"vz1/F");
    _tree->Branch("px1",&_px1,"px1/F");
    _tree->Branch("py1",&_py1,"py1/F");
    _tree->Branch("pz1",&_pz1,"pz1/F");
    _tree->Branch("vx2",&_vx2,"vx2/F");
    _tree->Branch("vy2",&_vy2,"vy2/F");
    _tree->Branch("vz2",&_vz2,"vz2/F");
    _tree->Branch("px2",&_px2,"px2/F");
    _tree->Branch("py2",&_py2,"py2/F");
    _tree->Branch("pz2",&_pz2,"pz2/F");

    return;
  }

  void Pi0Selection::resetTTree(TTree* _tree) {

    _nshower = std::numeric_limits<int>::min();
    _ntrack  = std::numeric_limits<int>::min();

    return;
  }

  std::pair<double,double> Pi0Selection::VtxCompatibility(const TVector3& nuvtx, const TVector3& shrvtx, const TVector3& shrdir) {

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
    
    
  }// end of vertex compatibility

  template <typename T> float Pi0Selection::PFPEnergy(const T& ass_clus_v) {

    float energy0 = 0;
    float energy1 = 0;
    float energy2 = 0;
    
    for (auto ass_clus : ass_clus_v) {

      //auto clus = clus_v[ass_clus.key()];
      
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
	
	
	if (ngamma == 0) {
	  _mcgamma0_e = mcp.E(0);
	  _mcgamma0_mom = mcp.Momentum(0).Vect().Unit();
	  _mcgamma0_px = _mcgamma0_mom.X();
	  _mcgamma0_py = _mcgamma0_mom.Y();
	  _mcgamma0_pz = _mcgamma0_mom.Z();
	}
	
	if (ngamma == 1) {
	  _mcgamma1_e = mcp.E(0);
	  _mcgamma1_mom = mcp.Momentum(0).Vect().Unit();
	  _mcgamma1_px = _mcgamma1_mom.X();
	  _mcgamma1_py = _mcgamma1_mom.Y();
	  _mcgamma1_pz = _mcgamma1_mom.Z();
	}
	
	ngamma += 1; // update photon count
	
	
      }// if a pi0 induced gamma
      
    }// for all mcparticles
    

  return;
  }

  void Pi0Selection::TrackFitdEdx(const searchingfornues::ProxyCaloElem_t& trk,
				  float& dedxU, float& dedxV ,float& dedxY)  {

    dedxU = -1;
    dedxV = -1;
    dedxY = -1;

    auto const trkcalos = trk.get<anab::Calorimetry>();
    
    for (const auto& tkcalo : trkcalos) {
      if (tkcalo->ResidualRange().size()==0) continue;
      std::vector<float> dedx4cm;
      for (size_t ic=0; ic<tkcalo->ResidualRange().size(); ++ic) {
	if ( (tkcalo->ResidualRange().back()-tkcalo->ResidualRange()[ic]) < 4.) {
	  //dedx4cm.push_back( tkcalo->dEdx()[ic] );
	  dedx4cm.push_back( tkcalo->dQdx()[ic] * 240. * (23.6/1e6) / 0.60 );
	}
      }
      float dedx4cm_med = -1.;
      if (dedx4cm.size()>0) {
	std::sort(dedx4cm.begin(), dedx4cm.end());
	if (dedx4cm.size()%2 == 1) dedx4cm_med = dedx4cm[dedx4cm.size()/2];
	else dedx4cm_med = 0.5*(dedx4cm[dedx4cm.size()/2] + dedx4cm[dedx4cm.size()/2 - 1]);
      }
      
      auto pl = tkcalo->PlaneID().Plane;
      if (pl == 0) { dedxU = dedx4cm_med; }
      if (pl == 1) { dedxV = dedx4cm_med; }
      if (pl == 2) { dedxY = dedx4cm_med; }
      
    }// for all calorimetry objects associated to the track
    
    return;
  }// TrackFitdEdx


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
    _energy1_Y = -1;
    _energy2_Y = -1;
    _dedx1_Y = -1;
    _dedx2_Y = -1;
    _energy1_V = -1;
    _energy2_V = -1;
    _dedx1_V = -1;
    _dedx2_V = -1;
    _energy1_U = -1;
    _energy2_U = -1;
    _dedx1_U = -1;
    _dedx2_U = -1;
    _ngamma = 0;
    _ntrack = 0;
    _nshower = 0;
    _shrscore1 = -1;
    _shrscore2 = -1;
    _gammadot = -1;
    _mass_Y = -1;
    _mass_V = -1;
    _mass_U = -1;

    _vx1 = _vx2 = _vy1 = _vy2 = _vz1 = _vz2 = 0;
    _px1 = _px2 = _py1 = _py2 = _pz1 = _pz2 = 0;
    
    return;
  }// end of reset function

  
  DEFINE_ART_CLASS_TOOL(Pi0Selection)
} // namespace selection

#endif
