#ifndef SELECTION_SELECTIONEXAMPLE_CXX
#define SELECTION_SELECTIONEXAMPLE_CXX

#include <iostream>
#include "SelectionToolBase.h"

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

    std::pair<double,double> VtxCompatibility(const TVector3& nuvtx, const ProxyPfpElem_t& pfp_pxy, TVector3& gammadir);

    /**
     * @brief calculate PFP energy based on hits associated to clusters
     */
    template <typename T> float PFPEnergy(const T& ass_clus_v);

    void Reset();

    // TTree variables
    int _nshower;
    int _ntrack;
    int _ngamma;
    float _radlen1, _radlen2;
    float _dot1, _dot2;
    float _energy1, _energy2;
    int _isshower1, _isshower2;
    float _gammadot;
    float _mass;
    float _rc_vtx_x, _rc_vtx_y, _rc_vtx_z; // reco neutrino vertex
    

    // module-specific settings
    bool _onlyshower; // should we use only showers to reconstruct pi0s?
    float _dmin;      // what is the minimum distance of the trk/shr vertex to the neutrino vertex?
    float _dotmin;    // maximum dot product between shower direction and vtx->start vector
    
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
    _onlyshower = pset.get< bool >  ("onlyshower");
    _dotmin     = pset.get< float > ("dotmin"    );
    _dmin       = pset.get< float > ("dmin"      );
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

    TVector3 gammadir1, gammadir2;

    // loop over all PFPs and get metadata, find the neutrino and save its vertex
    for (const auto& pfp_pxy : pfp_pxy_v) {

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

	// grab cluster associated
	auto ass_clus_v =  pfp_pxy.get<recob::Cluster>();// clus_pxy_v[pfp_pxy.key()];
	float energy = PFPEnergy(ass_clus_v);

	auto nshr = pfp_pxy.get<recob::Shower>().size();
	auto ntrk = pfp_pxy.get<recob::Track>().size();

	_nshower += nshr;
	_ntrack  += ntrk;
	
	if ( (nshr != 1) && (_onlyshower) ) continue;
	
	TVector3 gammadir;
	auto vtxcompat = VtxCompatibility(nuvtx,pfp_pxy, gammadir);
	
	// if blank result, continue
	if ( (vtxcompat.first == -1) && (vtxcompat.second == -1) ) continue;
	
	if ( (vtxcompat.second < _dmin) || (vtxcompat.first < _dotmin) ) continue;
	
	if (_ngamma == 0) {
	  if (nshr == 1) _isshower1 = 1;
	  _radlen1 = vtxcompat.second;
	  _dot1    = vtxcompat.first;
	  _energy1 = energy;
	  gammadir1 = gammadir;
	}
	
	if (_ngamma == 1) {
	  if (nshr == 1) _isshower2 = 1;
	  _radlen2 = vtxcompat.second;
	  _dot2    = vtxcompat.first;
	  _energy2 = energy;
	  gammadir2 = gammadir;
	}
	
	_ngamma += 1;

      }// if not neutrino
      
    }// loop over all PFParticles

    _gammadot = gammadir1.Dot(gammadir2);
    _mass = sqrt( 2 * _energy1 * _energy2 * (1 - _gammadot ) );

    return true;
  }
  
  
  void Pi0Selection::setBranches(TTree* _tree) {

    _tree->Branch("_nshower",&_nshower,"nshower/I");
    _tree->Branch("_ntrack" ,&_ntrack ,"ntrack/I" );
    _tree->Branch("_ngamma",&_ngamma,"ngamma/I");
    _tree->Branch("_radlen1",&_radlen1,"radlen1/F");
    _tree->Branch("_radlen2",&_radlen2,"radlen2/F");
    _tree->Branch("_dot1",&_dot1,"dot1/F");
    _tree->Branch("_dot2",&_dot2,"dot2/F");
    _tree->Branch("_energy1",&_energy1,"energy1/F");
    _tree->Branch("_energy2",&_energy2,"energy2/F");
    _tree->Branch("_isshower1",&_isshower1,"isshower1/I");
    _tree->Branch("_isshower2",&_isshower2,"isshower2/I");
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

  std::pair<double,double> Pi0Selection::VtxCompatibility(const TVector3& nuvtx, const ProxyPfpElem_t& pfp_pxy, TVector3& gammadir) {

    auto shr_v = pfp_pxy.get<recob::Shower>();
    auto trk_v = pfp_pxy.get<recob::Track>();

    if ( (shr_v.size() + trk_v.size()) != 1) {
      std::cout << "\t there are " << shr_v.size() << " showers associated to this PFP" << std::endl;
      std::cout << "\t there are " << trk_v.size() << " tracks  associated to this PFP" << std::endl;
      std::cout << "ERROR. PFP associated with != (shr+trk)." << std::endl;
      return std::make_pair(-1,-1);
    }

    if (shr_v.size() == 1) {
      
      auto shr = shr_v.at(0);
      
      // grab shower start point and direction
      auto shrvtx = shr->ShowerStart();
      auto shrdir = shr->Direction();
      gammadir = shrdir;
      
      // assess compatibility
      auto nuvtx2shrvtx = (shrvtx - nuvtx).Unit();
      auto shrdirnormed = shrdir.Unit();
      
      double dot  = nuvtx2shrvtx.Dot(shrdirnormed);
      double dist = (nuvtx-shrvtx).Mag();
      
      return std::make_pair(dot,dist);
      
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

  void Pi0Selection::Reset() {

    _dot1 = -1;
    _dot2 = -1;
    _radlen1 = -1;
    _radlen2 = -1;
    _energy1 = -1;
    _energy2 = -1;
    _ngamma = 0;
    _ntrack = 0;
    _nshower = 0;
    _isshower1 = 0;
    _isshower2 = 0;
    _gammadot = -1;
    _mass = -1;
    
    return;
  }// end of reset function

  
  DEFINE_ART_CLASS_TOOL(Pi0Selection)
} // namespace selection

#endif
