#ifndef SELECTION_SELECTIONEXAMPLE_CXX
#define SELECTION_SELECTIONEXAMPLE_CXX

#include <iostream>
#include "SelectionToolBase.h"

namespace selection
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       ShowerSelection
    // File:        ShowerSelection.cc
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
    
  class ShowerSelection : public SelectionToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    ShowerSelection(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~ShowerSelection(){};
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset);
    
    /**
     * @brief Selection function
     */
    bool selectEvent(art::Event const& e,
		     const std::vector<ProxyPfpElem_t>& pfp_pxy_v);

    /**
     * @brief set branches for TTree
     */
    void setBranches(TTree* _tree);

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree* _tree);
    
  private:

    /**
     * @brief calculate PFP energy based on hits associated to clusters
     */
    template <typename T> float PFPEnergy(const T& ass_clus_v);

    // TTree variables
    int _nshower;
    int _ntrack;
    float _emin, _emax;
    float _rc_vtx_x, _rc_vtx_y, _rc_vtx_z; // reco neutrino vertex
    
  };

  template <typename T> float ShowerSelection::PFPEnergy(const T& ass_clus_v) {

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
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  ShowerSelection::ShowerSelection(const fhicl::ParameterSet& pset)
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
  void ShowerSelection::configure(fhicl::ParameterSet const & pset)
  {}
  
  //----------------------------------------------------------------------------
  /// selectEvent
  ///
  /// Arguments:
  ///
  /// art::Event
  /// slice track pointer vector
  /// slice shower pointer vector
  ///
  bool ShowerSelection::selectEvent(art::Event const& e,
				    const std::vector<ProxyPfpElem_t>& pfp_pxy_v)
  {
    
    _nshower = 0;
    _ntrack  = 0;

    TVector3 nuvtx;
    Double_t xyz[3] = {};

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

      else {
	
	auto nshr = pfp_pxy.get<recob::Shower>().size();
	auto ntrk = pfp_pxy.get<recob::Track>().size();
	
	_nshower += nshr;
	_ntrack  += ntrk;
	
	// grab cluster associated
	if (nshr == 1) {
	  auto ass_clus_v =  pfp_pxy.get<recob::Cluster>();// clus_pxy_v[pfp_pxy.key()];
	  float energy = PFPEnergy(ass_clus_v);
	  if (energy > _emax) { _emax = energy; }
	  if (energy < _emin) { _emin = energy; }
	}// if pfp is shower-like
	
      }// if not neutrino

    }// for all PFP
      
    if ( _nshower >= 1)
      return true;
    
    return false;
  }
  
  void ShowerSelection::setBranches(TTree* _tree) {
    
    _tree->Branch("_nshower",&_nshower,"nshower/I");
    _tree->Branch("_ntrack" ,&_ntrack ,"ntrack/I" );
    _tree->Branch("_emin",&_emin,"emin/F");
    _tree->Branch("_emax",&_emax,"emax/F");
    _tree->Branch("_rc_vtx_x",&_rc_vtx_x,"rc_vtx_x/F");
    _tree->Branch("_rc_vtx_y",&_rc_vtx_y,"rc_vtx_y/F");
    _tree->Branch("_rc_vtx_z",&_rc_vtx_z,"rc_vtx_z/F");

    return;
  }

  void ShowerSelection::resetTTree(TTree* _tree) {

    _emin    = 1e6;
    _emax    = 0.;

    _nshower = std::numeric_limits<int>::min();
    _ntrack  = std::numeric_limits<int>::min();

    _rc_vtx_x = std::numeric_limits<int>::min();
    _rc_vtx_y = std::numeric_limits<int>::min();
    _rc_vtx_z = std::numeric_limits<int>::min();

    return;
  }
  
  
  DEFINE_ART_CLASS_TOOL(ShowerSelection)
} // namespace selection

#endif
