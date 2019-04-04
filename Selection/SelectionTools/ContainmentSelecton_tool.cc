#ifndef SELECTION_CONTAINMENTSELECTION_CXX
#define SELECTION_CONTAINMENTSELECTION_CXX

#include <iostream>
#include "SelectionToolBase.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

namespace selection
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       ContainmentSelection
    // File:        ContainmentSelection.cc
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
    
  class ContainmentSelection : public SelectionToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    ContainmentSelection(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~ContainmentSelection(){};
    
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

    bool IsFiducial(float x, float y, float z);

    float _FV; // FV boundary to apply

    // TTree variables
    float _rc_vtx_x, _rc_vtx_y, _rc_vtx_z;
    
  };

  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  ContainmentSelection::ContainmentSelection(const fhicl::ParameterSet& pset)
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
  void ContainmentSelection::configure(fhicl::ParameterSet const & pset)
  {

    _FV = pset.get<float>("FV");
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
  bool ContainmentSelection::selectEvent(art::Event const& e,
					 const std::vector<ProxyPfpElem_t>& pfp_pxy_v)
  {
    
    TVector3 nuvtx;
    Double_t xyz[3] = {};
    
    std::cout << "DAVIDC Containment" << std::endl;
    

    
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
	
	
	if (IsFiducial(nuvtx.X(), nuvtx.Y(), nuvtx.Z()) == false)
	  return false;
	
	    
      }// if neutrino PFP
	
      else { // if not the neutrino PFP
	
	auto ntrk = pfp_pxy.get<recob::Track>().size();
	
	if (ntrk == 1) {
	  
	  auto trk = pfp_pxy.get<recob::Track>().at(0);
	  
	  auto trkstart = trk->Vertex();
	  auto trkend   = trk->End();
	  
	  if (IsFiducial(trkstart.X(), trkstart.Y(), trkstart.Z() ) == false)
	    return false;
	  
	  if (IsFiducial(trkend.X(), trkend.Y(), trkend.Z() ) == false)
	    return false;
	  
	}// if associated to a track
	
      }// if not the neurino PFP
      
    }// for all PFP
    
    // made it this far, no vertex or track start/end point is out of the FV
    return true;
  }
  
  void ContainmentSelection::setBranches(TTree* _tree) {
    
    _tree->Branch("_rc_vtx_x",&_rc_vtx_x,"rc_vtx_x/F");
    _tree->Branch("_rc_vtx_y",&_rc_vtx_y,"rc_vtx_y/F");
    _tree->Branch("_rc_vtx_z",&_rc_vtx_z,"rc_vtx_z/F");

    return;
  }

  void ContainmentSelection::resetTTree(TTree* _tree) {

    _rc_vtx_x = std::numeric_limits<int>::min();
    _rc_vtx_y = std::numeric_limits<int>::min();
    _rc_vtx_z = std::numeric_limits<int>::min();

    return;
  }

  bool ContainmentSelection::IsFiducial(float x, float y, float z) {

    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    if (SCE->EnableCalSpatialSCE() == true) {
      
      auto offset = SCE->GetPosOffsets(geo::Point_t(x,y,z));
      x += offset.X();
      y -= offset.Y();
      z -= offset.Z();
      
    }// if spatial offset calibrations are enabled

    // are we within the FV?
    if (x < (0.    + _FV)) return false;
    if (x > (256.  - _FV)) return false;
    if (y < (-116. + _FV)) return false;
    if (y > (116.  - _FV)) return false;
    if (z < (0.    + _FV)) return false;
    if (z > (1036. - _FV)) return false;

    return true;
  }
  
  
  DEFINE_ART_CLASS_TOOL(ContainmentSelection)
} // namespace selection

#endif
