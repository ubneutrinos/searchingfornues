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
		     const std::vector<art::Ptr<recob::Track>  >& trkptr_v,
		     const std::vector<art::Ptr<recob::Shower> >& shrptr_v);

    /**
     * @brief set branches for TTree
     */
    void setBranches(TTree* _tree);

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree* _tree);
    
  private:

    // TTree variables
    int _nshower;
    int _ntrack;
    
  };
  
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
				     const std::vector<art::Ptr<recob::Track>  >& trkptr_v,
				     const std::vector<art::Ptr<recob::Shower> >& shrptr_v)
  {
    
    _nshower = shrptr_v.size();
    _ntrack  = trkptr_v.size();

    if ( _nshower >= 1)
      return true;
    
    return false;
  }

  void ShowerSelection::setBranches(TTree* _tree) {

    _tree->Branch("_nshower",&_nshower,"nshower/I");
    _tree->Branch("_ntrack" ,&_ntrack ,"ntrack/I" );

    return;
  }

  void ShowerSelection::resetTTree(TTree* _tree) {

    _nshower = std::numeric_limits<int>::min();
    _ntrack  = std::numeric_limits<int>::min();

    return;
  }
  
  
  DEFINE_ART_CLASS_TOOL(ShowerSelection)
} // namespace selection

#endif
