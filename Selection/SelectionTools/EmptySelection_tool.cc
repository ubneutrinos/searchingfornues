#ifndef SELECTION_EMPTYSELECTION_CXX
#define SELECTION_EMPTYSELECTION_CXX

#include <iostream>
#include "SelectionToolBase.h"

namespace selection
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       EmptySelection
    // File:        EmptySelection.cc
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
    
  class EmptySelection : public SelectionToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    EmptySelection(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~EmptySelection(){};
    
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
    void setBranches(TTree* _tree){};

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree* _tree){};
    
  private:
    
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  EmptySelection::EmptySelection(const fhicl::ParameterSet& pset)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void EmptySelection::configure(fhicl::ParameterSet const & pset)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  bool EmptySelection::selectEvent(art::Event const& e,
				     const std::vector<ProxyPfpElem_t>& pfp_pxy_v)
  {
    
    return true;
  }
  
  
  DEFINE_ART_CLASS_TOOL(EmptySelection)
} // namespace selection

#endif
