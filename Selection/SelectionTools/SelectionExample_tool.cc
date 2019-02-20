#ifndef SELECTION_SELECTIONEXAMPLE_CXX
#define SELECTION_SELECTIONEXAMPLE_CXX

#include <iostream>
#include "SelectionToolBase.h"

namespace selection
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       SelectionExample
    // File:        SelectionExample.cc
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
    
  class SelectionExample : public SelectionToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    SelectionExample(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~SelectionExample(){};
    
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

    // multiplicity requirement on slice
    unsigned int _multiplicity;
    
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  SelectionExample::SelectionExample(const fhicl::ParameterSet& pset)
  {
    _multiplicity = 0;
    configure(pset);
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void SelectionExample::configure(fhicl::ParameterSet const & pset)
  {
    _multiplicity = pset.get<unsigned int>("multiplicity");
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  bool SelectionExample::selectEvent(art::Event const& e,
				     const std::vector<ProxyPfpElem_t>& pfp_pxy_v)
  {
    
    if ( pfp_pxy_v.size() >= _multiplicity )
      return true;
    
    return false;
  }
  
  
  DEFINE_ART_CLASS_TOOL(SelectionExample)
} // namespace selection

#endif
