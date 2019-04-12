#ifndef ANALYSIS_SHOWERANALYSIS_CXX
#define ANALYSIS_SHOWERANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       ShowerAnalysis
    // File:        ShowerAnalysis.cc
    //
    //              A basic analysis example
    //
    // Configuration parameters:
    //
    // TBD
    //
    // Created by Giuseppe Cerati (cerati@fnal.gov) on 03/15/2019
    //
    ////////////////////////////////////////////////////////////////////////

  class ShowerAnalysis : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    ShowerAnalysis(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~ShowerAnalysis(){ };
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset);
    
    /**
     * @brief Analysis function
     */
    void analyzeEvent(art::Event const& e, bool fData) override;

    /**
     * @brief Analyze slice
     */
    void analyzeSlice(art::Event const& e, std::vector<ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) override;

    /**
     * @brief Save truth info for event associated to neutrino
     */
    void SaveTruth(art::Event const& e);
    
    /**
     * @brief set branches for TTree
     */
    void setBranches(TTree* _tree) override;

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree* _tree) override;
    
  private:


  };

  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  ShowerAnalysis::ShowerAnalysis(const fhicl::ParameterSet& p)
  {

  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void ShowerAnalysis::configure(fhicl::ParameterSet const & p)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void ShowerAnalysis::analyzeEvent(art::Event const& e, bool fData)
  {

  }

  void ShowerAnalysis::resetTTree(TTree* _tree)
  {

  }

  DEFINE_ART_CLASS_TOOL(ShowerAnalysis)
} // namespace analysis

#endif
