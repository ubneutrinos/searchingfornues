#ifndef ANALYSISTOOLBASE_H
#define ANALYSISTOOLBASE_H

// art TOOLS
#include "ubana/ubana/searchingfornues/Selection/SelectionTools/SelectionToolBase.h"

#include "TTree.h"
#include <limits>

namespace analysis {

  using ProxyPfpElem_t = selection::ProxyPfpElem_t;
  using ProxyClusColl_t = selection::ProxyClusColl_t;
  
class AnalysisToolBase {

public:

    /**
     *  @brief  Virtual Destructor
     */
    virtual ~AnalysisToolBase() noexcept = default;
    
    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    void configure(const fhicl::ParameterSet&){};

    /**
     * @brief Analysis function
     *
     * @param art::Event event record for analysis
     */
    virtual void analyzeEvent(art::Event const& e, bool fData) = 0;

    /**
     * @brief Analysis function
     *
     * @param art::Event event record for analysis
     */
    virtual void analyzeSlice(art::Event const& e, std::vector<ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) = 0;

    /**
     * @brief set branches for TTree
     */
    virtual void setBranches(TTree* _tree) = 0;

    
    /**
     * @brief resetset TTree branches
     */
    virtual void resetTTree(TTree* _tree) = 0;


};

} // analysis namespace

#endif
