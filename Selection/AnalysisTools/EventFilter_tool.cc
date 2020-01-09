#ifndef SELECTION_EVENTFILTER_CXX
#define SELECTION_EVENTFILTER_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/TrackFitterFunctions.h"
#include "../CommonDefs/SCECorrections.h"

#include "canvas/Persistency/Common/TriggerResults.h" 
#include "fhiclcpp/ParameterSetRegistry.h" 

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       EventFilter
    // File:        EventFilter.cc
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
    
  class EventFilter : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    EventFilter(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~EventFilter(){};
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset);
    
    /**
     * @brief Analysis function
     */
    void analyzeEvent(art::Event const& e, bool fData) override;

    /**
     * @brief Analyze slice
     */
    void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

    /**
     * @brief set branches for TTree
     */
    void setBranches(TTree* _tree);

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree* _tree);
    
  private:

    void Reset();

    // TTree variables
    
    int _filter_antibdt;
    int _filter_ccinclusive;
    int _filter_ncpi0;
    int _filter_pi0;

  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  EventFilter::EventFilter(const fhicl::ParameterSet& pset)
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
  void EventFilter::configure(fhicl::ParameterSet const & pset)
  {

}

  
void EventFilter::analyzeEvent(art::Event const &e, bool fData)
{

  Reset();

  // filters stored only for overlay events
  if (!fData) {
    
    art::InputTag trigResInputTag("TriggerResults","","OverlayFiltersPostStage2"); // the last is the name of process where the filters were run
    art::ValidHandle<art::TriggerResults> trigRes = e.getValidHandle<art::TriggerResults>(trigResInputTag);
    fhicl::ParameterSet pset;
    if (!fhicl::ParameterSetRegistry::get(trigRes->parameterSetID(), pset)) { throw cet::exception("PSet Not Found???"); }
    std::vector<std::string> trigger_path_names = pset.get<std::vector<std::string> >("trigger_paths", {});
    if (trigger_path_names.size()!=trigRes->size()) { throw cet::exception("Size mismatch???"); }
    for (size_t itp=0;itp<trigRes->size();itp++) {
      //
      if (trigger_path_names.at(itp)=="NuCC")    { _filter_ccinclusive = trigRes->at(itp).accept(); }
      if (trigger_path_names.at(itp)=="antibdt") { _filter_antibdt     = trigRes->at(itp).accept(); }
      if (trigger_path_names.at(itp)=="ncpi0")   { _filter_ncpi0       = trigRes->at(itp).accept(); }
      if (trigger_path_names.at(itp)=="pi0")     { _filter_pi0         = trigRes->at(itp).accept(); }
    }
  
  }// if not data
  
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
  void EventFilter::analyzeSlice(art::Event const& e, std::vector<ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected)
  {

    return;
  }
  
  
  void EventFilter::setBranches(TTree* _tree) {

    _tree->Branch("filter_antibdt",&_filter_antibdt,"filter_antibdt/I");
    _tree->Branch("filter_ncpi0",&_filter_ncpi0,"filter_ncpi0/I");
    _tree->Branch("filter_pi0",&_filter_pi0,"filter_pi0/I");
    _tree->Branch("filter_ccinclusive",&_filter_ccinclusive,"filter_ccinclusive/I");

    return;
  }

  void EventFilter::resetTTree(TTree* _tree) {

    return;
  }

  void EventFilter::Reset() {

    _filter_antibdt = 0;
    _filter_ccinclusive = 0;
    _filter_ncpi0 = 0;
    _filter_pi0 = 0;
    
    return;
  }// end of reset function

  
  DEFINE_ART_CLASS_TOOL(EventFilter)
} // namespace selection

#endif
