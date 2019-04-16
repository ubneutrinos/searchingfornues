#ifndef ANALYSIS_TRACKPIDANALYSIS_CXX
#define ANALYSIS_TRACKPIDANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

// backtracking tools
#include "ubana/ubana/searchingfornues/Selection/CommonDefs/Typedefs.h"
#include "ubana/ubana/searchingfornues/Selection/CommonDefs/BacktrackingFuncs.h"
#include "ubana/ubana/searchingfornues/Selection/CommonDefs/TrackShowerScoreFuncs.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       TrackPIDAnalysis
// File:        TrackPIDAnalysis.cc
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

class TrackPIDAnalysis : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  TrackPIDAnalysis(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~TrackPIDAnalysis(){};

  // provide for initialization
  void configure(fhicl::ParameterSet const &pset);

  /**
     * @brief Analysis function
     */
  void analyzeEvent(art::Event const &e, bool fData) override;

  /**
     * @brief Analyze slice
     */
  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

  /**
     * @brief Save truth info for event associated to neutrino
     */
  void SaveTruth(art::Event const &e);

  /**
     * @brief set branches for TTree
     */
  void setBranches(TTree *_tree) override;

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree) override;

private:

  art::InputTag fCALOproducer;
  art::InputTag fPIDproducer;

  std::vector<float> _trkscore_v;

};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
TrackPIDAnalysis::TrackPIDAnalysis(const fhicl::ParameterSet &p)
{

  fCALOproducer = p.get< art::InputTag > ("CALOproducer");
  fPIDproducer  = p.get< art::InputTag > ("PIDproducer" );

}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void TrackPIDAnalysis::configure(fhicl::ParameterSet const &p)
{
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void TrackPIDAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
}

void TrackPIDAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{

  searchingfornues::ProxyCaloColl_t const& calo_proxy = proxy::getCollection<std::vector<recob::Track> >(e, fCALOproducer,
													 proxy::withAssociated<anab::Calorimetry>(fCALOproducer));

  searchingfornues::ProxyPIDColl_t const& pid_proxy = proxy::getCollection<std::vector<recob::Track> >(e, fPIDproducer,
													 proxy::withAssociated<anab::ParticleID>(fPIDproducer));



  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
  {
    auto PDG = fabs(slice_pfp_v[i_pfp]->PdgCode());

    if (PDG == 12 || PDG == 14)
      continue;

    auto trk_v = slice_pfp_v[i_pfp].get<recob::Track>(); 

    auto ntrk = trk_v.size();

    if (ntrk != 1) continue;

    auto trk = trk_v.at(0);
    // get track proxy in order to fetch calorimtry
    auto trkpxy1 = calo_proxy[trk.key()];
    auto calopxy_v = trkpxy1.get<anab::Calorimetry>();
    std::cout << "There are "<< calopxy_v.size() << " associated calo objects " << std::endl;

    // get trk proxy in order to fetch PID
    auto trkpxy2 = pid_proxy[trk.key()];
    auto pidpxy_v = trkpxy2.get<anab::ParticleID>();
    std::cout << "There are "<< pidpxy_v.size() << " associated PID objects " << std::endl;

    _trkscore_v.push_back( searchingfornues::GetTrackShowerScore(slice_pfp_v[i_pfp]) );
    
  }// for all PFParticles
  
}
  
void TrackPIDAnalysis::setBranches(TTree *_tree)
{

  _tree->Branch("trkscore_v", "std::vector<float>",  &_trkscore_v);
}

void TrackPIDAnalysis::resetTTree(TTree *_tree)
{
  _trkscore_v.clear();
}

DEFINE_ART_CLASS_TOOL(TrackPIDAnalysis)
} // namespace analysis

#endif
