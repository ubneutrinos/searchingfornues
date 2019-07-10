#ifndef ANALYSIS_COSMICIP_CXX
#define ANALYSIS_COSMICIP_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       FlashMatching
    // File:        FlashMatching.cc
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

  class FlashMatching : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    FlashMatching(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~FlashMatching(){ };
    
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
     * @brief set branches for TTree
     */
    void setBranches(TTree* _tree) override;

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree* _tree) override;
    
  private:

    art::InputTag fPFPproducer;
    art::InputTag fT0producer;

    float _nu_flashmatch_score;
    float _best_cosmic_flashmatch_score;
    float _best_obviouscosmic_flashmatch_score;
    std::vector<float> _cosmic_flashmatch_score_v;
    
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  FlashMatching::FlashMatching(const fhicl::ParameterSet& p)
  {

    fPFPproducer = p.get< art::InputTag >("PFPproducer");
    fT0producer  = p.get< art::InputTag >("T0producer" );
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void FlashMatching::configure(fhicl::ParameterSet const & p)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void FlashMatching::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
  {

    // load tracks previously created for which T0 reconstruction is requested                                                                                                                                
    art::Handle<std::vector<anab::T0> > t0_h;
    e.getByLabel( fT0producer , t0_h );

    // make sure tracks look good                                                                                                                                                                             
    if(!t0_h.isValid()) {
      std::cerr<<"\033[93m[WARNING]\033[00m no anab::T0 for flash-matching. Skip flash-matching"<<std::endl;
      return;
    }



  art::ValidHandle<std::vector<recob::PFParticle>> pfp_h = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
  // grab PFP -> T0 flash-matching association for the event
  art::FindManyP< anab::T0 > pfp_t0_assn_v(pfp_h, e, fT0producer);
  // grab associated metadata
  art::FindManyP< larpandoraobj::PFParticleMetadata > pfp_meta_assn_v(pfp_h, e, fPFPproducer);

  // figure out which PFP is the neutrino
  size_t nupfp = 0;
  for (auto pfp : slice_pfp_v) {
    
    if (pfp->IsPrimary() == false) continue;
    
    if ( (pfp->PdgCode() == 12) || (pfp->PdgCode() == 14) ) {

      nupfp = pfp->Self();
      
    }// if neutrino
  }// for all PFPs in slice

  _nu_flashmatch_score = 1e6;
  _best_cosmic_flashmatch_score = 1e6;
  _best_obviouscosmic_flashmatch_score = 1e6;
  _cosmic_flashmatch_score_v.clear();

  // loop through all PFParticles
  for (size_t p=0; p < pfp_h->size(); p++) {

    auto const& pfp = pfp_h->at(p);

    // only primary PFPs have a flash-match score
    if (pfp.IsPrimary() == false) continue;

    // get flash-match score
    if (pfp_t0_assn_v.size() <= p) { std::cout << "NO T0!" << std::endl; continue; }
    if (pfp_t0_assn_v.at(p).size() != 1) { std::cout << "NO T0!" << std::endl; continue; }

    auto fmscore = pfp_t0_assn_v.at(p).at(0)->TriggerConfidence();

    // is this the neutrino?
    if (pfp.Self() == nupfp) { 
      _nu_flashmatch_score = fmscore;
      continue;
    }

    // if not the neutrino...
    _cosmic_flashmatch_score_v.push_back( fmscore );

    // get metadata
    if (pfp_meta_assn_v.size() <= p) { std::cout << "NO METADATA!" << std::endl; continue; }

    auto metadatalist = pfp_meta_assn_v.at(p);
    float nuscore = -1;

    bool clearcosmic = false;
    
    if (metadatalist.empty() == false) {
      
      for (unsigned int j=0; j<metadatalist.size(); ++j) {
	const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata(metadatalist.at(j));
	auto particleproperties = metadata->GetPropertiesMap();
	if (!particleproperties.empty())
	  for (std::map<std::string, float>::const_iterator it = particleproperties.begin(); it != particleproperties.end(); ++it) {
	    if ( it->first == "IsClearCosmic" ) 
	      clearcosmic = true;
	    if ( it->first == "NuScore" )
	      nuscore = it->second;
	  }// for all metadata items in the particle metadata
      }// for entries in list
      
    }// if there is metadata available

    if (fmscore < _best_cosmic_flashmatch_score) { _best_cosmic_flashmatch_score = fmscore; }
    if ( (clearcosmic == true) && (fmscore < _best_obviouscosmic_flashmatch_score) ) 
      _best_obviouscosmic_flashmatch_score = fmscore;
    
    std::cout << "Slice has flash-match score of " << fmscore << " with NuScore of " << nuscore << std::endl;


  }// for all PFPs
  
  return;
  }
  
  void FlashMatching::analyzeEvent(art::Event const &e, bool fData)
  {
    // std::cout << "analyze event" << std::endl;
  }

  void FlashMatching::setBranches(TTree* _tree)
  {
    _tree->Branch("nu_flashmatch_score",&_nu_flashmatch_score,"nu_flashmatch_score/F");
    _tree->Branch("best_cosmic_flashmatch_score",&_best_cosmic_flashmatch_score,"best_cosmic_flashmatch_score/F");
    _tree->Branch("best_obviouscosmic_flashmatch_score",&_best_obviouscosmic_flashmatch_score,"best_obviouscosmic_flashmatch_score/F");
    _tree->Branch("cosmic_flashmatch_score_v","std::vector<float>",&_cosmic_flashmatch_score_v);
  }
  
  void FlashMatching::resetTTree(TTree* _tree)
  {
    _nu_flashmatch_score = 1e6;
    _best_cosmic_flashmatch_score = 1e6;
    _best_obviouscosmic_flashmatch_score = 1e6;
    _cosmic_flashmatch_score_v.clear();
  }

  
  DEFINE_ART_CLASS_TOOL(FlashMatching)
} // namespace analysis

#endif
