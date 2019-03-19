#ifndef ANALYSIS_OBVIOUSCOSMICFLASHMATCHING_CXX
#define ANALYSIS_OBVIOUSCOSMICFLASHMATCHING_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"

#include "nusimdata/SimulationBase/MCTruth.h"

// backtracking tools
#include "ubana/ubana/searchingfornues/Selection/CommonDefs/BacktrackingFuncs.h"

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       ObviousCosmicFlashMatching
    // File:        ObviousCosmicFlashMatching.cc
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

  class ObviousCosmicFlashMatching : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    ObviousCosmicFlashMatching(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~ObviousCosmicFlashMatching(){ };
    
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

    art::InputTag fCRTVetoproducer; // producer for CRT veto ass tag [anab::T0 <-> recob::OpFlash]
    art::InputTag fCLSproducer; // cluster associated to PFP
    art::InputTag fMCTproducer;
    art::InputTag fBacktrackTag;

  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  ObviousCosmicFlashMatching::ObviousCosmicFlashMatching(const fhicl::ParameterSet& p)
  {
    fCRTVetoproducer = p.get< art::InputTag > ("CRTVetoproducer",""); // default is no CRT veto
    fCLSproducer = p.get< art::InputTag > ("CLSproducer");
    fMCTproducer = p.get< art::InputTag > ("MCTproducer");
    fBacktrackTag = p.get< art::InputTag > ("BacktrackTag");
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void ObviousCosmicFlashMatching::configure(fhicl::ParameterSet const & p)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void ObviousCosmicFlashMatching::analyzeEvent(art::Event const& e, bool fData)
  {

    /*
    // grab PFParticles in event
    auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);
    
    // grab spacepoints associated with PFParticles
    art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fPFPproducer);
    // grab tracks associated with PFParticles
    art::FindManyP<recob::Track> pfp_track_assn_v(pfp_h, e, fPFPproducer);
    
    // grab associated metadata
    art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfp_h, e, fPFPproducer);    
    
    auto const& spacepoint_h = e.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointproducer);
    
    art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, e, fSpacePointproducer);
    
    // ADDITION FROM PETRILLO
    //e.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointproducer);
    
    // grab the hits associated to the PFParticles
    //auto pfp_hit_assn_v = lar::FindManyInChainP<recob::Hit, recob::SpacePoint>::find(pfp_h, e, fPFPproducer);
    
    std::cout << "There are " << pfp_h->size() << " pfparticles in the event " << std::endl;
    
    _best_score = 1e5;
    _obvious = 0;
    
    // fill map: pfparticle Self() -> index/key
    _pfpmap.clear();
    for (unsigned int p=0; p < pfp_h->size(); p++)
      _pfpmap[pfp_h->at(p).Self()] = p;
    
    for (unsigned int p=0; p < pfp_h->size(); p++){
      
      auto const& pfp = pfp_h->at(p);
      
      // start from primary PFParticles
      if (pfp.IsPrimary() == false) continue;
      
      const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
      
      // now build vectors of PFParticles, space-points, and hits for this slice
      std::vector<recob::PFParticle> pfp_v;
      std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
      std::vector<std::vector<art::Ptr<recob::SpacePoint>>> spacepoint_v_v;
      std::vector<std::vector<art::Ptr<recob::Hit>>> hit_v_v;
      
      std::cout << "creating PFP hierarchy." << std::endl;
      AddDaughters(pfp_ptr, pfp_h, pfp_ptr_v);
      std::cout << "There are " << pfp_ptr_v.size() << " PFParticles in this hierarchy " << std::endl << std::endl;
      
      // go through these pfparticles and fill info needed for matching
      for (size_t i=0; i < pfp_ptr_v.size(); i++) {
	
	auto key = pfp_ptr_v.at(i).key();
	recob::PFParticle pfp = *pfp_ptr_v.at(i);
	
	pfp_v.push_back(pfp);
	
	//auto const& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
	const std::vector< art::Ptr<recob::SpacePoint> >& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
	
	std::vector< art::Ptr<recob::Hit> > hit_ptr_v;
	
	for (size_t sp=0; sp < spacepoint_ptr_v.size(); sp++) {
	  auto const& spkey = spacepoint_ptr_v.at(sp).key();
	  const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = spacepoint_hit_assn_v.at( spkey );
	  for (size_t h=0; h < this_hit_ptr_v.size(); h++) {
	    hit_ptr_v.push_back( this_hit_ptr_v.at( h ) );
	  }// for all hits associated to this spacepoint
	}// fpr all spacepoints
	
	spacepoint_v_v.push_back( spacepoint_ptr_v );
	hit_v_v.push_back( hit_ptr_v );
	
      }// for all pfp pointers
      
      // ignore events where the primary is not a muon [we are looking for obvious cosmics!]
      if (pfp_track_assn_v.at(p).size() != 1) continue;
      
      auto const trk = pfp_track_assn_v.at(p).at(0);
      
      // tracks too short should not be considered...
      if (trk->Length() < 10) continue;
      
      // do not attempt flash-matching with tracks not in time with the beam drift window
      if ( TrackInTime(trk) == false ) continue;
      
      // ready to call flash-matching
      _score = _flashmatchTool->ClassifySlice(e, pfp_ptr_v, spacepoint_v_v, hit_v_v, _peSpectrum, _peHypothesis);
      if (_score <= 0) continue;
      
      if (_score < _best_score) { 
	_best_score = _score; 
	// is this hierarchy an obvious cosmic?
	_obvious = 0;
	if ( ObviousCosmic(pfp_track_assn_v.at(p)) ) _obvious = 1;
      }// if the best score yet
      
    }// for all PFParticles
    */
    
    return;
  }

  void ObviousCosmicFlashMatching::analyzeSlice(art::Event const& e, std::vector<ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected)
  {

    return;
  }
  
  void ObviousCosmicFlashMatching::setBranches(TTree* _tree) 
  {
    //_tree->Branch("_evt"   ,&_evt   ,"evt/I"   );

  }

  void ObviousCosmicFlashMatching::resetTTree(TTree* _tree)
  {

  }

  DEFINE_ART_CLASS_TOOL(ObviousCosmicFlashMatching)
} // namespace analysis

#endif
