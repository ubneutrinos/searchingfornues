#ifndef ANALYSIS_OBVIOUSCOSMICFLASHMATCHING_CXX
#define ANALYSIS_OBVIOUSCOSMICFLASHMATCHING_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

// backtracking tools
#include "ubana/ubana/searchingfornues/Selection/CommonDefs/BacktrackingFuncs.h"
// flash-matching tools
#include "ubana/ubana/searchingfornues/FlashMatching/FlashMatchingToolBase_tool.h"

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

    void AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v);

    bool TrackInTime(const art::Ptr<recob::Track>& pfp_track_assn_v);

    art::InputTag fPFPproducer;
    art::InputTag fSpacePointproducer;

    float _obvious_flashmatch_score;
    float _neutrino_score;
    float _score;
    int   _obvious; // is the best score an obvious?
    int   _obvious_cosmics;
    float _obvious_trklen;
    float _obvious_startx, _obvious_starty, _obvious_startz;
    float _obvious_endx, _obvious_endy, _obvious_endz;

    std::vector<float> _peSpectrum, _peHypothesis, _peHypothesisNu, _peHypothesisCosmic;
    
    // PFP map
    std::map<unsigned int, unsigned int> _pfpmap;
    
    std::unique_ptr<flashmatch::FlashMatchingToolBase> _flashmatchTool;             ///< The slice id tool

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

    fPFPproducer         = p.get< art::InputTag >("PFPproducer");
    fSpacePointproducer  = p.get< art::InputTag >("SpacePointproducer");

    const fhicl::ParameterSet& flashconfig = p.get<fhicl::ParameterSet>("SliceTool");  
    
    _flashmatchTool = art::make_tool<flashmatch::FlashMatchingToolBase>(flashconfig);

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
    
    _obvious_flashmatch_score = 1e5;
    _obvious_cosmics = 0;
    _obvious_trklen  = 0;
    
    // fill map: pfparticle Self() -> index/key
    _pfpmap.clear();
    for (unsigned int p=0; p < pfp_h->size(); p++)
      _pfpmap[pfp_h->at(p).Self()] = p;
    
    for (unsigned int p=0; p < pfp_h->size(); p++){
      
      auto const& pfp = pfp_h->at(p);

      // start from primary PFParticles
      if (pfp.IsPrimary() == false) continue;

      bool clearCosmic = false;

      // get metadata for this PFP
      const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(p));

      if (pfParticleMetadataList.empty()) continue;
	
      for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j) {
	const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
	auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
	if (!pfParticlePropertiesMap.empty())
	  for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	    if ( it->first == "IsClearCosmic" ) {
	      clearCosmic = true;
	      break;
	    }
	  }// for all metadata items in the particle metadata
      }// for entries in list

      //if (clearCosmic == false) continue;
      
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

      // ready to call flash-matching
      _score = _flashmatchTool->ClassifySlice(e, pfp_ptr_v, spacepoint_v_v, hit_v_v, _peSpectrum, _peHypothesis);

      // now decide if neutrino or cosmic
      if (pfp.PdgCode() == 12 || pfp.PdgCode() == 14) {
	_neutrino_score = _score;
	_peHypothesisNu = _peHypothesis;
      }
      
      // ignore events where the primary is not a muon [we are looking for obvious cosmics!]
      if (pfp_track_assn_v.at(p).size() != 1) continue;
      
      auto const trk = pfp_track_assn_v.at(p).at(0);
      
      // tracks too short should not be considered...
      if (trk->Length() < 10) continue;
      
      // do not attempt flash-matching with tracks not in time with the beam drift window
      if ( TrackInTime(trk) == false ) continue;

      std::cout << "match score = " << _score << std::endl;
      if (_score <= 0) continue;
      
      if (_score < _obvious_flashmatch_score) { 
	if (clearCosmic) { 
	  _obvious = 1; 
	  _obvious_cosmics += 1;
	}
	else { _obvious  = 0; }
	_obvious_trklen = trk->Length();
	_obvious_flashmatch_score = _score;
	_peHypothesisCosmic = _peHypothesis;
	auto vtx = trk->Vertex();
	auto end = trk->End();
	_obvious_startx = vtx.X();
	_obvious_starty = vtx.Y();
	_obvious_startz = vtx.Z();
	_obvious_endx = end.X();
	_obvious_endy = end.Y();
	_obvious_endz = end.Z();
	// is this hierarchy an obvious cosmic?
      }// if the best score yet
      
    }// for all PFParticles
    
    return;
  }

  void ObviousCosmicFlashMatching::analyzeSlice(art::Event const& e, std::vector<ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected)
  {

    return;
  }
  
  void ObviousCosmicFlashMatching::setBranches(TTree* _tree) 
  {
    _tree->Branch("_obvious_flashmatch_score"   ,&_obvious_flashmatch_score   ,"obvious_flashmatch_score/F");
    _tree->Branch("_neutrino_score"             ,&_neutrino_score              ,"neutrino_score/F"         );
    _tree->Branch("_obvious_cosmics"            ,&_obvious_cosmics            ,"obvious_cosmics/I"         );
    _tree->Branch("_obvious_trklen"             ,&_obvious_trklen             ,"obvious_trklen/F"          );
    _tree->Branch("_obvious_startx"             ,&_obvious_startx             ,"obvious_startx/F"          );
    _tree->Branch("_obvious_starty"             ,&_obvious_starty             ,"obvious_starty/F"          );
    _tree->Branch("_obvious_startz"             ,&_obvious_startz             ,"obvious_startz/F"          );
    _tree->Branch("_obvious_endx"               ,&_obvious_endx               ,"obvious_endx/F"            );
    _tree->Branch("_obvious_endy"               ,&_obvious_endy               ,"obvious_endy/F"            );
    _tree->Branch("_obvious_endz"               ,&_obvious_endz               ,"obvious_endz/F"            );
    _tree->Branch("_obvious"                    ,&_obvious                    ,"obvious/I"                 );
  _tree->Branch("peSpectrum"  ,"std::vector<float>",&_peSpectrum);
  _tree->Branch("peHypothesisNu","std::vector<float>",&_peHypothesisNu);
  _tree->Branch("peHypothesisCosmic","std::vector<float>",&_peHypothesisCosmic);
  }

  void ObviousCosmicFlashMatching::resetTTree(TTree* _tree)
  {
    _obvious_flashmatch_score = 1e6;
    _neutrino_score           = 1e6;
    _obvious_cosmics = 0;
    _peHypothesisNu.clear();
    _peHypothesisCosmic.clear();
  }


  bool ObviousCosmicFlashMatching::TrackInTime(const art::Ptr<recob::Track>& trk) {
    
    auto vtx = trk->Vertex();
    auto end = trk->End();
    
    if (vtx.X() < -10 || vtx.X() > 260.) return false;
    if (end.X() < -10 || end.X() > 260.) return false;
    
    return true;
  }// is the track in time with the beam?

  void ObviousCosmicFlashMatching::AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v) {
    
    auto daughters = pfp_ptr->Daughters();
    
    pfp_v.push_back(pfp_ptr);
    
    std::cout << "\t PFP w/ PdgCode " << pfp_ptr->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;
    
    for(auto const& daughterid : daughters) {
      
      if (_pfpmap.find(daughterid) == _pfpmap.end()) {
	std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
	continue;
      }
      
      const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, _pfpmap.at(daughterid) );
      
      AddDaughters(pfp_ptr, pfp_h, pfp_v);
      
    }// for all daughters
    
    return;
  }

  DEFINE_ART_CLASS_TOOL(ObviousCosmicFlashMatching)
} // namespace analysis

#endif
