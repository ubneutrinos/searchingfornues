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

#include "../../FlashMatching/SliceCandidate.h"

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

    void AddDaughters(const art::Ptr<recob::PFParticle>& pfp,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v);
    
  private:

    art::InputTag fPFPproducer;
    art::InputTag fT0producer;
    art::InputTag fFLASHproducer;
    art::InputTag fSpacePointproducer;
    art::InputTag fAllPFPproducer;
    art::InputTag fAllSpacePointproducer;
    // art::InputTag fAllT0producer;
    art::InputTag fHproducer;
    art::InputTag fHTproducer;

    float _flash_pe;
    float _flash_time, _flash_timewidth;
    float _flash_y, _flash_z, _flash_ywidth, _flash_zwidth;
    std::vector<float> _flash_pe_v; //* reconstructed PE spectrum across PMT array *//
    std::vector<float> _slice_pe_v; //* hypothesized PE spectrum based on TPC charge across PMT array *//

    float _nu_flashmatch_score;
    float _nu_centerX, _nu_centerY, _nu_centerZ, _nu_totalCharge;

    float _best_cosmic_flashmatch_score;
    float _best_obviouscosmic_flashmatch_score;
    std::vector<float> _cosmic_flashmatch_score_v;
    std::vector<float> _cosmic_topological_score_v;
    std::vector<float> _cosmic_centerX_v, _cosmic_centerY_v, _cosmic_centerZ_v, _cosmic_totalCharge_v;
    std::vector<int> _cosmic_nhits_v, _cosmic_nunhits_v, _cosmic_isclear_v;

    float  m_chargeToNPhotonsTrack;   ///< The conversion factor between charge and number of photons for tracks
    float  m_chargeToNPhotonsShower;  ///< The conversion factor between charge and number of photons for showers

    std::map<unsigned int, unsigned int> _pfpmap;

    flashana::FlashMatchManager    m_flashMatchManager;       ///< The flash match manager
    bool m_applyLifetimeCorr;
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
    fFLASHproducer  = p.get< art::InputTag >("FLASHproducer" );
    fSpacePointproducer = p.get< art::InputTag >("SpacePointproducer");
    fAllPFPproducer = p.get< art::InputTag >("AllPFPproducer");
    fAllSpacePointproducer = p.get< art::InputTag >("AllSpacePointproducer");
    // fAllT0producer  = p.get< art::InputTag >("AllT0producer" );
    fHproducer = p.get<art::InputTag>("Hproducer");
    fHTproducer = p.get<art::InputTag>("HTproducer");
    m_chargeToNPhotonsTrack = p.get< float >("ChargeToNPhotonsTrack");
    m_chargeToNPhotonsShower = p.get< float >("ChargeToNPhotonsShower");
    m_flashMatchManager.Configure(p.get<flashana::Config_t>("FlashMatchConfig"));
    m_applyLifetimeCorr = p.get<bool>("ApplyLifetimeCorr");
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

    art::ValidHandle<std::vector<recob::PFParticle>> pfp_h = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    // grab PFP -> T0 flash-matching association for the event
    art::FindManyP< anab::T0 > pfp_t0_assn_v(pfp_h, e, fT0producer);
    // grab associated metadata
    art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fPFPproducer);
    art::FindManyP<recob::Hit> spacepoint_hit_assn_v(e.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointproducer), e, fSpacePointproducer);

    std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
    std::vector<std::vector<art::Ptr<recob::SpacePoint> > > spacepoint_v_v;
    std::vector<std::vector<art::Ptr<recob::Hit> > > hit_v_v;

    // figure out which PFP is the neutrino
    for (auto pfp : slice_pfp_v) {

      if ( (pfp->PdgCode() == 12) || (pfp->PdgCode() == 14) ) {
	size_t key = pfp_h->size();
	for (size_t p=0;p<pfp_h->size();p++) {
	  if (pfp->Self() == pfp_h->at(p).Self()) {
	      key = p;
	      break;
	    }
	}

	// get flash-match score
	if (pfp_t0_assn_v.size() <= key) { std::cout << "NO T0!" << std::endl; continue; }
	if (pfp_t0_assn_v.at(key).size() != 1) { std::cout << "NO T0!" << std::endl; continue; }

	auto fmscore = pfp_t0_assn_v.at(key).at(0)->TriggerConfidence();
	_nu_flashmatch_score = fmscore;
	//std::cout << "nu fmscore=" << fmscore << std::endl;

	continue;
      }

      size_t p=0;
      for (;p<pfp_h->size();p++) if (pfp->Self() == pfp_h->at(p).Self()) break;
      const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
      const std::vector< art::Ptr<recob::SpacePoint> >& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(pfp_ptr.key());
      std::vector< art::Ptr<recob::Hit> > hit_ptr_v;
      for (size_t sp=0; sp < spacepoint_ptr_v.size(); sp++) {
	auto const& spkey = spacepoint_ptr_v.at(sp).key();
	const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = spacepoint_hit_assn_v.at( spkey );
	for (size_t h=0; h < this_hit_ptr_v.size(); h++) {
	  hit_ptr_v.push_back( this_hit_ptr_v.at( h ) );
	}// for all hits associated to this spacepoint
      }// fpr all spacepoints
      pfp_ptr_v.push_back(pfp_ptr);
      spacepoint_v_v.push_back( spacepoint_ptr_v );
      hit_v_v.push_back( hit_ptr_v );
    }
    flashmatch::SliceCandidate slice(pfp_ptr_v, spacepoint_v_v, hit_v_v, m_chargeToNPhotonsTrack, m_chargeToNPhotonsShower);
    _nu_centerX = slice.m_centerX;
    _nu_centerY = slice.m_centerY;
    _nu_centerZ = slice.m_centerZ;
    _nu_totalCharge = slice.m_totalCharge;
    //std::cout << "nu center Z=" << _nu_centerZ << std::endl;

    return;
  }

  void FlashMatching::analyzeEvent(art::Event const &e, bool fData)
  {

    art::Handle<std::vector<recob::OpFlash> > flash_h;
    e.getByLabel( fFLASHproducer , flash_h );

    size_t ibeamFlash = flash_h->size();
    for (size_t f=0; f < flash_h->size(); f++) {

      auto flash = flash_h->at(f);

      if (flash.TotalPE() > _flash_pe) {
	_flash_pe    = flash.TotalPE();
	_flash_time  = flash.Time();
	_flash_y = flash.YCenter();
	_flash_z = flash.ZCenter();
	_flash_timewidth = flash.TimeWidth();
	_flash_ywidth = flash.YWidth();
	_flash_zwidth = flash.ZWidth();
	ibeamFlash = f;
      }// if larger then other flashes

    }// for all flashes

    art::ValidHandle<std::vector<recob::PFParticle>> pfp_h = e.getValidHandle<std::vector<recob::PFParticle>>(fAllPFPproducer);
    // grab associated metadata
    // art::FindManyP< anab::T0 > pfp_t0_assn_v(pfp_h, e, fAllT0producer);
    art::FindManyP< larpandoraobj::PFParticleMetadata > pfp_meta_assn_v(pfp_h, e, fAllPFPproducer);

    art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fAllPFPproducer);
    art::FindManyP<recob::Hit> spacepoint_hit_assn_v(e.getValidHandle<std::vector<recob::SpacePoint> >(fAllSpacePointproducer), e, fAllSpacePointproducer);

    auto assocPfpSlice = std::unique_ptr<art::FindManyP<recob::Slice> >(new art::FindManyP<recob::Slice>(pfp_h, e, fAllPFPproducer));
    auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit> >(new art::FindManyP<recob::Hit>(e.getValidHandle<std::vector<recob::Slice> >(fAllPFPproducer), e, fAllPFPproducer));
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
    if (!fData) assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(e.getValidHandle<std::vector<recob::Hit>>(fHproducer), e, fHTproducer));

    // fill map: pfparticle Self() -> index/key
    _pfpmap.clear();
    for (unsigned int p=0; p < pfp_h->size(); p++) _pfpmap[pfp_h->at(p).Self()] = p;

    // loop through all PFParticles
    for (size_t p=0; p < pfp_h->size(); p++) {

      auto const& pfp = pfp_h->at(p);
      //std::cout << "pfp id=" << p << " pdg=" << pfp.PdgCode() << " primary=" << pfp.IsPrimary() << std::endl;
      const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);

      // only primary PFPs have a flash-match score
      if (pfp.IsPrimary() == false) continue;

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

      if ( !(pfp.PdgCode() == 12 || pfp.PdgCode() == 14 || clearcosmic) ) continue;

      // now build vectors of PFParticles, space-points, and hits for this slice
      std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
      std::vector<std::vector<art::Ptr<recob::SpacePoint> > > spacepoint_v_v;
      std::vector<std::vector<art::Ptr<recob::Hit> > > hit_v_v;

      AddDaughters(pfp_ptr, pfp_h, pfp_ptr_v);

      // go through these pfparticles and fill info needed for matching
      for (size_t i=0; i < pfp_ptr_v.size(); i++) {
	auto key = pfp_ptr_v.at(i).key();
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
	//std::cout << "size sps=" << spacepoint_ptr_v.size() << " hits=" << hit_ptr_v.size() << std::endl;
      }// for all pfp pointers
      flashmatch::SliceCandidate slice(pfp_ptr_v, spacepoint_v_v, hit_v_v, m_chargeToNPhotonsTrack, m_chargeToNPhotonsShower, m_applyLifetimeCorr);
      float fmscore = 1e6;
      if (ibeamFlash != flash_h->size()) {
	flashmatch::FlashCandidate beamFlash(e,flash_h->at(ibeamFlash));
	fmscore = slice.GetFlashMatchScore(beamFlash, m_flashMatchManager);
	// store flash PE spectrum in the ordering as used by Flash-Matching 
	_flash_pe_v = beamFlash.m_peSpectrum;
	_slice_pe_v = beamFlash.m_peHypSpectrum;
      }
      // get flash-match score (in case there is a valid T0 for all outcomes...)
      // size_t key = pfp_ptr.key();
      // if (pfp_t0_assn_v.size() <= key) { std::cout << "NO T0!" << std::endl; continue; }
      // if (pfp_t0_assn_v.at(key).size() != 1) { std::cout << "NO T0!" << std::endl; continue; }
      // auto fmscore = pfp_t0_assn_v.at(key).at(0)->TriggerConfidence();
      //std::cout << "cosmic fmscore=" << fmscore << std::endl;

      //std::cout << "center Z=" << slice.m_centerZ << std::endl;
      // if not the neutrino...

      if (clearcosmic) {
	if ( fmscore>-1. && fmscore < _best_obviouscosmic_flashmatch_score ) _best_obviouscosmic_flashmatch_score = fmscore;
      } else {

	_cosmic_flashmatch_score_v.push_back( fmscore );
	_cosmic_centerX_v.push_back( slice.m_centerX );
	_cosmic_centerY_v.push_back( slice.m_centerY );
	_cosmic_centerZ_v.push_back( slice.m_centerZ );
	_cosmic_totalCharge_v.push_back( slice.m_totalCharge );

	auto aslice = assocPfpSlice->at(pfp_ptr.key());
	int nhits = 0, nunhits = 0;
	if (aslice.size()>0) {
	  auto slhits = assocSliceHit->at(aslice[0].key());
	  nhits = slhits.size();
	  if (!fData) {
	    for (size_t ih=0;ih<slhits.size();ih++) {
	      if (assocMCPart->at(slhits[ih].key()).size()) nunhits++;
	    }
	  }
	}
	_cosmic_nhits_v.push_back( nhits );
	_cosmic_nunhits_v.push_back( nunhits );

	if (fmscore < _best_cosmic_flashmatch_score) { _best_cosmic_flashmatch_score = fmscore; }

	_cosmic_isclear_v.push_back( clearcosmic );
	_cosmic_topological_score_v.push_back( nuscore );
	std::cout << "Slice has flash-match score of " << fmscore << " with NuScore of " << nuscore << std::endl;
      }

    }// for all PFPs


    return;
  }

  void FlashMatching::setBranches(TTree* _tree)
  {

    _tree->Branch("flash_pe",&_flash_pe,"flash_pe/F");
    _tree->Branch("flash_pe_v","std::vector<float>",&_flash_pe_v);
    _tree->Branch("slice_pe_v","std::vector<float>",&_slice_pe_v);
    _tree->Branch("flash_time",&_flash_time,"flash_time/F");
    _tree->Branch("flash_y",&_flash_y,"flash_y/F");
    _tree->Branch("flash_z",&_flash_z,"flash_z/F");
    _tree->Branch("flash_timewidth",&_flash_timewidth,"flash_timewidth/F");
    _tree->Branch("flash_ywidth",&_flash_ywidth,"flash_ywidth/F");
    _tree->Branch("flash_zwidth",&_flash_zwidth,"flash_zwidth/F");

    _tree->Branch("nu_flashmatch_score",&_nu_flashmatch_score,"nu_flashmatch_score/F");
    _tree->Branch("nu_centerX",&_nu_centerX,"nu_centerX/F");
    _tree->Branch("nu_centerY",&_nu_centerY,"nu_centerY/F");
    _tree->Branch("nu_centerZ",&_nu_centerZ,"nu_centerZ/F");
    _tree->Branch("nu_totalCharge",&_nu_totalCharge,"nu_totalCharge/F");
    _tree->Branch("best_cosmic_flashmatch_score",&_best_cosmic_flashmatch_score,"best_cosmic_flashmatch_score/F");
    _tree->Branch("best_obviouscosmic_flashmatch_score",&_best_obviouscosmic_flashmatch_score,"best_obviouscosmic_flashmatch_score/F");
    _tree->Branch("cosmic_flashmatch_score_v","std::vector<float>",&_cosmic_flashmatch_score_v);
    _tree->Branch("cosmic_topological_score_v","std::vector<float>",&_cosmic_topological_score_v);
    _tree->Branch("cosmic_centerX_v","std::vector<float>",&_cosmic_centerX_v);
    _tree->Branch("cosmic_centerY_v","std::vector<float>",&_cosmic_centerY_v);
    _tree->Branch("cosmic_centerZ_v","std::vector<float>",&_cosmic_centerZ_v);
    _tree->Branch("cosmic_totalCharge_v","std::vector<float>",&_cosmic_totalCharge_v);
    _tree->Branch("cosmic_nhits_v","std::vector<int>",&_cosmic_nhits_v);
    _tree->Branch("cosmic_nunhits_v","std::vector<int>",&_cosmic_nunhits_v);
    _tree->Branch("cosmic_isclear_v","std::vector<int>",&_cosmic_isclear_v);
  }
  
  void FlashMatching::resetTTree(TTree* _tree)
  {
    _nu_flashmatch_score = 1e6;
    _nu_centerX = std::numeric_limits<float>::lowest();
    _nu_centerY = std::numeric_limits<float>::lowest();
    _nu_centerZ = std::numeric_limits<float>::lowest();
    _nu_totalCharge = std::numeric_limits<float>::lowest();
    _best_cosmic_flashmatch_score = 1e6;
    _best_obviouscosmic_flashmatch_score = 1e6;
    _cosmic_flashmatch_score_v.clear();
    _cosmic_topological_score_v.clear();
    _cosmic_centerX_v.clear();
    _cosmic_centerY_v.clear();
    _cosmic_centerZ_v.clear();
    _cosmic_totalCharge_v.clear();
    _cosmic_nhits_v.clear();
    _cosmic_nunhits_v.clear();
    _cosmic_isclear_v.clear();
    _flash_pe   = std::numeric_limits<float>::lowest();
    _flash_time = std::numeric_limits<float>::lowest();
    _flash_y = std::numeric_limits<float>::lowest();
    _flash_z = std::numeric_limits<float>::lowest();
    _flash_timewidth = std::numeric_limits<float>::lowest();
    _flash_ywidth = std::numeric_limits<float>::lowest();
    _flash_zwidth = std::numeric_limits<float>::lowest();
  }

  void FlashMatching::AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v) {
  
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
  
  DEFINE_ART_CLASS_TOOL(FlashMatching)
} // namespace analysis

#endif
