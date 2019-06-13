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
    // Class:       CosmicIP
    // File:        CosmicIP.cc
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

  class CosmicIP : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    CosmicIP(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~CosmicIP(){ };
    
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

    //void AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v);

    //bool TrackInTime(const art::Ptr<recob::Track>& pfp_track_assn_v);

    art::InputTag fPFPproducer;
    art::InputTag fSpacePointproducer;

    float fTrkShrScore;
    float _CosmicIP;
    
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  CosmicIP::CosmicIP(const fhicl::ParameterSet& p)
  {

    fPFPproducer         = p.get< art::InputTag >("PFPproducer");
    fSpacePointproducer  = p.get< art::InputTag >("SpacePointproducer");
    fTrkShrScore         = p.get< float >("TrkShrScore");
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void CosmicIP::configure(fhicl::ParameterSet const & p)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void CosmicIP::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
  {

    // std::cout << "[NEW EVENT]" << e.event() << std::endl;

    // set defaults
    _CosmicIP = 9999.;

    // first get the candidate shower's start point
    float shr_candidate_energy = 0.;
    TVector3 shrStart;
    // loop over all PFPs and get metadata, find the neutrino and save its vertex
    for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++) {
      
      auto const &pfp_pxy = slice_pfp_v.at(i_pfp);
      
      auto PDG = fabs(pfp_pxy->PdgCode());
      
      // skip neutrino PFP
      if ((PDG == 12) || (PDG == 14))
	continue;
      
      // grab shower/track score
      auto trkshrscore = searchingfornues::GetTrackShowerScore(pfp_pxy);
      
      // 1 -> track-like
      if (trkshrscore > fTrkShrScore)
	continue;
      
      if ( pfp_pxy.get<recob::Shower>().size() == 0 )
	continue;
      
      auto const &shr = pfp_pxy.get<recob::Shower>().at(0);
      
      if (shr->Energy()[2] > shr_candidate_energy) {
	shr_candidate_energy = shr->Energy()[2];
	shrStart = shr->ShowerStart();
      }
    }// for all PFParticles

    if (shr_candidate_energy == 0)
      return;
    
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

    for (unsigned int p=0; p < pfp_h->size(); p++){
      
      auto const& pfp = pfp_h->at(p);

      // start from primary PFParticles
      if (pfp.IsPrimary() == false) continue;

      // cut on PDG code. We do not want the neutrino candidate!
      if (pfp.PdgCode() != 13) continue;

      // bool clearCosmic = false;

      // get metadata for this PFP
      const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(p));

      if (pfParticleMetadataList.empty()) continue;
	
      for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j) {
	const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
	auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
	if (!pfParticlePropertiesMap.empty())
	  for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	    if ( it->first == "IsClearCosmic" ) {
	      // clearCosmic = true;
	      break;
	    }
	  }// for all metadata items in the particle metadata
      }// for entries in list

      // if (clearCosmic) { std::cout << "ClearCosmic" << std::endl; }

      // get spacepoints associated to PFParticle
      const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
      const std::vector< art::Ptr<recob::SpacePoint> >& spacepoint_ptr_v = pfp_spacepoint_assn_v.at( pfp_ptr.key() );

      float dmax = 9999.;
      
      // loop through spacepoints, find closest ditance to shower start point
      for (size_t s=0; s < spacepoint_ptr_v.size(); s++) {
	auto sps = spacepoint_ptr_v[s];
	auto spsx = sps->XYZ()[0];
	auto spsy = sps->XYZ()[1];
	auto spsz = sps->XYZ()[2];
	double d = sqrt(  ((spsx-shrStart.X())*(spsx-shrStart.X())) + 
			  ((spsy-shrStart.Y())*(spsy-shrStart.Y())) +
			  ((spsz-shrStart.Z())*(spsz-shrStart.Z())) );
	if (d < dmax) {
	  dmax = d;
	}
      }// for all spacepoints

      if (dmax < _CosmicIP) 
	_CosmicIP = dmax;

      // std::cout << "DMAX is      : " << dmax << std::endl;
      // std::cout << "Cosmic IP is : " << _CosmicIP << std::endl;
      // std::cout << std::endl;

    }// for all pfparticles

    // std::cout << "COSMIC IP is " << _CosmicIP << std::endl;
    
    return;
  }

  void CosmicIP::analyzeEvent(art::Event const &e, bool fData)
  {
    // std::cout << "analyze event" << std::endl;
  }

  void CosmicIP::setBranches(TTree* _tree)
  {
    _tree->Branch("CosmicIP",&_CosmicIP,"CosmicIP/F");
  }
  
  void CosmicIP::resetTTree(TTree* _tree)
  {
    _CosmicIP = 9999.;
  }

  
  DEFINE_ART_CLASS_TOOL(CosmicIP)
} // namespace analysis

#endif
