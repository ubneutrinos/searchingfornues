#ifndef SELECTION_SELECTIONEXAMPLE_CXX
#define SELECTION_SELECTIONEXAMPLE_CXX

#include <iostream>
#include "SelectionToolBase.h"

namespace selection
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       ShowerSelection
    // File:        ShowerSelection.cc
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
    
  class ShowerSelection : public SelectionToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    ShowerSelection(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~ShowerSelection(){};
    
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
    void setBranches(TTree* _tree);

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree* _tree);

    /**
     * @brief reset module
     */
    void Reset();
    
  private:

    /**
     * @brief calculate PFP energy based on hits associated to clusters
     */
    template <typename T> float PFPEnergy(const T& ass_clus_v);
    
    /**
     * @brief given a PFP (the neutrino one) grab the associated vertex and store vertex by reference
     */
    void StoreVertex(const ProxyPfpElem_t& pfp_pxy, TVector3& nuvtx);

    /*
    * @brief get shower score
    */
    float GetTrackShowerScore(const ProxyPfpElem_t& pfp_pxy);

    // TTree variables
    int _nshower;
    float _shr_score;  // shower classification score (1 = track-like)
    float _shr_energy; // shower energy
    float _shr_dedx;   // shower dEdx
    float _shr_dist;   // shower start point distance to vertex

    // input parameters
    float _trkshrscore;
    
  };

  template <typename T> float ShowerSelection::PFPEnergy(const T& ass_clus_v) {

    float energy0 = 0;
    float energy1 = 0;
    float energy2 = 0;
    
    for (auto ass_clus : ass_clus_v) {

      //auto clus = clus_v[ass_clus.key()];
      
      std::cout << "cluster integral : " << ass_clus->Integral() << std::endl;

      if (ass_clus->View() == 0) { energy0 = ass_clus->Integral(); }
      if (ass_clus->View() == 1) { energy1 = ass_clus->Integral(); }
      if (ass_clus->View() == 2) { energy2 = ass_clus->Integral(); }

    }// for all clusters

    if (energy2 != 0) return energy2;
    return (energy1+energy0)/2.;
  }// calculte PFP energy based on associated hit charge
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  ShowerSelection::ShowerSelection(const fhicl::ParameterSet& pset)
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
  void ShowerSelection::configure(fhicl::ParameterSet const & pset)
  {
    _trkshrscore = pset.get< float > ("trkshrscore");
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
  bool ShowerSelection::selectEvent(art::Event const& e,
				    const std::vector<ProxyPfpElem_t>& pfp_pxy_v)
  {

    Reset();

    // container to store vertex
    TVector3 nuvtx;

    // loop over all PFPs and get metadata, find the neutrino and save its vertex
    for (const auto& pfp_pxy : pfp_pxy_v) {
      
      auto PDG = fabs(pfp_pxy->PdgCode());

      // skip neutrino PFP
      if ( (PDG == 12) || (PDG == 14) ) 
	StoreVertex(pfp_pxy,nuvtx);
      
      // if non-neutrino PFP
      else {
	
	// grab shower/track score
	auto trkshrscore = GetTrackShowerScore(pfp_pxy);

	std::cout << "DAVIDC PFP has trkscore : " << trkshrscore << std::endl;

	// 1 -> track-like
	if (trkshrscore > _trkshrscore)  continue;

	_shr_score = trkshrscore;
	
	auto nshr = pfp_pxy.get<recob::Shower>().size();
	
	_nshower += nshr;
	
	// 1 -> track-like
	if (trkshrscore > _trkshrscore)  continue;
	
	if (nshr != 1) continue;
	
	auto const& shr = pfp_pxy.get<recob::Shower>().at(0);
	
	// if this is the highest energy shower, save as shower candidate
	if (shr->Energy()[2] > _shr_energy) {
	  _shr_energy = shr->Energy()[2];
	  _shr_dedx   = shr->dEdx()[2];
	  _shr_dist   = (shr->Direction() - nuvtx).Mag();
	}// if highest energy shower so far

      }// if non-neutrino PFP
    }// for all PFParticles

    
    if ( _nshower >= 1)
      return true;
    
    return false;
  }

  float ShowerSelection::GetTrackShowerScore(const ProxyPfpElem_t& pfp_pxy) {

    const auto& pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();

    if (pfParticleMetadataList.size() == 0) 
      return 1;

    for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j) {
      
      const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
      auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
      if (!pfParticlePropertiesMap.empty()) {
	for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	  if (it->first == "TrackScore")
	    return it->second;
	}// for map elements
      }// if pfp metadata map not empty
    }// for list

    return 1;
  }
  
  void ShowerSelection::StoreVertex(const ProxyPfpElem_t& pfp_pxy, TVector3& nuvtx) {
    
    Double_t xyz[3] = {};

    // grab vertex
    auto vtx = pfp_pxy.get<recob::Vertex>();
    if (vtx.size() != 1) {
      std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
      return;
    }
    
    // save vertex to array
    vtx.at(0)->XYZ(xyz);
    nuvtx = TVector3(xyz[0],xyz[1],xyz[2]);
  
    return;
  }// end store vertex
  
  void ShowerSelection::setBranches(TTree* _tree) {
    
    _tree->Branch("_nshower",&_nshower,"nshower/I");
    _tree->Branch("_shr_score" ,&_shr_score ,"shr_score/F" );
    _tree->Branch("_shr_energy",&_shr_energy,"shr_energy/F");
    _tree->Branch("_shr_dedx"  ,&_shr_dedx  ,"shr_dedx/F"  );
    _tree->Branch("_shr_dist"  ,&_shr_dist  ,"shr_dist/F"  );
    
    return;
  }

  void ShowerSelection::Reset() {

    _shr_score  = -1;    
    _shr_energy = 0;
    _shr_dedx   = 0;
    _shr_dist   = -1;

    _nshower = 0;

    return;
  }

  void ShowerSelection::resetTTree(TTree* _tree) {
    return;
  }
  
  
  DEFINE_ART_CLASS_TOOL(ShowerSelection)
} // namespace selection

#endif
