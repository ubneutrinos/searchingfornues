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
     * @brief given a PFP (the neutrino one) grab the associated vertex and store vertex by reference
     */
    void StoreVertex(const ProxyPfpElem_t& pfp_pxy, TVector3& nuvtx);

    /*
    * @brief get shower score
    */
    float GetTrackShowerScore(const ProxyPfpElem_t& pfp_pxy);

    // TTree variables
    int _nshower;
    int _ntrack;
    float _maxtrklen; // maximum track length for any track associated in the slice
    float _shr_score;  // shower classification score (1 = track-like)
    float _shr_energy_Y, _shr_energy_V, _shr_energy_U; // shower energy
    float _shr_dedx_Y, _shr_dedx_V, _shr_dedx_U;       // shower dEdx
    float _shr_dist;   // shower start point distance to vertex
    float _shr_x, _shr_y, _shr_z; // shower start poisition
    float _shr_px, _shr_py, _shr_pz; // shower momentum vector

    // input parameters
    float fTrkShrscore;
    float fShrdedxmax;
    float fShrRadlen;
    float fShrEnergy;
    float fMaxTrklen;
    
  };

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
    fTrkShrscore = pset.get< float > ("TrkShrscore");
    fShrdedxmax  = pset.get< float > ("Shrdedxmax");
    fShrRadlen   = pset.get< float > ("ShrRadlen");
    fShrEnergy   = pset.get< float > ("ShrEnergy");
    fMaxTrklen   = pset.get< float > ("MaxTrkLen");
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

	auto ntrk = pfp_pxy.get<recob::Track>().size();
	_ntrack += ntrk;

	// save track-length of longest track
	if (ntrk == 1) {
	  if (pfp_pxy.get<recob::Track>()[0]->Length() > _maxtrklen)
	    _maxtrklen = pfp_pxy.get<recob::Track>()[0]->Length();
	}// if there is a track associated

	// 1 -> track-like
	if (trkshrscore > fTrkShrscore)  continue;

	_shr_score = trkshrscore;
	
	auto nshr = pfp_pxy.get<recob::Shower>().size();
	_nshower += nshr;

	// 1 -> track-like
	if (trkshrscore > fTrkShrscore)  continue;
	
	if (nshr != 1) continue;
	
	auto const& shr = pfp_pxy.get<recob::Shower>().at(0);
	
	// if this is the highest energy shower, save as shower candidate
	if (shr->Energy()[2] > _shr_energy_Y) {
	  _shr_energy_Y = shr->Energy()[2];
	  _shr_dedx_Y   = shr->dEdx()[2];
	  _shr_energy_V = shr->Energy()[1];
	  _shr_dedx_V   = shr->dEdx()[1];
	  _shr_energy_U = shr->Energy()[0];
	  _shr_dedx_U   = shr->dEdx()[0];
	  _shr_x        = shr->ShowerStart().X();
	  _shr_y        = shr->ShowerStart().Y();
	  _shr_z        = shr->ShowerStart().Z();
	  _shr_px        = shr->Direction().X();
	  _shr_py        = shr->Direction().Y();
	  _shr_pz        = shr->Direction().Z();
	  _shr_dist   = (shr->ShowerStart() - nuvtx).Mag();
	}// if highest energy shower so far

      }// if non-neutrino PFP
    }// for all PFParticles


    if (_nshower < 1)
      return false;
    if (_shr_dist > fShrRadlen)
      return false;
    if (_shr_dedx_Y > fShrdedxmax)
      return false;
    if (_shr_energy_Y < fShrEnergy)
      return false;
    if (_shr_score > fTrkShrscore)
      return false;
    if (_maxtrklen > fMaxTrklen)
      return false;
    
    return true;
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
    _tree->Branch("_ntrack" ,&_ntrack ,"ntrack/I" );
    _tree->Branch("_maxtrklen",&_maxtrklen,"maxtrklen/F");
    _tree->Branch("_shr_score" ,&_shr_score ,"shr_score/F" );
    _tree->Branch("_shr_energy_Y",&_shr_energy_Y,"shr_energy_Y/F");
    _tree->Branch("_shr_dedx_Y"  ,&_shr_dedx_Y  ,"shr_dedx_Y/F"  );
    _tree->Branch("_shr_energy_V",&_shr_energy_V,"shr_energy_V/F");
    _tree->Branch("_shr_dedx_V"  ,&_shr_dedx_V  ,"shr_dedx_V/F"  );
    _tree->Branch("_shr_energy_U",&_shr_energy_U,"shr_energy_U/F");
    _tree->Branch("_shr_dedx_U"  ,&_shr_dedx_U  ,"shr_dedx_U/F"  );
    _tree->Branch("_shr_dist"  ,&_shr_dist  ,"shr_dist/F"  );
    _tree->Branch("_shr_x"  ,&_shr_x  ,"shr_x/F"  );
    _tree->Branch("_shr_y"  ,&_shr_y  ,"shr_y/F"  );
    _tree->Branch("_shr_z"  ,&_shr_z  ,"shr_z/F"  );
    _tree->Branch("_shr_px"  ,&_shr_px  ,"shr_px/F"  );
    _tree->Branch("_shr_py"  ,&_shr_py  ,"shr_py/F"  );
    _tree->Branch("_shr_pz"  ,&_shr_pz  ,"shr_pz/F"  );
    
    return;
  }

  void ShowerSelection::Reset() {

    _shr_score  = -1;    
    _shr_energy_Y = 0;
    _shr_dedx_Y   = 0;
    _shr_energy_V = 0;
    _shr_dedx_V   = 0;
    _shr_energy_U = 0;
    _shr_dedx_U   = 0;
    _shr_dist   = -1;
    _shr_x        = 0;
    _shr_y        = 0;
    _shr_z        = 0;
    _shr_px       = 0;
    _shr_py       = 0;
    _shr_pz       = 0;

    _maxtrklen = 0;

    _nshower = 0;
    _ntrack  = 0;

    return;
  }

  void ShowerSelection::resetTTree(TTree* _tree) {
    return;
  }
  
  
  DEFINE_ART_CLASS_TOOL(ShowerSelection)
} // namespace selection

#endif
