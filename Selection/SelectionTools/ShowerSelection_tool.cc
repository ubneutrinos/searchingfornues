#ifndef SELECTION_SELECTIONEXAMPLE_CXX
#define SELECTION_SELECTIONEXAMPLE_CXX

#include <iostream>
#include "SelectionToolBase.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"

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

class ShowerSelection : public SelectionToolBase
{
  
public:
  /**
   *  @brief  Constructor
   *
   *  @param  pset
   */
  ShowerSelection(const fhicl::ParameterSet &pset);
  
  /**
   *  @brief  Destructor
   */
  ~ShowerSelection(){};
  
  // provide for initialization
  void configure(fhicl::ParameterSet const &pset);
  
  /**
   * @brief Selection function
   */
  bool selectEvent(art::Event const &e,
                   const std::vector<ProxyPfpElem_t> &pfp_pxy_v);
  
  /**
   * @brief set branches for TTree
   */
  void setBranches(TTree *_tree);
  
  /**
   * @brief reset ttree branches
   */
  void resetTTree(TTree *_tree);
  
  /**
   * @brief reset module
   */
  void Reset();
  
private:
  
  // obtain track fit dedx in first 4 cm from track calo
  void TrackFitdEdx(const searchingfornues::ProxyCaloElem_t& trk,
		    float& dedxU, float& dedxV ,float& dedxY);
  
  /**
   * @brief given a PFP (the neutrino one) grab the associated vertex and store vertex by reference
   */
  void StoreVertex(const ProxyPfpElem_t &pfp_pxy, TVector3 &nuvtx);
  
  /*
   * @brief get shower score
   */
  float GetTrackShowerScore(const ProxyPfpElem_t& pfp_pxy);
  
  // TTree variables
  int _nshower;
  int _ntrack;
  int _nupdgreco;
  float _maxtrklen; // maximum track length for any track associated in the slice
  float _shr_score;  // shower classification score (1 = track-like)
  float _shr_energy_Y, _shr_energy_V, _shr_energy_U; // shower energy
  float _shr_dedx_Y, _shr_dedx_V, _shr_dedx_U;       // shower dEdx
  float _shr_dedx_fit_Y, _shr_dedx_fit_V, _shr_dedx_fit_U; // shower dEdx from track-fitter
  float _shr_dist;   // shower start point distance to vertex
  float _shr_x, _shr_y, _shr_z; // shower start poisition
  float _shr_px, _shr_py, _shr_pz; // shower momentum vector
  
  size_t _shr_maxe_pfp_i;
  
  // input parameters
  float fTrkShrscore;
  float fShrdedxmax;
  float fShrRadlen;
  float fShrEnergy;
  float fMaxTrklen;
  art::InputTag fTRKproducer;
  art::InputTag fCALproducer;
  
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

    fTRKproducer = pset.get< art::InputTag > ("TRKproducer", "");
    fCALproducer = pset.get< art::InputTag > ("CALproducer", "");

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

    searchingfornues::ProxyCaloColl_t const* tkcalo_proxy = NULL;
    if (fTRKproducer!="") {
      tkcalo_proxy = new searchingfornues::ProxyCaloColl_t( proxy::getCollection<std::vector<recob::Track> >(e,fTRKproducer,proxy::withAssociated<anab::Calorimetry>(fCALproducer)) );
    }

    // container to store vertex
    TVector3 nuvtx;

    // loop over all PFPs and get metadata, find the neutrino and save its vertex
    for (size_t i_pfp=0; i_pfp < pfp_pxy_v.size(); i_pfp++) {

      auto const& pfp_pxy = pfp_pxy_v.at(i_pfp);
      
      auto PDG = fabs(pfp_pxy->PdgCode());
      
      // skip neutrino PFP
      if ( (PDG == 12) || (PDG == 14) ) {
	StoreVertex(pfp_pxy,nuvtx);
	_nupdgreco = PDG;
      }
      
      // if non-neutrino PFP
      else {
	
	// grab shower/track score
	auto trkshrscore = searchingfornues::GetTrackShowerScore(pfp_pxy);
	
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
	  _shr_maxe_pfp_i = i_pfp; // Index of the PFParticle corresponding to the most energetic shower
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

	  if (tkcalo_proxy!=NULL) {
	    
	    for (const searchingfornues::ProxyCaloElem_t& tk : *tkcalo_proxy) {
	      
	      // find track with ID matching the pfp index (this convention apparently works only for shower fits...)
	      if (tk->ID()==int(pfp_pxy_v[i_pfp].index())) {

		TrackFitdEdx(tk, _shr_dedx_fit_U, _shr_dedx_fit_V, _shr_dedx_fit_Y);
		
	      }// if track matches shower index -> this is the track-fitted to the shower
	    }// for all track fits to showers
	  }// if track-fits to showers exist


	}// if highest energy shower so far

      }// if non-neutrino PFP
    }// for all PFParticles


    if (_nupdgreco != 12)
      return false;
    /*
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
    */
    
    return true;
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
  
  
  void ShowerSelection::TrackFitdEdx(const searchingfornues::ProxyCaloElem_t& trk,
				     float& dedxU, float& dedxV ,float& dedxY)  {
    
    dedxU = -1;
    dedxV = -1;
    dedxY = -1;
    
    auto const trkcalos = trk.get<anab::Calorimetry>();
    
    for (const auto& tkcalo : trkcalos) {
      if (tkcalo->ResidualRange().size()==0) continue;
      std::vector<float> dedx4cm;
      for (size_t ic=0; ic<tkcalo->ResidualRange().size(); ++ic) {
	if ( (tkcalo->ResidualRange().back()-tkcalo->ResidualRange()[ic]) < 4.) {
	  dedx4cm.push_back( tkcalo->dEdx()[ic] );
	}
      }
      float dedx4cm_med = -1.;
      if (dedx4cm.size()>0) {
	std::sort(dedx4cm.begin(), dedx4cm.end());
	if (dedx4cm.size()%2 == 1) dedx4cm_med = dedx4cm[dedx4cm.size()/2];
	else dedx4cm_med = 0.5*(dedx4cm[dedx4cm.size()/2] + dedx4cm[dedx4cm.size()/2 - 1]);
      }
      
      auto pl = tkcalo->PlaneID().Plane;
      if (pl == 0) { dedxU = dedx4cm_med; }
      if (pl == 1) { dedxV = dedx4cm_med; }
      if (pl == 2) { dedxY = dedx4cm_med; }
      
    }// for all calorimetry objects associated to the track
    
    return;
  }// TrackFitdEdx
  
  void ShowerSelection::setBranches(TTree* _tree) {
    
    _tree->Branch("nshower",&_nshower,"nshower/I");
    _tree->Branch("ntrack" ,&_ntrack ,"ntrack/I" );
    _tree->Branch("maxtrklen",&_maxtrklen,"maxtrklen/F");
    _tree->Branch("shr_score" ,&_shr_score ,"shr_score/F" );
    _tree->Branch("shr_energy_Y",&_shr_energy_Y,"shr_energy_Y/F");
    _tree->Branch("shr_dedx_Y"  ,&_shr_dedx_Y  ,"shr_dedx_Y/F"  );
    _tree->Branch("shr_dedx_fit_Y"  ,&_shr_dedx_fit_Y  ,"shr_dedx_fit_Y/F"  );
    _tree->Branch("shr_energy_V",&_shr_energy_V,"shr_energy_V/F");
    _tree->Branch("shr_dedx_V"  ,&_shr_dedx_V  ,"shr_dedx_V/F"  );
    _tree->Branch("shr_dedx_fit_V"  ,&_shr_dedx_fit_V  ,"shr_dedx_fit_V/F"  );
    _tree->Branch("shr_energy_U",&_shr_energy_U,"shr_energy_U/F");
    _tree->Branch("shr_dedx_U"  ,&_shr_dedx_U  ,"shr_dedx_U/F"  );
    _tree->Branch("shr_dedx_fit_U"  ,&_shr_dedx_fit_U  ,"shr_dedx_fit_U/F"  );
    _tree->Branch("shr_dist"  ,&_shr_dist  ,"shr_dist/F"  );
    _tree->Branch("shr_x"  ,&_shr_x  ,"shr_x/F"  );
    _tree->Branch("shr_y"  ,&_shr_y  ,"shr_y/F"  );
    _tree->Branch("shr_z"  ,&_shr_z  ,"shr_z/F"  );
    _tree->Branch("shr_px"  ,&_shr_px  ,"shr_px/F"  );
    _tree->Branch("shr_py"  ,&_shr_py  ,"shr_py/F"  );
    _tree->Branch("shr_pz"  ,&_shr_pz  ,"shr_pz/F"  );
  _tree->Branch("shr_maxe_pfp_i", &_shr_maxe_pfp_i, "shr_maxe_pfp_i/i");
    
    return;
  }
  
  
  void ShowerSelection::Reset()
  {

    _nupdgreco = 0;
    _shr_score = -1;
    _shr_energy_Y = 0;
    _shr_dedx_Y = 0;
    _shr_energy_V = 0;
    _shr_dedx_V = 0;
    _shr_energy_U = 0;
    _shr_dedx_U = 0;
    _shr_dist = -1;
    _shr_x = 0;
    _shr_y = 0;
    _shr_z = 0;
    _shr_px = 0;
    _shr_py = 0;
    _shr_pz = 0;
    _shr_maxe_pfp_i = 0;
    _maxtrklen = 0;
    
    _nshower = 0;
    _ntrack = 0;
    
    return;
  }
  
  void ShowerSelection::resetTTree(TTree *_tree)
{
  return;
}
  
  DEFINE_ART_CLASS_TOOL(ShowerSelection)
} // namespace selection

#endif
