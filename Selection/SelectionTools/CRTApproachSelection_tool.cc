#ifndef SELECTION_CRTAPPROACHSEL_CXX
#define SELECTION_CRTAPPROACHSEL_CXX

#include <iostream>

#include "SelectionToolBase.h"
//#include "Analysis/anabase.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "canvas/Persistency/Common/FindMany.h"
//#include "RecoBase.h"
#include "TRandom3.h"


const int kMaxCosm = 200;
int _kMaxCosm=200;
namespace selection
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       CRTApproachSelection
    // File:        CRTApproachSelection.cc
    //
    // Select (or reject) neutrino verteces close to tagged CRT track
    //
    // Configuration parameters:
    //
    // TBD
    //
    // Created by Elena Gramellini (elenag@fnal.gov) on 03/07/2019
    // Based off of David's Caratelli basic example
    // 
    // 
    // TO DO:
    // [  ] assigne a sensible track ID to closest tagged cosmic track
    // [  ] assign a sensible t0 to closest tagged cosmic track
    // 
    ////////////////////////////////////////////////////////////////////////

    
  class CRTApproachSelection : public SelectionToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    CRTApproachSelection(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~CRTApproachSelection(){};
    
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
    void resetTTree(TTree* _tree){};


    
  private:

    double Point2PointDistance(const TVector3& nuvtx, TVector3& taggedTrackSpt);
    void Reset();

    //LArSoft Config
    std::string  fTrackAssnModuleLabel;
    
    // Random variable
    TRandom3 rand;    
    //Fiducial Volume <-- hardcoding makes Jesus cry
    Double_t fidVolMinX =    0;
    Double_t fidVolMaxX =  256;
    Double_t fidVolMinY = -116;
    Double_t fidVolMaxY =  116;
    Double_t fidVolMinZ =    0;
    Double_t fidVolMaxZ = 1030;


    
    // TTree variables associated with cosmic tracks
    int _nAssnCosmics;
    Double_t _t0Cosm_v[kMaxCosm];
    Double_t _t0timesVCosm_v[kMaxCosm];
    Double_t _xStartCosm_v[kMaxCosm];
    Double_t _yStartCosm_v[kMaxCosm];
    Double_t _zStartCosm_v[kMaxCosm];
    Double_t _xEndCosm_v[kMaxCosm];
    Double_t _yEndCosm_v[kMaxCosm];
    Double_t _zEndCosm_v[kMaxCosm];

    // TTree variables associated with reconstructed nu
    Double_t _recoNu_vtx_x, _recoNu_vtx_y, _recoNu_vtx_z; // reco neutrino vertex
    Double_t _t0_nu_cosmic;                               // Time of the associated CRT hit 
    Double_t _nu_cosmic_x, _nu_cosmic_y, _nu_cosmic_z;    // position of closest trjPoint of the cosmic to the neutrino 
    Double_t _nu_cosmic_Length;                           // length of the cosmic track
    Double_t _nu_cosmic_Start_x,_nu_cosmic_Start_y,_nu_cosmic_Start_z; // Start point of the cosmics track
    Double_t _nu_cosmic_End_x  ,_nu_cosmic_End_y  ,_nu_cosmic_End_z  ; // End point of the cosmic track
    Double_t _nu_cosmic_TrackID;
    Double_t _closestNuCosmicDist  =  999999999.;

    // TTree variables associated with random point
    Double_t _rand_vtx_x, _rand_vtx_y, _rand_vtx_z; // randomly generated neutrino vertex    
    Double_t _t0_rand_cosmic;                               // Time of the associated CRT hit 
    Double_t _rand_cosmic_x, _rand_cosmic_y, _rand_cosmic_z;    // position of closest trjPoint of the cosmic to the neutrino 
    Double_t _rand_cosmic_Length;                           // length of the cosmic track
    Double_t _rand_cosmic_Start_x,_rand_cosmic_Start_y,_rand_cosmic_Start_z; // Start point of the cosmics track
    Double_t _rand_cosmic_End_x  ,_rand_cosmic_End_y  ,_rand_cosmic_End_z  ; // End point of the cosmic track
    Double_t _rand_cosmic_TrackID;
    Double_t _closestRandCosmicDist=  999999999.;
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  CRTApproachSelection::CRTApproachSelection(const fhicl::ParameterSet& pset)
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
  void CRTApproachSelection::configure(fhicl::ParameterSet const & pset)
  {
    fTrackAssnModuleLabel = pset.get< std::string >  ("TrackAssnModuleLabel","trackmatch");

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
  bool CRTApproachSelection::selectEvent(art::Event const& e,
					 const std::vector<ProxyPfpElem_t>& pfp_pxy_v)
  {
    std::cout<<"--------------------------------------------------- \n";
    Reset();
    
    // Let's load the tracks that are tagged by the CRT
    // via their T0 association
    
    // Containers declaration
    std::vector<anab::T0>        t0tags;
    std::vector<recob::Track>    tracks;
    t0tags.clear();
    tracks.clear();
    std::cout<<"1\n";
    // get handle to the association             
    auto const& assoc_handle   = e.getValidHandle<art::Assns<recob::Track, anab::T0       >>(fTrackAssnModuleLabel);
    if (assoc_handle->size() == 0) {
      std::cout<<"Before I return from Assn -----------------------------------------------------------> \n";
      return false; // no t0 tag
    }
    // Fill the track vector 
    // Silly Association counter
    size_t sillyCounter = 0; 
    for (auto &ass : *assoc_handle) {  
      art::Ptr<recob::Track> t     = ass.first;
      art::Ptr<anab::T0>     time0 = ass.second;
      tracks.emplace_back(*t);
      t0tags.emplace_back(*time0);
      if (sillyCounter < kMaxCosm)
	{
	  _t0Cosm_v[sillyCounter]       = (*time0).Time();
	  _t0timesVCosm_v[sillyCounter] = _t0Cosm_v[sillyCounter]*0.114; // Just multiplying for the drift velocity
	  _xStartCosm_v[sillyCounter]   = (*t).Start().X();
	  _yStartCosm_v[sillyCounter]   = (*t).Start().Y();
	  _zStartCosm_v[sillyCounter]   = (*t).Start().Z();
	  _xEndCosm_v[sillyCounter]     = (*t).End().X();
	  _yEndCosm_v[sillyCounter]     = (*t).End().Y();
	  _zEndCosm_v[sillyCounter]     = (*t).End().Z();
	  sillyCounter++;
	}
    }

    _nAssnCosmics = (int)t0tags.size();


    std::cout<<"3\n";
    // Containers for reconstructed neutrino vertex
    TVector3 nuvtx;
    Double_t xyz[3] = {};

    // Containers for random neutrino vertex
    TVector3 rndvtx;

    // Containers for crt tagged cosmic track
    TVector3 nuCosmicTaggedPos;
    Double_t xyzNuCosm[3] = {};

    // loop over all PFPs and get metadata, find the neutrino and save its vertex
    for (const auto& pfp_pxy : pfp_pxy_v) {
      auto PDG = fabs(pfp_pxy->PdgCode());
      if ( (PDG == 12) || (PDG == 14) ) {
	// grab vertex
	auto vtx = pfp_pxy.get<recob::Vertex>();
	if (vtx.size() != 1) {
	  std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
	  return false;
	}
	std::cout<<"Neutrino! \n";
	// save vertex to array
	vtx.at(0)->XYZ(xyz);
	nuvtx = TVector3(xyz[0],xyz[1],xyz[2]);
	_recoNu_vtx_x = nuvtx.X();
	_recoNu_vtx_y = nuvtx.Y();
	_recoNu_vtx_z = nuvtx.Z();

	// random generated "neutrino" vtx in the fiducial volume
	_rand_vtx_x = rand.Uniform(fidVolMinX,fidVolMaxX);
	_rand_vtx_y = rand.Uniform(fidVolMinY,fidVolMaxY);
	_rand_vtx_z = rand.Uniform(fidVolMinZ,fidVolMaxZ);
	
	rndvtx = TVector3(_rand_vtx_x,_rand_vtx_y,_rand_vtx_z);
	// Loop on the tracks to get the tracking points
	//for (auto &track : tracks) {
	std::cout<<"Track Size "<<tracks.size()<<"\n";
	for(size_t iT = 0; iT<tracks.size(); iT++){
	    auto track = tracks[iT];
	    for (size_t iPt = 0; iPt < track.NumberTrajectoryPoints(); iPt++ )
	    {
	      auto thisTrajPointPosition = track.TrajectoryPoint (iPt).position;
	      xyzNuCosm[0] = thisTrajPointPosition.X();
	      xyzNuCosm[1] = thisTrajPointPosition.Y();
	      xyzNuCosm[2] = thisTrajPointPosition.Z();
	      nuCosmicTaggedPos = TVector3(xyzNuCosm[0],xyzNuCosm[1],xyzNuCosm[2]);

	      // Calculate distance between Nu Vertex and Cosmic Trajectory
	      double currentDistNu = Point2PointDistance(nuvtx, nuCosmicTaggedPos);	     
	      if (currentDistNu < _closestNuCosmicDist ) {
		_closestNuCosmicDist = currentDistNu;	
		_t0_nu_cosmic      = t0tags[iT].Time();
		_nu_cosmic_x       = xyzNuCosm[0]; 
		_nu_cosmic_y       = xyzNuCosm[1]; 
		_nu_cosmic_z       = xyzNuCosm[2];
		_nu_cosmic_Length  = track.Length();                       
		_nu_cosmic_Start_x = track.Start().X();
		_nu_cosmic_Start_y = track.Start().Y();
		_nu_cosmic_Start_z = track.Start().Z();
		_nu_cosmic_End_x   = track.End().X();
		_nu_cosmic_End_y   = track.End().Y();
		_nu_cosmic_End_z   = track.End().Z();
		_nu_cosmic_TrackID = 0.;// Needs to be something like track.key();  CHANGE CHANGE CHANGE

	      }

	      // Calculate dicstance between  Random Vertex and Cosmic Trajectory Point
	      double currentDistRand = Point2PointDistance(rndvtx, nuCosmicTaggedPos);	     
	      if (currentDistRand < _closestRandCosmicDist) {
		_closestRandCosmicDist = currentDistRand;	
		_t0_rand_cosmic      = t0tags[iT].Time();
		_rand_cosmic_x       = xyzNuCosm[0]; 
		_rand_cosmic_y       = xyzNuCosm[1]; 
		_rand_cosmic_z       = xyzNuCosm[2];
		_rand_cosmic_Length  = track.Length();                       
		_rand_cosmic_Start_x = track.Start().X();
		_rand_cosmic_Start_y = track.Start().Y();
		_rand_cosmic_Start_z = track.Start().Z();
		_rand_cosmic_End_x   = track.End().X();
		_rand_cosmic_End_y   = track.End().Y();
		_rand_cosmic_End_z   = track.End().Z();
		_rand_cosmic_TrackID = 0.;// Needs to be something like track.key();  CHANGE CHANGE CHANGE
	      }
	      
	    }// Loop on track points
	}// Loop on tracks
      }// If neutrinos
    }// Paticle loop 
    return true;
  }

  
  void CRTApproachSelection::setBranches(TTree* _tree) {



    _tree->Branch("_nAssnCosmics"     ,&_nAssnCosmics     ,"_nAssnCosmics/I");
    _tree->Branch("_kMaxCosm"         ,&_kMaxCosm          ,"_kMaxCosm/I");
    _tree->Branch("_t0Cosm_v"         ,_t0Cosm_v          ,"_t0Cosm_v[_kMaxCosm]/D");
    _tree->Branch("_t0timesVCosm_v"   ,_t0timesVCosm_v    ,"_t0timesVCosm_v[_kMaxCosm]/D");
    _tree->Branch("_xStartCosm_v"     ,_xStartCosm_v      ,"_xStartCosm_v[_kMaxCosm]/D");
    _tree->Branch("_yStartCosm_v"     ,_yStartCosm_v      ,"_yStartCosm_v[_kMaxCosm]/D");
    _tree->Branch("_zStartCosm_v"     ,_zStartCosm_v      ,"_zStartCosm_v[_kMaxCosm]/D");
    _tree->Branch("_xEndCosm_v"       ,_xEndCosm_v        ,"_xEndCosm_v[_kMaxCosm]/D");
    _tree->Branch("_yEndCosm_v"       ,_yEndCosm_v        ,"_yEndCosm_v[_kMaxCosm]/D");
    _tree->Branch("_zEndCosm_v"       ,_zEndCosm_v        ,"_zEndCosm_v[_kMaxCosm]/D");


    _tree->Branch("_recoNu_vtx_x"     ,&_recoNu_vtx_x     ,"reco_vtx_x/D");
    _tree->Branch("_recoNu_vtx_y"     ,&_recoNu_vtx_y     ,"reco_vtx_y/D");
    _tree->Branch("_recoNu_vtx_z"     ,&_recoNu_vtx_z     ,"reco_vtx_z/D");
    _tree->Branch("_t0_nu_cosmic"     ,&_t0_nu_cosmic     ,"_t0_nu_cosmic/D");
    _tree->Branch("_nu_cosmic_x"      ,&_nu_cosmic_x      ,"_nu_cosmic_x/D");
    _tree->Branch("_nu_cosmic_y"      ,&_nu_cosmic_y      ,"_nu_cosmic_y/D");
    _tree->Branch("_nu_cosmic_z"      ,&_nu_cosmic_z      ,"_nu_cosmic_z/D");
    _tree->Branch("_nu_cosmic_Length" ,&_nu_cosmic_Length ,"_nu_cosmic_Length/D");
    _tree->Branch("_nu_cosmic_Start_x",&_nu_cosmic_Start_x,"_nu_cosmic_Start_x/D");
    _tree->Branch("_nu_cosmic_Start_y",&_nu_cosmic_Start_y,"_nu_cosmic_Start_y/D");
    _tree->Branch("_nu_cosmic_Start_z",&_nu_cosmic_Start_z,"_nu_cosmic_Start_z/D");
    _tree->Branch("_nu_cosmic_End_x"  ,&_nu_cosmic_End_x  ,"_nu_cosmic_End_x/D");
    _tree->Branch("_nu_cosmic_End_y"  ,&_nu_cosmic_End_y  ,"_nu_cosmic_End_y/D");
    _tree->Branch("_nu_cosmic_End_z"  ,&_nu_cosmic_End_z  ,"_nu_cosmic_End_z/D");
    _tree->Branch("_nu_cosmic_TrackID",&_nu_cosmic_TrackID,"_nu_cosmic_TrackID/D");
    _tree->Branch("_closestNuCosmicDist",&_closestNuCosmicDist,"_closestNuCosmicDist/D");

    _tree->Branch("_rand_vtx_x"         ,&_rand_vtx_x         ,"rand_vtx_x/D");
    _tree->Branch("_rand_vtx_y"         ,&_rand_vtx_y         ,"rand_vtx_y/D");
    _tree->Branch("_rand_vtx_z"         ,&_rand_vtx_z         ,"rand_vtx_z/D");
    _tree->Branch("_t0_rand_cosmic"     ,&_t0_rand_cosmic     ,"_t0_rand_cosmic/D");
    _tree->Branch("_rand_cosmic_x"      ,&_rand_cosmic_x      ,"_rand_cosmic_x/D");
    _tree->Branch("_rand_cosmic_y"      ,&_rand_cosmic_y      ,"_rand_cosmic_y/D");
    _tree->Branch("_rand_cosmic_z"      ,&_rand_cosmic_z      ,"_rand_cosmic_z/D");
    _tree->Branch("_rand_cosmic_Length" ,&_rand_cosmic_Length ,"_rand_cosmic_Length/D");
    _tree->Branch("_rand_cosmic_Start_x",&_rand_cosmic_Start_x,"_rand_cosmic_Start_x/D");
    _tree->Branch("_rand_cosmic_Start_y",&_rand_cosmic_Start_y,"_rand_cosmic_Start_y/D");
    _tree->Branch("_rand_cosmic_Start_z",&_rand_cosmic_Start_z,"_rand_cosmic_Start_z/D");
    _tree->Branch("_rand_cosmic_End_x"  ,&_rand_cosmic_End_x  ,"_rand_cosmic_End_x/D");
    _tree->Branch("_rand_cosmic_End_y"  ,&_rand_cosmic_End_y  ,"_rand_cosmic_End_y/D");
    _tree->Branch("_rand_cosmic_End_z"  ,&_rand_cosmic_End_z  ,"_rand_cosmic_End_z/D");
    _tree->Branch("_rand_cosmic_TrackID",&_rand_cosmic_TrackID,"_rand_cosmic_TrackID/D");
    _tree->Branch("_closestRandCosmicDist",&_closestRandCosmicDist,"_closestRandCosmicDist/D");

    return;
  }




  double CRTApproachSelection::Point2PointDistance(const TVector3& nuvtx, TVector3& taggedTrackSpt){
    return  (nuvtx-taggedTrackSpt).Mag();
  }
  


  void CRTApproachSelection::Reset() {
    _nAssnCosmics = -1;
    for (size_t i = 0; i < kMaxCosm; i++){
      _t0Cosm_v[i]      = -9999.;
      _t0timesVCosm_v[i]= -9999.;
      _xStartCosm_v[i]=-9999.;
      _yStartCosm_v[i]=-9999.;
      _zStartCosm_v[i]=-9999.;
      _xEndCosm_v[i]=-9999.;
      _yEndCosm_v[i]=-9999.;
      _zEndCosm_v[i]=-9999.;
    }

    _recoNu_vtx_x      = -9999.;
    _recoNu_vtx_y      = -9999.;
    _recoNu_vtx_z      = -9999.;
    _t0_nu_cosmic      = -9999.;
    _nu_cosmic_x       = -9999.;
    _nu_cosmic_y       = -9999.;
    _nu_cosmic_z       = -9999.;
    _nu_cosmic_Length  = -9999.;
    _nu_cosmic_Start_x = -9999.;
    _nu_cosmic_Start_y = -9999.;
    _nu_cosmic_Start_z = -9999.;
    _nu_cosmic_End_x   = -9999.;
    _nu_cosmic_End_y   = -9999.;
    _nu_cosmic_End_z   = -9999.;
    _nu_cosmic_TrackID = -9999.;
    _closestNuCosmicDist  =  999999999.;

    _rand_vtx_x            = -9999.;
    _rand_vtx_y            = -9999.;
    _rand_vtx_z            = -9999.;
    _t0_rand_cosmic        = -9999.;
    _rand_cosmic_x         = -9999.;
    _rand_cosmic_y         = -9999.;
    _rand_cosmic_z         = -9999.;
    _rand_cosmic_Length    = -9999.;
    _rand_cosmic_Start_x   = -9999.;
    _rand_cosmic_Start_y   = -9999.;
    _rand_cosmic_Start_z   = -9999.;
    _rand_cosmic_End_x     = -9999.;
    _rand_cosmic_End_y     = -9999.;
    _rand_cosmic_End_z     = -9999.;
    _rand_cosmic_TrackID   = -9999.;
    _closestRandCosmicDist =  999999999.;

    return;
  }// end of reset function

  
  DEFINE_ART_CLASS_TOOL(CRTApproachSelection)
} // namespace selection

#endif
