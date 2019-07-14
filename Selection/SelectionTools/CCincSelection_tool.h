#ifndef SELECTION_CCINC_CXX
#define SELECTION_CCINC_CXX

#include <iostream>
#include "SelectionToolBase.h"
#include "../CommonDefs/Typedefs.h"
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/SCECorrections.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

namespace selection
{
////////////////////////////////////////////////////////////////////////
//
// Class:       CCincSelection
// File:        CCincSelection.cc
//
//              Electron neutrino charged current selection
//
// Created by Wouter Van De Pontseele July 2019
//
////////////////////////////////////////////////////////////////////////

class CCincSelection : public SelectionToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  CCincSelection(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~CCincSelection(){};

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
  void setBranches(TTree *_tree){};

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree){};

  /**
     * @brief Check if point is inside fiducial volume
     *
     * @param x array of coordinates
     */
  bool isFiducial(const double x[3]) const;

  /**
     *  @brief  Fill the tree for the electron candidate.
     *
     *  @param  pfparticle ptr The Pfparticle corresponding to the daughter.
     *  @return 1 if succesful, 0 if failure.
     */
  bool fillElectron(const art::Ptr<recob::PFParticle> &pfp,
                    const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle,
                    const art::FindManyP<anab::ParticleID> &trackPIDAssn);

private:
  // Fields for electron candidate
  size_t m_shrPfpId; /**< Index of the leading shower in the PFParticle vector */
  float m_shrHits;   /**< Number of hits of the candidate electron shower */
  float m_shrHitsU;  /**< Number of hits of the candidate electron shower, U plane */
  float m_shrHitsV;  /**< Number of hits of the candidate electron shower, V plane */
  float m_shrHitsY;  /**< Number of hits of the candidate electron shower, Y plane */
  float m_shrSps;    /**< Number of spacepoints of the candidate electron shower */

  float m_shrEnergy;     /**< Energy of the candidate electron shower (in GeV) */
  float m_shrEnergyCali; /**< Energy of the calibrated candidate electron shower (in GeV) */

  float m_shrStartX;     /**< Start x coordinate of the leading shower */
  float m_shrStartY;     /**< Start y coordinate of the leading shower */
  float m_shrStartZ;     /**< Start z coordinate of the leading shower */
  float m_shrDedxY;      /**< dE/dx of the leading shower on the Y plane with the 1x4 cm box method */
  float m_shrDedxV;      /**< dE/dx of the leading shower on the V plane with the 1x4 cm box method */
  float m_shrDedxU;      /**< dE/dx of the leading shower on the U plane with the 1x4 cm box method */
  float m_shrDedxYCali;  /**< Calibrated dE/dx of the leading shower on the Y plane with the 1x4 cm box method */
  float m_shrDedxVCali;  /**< Calibrated dE/dx of the leading shower on the V plane with the 1x4 cm box method */
  float m_shrDedxUCali;  /**< Calibrated dE/dx of the leading shower on the U plane with the 1x4 cm box method */
  float m_shrDistance;   /**< Distance between leading shower vertex and reconstructed neutrino vertex */
  float m_shrScore;      /**< Pandora track score for the leading shower */
  float m_shrTheta;      /**< Reconstructed theta angle for the leading shower */
  float m_shrPhi;        /**< Reconstructed phi angle for the leading shower */
  float m_shrPx;         /**< X component of the reconstructed momentum of the leading shower (in GeV/c) */
  float m_shrPy;         /**< Y component of the reconstructed momentum of the leading shower (in GeV/c) */
  float m_shrPz;         /**< Z component of the reconstructed momentum of the leading shower (in GeV/c) */
  float m_shrPca0;       /**< First eigenvalue of the PCAxis of the leading shower */
  float m_shrPca1;       /**< Second eigenvalue of the PCAxis of the leading shower */
  float m_shrPca2;       /**< Third eigenvalue of the PCAxis of the leading shower */
  float m_shrOpenangle;  /**< Opening angle of the shower */
  float m_shrPidchipr;   /**< Chi2 proton score for the leading shower (with the shower reconstructed as track) */
  float m_shrPidchimu;   /**< Chi2 muon score for the leading shower (with the shower reconstructed as track) */
  float m_shrBragg_p;    /**< Proton Bragg likelihood score for the leading shower (with the shower reconstructed as track) */
  float m_shrBragg_mu;   /**< Muon Bragg likelihood score for the leading shower (with the shower reconstructed as track) */
  float m_shrBragg_mip;  /**< MIP Bragg likelihood score for the leading shower (with the shower reconstructed as track) */
  float m_shrBragg_pion; /**< Pion Bragg likelihood for the leading shower (with the shower reconstructed as track) */
  float m_shrBragg_kaon; /**< Kaon Bragg likelihood for the leading shower (with the shower reconstructed as track) */

  //To Do: track length, range momentum, mcs momentum

  float m_shrBktPurity;       /**< Purity of the leading shower */
  float m_shrBktCompleteness; /**< Completeness of the leading shower */
  float m_shrBktE;            /**< Energy of the MCParticle matched to the leading shower */
  int m_shrBktPdg;            /**< PDG code of the MCParticle matched to the leading shower */

  float m_shrTkfitStartX; /**< Start x coordinate of the leading shower obtained with the track fitting */
  float m_shrTkfitStartY; /**< Start y coordinate of the leading shower obtained with the track fitting */
  float m_shrTkfitStartZ; /**< Start z coordinate of the leading shower obtained with the track fitting */
  float m_shrTkfitphi;    /**< Phi angle of the leading shower obtained with the track fitting */
  float m_shrTkfittheta;  /**< Track angle of the leading shower obtained with the track fitting */
  float m_shrTkfitdedxY;  /**< dE/dx of the leading shower on the Y plane with the track fitting */
  float m_shrTkfitdedxV;  /**< dE/dx of the leading shower on the V plane with the track fitting */
  float m_shrTkfitdedxU;  /**< dE/dx of the leading shower on the U plane with the track fitting */
  uint m_shrTkfitnhitsY;  /**< Number of hits in the 1x4 cm box on the Y plane with the track fitting */
  uint m_shrTkfitnhitsV;  /**< Number of hits in the 1x4 cm box on the V plane with the track fitting */
  uint m_shrTkfitnhitsU;  /**< Number of hits in the 1x4 cm box on the U plane with the track fitting */

  // Fields for the event
  bool m_vtxFiducial;        /**< The vertex is in the fiducial volume */
  bool m_trackFiducial;      /**< Start and end points of all pfp reconstructed as trackare in the fiducial volume */
  uint m_hitsOutfv;          /**< Number of hits of PFParticles outside the fiducial volume */
  float m_containedFraction; /**< Fraction of hits of the PFParticles contained in the fiducial volume */

  // Fields for fiducial volume
  float m_fidvolZstart;    /**< Fiducial volume distance from the start of the TPC on the z axis (default 10 cm) */
  float m_fidvolZend;      /**< Fiducial volume distance from the end of the TPC on the z axis (default 50 cm) */
  float m_fidvolYstart;    /**< Fiducial volume distance from the bottom of the TPC on the y axis (default 10 cm) */
  float m_fidvolYend;      /**< Fiducial volume distance from the top of the TPC on the y axis (default 10 cm) */
  float m_fidvolXstart;    /**< Fiducial volume distance from the start of the TPC on the x axis (default 10 cm) */
  float m_fidvolXend;      /**< Fiducial volume distance from the end of the TPC on the x axis (default 10 cm) */
  float m_fidVolContained; /**< Distance from the edge to be considered contained (default 10 cm) */

  // Fcl parameters
  float m_trkScore; /**< Threshold on the Pandora track score for electron candidate (default 0.9) */
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
CCincSelection::CCincSelection(const fhicl::ParameterSet &pset)
{
  m_trkScore = 0.9;
  configure(pset);
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void CCincSelection::configure(fhicl::ParameterSet const &pset)
{
  m_trkScore = pset.get<unsigned int>("ElectronMaxTrackScore");
  /*
  fCLSproducer = pset.get<art::InputTag>("CLSproducer", "pandora");
  fTRKproducer = pset.get<art::InputTag>("TRKproducer", "pandoraTrack");
  fTRKproducerTrkFit = pset.get<art::InputTag>("TRKproducerTrkFit", "shrreco3dKalmanShower");

  fPIDproducer = pset.get<art::InputTag>("PIDproducer", "pandoraTrackcalipid");
  fCALproducer = pset.get<art::InputTag>("CALproducer", "pandoraTrackcali");
  fCALproducerTrkFit = pset.get<art::InputTag>("CALproducerTrkFit", "shrreco3dKalmanShowercali");
  fHproducer = pset.get<art::InputTag>("Hproducer", "gaushit");
  fMCRproducer = pset.get<art::InputTag>("MCRproducer", "mcreco");
  fBacktrackTag = pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
  fMCPproducer = pset.get<art::InputTag>("MCPproducer", "largeant");
  */
  m_fidvolZstart = pset.get<float>("FidvolZstart", 10);
  m_fidvolZend = pset.get<float>("FidvolZend", 50);
  m_fidvolYstart = pset.get<float>("FidvolYstart", 10);
  m_fidvolYend = pset.get<float>("FidvolYend", 10);
  m_fidvolXstart = pset.get<float>("FidvolXstart", 10);
  m_fidvolXend = pset.get<float>("FidvolXend", 10);
  m_fidVolContained = pset.get<float>("ContainmentBorder", 10);
}

DEFINE_ART_CLASS_TOOL(CCincSelection)
} // namespace selection

#endif
