#ifndef SELECTION_CCINC_CXX
#define SELECTION_CCINC_CXX

#include <iostream>
#include "SelectionToolBase.h"

#include "../CommonDefs/Typedefs.h"
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/CalibrationFuncs.h"
#include "../CommonDefs/PIDFuncs.h"

#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

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
   void setBranches(TTree *_tree);

   /**
     * @brief reset ttree branches
     */
   void resetTTree(TTree *_tree);

   /**
     * @brief Check if point is inside fiducial volume
     *
     * @param x array of coordinates
     */
   bool isFiducial(const double x[3]) const;
   bool isContained(const double x[3], float tolerance) const;

   /**
     *  @brief  Fill the tree for the electron candidate.
     *
     *  @param  pfp_pxy ProxyPfpElem_t corresponding to the electron candidate.
     *  @return 1 if succesful, 0 if failure.
     */
   bool FillElectronCandidate(art::Event const &e,
                              const ProxyPfpElem_t &pfp_pxy);

private:
   // Fields for electron candidate
   bool m_electron_candidate; /**< Is there a candidate for the electron */
   size_t m_shrPfpId;         /**< Index of the leading shower in the PFParticle vector */
   float m_shrDistance;       /**< Distance between leading shower vertex and reconstructed neutrino vertex */
   size_t m_shrHits;          /**< Number of hits of the candidate electron shower */
   size_t m_shrHitsU;         /**< Number of hits of the candidate electron shower, U plane */
   size_t m_shrHitsV;         /**< Number of hits of the candidate electron shower, V plane */
   size_t m_shrHitsY;         /**< Number of hits of the candidate electron shower, Y plane */
   size_t m_shrSps;           /**< Number of spacepoints of the candidate electron shower */

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

   float m_shrPidchipr;     /**< Chi2 proton score for the leading shower (with the shower reconstructed as track) */
   float m_shrPidchimu;     /**< Chi2 muon score for the leading shower (with the shower reconstructed as track) */
   float m_shrBragg_p;      /**< Proton Bragg likelihood score for the leading shower (with the shower reconstructed as track) */
   float m_shrBragg_mu;     /**< Muon Bragg likelihood score for the leading shower (with the shower reconstructed as track) */
   float m_shrBragg_mip;    /**< MIP Bragg likelihood score for the leading shower (with the shower reconstructed as track) */
   float m_shrTrkLength;    /**< Length (with the shower reconstructed as track) */
   bool m_shrTrkContained;  /**< End point is contained (with the shower reconstructed as track) */
   float m_shrRangeMomMuon; /**< Muon momentum from range (with the shower reconstructed as track) */
   float m_shrMCSMomMuon;   /**< Muon momentum from MCS (with the shower reconstructed as track) */

   //float m_shrBktPurity;       /**< Purity of the leading shower */
   //float m_shrBktCompleteness; /**< Completeness of the leading shower */
   //float m_shrBktE;            /**< Energy of the MCParticle matched to the leading shower */
   //int m_shrBktPdg;            /**< PDG code of the MCParticle matched to the leading shower */

   //float m_shrTkfitStartX; /**< Start x coordinate of the leading shower obtained with the track fitting */
   //float m_shrTkfitStartY; /**< Start y coordinate of the leading shower obtained with the track fitting */
   //float m_shrTkfitStartZ; /**< Start z coordinate of the leading shower obtained with the track fitting */
   //float m_shrTkfitphi;    /**< Phi angle of the leading shower obtained with the track fitting */
   //float m_shrTkfittheta;  /**< Track angle of the leading shower obtained with the track fitting */
   //float m_shrTkfitdedxY;  /**< dE/dx of the leading shower on the Y plane with the track fitting */
   //float m_shrTkfitdedxV;  /**< dE/dx of the leading shower on the V plane with the track fitting */
   //float m_shrTkfitdedxU;  /**< dE/dx of the leading shower on the U plane with the track fitting */
   //uint m_shrTkfitnhitsY;  /**< Number of hits in the 1x4 cm box on the Y plane with the track fitting */
   //uint m_shrTkfitnhitsV;  /**< Number of hits in the 1x4 cm box on the V plane with the track fitting */
   //uint m_shrTkfitnhitsU;  /**< Number of hits in the 1x4 cm box on the U plane with the track fitting */

   // Fields for the event
   //bool m_vtxFiducial;        /**< The vertex is in the fiducial volume */
   //bool m_trackFiducial;      /**< Start and end points of all pfp reconstructed as trackare in the fiducial volume */
   //uint m_hitsOutfv;          /**< Number of hits of PFParticles outside the fiducial volume */
   //float m_containedFraction; /**< Fraction of hits of the PFParticles contained in the fiducial volume */

   // Fcl parameters
   float m_trkScore; /**< Threshold on the Pandora track score for electron candidate (default 0.9) */
   // Fields for fiducial volume
   float m_fidvolZstart;    /**< Fiducial volume distance from the start of the TPC on the z axis (default 10 cm) */
   float m_fidvolZend;      /**< Fiducial volume distance from the end of the TPC on the z axis (default 50 cm) */
   float m_fidvolYstart;    /**< Fiducial volume distance from the bottom of the TPC on the y axis (default 10 cm) */
   float m_fidvolYend;      /**< Fiducial volume distance from the top of the TPC on the y axis (default 10 cm) */
   float m_fidvolXstart;    /**< Fiducial volume distance from the start of the TPC on the x axis (default 10 cm) */
   float m_fidvolXend;      /**< Fiducial volume distance from the end of the TPC on the x axis (default 10 cm) */
   float m_fidVolContained; /**< Distance from the edge to be considered contained (default 10 cm) */

   art::InputTag fCLSproducer;
   art::InputTag fTRKproducer;
   art::InputTag fTRKproducerTrkFit;
   art::InputTag fPIDproducer;

   trkf::TrackMomentumCalculator m_trkmom;
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
   m_trkScore = pset.get<float>("ElectronMaxTrackScore", 0.9);

   fCLSproducer = pset.get<art::InputTag>("CLSproducer", "pandora");
   fTRKproducer = pset.get<art::InputTag>("TRKproducer", "pandoraTrack");
   fTRKproducerTrkFit = pset.get<art::InputTag>("TRKproducerTrkFit", "shrreco3dKalmanShower");
   fPIDproducer = pset.get<art::InputTag>("PIDproducer", "pandoraTrackcalipid");

   /*
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

void CCincSelection::resetTTree(TTree *_tree)
{
   m_electron_candidate = false;
   m_shrPfpId = 0;
   m_shrDistance = std::numeric_limits<float>::lowest();

   m_shrHits = 0;
   m_shrHitsU = 0;
   m_shrHitsV = 0;
   m_shrHitsY = 0;
   m_shrSps = 0;

   m_shrEnergy = 0;
   m_shrEnergyCali = 0;

   m_shrStartX = std::numeric_limits<float>::lowest();
   m_shrStartY = std::numeric_limits<float>::lowest();
   m_shrStartZ = std::numeric_limits<float>::lowest();

   m_shrDedxY = std::numeric_limits<float>::lowest();
   m_shrDedxV = std::numeric_limits<float>::lowest();
   m_shrDedxU = std::numeric_limits<float>::lowest();
   m_shrDedxYCali = std::numeric_limits<float>::lowest();
   m_shrDedxVCali = std::numeric_limits<float>::lowest();
   m_shrDedxUCali = std::numeric_limits<float>::lowest();

   m_shrScore = -1;
   m_shrTheta = std::numeric_limits<float>::lowest();
   m_shrPhi = std::numeric_limits<float>::lowest();

   m_shrPx = std::numeric_limits<float>::lowest();
   m_shrPy = std::numeric_limits<float>::lowest();
   m_shrPz = std::numeric_limits<float>::lowest();
   m_shrPca0 = std::numeric_limits<float>::lowest();
   m_shrPca1 = std::numeric_limits<float>::lowest();
   m_shrPca2 = std::numeric_limits<float>::lowest();
   m_shrOpenangle = std::numeric_limits<float>::lowest();

   m_shrPidchipr = std::numeric_limits<float>::lowest();
   m_shrPidchimu = std::numeric_limits<float>::lowest();
   m_shrBragg_p = std::numeric_limits<float>::lowest();
   m_shrBragg_mu = std::numeric_limits<float>::lowest();
   m_shrBragg_mip = std::numeric_limits<float>::lowest();

   m_shrTrkLength = std::numeric_limits<float>::lowest();
   m_shrTrkContained = false;
   m_shrRangeMomMuon = std::numeric_limits<float>::lowest();
   m_shrMCSMomMuon = std::numeric_limits<float>::lowest();
}

void CCincSelection::setBranches(TTree *_tree)
{
   _tree->Branch("e_candidate", &m_electron_candidate, "e_candidate/O");
   _tree->Branch("e_candidate_id", &m_shrPfpId, "e_candidate_id/i");
   _tree->Branch("e_nu_distance", &m_shrDistance, "e_nu_distance/F");

   _tree->Branch("e_hits_total", &m_shrHits, "e_hits_total/i");
   _tree->Branch("e_hits_u", &m_shrHitsU, "e_hits_u/i");
   _tree->Branch("e_hits_v", &m_shrHitsV, "e_hits_v/i");
   _tree->Branch("e_hits_y", &m_shrHitsY, "e_hits_y/i");
   _tree->Branch("e_sps_total", &m_shrSps, "e_sps_total/i");

   _tree->Branch("e_energy", &m_shrEnergy, "e_energy/F");
   _tree->Branch("e_energy_cali", &m_shrEnergyCali, "e_energy_cali/F");

   _tree->Branch("e_start_x", &m_shrStartX, "e_start_x/F");
   _tree->Branch("e_start_y", &m_shrStartY, "e_start_y/F");
   _tree->Branch("e_start_z", &m_shrStartZ, "e_start_z/F");

   _tree->Branch("e_dedx_y", &m_shrDedxY, "e_dedx_y/F");
   _tree->Branch("e_dedx_v", &m_shrDedxV, "e_dedx_v/F");
   _tree->Branch("e_dedx_u", &m_shrDedxU, "e_dedx_u/F");
   _tree->Branch("e_dedx_y_cali", &m_shrDedxYCali, "e_dedx_y_cali/F");
   _tree->Branch("e_dedx_v_cali", &m_shrDedxVCali, "e_dedx_v_cali/F");
   _tree->Branch("e_dedx_u_cali", &m_shrDedxUCali, "e_dedx_u_cali/F");

   _tree->Branch("e_trackscore", &m_shrScore, "e_dedx_y/F");
   _tree->Branch("e_phi", &m_shrTheta, "e_dedx_y/F");
   _tree->Branch("e_theta", &m_shrPhi, "e_dedx_y/F");

   _tree->Branch("e_px", &m_shrPx, "e_px/F");
   _tree->Branch("e_py", &m_shrPy, "e_py/F");
   _tree->Branch("e_pz", &m_shrPz, "e_pz/F");

   _tree->Branch("e_pca0", &m_shrPca0, "e_pca0/F");
   _tree->Branch("e_pca1", &m_shrPca1, "e_pca1/F");
   _tree->Branch("e_pca2", &m_shrPca2, "e_pca2/F");
   _tree->Branch("e_openangle", &m_shrOpenangle, "e_openangle/F");

   _tree->Branch("e_pid_chi_proton", &m_shrPidchipr, "e_pid_chi_proton/F");
   _tree->Branch("e_pid_chi_muon", &m_shrPidchimu, "e_pid_chi_muon/F");
   _tree->Branch("e_pid_bragg_proton", &m_shrBragg_p, "e_pid_bragg_proton/F");
   _tree->Branch("e_pid_bragg_muon", &m_shrBragg_mu, "e_pid_bragg_muon/F");
   _tree->Branch("e_pid_bragg_mip", &m_shrBragg_mip, "e_pid_bragg_mip/F");

   _tree->Branch("e_track_length", &m_shrTrkLength, "e_track_length/F");
   _tree->Branch("e_track_contained", &m_shrTrkContained, "e_track_contained/O");
   _tree->Branch("e_track_range_momentum", &m_shrRangeMomMuon, "e_track_range_momentum/F");
   _tree->Branch("e_track_mcs_momentum", &m_shrRangeMomMuon, "e_track_mcs_momentum/F");
}

DEFINE_ART_CLASS_TOOL(CCincSelection)
} // namespace selection

#endif
