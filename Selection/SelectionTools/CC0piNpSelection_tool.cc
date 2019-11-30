#ifndef SELECTION_CC0PINPSELECTION_CXX
#define SELECTION_CC0PINPSELECTION_CXX

#include <iostream>
#include "SelectionToolBase.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "../CommonDefs/Typedefs.h"
#include "../CommonDefs/PIDFuncs.h"
#include "../CommonDefs/Containment.h"
#include "../CommonDefs/TrackFitterFunctions.h"
#include "../CommonDefs/CalibrationFuncs.h"
#include "../CommonDefs/PFPHitDistance.h"
#include "../CommonDefs/ProximityClustering.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/ShowerBranchTagger.h"
//#include "../CommonDefs/PIDFuncs.h"
//#include "../CommonDefs/LLR_PID.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/SpacePoint.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibService.h"

namespace selection
{

/**
 *  Selection of electron neutrinos with 0 pions and at least on proton in the final state
 *  Author: Stefano Roberto Soleti (srsoleti@fnal.gov)
 */

class CC0piNpSelection : public SelectionToolBase
{

public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset List of parameters in the FCL file
     */
    CC0piNpSelection(const fhicl::ParameterSet &pset);

    /**
     *  @brief  Destructor
     */
    ~CC0piNpSelection(){};

    // provide for initialization
    void configure(fhicl::ParameterSet const &pset);

    /**
     * @brief Selection function
     *
     * @param e art Event
     * @param pfp_pxy_v Proxy of PFParticle vector in the slice
     */
    bool selectEvent(art::Event const &e,
                     const std::vector<ProxyPfpElem_t> &pfp_pxy_v);

    /**
     * @brief Set branches for TTree
     *
     * @param _tree ROOT TTree with the selection information
     */
    void setBranches(TTree *_tree);

    /**
     * @brief Reset TTree branches
     *
     * @param _tree ROOT TTree with the selection information
     */
    void resetTTree(TTree *_tree);

private:

    const trkf::TrackMomentumCalculator _trkmom;
    const trkf::TrajectoryMCSFitter mcsfitter;

  //searchingfornues::LLRPID LLR_PID;

    TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
    TParticlePDG *electron = TDatabasePDG::Instance()->GetParticle(11);
    TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);
    TParticlePDG *pion = TDatabasePDG::Instance()->GetParticle(211);
    TParticlePDG *kaon = TDatabasePDG::Instance()->GetParticle(321);

    float fTrkShrscore;  /**< Threshold on the Pandora track score (default 0.5) */
    float fFidvolZstart; /**< Fiducial volume distance from the start of the TPC on the z axis (default 10 cm) */
    float fFidvolZend;   /**< Fiducial volume distance from the end of the TPC on the z axis (default 50 cm) */
    float fFidvolYstart; /**< Fiducial volume distance from the bottom of the TPC on the y axis (default 15 cm) */
    float fFidvolYend;   /**< Fiducial volume distance from the top of the TPC on the y axis (default 15 cm) */
    float fFidvolXstart; /**< Fiducial volume distance from the start of the TPC on the x axis (default 10 cm) */
    float fFidvolXend;   /**< Fiducial volume distance from the end of the TPC on the x axis (default 10 cm) */

    float fdEdxcmSkip, fdEdxcmLen; // cm from vtx to skip for dE/dx and dE/dx segment length, respectively
    bool fSaveMoreDedx;            // save more track fit dedx definitions (currently with 0.5 and 1.0 cm gap at start)
    bool fLocaldEdx;               // use local dE/dx from calorimetry (true) or globally convert Q -> MeV (false)

    art::InputTag fCLSproducer;
    art::InputTag fPIDproducer;
    art::InputTag fTRKproducer;
    art::InputTag fTRKproducerTrkFit;
    art::InputTag fBacktrackTag;
    art::InputTag fHproducer;
    art::InputTag fMCRproducer;
    art::InputTag fMCPproducer;
    art::InputTag fCALproducer;
    art::InputTag fCALproducerTrkFit;

  float _wire2cm, _time2cm;


    unsigned int _n_showers_contained; /**< Number of showers with a starting point within the fiducial volume */
    unsigned int _n_tracks_contained;  /**< Number of tracks fully contained in the fiducial volume */
    unsigned int _shr_hits_max;        /**< Number of hits of the leading shower */
    unsigned int _trk_hits_max;        /**< Number of hits of the leading track */
    unsigned int _shr_hits_tot;        /**< Total number of shower hits */
    unsigned int _trk_hits_tot;        /**< Total number of track hits */
    unsigned int _trk_hits_y_tot;      /**< Total number of track hits on the Y plane */
    unsigned int _trk_hits_v_tot;      /**< Total number of track hits on the V plane */
    unsigned int _trk_hits_u_tot;      /**< Total number of track hits on the U plane */
    unsigned int _shr_hits_y_tot;      /**< Total number of shower hits on the Y plane */
    unsigned int _shr_hits_v_tot;      /**< Total number of shower hits on the V plane */
    unsigned int _shr_hits_u_tot;      /**< Total number of shower hits on the U plane */
    float _shr_energy;                 /**< Energy of the shower with the largest number of hits (in GeV) */
    float _shr_energy_tot;             /**< Sum of the energy of the showers (in GeV) */
    float _shr_energy_cali;            /**< Energy of the calibrated shower with the largest number of hits (in GeV) */
    float _shr_energy_tot_cali;        /**< Sum of the energy of the calibrated showers (in GeV) */
    float _shr_dedx_Y;                 /**< dE/dx of the leading shower on the Y plane with the 1x4 cm box method */
    float _shr_dedx_V;                 /**< dE/dx of the leading shower on the V plane with the 1x4 cm box method */
    float _shr_dedx_U;                 /**< dE/dx of the leading shower on the U plane with the 1x4 cm box method */
    float _shr_dedx_Y_cali;            /**< Calibrated dE/dx of the leading shower on the Y plane with the 1x4 cm box method */
    float _shr_dedx_V_cali;            /**< Calibrated dE/dx of the leading shower on the V plane with the 1x4 cm box method */
    float _shr_dedx_U_cali;            /**< Calibrated dE/dx of the leading shower on the U plane with the 1x4 cm box method */
    float _shr_distance;               /**< Distance between leading shower vertex and reconstructed neutrino vertex */
    float _tksh_distance;              /**< Distance between leading shower vertex and longest track vertex */
    float _tksh_angle;                 /**< Angle between leading shower vertex and longest track vertex */
    float _tksh_angle_muon;            /**< Angle between leading shower vertex and longest track vertex */
    float _shr_score;                  /**< Pandora track score for the leading shower */
    float _shr_theta;                  /**< Reconstructed theta angle for the leading shower */
    float _shr_phi;                    /**< Reconstructed phi angle for the leading shower */
    float _shr_px;                     /**< X component of the reconstructed momentum of the leading shower (in GeV/c) */
    float _shr_py;                     /**< Y component of the reconstructed momentum of the leading shower (in GeV/c) */
    float _shr_pz;                     /**< Z component of the reconstructed momentum of the leading shower (in GeV/c) */
    float _shr_pca_0;                  /**< First eigenvalue of the PCAxis of the leading shower */
    float _shr_pca_1;                  /**< Second eigenvalue of the PCAxis of the leading shower */
    float _shr_pca_2;                  /**< Third eigenvalue of the PCAxis of the leading shower */
    float _shr_openangle;              /**< Opening angle of the shower */
    float _shr_pidchipr;               /**< Chi2 proton score for the leading shower (with the shower reconstructed as track) */
    float _shr_pidchimu;               /**< Chi2 muon score for the leading shower (with the shower reconstructed as track) */
    float _shr_bragg_p;                /**< Proton Bragg likelihood score for the leading shower (with the shower reconstructed as track) */
    float _shr_bragg_mu;               /**< Muon Bragg likelihood score for the leading shower (with the shower reconstructed as track) */
    float _shr_bragg_mip;              /**< MIP Bragg likelihood score for the leading shower (with the shower reconstructed as track) */
    float _shr_bragg_pion;             /**< Pion Bragg likelihood for the leading shower (with the shower reconstructed as track) */
    float _shr_bragg_kaon;             /**< Kaon Bragg likelihood for the leading shower (with the shower reconstructed as track) */

    size_t _shr_pfp_id; /**< Index of the leading shower in the PFParticle vector */

    float _trk_len;             /**< Length of the longest track */
    float _trk_energy;          /**< Energy of the longest track assuming it's a proton and using stopping power in LAr */
    float _trk_energy_muon;     /**< Energy of the longest track assuming it's a muon and using stopping power in LAr */
    float _trk_energy_muon_mcs; /**< Energy of the longest track assuming it's a muon and using MCS */
    float _trk_energy_tot;      /**< Sum of the track energies assuming they are protons and using stopping power in LAr */
    float _trk_energy_muon_tot; /**< Sum of the track energies assuming they are muons and using stopping power in LAr */
    float _trk_distance;        /**< Distance between longest track and reconstructed neutrino vertex */
    float _trk_theta;           /**< Reconstructed theta angle for the longest track */
    float _trk_phi;             /**< Reconstructed phi angle for the longest track */
    size_t _trk_pfp_id;         /**< Index of the longest track in the PFParticle vector */

    float _hits_ratio;     /**< Ratio between hits from showers and total number of hits */
    float _trk_bragg_p;    /**< Proton Bragg likelihood score for the longest track */
    float _trk_bragg_mu;   /**< Muon Bragg likelihood score for the longest track */
    float _trk_bragg_mip;  /**< MIP Bragg likelihood score for the longest track */
    float _trk_bragg_pion; /**< Pion Bragg likelihood score for the longest track */
    float _trk_bragg_kaon; /**< Kaon Bragg likelihood score for the longest track */

    float _trk_pidchipr;       /**< Chi2 proton score for the longest track */
    float _trk_pidchipr_best;  /**< Best Chi2 proton score amongst all the tracks */
    float _trk_pidchipr_worst; /**< Worst Chi2 proton score amongst all the tracks */

    float _trk_pidchimu;       /**< Chi2 muon score for the longest track */
    float _trk_pidchimu_best;  /**< Best Chi2 muon score amongst all the tracks */
    float _trk_pidchimu_worst; /**< Worst Chi2 muon score amongst all the tracks */

    float _trk_pida;  /**< PIDA score for the longest track */
    float _trk_score; /**< Pandora track score for the longest track */

    float _pt;             /**< Total reconstructed transverse momentum, assuming all the tracks are protons and all the showers are electrons */
    float _pt_assume_muon; /**< Total reconstructed transverse momentum, assuming all the tracks are muons and all the showers are electrons */
    float _p;              /**< Total reconstructed momentum, assuming all the tracks are protons and all the showers are electrons */
    float _p_assume_muon;  /**< Total reconstructed momentum, assuming all the tracks are muons and all the showers are electrons */

    float _shr_bkt_purity;       /**< Purity of the leading shower */
    float _shr_bkt_completeness; /**< Completeness of the leading shower */
    float _shr_bkt_E;            /**< Energy of the MCParticle matched to the leading shower */
    int _shr_bkt_pdg;            /**< PDG code of the MCParticle matched to the leading shower */
    float _trk_bkt_purity;       /**< Purity of the longest track */
    float _trk_bkt_completeness; /**< Completeness of the longest track */
    float _trk_bkt_E;            /**< Energy of the MCParticle matched to the longest track */
    int _trk_bkt_pdg;            /**< PDG code of the MCParticle matched to the longest track */

    float _matched_E; /**< Total kinetic energy of the MCParticles matched to PFParticles */

    int _shrsubclusters0, _shrsubclusters1, _shrsubclusters2; /**< in how many sub-clusters can the shower be broken based on proximity clustering? */
    float _shrclusfrac0, _shrclusfrac1, _shrclusfrac2;        /**< what fraction of the total charge does the dominant shower sub-cluster carry? */
    float _shrclusdir0, _shrclusdir1, _shrclusdir2;        /**< 2D charge-weighted direction of shower hits calculated from neutrino vertex w.r.t. vertical in plane */
    float _trkshrhitdist0, _trkshrhitdist1, _trkshrhitdist2;  /**< distance between hits of shower and track in 2D on each palne based on hit-hit distances */

  float _shrmoliereavg; /**< avg of moliere angle */
  float _shrmoliererms; /**< rms of moliere angle */

    int _shr_tkfit_npoints;      // number of points associated to shower fitted track
    int _shr_tkfit_npointsvalid; // number of VALID points associated to shower fitted track

    float _shr_trkfitmedangle; // median angle for first N cm (default 10) for track-fitter track

    float _shr_tkfit_start_x; /**< Start x coordinate of the leading shower obtained with the track fitting */
    float _shr_tkfit_start_y; /**< Start y coordinate of the leading shower obtained with the track fitting */
    float _shr_tkfit_start_z; /**< Start z coordinate of the leading shower obtained with the track fitting */
    float _shr_start_x;       /**< Start x coordinate of the leading shower */
    float _shr_start_y;       /**< Start y coordinate of the leading shower */
    float _shr_start_z;       /**< Start z coordinate of the leading shower */
    float _shr_tkfit_phi;     /**< Phi angle of the leading shower obtained with the track fitting */
    float _shr_tkfit_theta;   /**< Track angle of the leading shower obtained with the track fitting */
    float _shr_tkfit_dedx_Y;  /**< dE/dx of the leading shower on the Y plane with the track fitting */
    float _shr_tkfit_dedx_V;  /**< dE/dx of the leading shower on the V plane with the track fitting */
    float _shr_tkfit_dedx_U;  /**< dE/dx of the leading shower on the U plane with the track fitting */

    float _shr_tkfit_dedx_Y_alt;  /**< dE/dx of the leading shower on the Y plane with the track fitting  [calculated using XYZ instead of RR]  */
    float _shr_tkfit_dedx_V_alt;  /**< dE/dx of the leading shower on the V plane with the track fitting  [calculated using XYZ instead of RR]  */
    float _shr_tkfit_dedx_U_alt;  /**< dE/dx of the leading shower on the U plane with the track fitting  [calculated using XYZ instead of RR]  */

    unsigned int _shr_tkfit_nhits_Y; /**< Number of hits in the 1x4 cm box on the Y plane with the track fitting */
    unsigned int _shr_tkfit_nhits_V; /**< Number of hits in the 1x4 cm box on the V plane with the track fitting */
    unsigned int _shr_tkfit_nhits_U; /**< Number of hits in the 1x4 cm box on the U plane with the track fitting */

    unsigned int _shr_tkfit_nhits_Y_alt; /**< Number of hits in the 1x4 cm box on the Y plane with the track fitting [calculated using XYZ instead of RR] */
    unsigned int _shr_tkfit_nhits_V_alt; /**< Number of hits in the 1x4 cm box on the V plane with the track fitting  [calculated using XYZ instead of RR] */
    unsigned int _shr_tkfit_nhits_U_alt; /**< Number of hits in the 1x4 cm box on the U plane with the track fitting [calculated using XYZ instead of RR] */

    float _shr_tkfit_2cm_dedx_Y;         /**< dE/dx of the leading shower on the Y plane with the track fitting, use first 2 cm */
    float _shr_tkfit_2cm_dedx_V;         /**< dE/dx of the leading shower on the V plane with the track fitting, use first 2 cm */
    float _shr_tkfit_2cm_dedx_U;         /**< dE/dx of the leading shower on the U plane with the track fitting, use first 2 cm */
    unsigned int _shr_tkfit_2cm_nhits_Y; /**< Number of hits in the 1x4 cm box on the Y plane with the track fitting, use first 2 cm */
    unsigned int _shr_tkfit_2cm_nhits_V; /**< Number of hits in the 1x4 cm box on the V plane with the track fitting, use first 2 cm */
    unsigned int _shr_tkfit_2cm_nhits_U; /**< Number of hits in the 1x4 cm box on the U plane with the track fitting, use first 2 cm */

    float _shr_tkfit_gap05_dedx_Y;         /**< dE/dx of the leading shower on the Y plane with the track fitting, skip first 5 mm */
    float _shr_tkfit_gap05_dedx_V;         /**< dE/dx of the leading shower on the V plane with the track fitting, skip first 5 mm */
    float _shr_tkfit_gap05_dedx_U;         /**< dE/dx of the leading shower on the U plane with the track fitting, skip first 5 mm */
    unsigned int _shr_tkfit_gap05_nhits_Y; /**< Number of hits in the 1x4 cm box on the Y plane with the track fitting, skip first 5 mm */
    unsigned int _shr_tkfit_gap05_nhits_V; /**< Number of hits in the 1x4 cm box on the V plane with the track fitting, skip first 5 mm */
    unsigned int _shr_tkfit_gap05_nhits_U; /**< Number of hits in the 1x4 cm box on the U plane with the track fitting, skip first 5 mm */

    float _shr_tkfit_gap10_dedx_Y;         /**< dE/dx of the leading shower on the Y plane with the track fitting, skip first 10 mm */
    float _shr_tkfit_gap10_dedx_V;         /**< dE/dx of the leading shower on the V plane with the track fitting, skip first 10 mm */
    float _shr_tkfit_gap10_dedx_U;         /**< dE/dx of the leading shower on the U plane with the track fitting, skip first 10 mm */
    unsigned int _shr_tkfit_gap10_nhits_Y; /**< Number of hits in the 1x4 cm box on the Y plane with the track fitting, skip first 10 mm */
    unsigned int _shr_tkfit_gap10_nhits_V; /**< Number of hits in the 1x4 cm box on the V plane with the track fitting, skip first 10 mm */
    unsigned int _shr_tkfit_gap10_nhits_U; /**< Number of hits in the 1x4 cm box on the U plane with the track fitting, skip first 10 mm */

    unsigned int _hits_outfv;      /**< Number of hits of PFParticles outside the fiducial volume */
    float _contained_fraction;     /**< Fraction of hits of the PFParticles contained in the fiducial volume */
    float _sps_contained_fraction; /**< Fraction of SpacePoints of the PFParticles contained in the fiducial volume */

    float _trk_energy_hits_tot; /**< Sum of the energy of the tracks obtained with the deposited charge */

    unsigned int _total_hits_y; /**< Total number of hits on the Y plane */
    float _extra_energy_y;      /**< Total energy of the unclustered hits on the Y plane */
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
CC0piNpSelection::CC0piNpSelection(const fhicl::ParameterSet &pset)
    : mcsfitter(fhicl::Table<trkf::TrajectoryMCSFitter::Config>(pset.get<fhicl::ParameterSet>("mcsfitmu")))
{
    configure(pset);

    // get detector specific properties
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    _wire2cm = geom->WirePitch(0,0,0);
    _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );

    /*
  LLR_PID.set_dedx_binning(0, dedx_num_bins_pl_0, dedx_edges_pl_0);
  std::vector<size_t> parameters_num_bins_0 = {parameter_0_num_bins_pl_0, parameter_1_num_bins_pl_0};
  std::vector<std::vector<float>> parameters_bin_edges_0 = {parameter_0_edges_pl_0, parameter_1_edges_pl_0};
  LLR_PID.set_par_binning(0, parameters_num_bins_0, parameters_bin_edges_0);
  LLR_PID.set_lookup_tables(0, dedx_pdf_pl_0);

  LLR_PID.set_dedx_binning(1, dedx_num_bins_pl_1, dedx_edges_pl_1);
  std::vector<size_t> parameters_num_bins_1 = {parameter_0_num_bins_pl_1, parameter_1_num_bins_pl_1};
  std::vector<std::vector<float>> parameters_bin_edges_1 = {parameter_0_edges_pl_1, parameter_1_edges_pl_1};
  LLR_PID.set_par_binning(1, parameters_num_bins_1, parameters_bin_edges_1);
  LLR_PID.set_lookup_tables(1, dedx_pdf_pl_1);

  LLR_PID.set_dedx_binning(2, dedx_num_bins_pl_2, dedx_edges_pl_2);
  std::vector<size_t> parameters_num_bins_2 = {parameter_0_num_bins_pl_2, parameter_1_num_bins_pl_2};
  std::vector<std::vector<float>> parameters_bin_edges_2 = {parameter_0_edges_pl_2, parameter_1_edges_pl_2};
  LLR_PID.set_par_binning(2, parameters_num_bins_2, parameters_bin_edges_2);
  LLR_PID.set_lookup_tables(2, dedx_pdf_pl_2);
    */
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void CC0piNpSelection::configure(fhicl::ParameterSet const &pset)
{
    fTrkShrscore = pset.get<float>("TrkShrscore", 0.5);
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

    fFidvolZstart = pset.get<float>("FidvolZstart", 10);
    fFidvolZend = pset.get<float>("FidvolZend", 50);
    fFidvolYstart = pset.get<float>("FidvolYstart", 15);
    fFidvolYend = pset.get<float>("FidvolYend", 15);
    fFidvolXstart = pset.get<float>("FidvolXstart", 10);
    fFidvolXend = pset.get<float>("FidvolXend", 10);

    fdEdxcmSkip = pset.get<float>("dEdxcmSkip", 0.0);     // how many cm to skip @ vtx for dE/dx calculation
    fdEdxcmLen = pset.get<float>("dEdxcmLen", 4.0);       // how long the dE/dx segment should be
    fSaveMoreDedx = pset.get<bool>("SaveMoreDedx", true); // save additional track fit dedx definitions
    fLocaldEdx = pset.get<bool>("LocaldEdx", true);       // use dE/dx from calo?
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
bool CC0piNpSelection::selectEvent(art::Event const &e,
                                   const std::vector<ProxyPfpElem_t> &pfp_pxy_v)
{
    ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
                                                                                          proxy::withAssociated<recob::Hit>(fCLSproducer));
    searchingfornues::ProxyPIDColl_t const &pid_proxy = proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer,
                                                                                                        proxy::withAssociated<anab::ParticleID>(fPIDproducer));

    art::ValidHandle<std::vector<recob::PFParticle>> inputPfParticle = e.getValidHandle<std::vector<recob::PFParticle>>(fCLSproducer);
    art::ValidHandle<std::vector<recob::Slice>> inputSlice = e.getValidHandle<std::vector<recob::Slice>>(fCLSproducer);
    art::ValidHandle<std::vector<recob::SpacePoint>> inputSpacePoint = e.getValidHandle<std::vector<recob::SpacePoint>>(fCLSproducer);

    auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSlice, e, fCLSproducer));
    auto assocSlice = std::unique_ptr<art::FindManyP<recob::Slice>>(new art::FindManyP<recob::Slice>(inputPfParticle, e, fCLSproducer));
    auto assocSpacePoint = std::unique_ptr<art::FindManyP<recob::SpacePoint>>(new art::FindManyP<recob::SpacePoint>(inputPfParticle, e, fCLSproducer));

    auto assocSpacePointHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSpacePoint, e, fCLSproducer));

    art::FindManyP<larpandoraobj::PFParticleMetadata> pfPartToMetadataAssoc(inputPfParticle, e, fCLSproducer);
    searchingfornues::ProxyCaloColl_t const *tkcalo_proxy = NULL;
    if (fTRKproducerTrkFit != "")
    {
        tkcalo_proxy = new searchingfornues::ProxyCaloColl_t(proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducerTrkFit, proxy::withAssociated<anab::Calorimetry>(fCALproducerTrkFit)));
    }
    // load backtrack information
    std::vector<searchingfornues::BtPart> btparts_v;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
    if (!fData)
    {
        const std::vector<sim::MCShower> &inputMCShower = *(e.getValidHandle<std::vector<sim::MCShower>>(fMCRproducer));
        const std::vector<sim::MCTrack> &inputMCTrack = *(e.getValidHandle<std::vector<sim::MCTrack>>(fMCRproducer));
        art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
        assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
        btparts_v = searchingfornues::initBacktrackingParticleVec(inputMCShower, inputMCTrack, *inputHits, assocMCPart);
    }

    for (unsigned int inpf = 0; inpf < inputPfParticle->size(); ++inpf)
    {
        art::Ptr<recob::PFParticle> npfp(inputPfParticle, inpf);
        bool isTheNeutrino = false;

        const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> &pfParticleMetadataList(pfPartToMetadataAssoc.at(inpf));
        if (!pfParticleMetadataList.empty())
        {
            for (unsigned int j = 0; j < pfParticleMetadataList.size(); ++j)
            {
                const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));

                const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
                for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
                {
                    if (it->first == "IsNeutrino" && it->second == 1)
                        isTheNeutrino = true;
                }
            }
        }
        if (npfp->IsPrimary() == false)
            continue;
        //  recob::Slices are unique for the neutrino slice and for clear cosmics. Ambiguous cosmics may have >1 primary PFParticle per slice
        //(i.e. the slice is not unique in terms of primary PFParticles)
        auto slices = assocSlice->at(npfp.key());
        if (slices.size() != 1)
        {
            std::cout << "WRONG!!! n slices = " << slices.size() << std::endl;
        }
        if (isTheNeutrino == true)
        {
            auto slicehit = assocSliceHit->at(slices[0].key());
            for (unsigned int isl = 0; isl < slicehit.size(); isl++)
            {
                const auto &slhit = slicehit.at(isl);
                int slhit_pl = slhit->WireID().Plane;
                if (slhit_pl == 2)
                {
                    _total_hits_y++;
                    _extra_energy_y += 1.01 * slhit->Integral() * 238 * 23.6e-6 / 0.6 / 1000;
                }
            }
        }
    }

    // START checking if vertex is in the fiducial volume
    double nu_vtx[3] = {};
    TVector3 nuvtx;
    for (size_t i_pfp = 0; i_pfp < pfp_pxy_v.size(); i_pfp++)
    {
        auto PDG = fabs(pfp_pxy_v[i_pfp]->PdgCode());

        if (PDG == 12 || PDG == 14)
        {

            auto vtx = pfp_pxy_v[i_pfp].get<recob::Vertex>();
            if (vtx.size() != 1)
            {
                std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
            }
            else
            {
                vtx.at(0)->XYZ(nu_vtx);
                if (!searchingfornues::isFiducial(nu_vtx,fFidvolXstart,fFidvolYstart,fFidvolZstart,fFidvolXend,fFidvolYend,fFidvolZend))
                {
                    return false;
                }
                nuvtx.SetXYZ(nu_vtx[0], nu_vtx[1], nu_vtx[2]);
            }

            break;
        }
    }// for all PFParticles
    // DONE checking if vertex is in the fiducial volume

    TVector3 total_p;
    TVector3 total_p_mu;
    TVector3 trk_vtx;
    TVector3 shr_vtx;
    TVector3 trk_p;
    TVector3 trk_p_mu;
    TVector3 shr_p;
    std::vector<float> matched_energies;

    // NEXT STEP : find showers
    // all showers will be merged in order to determine the shower energy
    // the largest shower (by number of hits) will be the primry electron
    // used to compute vertex, dE/dx, direction, ...
    size_t sps_fv = 0, sps_all = 0;
    for (size_t i_pfp = 0; i_pfp < pfp_pxy_v.size(); i_pfp++)
    {
        auto const &pfp_pxy = pfp_pxy_v.at(i_pfp);

        auto PDG = fabs(pfp_pxy->PdgCode());

        if (PDG == 12 || PDG == 14)
            continue;

        std::vector<art::Ptr<recob::SpacePoint>> spcpnts = assocSpacePoint->at(i_pfp);
        sps_all += spcpnts.size();
        for (auto &sp : spcpnts)
        {
	  if (searchingfornues::isFiducial(sp->XYZ(),fFidvolXstart,fFidvolYstart,fFidvolZstart,fFidvolXend,fFidvolYend,fFidvolZend))
                sps_fv++;
        }

        auto trkshrscore = searchingfornues::GetTrackShowerScore(pfp_pxy);
        if (trkshrscore < fTrkShrscore)
        {
	  unsigned int shr_hits = 0;

            for (const auto &shr : pfp_pxy.get<recob::Shower>())
            {

                double shr_vertex[3] = {shr->ShowerStart().X(), shr->ShowerStart().Y(), shr->ShowerStart().Z()};
                auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();
                std::vector<art::Ptr<recob::Hit>> hit_v;

                if (!searchingfornues::isFiducial(shr_vertex,fFidvolXstart,fFidvolYstart,fFidvolZstart,fFidvolXend,fFidvolYend,fFidvolZend))
                {
                    for (auto ass_clus : clus_pxy_v)
                    {
                        const auto &clus = clus_proxy[ass_clus.key()];
                        auto clus_hit_v = clus.get<recob::Hit>();
                        _hits_outfv += clus_hit_v.size();
                    }
                    continue;
                }

                _n_showers_contained++;

                auto &pca_pxy_v = pfp_pxy.get<recob::PCAxis>();
                if (pca_pxy_v.size() > 0)
                {
                    _shr_pca_0 = pca_pxy_v[0]->getEigenValues()[0];
                    _shr_pca_1 = pca_pxy_v[0]->getEigenValues()[1];
                    _shr_pca_2 = pca_pxy_v[0]->getEigenValues()[2];
                }

                // store hits for each plane
                std::vector<std::vector<art::Ptr<recob::Hit>>> cluster_hits_v_v(3, std::vector<art::Ptr<recob::Hit>>());

                for (auto ass_clus : clus_pxy_v)
                {
                    // get cluster proxy
                    const auto &clus = clus_proxy[ass_clus.key()];
                    auto clus_hit_v = clus.get<recob::Hit>();
                    shr_hits += clus_hit_v.size();
                    for (const auto &hit : clus_hit_v)
                    {
                        hit_v.push_back(hit);
                    }
                    if (clus->Plane().Plane == 0)
                    {
                        _shr_hits_u_tot += clus_hit_v.size();
                        // gather hits from this plane's cluster
                        for (size_t h = 0; h < clus_hit_v.size(); h++)
                            cluster_hits_v_v[0].push_back(clus_hit_v[h]);
                    }
                    else if (clus->Plane().Plane == 1)
                    {
                        _shr_hits_v_tot += clus_hit_v.size();
                        // gather hits from this plane's cluster
                        for (size_t h = 0; h < clus_hit_v.size(); h++)
                            cluster_hits_v_v[1].push_back(clus_hit_v[h]);
                    }
                    else if (clus->Plane().Plane == 2)
                    {
                        _shr_hits_y_tot += clus_hit_v.size();
                        // gather hits from this plane's cluster
                        for (size_t h = 0; h < clus_hit_v.size(); h++)
                            cluster_hits_v_v[2].push_back(clus_hit_v[h]);
                    }
                } // for all clusters associated to this shower

                for (size_t pl = 0; pl < 3; pl++)
                {
                    int nclus = 0;
                    float hitfracmax = 0.;
		    // store direction of 2Dhits in cluster w.r.t. vertical on plane
		    float clusterdir = 0;
                    std::vector<std::vector<unsigned int>> out_cluster_v;
                    if (cluster_hits_v_v[pl].size())
                    {
		      searchingfornues::cluster(cluster_hits_v_v[pl], out_cluster_v, 2.0, 1.0);
		      clusterdir = searchingfornues::ClusterVtxDirection(nuvtx,cluster_hits_v_v[pl],_wire2cm,_time2cm);
                        // find how many clusters above some # of hit threshold there are
                        // find cluste with largest fraction of all hits
                        float tothits = cluster_hits_v_v[pl].size();
                        for (size_t nc = 0; nc < out_cluster_v.size(); nc++)
                        {
                            auto clus_hit_idx_v = out_cluster_v.at(nc);
                            int nhitclus = clus_hit_idx_v.size();
                            if (nhitclus > 3.)
                                nclus += 1;
                            float hitfrac = nhitclus / tothits;
                            if (hitfrac > hitfracmax)
                                hitfracmax = hitfrac;
                        } // for all sub-clusters
                    }     // if there are any hits on this plane
                    if (pl == 0)
                    {
                        _shrsubclusters0 = nclus;
                        _shrclusfrac0 = hitfracmax;
			_shrclusdir0  = clusterdir;
                    }
                    if (pl == 1)
                    {
                        _shrsubclusters1 = nclus;
                        _shrclusfrac1 = hitfracmax;
			_shrclusdir1  = clusterdir;
                    }
                    if (pl == 2)
                    {
                        _shrsubclusters2 = nclus;
                        _shrclusfrac2 = hitfracmax;
			_shrclusdir2  = clusterdir;
                    }
                } // for all planes

                _shr_energy_tot += shr->Energy()[2] / 1000;

                std::vector<float> cali_corr(3);
                searchingfornues::getCali(spcpnts, *assocSpacePointHit, cali_corr);
                _shr_energy_tot_cali += shr->Energy()[2] / 1000 * cali_corr[2];

                _shr_hits_tot += shr_hits;
		
		// if this is the shower with most hits, take as the main shower
                if (shr_hits > _shr_hits_max)
                {
                    if (!fData)
                    {
                        if (clus_pxy_v.size() != 0)
                        {
                            float purity = 0., completeness = 0.;
                            int ibt = searchingfornues::getAssocBtPart(hit_v, assocMCPart, btparts_v, purity, completeness);
                            if (ibt >= 0)
                            {
                                auto &mcp = btparts_v[ibt];
                                auto PDG = mcp.pdg;
                                auto E = mcp.e;
                                _shr_bkt_pdg = PDG;
                                _shr_bkt_purity = purity;
                                _shr_bkt_completeness = completeness;
                                _shr_bkt_E = E;
                                bool already_matched = (std::find(matched_energies.begin(), matched_energies.end(), E) != matched_energies.end());
                                if (!already_matched)
                                {
                                    TParticlePDG *particle_pdg = TDatabasePDG::Instance()->GetParticle(PDG);
                                    if (particle_pdg != NULL)
                                    { // PDG codes corresponding to ions e.g. 2000000101 are not in the database
                                        float ke = E - particle_pdg->Mass();
                                        _matched_E += ke;
                                        matched_energies.push_back(E);
                                    }
                                }
                            }
                        }
                    }

                    _shr_distance = (shr->ShowerStart() - nuvtx).Mag();
                    shr_vtx = shr->ShowerStart();
                    shr_p.SetXYZ(shr->Direction().X(), shr->Direction().Y(), shr->Direction().Z());
                    shr_p.SetMag(sqrt(pow(shr->Energy()[2] / 1000. + electron->Mass(), 2) - pow(electron->Mass(), 2)));
                    total_p += shr_p;
                    total_p_mu += shr_p;
                    std::vector<float> dqdx_cali(3);
                    searchingfornues::getDQdxCali(shr, dqdx_cali);
                    _shr_dedx_Y = shr->dEdx()[2];
                    _shr_dedx_V = shr->dEdx()[1];
                    _shr_dedx_U = shr->dEdx()[0];
                    _shr_dedx_Y_cali = shr->dEdx()[2] * dqdx_cali[2];
                    _shr_dedx_V_cali = shr->dEdx()[1] * dqdx_cali[1];
                    _shr_dedx_U_cali = shr->dEdx()[0] * dqdx_cali[0];
                    _shr_start_x = shr->ShowerStart().X();
                    _shr_start_y = shr->ShowerStart().Y();
                    _shr_start_z = shr->ShowerStart().Z();

                    _shr_energy = shr->Energy()[2] / 1000; // GeV
                    _shr_energy_cali = _shr_energy * cali_corr[2];

                    _shr_pfp_id = i_pfp;
                    _shr_hits_max = shr_hits;
                    _shr_score = trkshrscore;
                    _shr_theta = shr->Direction().Theta();
                    _shr_phi = shr->Direction().Phi();
                    _shr_px = shr_p.X();
                    _shr_py = shr_p.Y();
                    _shr_pz = shr_p.Z();
                    _shr_openangle = shr->OpenAngle();

                    if (tkcalo_proxy == NULL)
                    {
                        _shr_tkfit_start_x = std::numeric_limits<float>::lowest();
                        _shr_tkfit_start_y = std::numeric_limits<float>::lowest();
                        _shr_tkfit_start_z = std::numeric_limits<float>::lowest();

                        _shr_tkfit_phi = std::numeric_limits<float>::lowest();
                        _shr_tkfit_theta = std::numeric_limits<float>::lowest();
                        _shr_tkfit_dedx_Y = std::numeric_limits<float>::lowest();
                        _shr_tkfit_dedx_V = std::numeric_limits<float>::lowest();
                        _shr_tkfit_dedx_U = std::numeric_limits<float>::lowest();

                        _shr_tkfit_2cm_dedx_Y = std::numeric_limits<float>::lowest();
                        _shr_tkfit_2cm_dedx_V = std::numeric_limits<float>::lowest();
                        _shr_tkfit_2cm_dedx_U = std::numeric_limits<float>::lowest();

                        _shr_tkfit_gap05_dedx_Y = std::numeric_limits<float>::lowest();
                        _shr_tkfit_gap05_dedx_V = std::numeric_limits<float>::lowest();
                        _shr_tkfit_gap05_dedx_U = std::numeric_limits<float>::lowest();
                        _shr_tkfit_gap10_dedx_Y = std::numeric_limits<float>::lowest();
                        _shr_tkfit_gap10_dedx_V = std::numeric_limits<float>::lowest();
                        _shr_tkfit_gap10_dedx_U = std::numeric_limits<float>::lowest();

                        continue;
                    }

                    for (const searchingfornues::ProxyCaloElem_t &tk : *tkcalo_proxy)
                    {

                        // find track with ID matching the pfp index (this convention apparently works only for shower fits...)
                        if (tk->ID() != int(pfp_pxy_v[i_pfp].index()))
                            continue;

                        _shr_tkfit_npoints = tk->NumberTrajectoryPoints();
                        _shr_tkfit_npointsvalid = tk->CountValidPoints();

                        _shr_trkfitmedangle = searchingfornues::GetTrackRMSDeflection(tk, 10.);

			searchingfornues::GetMoliereRadius(pfp_pxy,_shrmoliereavg,_shrmoliererms);

                        _shr_tkfit_start_x = tk->Start().X();
                        _shr_tkfit_start_y = tk->Start().Y();
                        _shr_tkfit_start_z = tk->Start().Z();

			// SCE corrected shower start point
			float shr_tkfit_start_sce[3];
			searchingfornues::ApplySCECorrectionXYZ(_shr_tkfit_start_x,
								_shr_tkfit_start_y,
								_shr_tkfit_start_z,
								shr_tkfit_start_sce);

                        _shr_tkfit_phi = tk->StartDirection().Phi();
                        _shr_tkfit_theta = tk->StartDirection().Theta();

                        auto const tkcalos = tk.get<anab::Calorimetry>();

                        float calodEdx; // dEdx computed for track-fitter
                        int caloNpts;   // number of track-fitter dE/dx hits
                        float calodEdxXYZ; // dEdx computed for track-fitter [with XYZ and not RR]
                        int caloNptsXYZ;   // number of track-fitter dE/dx hits [with XYZ and not RR]

                        for (const auto &tkcalo : tkcalos)
                        {

                            if (tkcalo->ResidualRange().size() == 0)
                                continue;

                            // using function from CommonDefs/TrackFitterFunctions.h
                            searchingfornues::GetTrackFitdEdx(tkcalo, fdEdxcmSkip, fdEdxcmLen, fLocaldEdx, calodEdx, caloNpts);
			    // use xyz coordinates instead of RR to calculate distance
                            searchingfornues::GetTrackFitdEdx(tkcalo, fdEdxcmSkip, fdEdxcmLen, fLocaldEdx,
							      shr_tkfit_start_sce[0],shr_tkfit_start_sce[1],shr_tkfit_start_sce[2],
							      calodEdxXYZ, caloNptsXYZ);
			    
                            if (tkcalo->PlaneID().Plane == 2)
                            {
                                _shr_tkfit_dedx_Y = calodEdx;
                                _shr_tkfit_nhits_Y = caloNpts;
                                _shr_tkfit_dedx_Y_alt = calodEdxXYZ;
                                _shr_tkfit_nhits_Y_alt = caloNptsXYZ;
                            }
                            else if (tkcalo->PlaneID().Plane == 1)
                            {
                                _shr_tkfit_dedx_V = calodEdx;
                                _shr_tkfit_nhits_V = caloNpts;
                                _shr_tkfit_dedx_V_alt = calodEdxXYZ;
                                _shr_tkfit_nhits_V_alt = caloNptsXYZ;
                            }
                            else if (tkcalo->PlaneID().Plane == 0)
                            {
                                _shr_tkfit_dedx_U = calodEdx;
                                _shr_tkfit_nhits_U = caloNpts;
                                _shr_tkfit_dedx_U_alt = calodEdxXYZ;
                                _shr_tkfit_nhits_U_alt = caloNptsXYZ;
                            }

                            if (fSaveMoreDedx == false)
                                continue;

                            // use only first 2 cm
                            searchingfornues::GetTrackFitdEdx(tkcalo, fdEdxcmSkip, 2.0, fLocaldEdx, calodEdx, caloNpts);
                            if (tkcalo->PlaneID().Plane == 2)
                            {
                                _shr_tkfit_2cm_dedx_Y = calodEdx;
                                _shr_tkfit_2cm_nhits_Y = caloNpts;
                            }
                            else if (tkcalo->PlaneID().Plane == 1)
                            {
                                _shr_tkfit_2cm_dedx_V = calodEdx;
                                _shr_tkfit_2cm_nhits_V = caloNpts;
                            }
                            else if (tkcalo->PlaneID().Plane == 0)
                            {
                                _shr_tkfit_2cm_dedx_U = calodEdx;
                                _shr_tkfit_2cm_nhits_U = caloNpts;
                            }

                            // Gap 0.5 cm
                            searchingfornues::GetTrackFitdEdx(tkcalo, 0.5, fdEdxcmLen, fLocaldEdx, calodEdx, caloNpts);
                            if (tkcalo->PlaneID().Plane == 2)
                            {
                                _shr_tkfit_gap05_dedx_Y = calodEdx;
                                _shr_tkfit_gap05_nhits_Y = caloNpts;
                            }
                            else if (tkcalo->PlaneID().Plane == 1)
                            {
                                _shr_tkfit_gap05_dedx_V = calodEdx;
                                _shr_tkfit_gap05_nhits_V = caloNpts;
                            }
                            else if (tkcalo->PlaneID().Plane == 0)
                            {
                                _shr_tkfit_gap05_dedx_U = calodEdx;
                                _shr_tkfit_gap05_nhits_U = caloNpts;
                            }

                            // Gap 1.0 cm
                            searchingfornues::GetTrackFitdEdx(tkcalo, 1.0, fdEdxcmLen, fLocaldEdx, calodEdx, caloNpts);
                            if (tkcalo->PlaneID().Plane == 2)
                            {
                                _shr_tkfit_gap10_dedx_Y = calodEdx;
                                _shr_tkfit_gap10_nhits_Y = caloNpts;
                            }
                            else if (tkcalo->PlaneID().Plane == 1)
                            {
                                _shr_tkfit_gap10_dedx_V = calodEdx;
                                _shr_tkfit_gap10_nhits_V = caloNpts;
                            }
                            else if (tkcalo->PlaneID().Plane == 0)
                            {
                                _shr_tkfit_gap10_dedx_U = calodEdx;
                                _shr_tkfit_gap10_nhits_U = caloNpts;
                            }

                        } // for all calo objects
                    }
                }
            }
            for (const auto &trk : pfp_pxy.get<recob::Track>())
            {
                if (shr_hits == _shr_hits_max)
                {
                    auto trkpxy2 = pid_proxy[trk.key()];
                    auto pidpxy_v = trkpxy2.get<anab::ParticleID>();
                    if (pidpxy_v.size() > 0)
                    {
                        _shr_pidchipr = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, proton->PdgCode(), 2);
                        _shr_pidchimu = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, muon->PdgCode(), 2);
                        _shr_bragg_p = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, proton->PdgCode(), 2),
                                                searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, proton->PdgCode(), 2));
                        _shr_bragg_mu = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, muon->PdgCode(), 2),
                                                 searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, muon->PdgCode(), 2));
                        _shr_bragg_pion = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, pion->PdgCode(), 2),
                                                   searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, pion->PdgCode(), 2));
                        _shr_bragg_kaon = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, kaon->PdgCode(), 2),
                                                   searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, kaon->PdgCode(), 2));
                        _shr_bragg_mip = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 2);
                    }
                }
            }
        }// if shower-like PFParticle

	// if track-like PFParticle
        for (const auto &trk : pfp_pxy.get<recob::Track>())
        {
            if (trkshrscore > fTrkShrscore)
            {
                double trk_start[3] = {trk->Start().X(), trk->Start().Y(), trk->Start().Z()};
                double trk_end[3] = {trk->End().X(), trk->End().Y(), trk->End().Z()};
		
                auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();

                if (!searchingfornues::isFiducial(trk_start,fFidvolXstart,fFidvolYstart,fFidvolZstart,fFidvolXend,fFidvolYend,fFidvolZend) || 
		    !searchingfornues::isFiducial(trk_end,fFidvolXstart,fFidvolYstart,fFidvolZstart,fFidvolXend,fFidvolYend,fFidvolZend) )
                {
                    for (auto ass_clus : clus_pxy_v)
                    {
                        const auto &clus = clus_proxy[ass_clus.key()];
                        auto clus_hit_v = clus.get<recob::Hit>();
                        _hits_outfv += clus_hit_v.size();
                    }
                    continue;
                }
                _n_tracks_contained++;
                unsigned int trk_hits = 0;

                std::vector<art::Ptr<recob::Hit>> hit_v;

                for (auto ass_clus : clus_pxy_v)
                {
                    // get cluster proxy
                    const auto &clus = clus_proxy[ass_clus.key()];
                    auto clus_hit_v = clus.get<recob::Hit>();
                    trk_hits += clus_hit_v.size();
                    for (const auto &hit : clus_hit_v)
                    {
                        hit_v.push_back(hit);
                        if (clus->Plane().Plane == 2)
                            _trk_energy_hits_tot += 1.01 * hit->Integral() * 238 * 23.6e-6 / 0.6 / 1000;
                    }
                    if (clus->Plane().Plane == 0)
                    {
                        _trk_hits_u_tot += clus_hit_v.size();
                    }
                    else if (clus->Plane().Plane == 1)
                    {
                        _trk_hits_v_tot += clus_hit_v.size();
                    }
                    else if (clus->Plane().Plane == 2)
                    {
                        _trk_hits_y_tot += clus_hit_v.size();
                    }
                }

                // Kinetic energy from stopping power of proton in LAr
                float energy_proton = std::sqrt(std::pow(_trkmom.GetTrackMomentum(trk->Length(), proton->PdgCode()), 2) + std::pow(proton->Mass(), 2)) - proton->Mass();
                _trk_energy_tot += energy_proton;
                _trk_hits_tot += trk_hits;

                float energy_muon = std::sqrt(std::pow(_trkmom.GetTrackMomentum(trk->Length(), muon->PdgCode()), 2) + std::pow(muon->Mass(), 2)) - muon->Mass();
                _trk_energy_muon_tot += energy_muon;

                auto trkpxy2 = pid_proxy[trk.key()];
                auto pidpxy_v = trkpxy2.get<anab::ParticleID>();
                float chipr = std::numeric_limits<float>::lowest();
                float chimu = std::numeric_limits<float>::lowest();
                if (pidpxy_v.size() != 0)
                {
                    chipr = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, proton->PdgCode(), 2);
                    if (chipr > 0 && chipr < _trk_pidchipr_best)
                        _trk_pidchipr_best = chipr;
                    if (chipr > 0 && chipr > _trk_pidchipr_worst)
                        _trk_pidchipr_worst = chipr;

                    chimu = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, muon->PdgCode(), 2);
                    if (chimu > 0 && chimu < _trk_pidchimu_best)
                        _trk_pidchimu_best = chimu;
                    if (chimu > 0 && chimu > _trk_pidchimu_worst)
                        _trk_pidchimu_worst = chimu;
                }

                if (trk_hits > _trk_hits_max)
                {
                    if (!fData)
                    {
                        if (clus_pxy_v.size() != 0)
                        {
                            float purity = 0., completeness = 0.;
                            int ibt = searchingfornues::getAssocBtPart(hit_v, assocMCPart, btparts_v, purity, completeness);
                            if (ibt >= 0)
                            {
                                auto &mcp = btparts_v[ibt];
                                auto PDG = mcp.pdg;
                                auto E = mcp.e;
                                _trk_bkt_pdg = PDG;
                                _trk_bkt_purity = purity;
                                _trk_bkt_completeness = completeness;
                                _trk_bkt_E = E;
                                bool already_matched = (std::find(matched_energies.begin(), matched_energies.end(), E) != matched_energies.end());
                                if (!already_matched)
                                {
                                    TParticlePDG *particle_pdg = TDatabasePDG::Instance()->GetParticle(PDG);
                                    if (particle_pdg != NULL)
                                    { // PDG codes corresponding to ions e.g. 2000000101 are not in the database
                                        float ke = E - particle_pdg->Mass();
                                        _matched_E += ke;
                                        matched_energies.push_back(E);
                                    }
                                }
                            }
                        }
                    }

                    trk_vtx.SetXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
                    _trk_distance = (trk_vtx - nu_vtx).Mag();
                    _trk_len = trk->Length();
                    _trk_energy = energy_proton;
                    _trk_energy_muon = energy_muon;
                    _trk_energy_muon_mcs = std::sqrt(std::pow(mcsfitter.fitMcs(trk->Trajectory(), muon->PdgCode()).bestMomentum(), 2) + std::pow(muon->Mass(), 2)) - muon->Mass();
                    _trk_pfp_id = i_pfp;
                    _trk_hits_max = trk_hits;
                    _trk_theta = trk->Theta();
                    _trk_phi = trk->Phi();
                    trk_p.SetXYZ(trk->StartDirection().X(), trk->StartDirection().Y(), trk->StartDirection().Z());
                    trk_p.SetMag(sqrt(pow(energy_proton + proton->Mass(), 2) - pow(proton->Mass(), 2)));
                    trk_p_mu.SetXYZ(trk->StartDirection().X(), trk->StartDirection().Y(), trk->StartDirection().Z());
                    trk_p_mu.SetMag(sqrt(pow(energy_muon + muon->Mass(), 2) - pow(muon->Mass(), 2)));
                    total_p += trk_p;
                    total_p_mu += trk_p_mu;
                    if (pidpxy_v.size() != 0)
                    {
                        _trk_bragg_p = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, proton->PdgCode(), 2),
                                                searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, proton->PdgCode(), 2));
                        _trk_bragg_mu = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, muon->PdgCode(), 2),
                                                 searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, muon->PdgCode(), 2));
                        _trk_bragg_pion = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, pion->PdgCode(), 2),
                                                   searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, pion->PdgCode(), 2));
                        _trk_bragg_kaon = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, kaon->PdgCode(), 2),
                                                   searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, kaon->PdgCode(), 2));
                        _trk_bragg_mip = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 2);
                        _trk_pidchipr = chipr;
                        _trk_pidchimu = chimu;
                        _trk_pida = searchingfornues::PID(pidpxy_v[0], "PIDA_median", anab::kPIDA, anab::kForward, 0, 2);
                    }
                    _trk_score = trkshrscore;
                }
            }
        }
    }
    
    auto trkshrhitdist_v = searchingfornues::GetPFPHitDistance(pfp_pxy_v[_trk_pfp_id], pfp_pxy_v[_shr_pfp_id], clus_proxy);
    _trkshrhitdist0 = trkshrhitdist_v[0];
    _trkshrhitdist1 = trkshrhitdist_v[1];
    _trkshrhitdist2 = trkshrhitdist_v[2];

    _extra_energy_y -= (_trk_energy_hits_tot + _shr_energy_tot);
    _pt = total_p.Perp();
    _p = total_p.Mag();
    _pt_assume_muon = total_p_mu.Perp();
    _p_assume_muon = total_p_mu.Mag();
    if (_trk_hits_tot + _shr_hits_tot > 0)
        _hits_ratio = (float)_shr_hits_tot / (_trk_hits_tot + _shr_hits_tot);
    _tksh_distance = (trk_vtx - shr_vtx).Mag();
    if (trk_p.Mag() * shr_p.Mag() > 0)
        _tksh_angle = trk_p.Dot(shr_p) / (trk_p.Mag() * shr_p.Mag());
    if (trk_p_mu.Mag() * shr_p.Mag() > 0)
        _tksh_angle_muon = trk_p_mu.Dot(shr_p) / (trk_p_mu.Mag() * shr_p.Mag());

    _contained_fraction = ((float)(_trk_hits_tot + _shr_hits_tot)) / (_trk_hits_tot + _shr_hits_tot + _hits_outfv);
    _sps_contained_fraction = float(sps_fv) / float(sps_all);

    if (_contained_fraction < 0.9)
        return false;

    if (!(_n_showers_contained > 0))
        return false;
    return true;
}



void CC0piNpSelection::resetTTree(TTree *_tree)
{

    _shr_pfp_id = 0;
    _trk_pfp_id = 0;

    _trk_hits_max = 0;
    _shr_hits_max = 0;

    _shr_dedx_Y = std::numeric_limits<float>::lowest();
    _shr_dedx_V = std::numeric_limits<float>::lowest();
    _shr_dedx_U = std::numeric_limits<float>::lowest();
    _shr_dedx_Y_cali = std::numeric_limits<float>::lowest();
    _shr_dedx_V_cali = std::numeric_limits<float>::lowest();
    _shr_dedx_U_cali = std::numeric_limits<float>::lowest();
    _shr_score = std::numeric_limits<float>::lowest();
    _shr_energy = 0;
    _shr_energy_cali = 0;
    _shr_energy_tot_cali = 0;
    _shr_energy_tot = 0;
    _shr_distance = std::numeric_limits<float>::lowest();
    _tksh_distance = std::numeric_limits<float>::lowest();
    _tksh_angle = std::numeric_limits<float>::lowest();

    _n_showers_contained = 0;
    _n_tracks_contained = 0;

    _trk_distance = std::numeric_limits<float>::lowest();
    _trk_len = std::numeric_limits<float>::lowest();
    _trk_energy = 0;
    _trk_energy_muon = 0;
    _trk_energy_muon_mcs = 0;
    _trk_energy_tot = 0;
    _trk_energy_muon_tot = 0;
    _hits_ratio = 0;
    _trk_hits_tot = 0;
    _trk_hits_y_tot = 0;
    _trk_hits_u_tot = 0;
    _trk_hits_v_tot = 0;

    _shr_hits_tot = 0;
    _shr_hits_y_tot = 0;
    _shr_hits_v_tot = 0;
    _shr_hits_u_tot = 0;

    _trk_bragg_p = std::numeric_limits<float>::lowest();
    _trk_bragg_mu = std::numeric_limits<float>::lowest();
    _trk_bragg_mip = std::numeric_limits<float>::lowest();
    _trk_bragg_pion = std::numeric_limits<float>::lowest();
    _trk_bragg_kaon = std::numeric_limits<float>::lowest();
    _trk_pidchipr = std::numeric_limits<float>::lowest();
    _trk_pidchipr_best = std::numeric_limits<float>::max();
    _trk_pidchipr_worst = std::numeric_limits<float>::lowest();

    _trk_pidchimu = std::numeric_limits<float>::lowest();
    _trk_pidchimu_best = std::numeric_limits<float>::max();
    _trk_pidchimu_worst = std::numeric_limits<float>::lowest();

    _trk_pida = std::numeric_limits<float>::lowest();
    _trk_score = std::numeric_limits<float>::lowest();
    _trk_bkt_pdg = 0;
    _trk_bkt_purity = std::numeric_limits<float>::lowest();
    _trk_bkt_completeness = std::numeric_limits<float>::lowest();
    _trk_bkt_E = std::numeric_limits<float>::lowest();

    _pt = 0;
    _p = 0;
    _pt_assume_muon = 0;
    _p_assume_muon = 0;

    _trkshrhitdist0 = std::numeric_limits<float>::lowest();
    _trkshrhitdist1 = std::numeric_limits<float>::lowest();
    _trkshrhitdist2 = std::numeric_limits<float>::lowest();

    _trk_theta = std::numeric_limits<float>::lowest();
    _trk_phi = std::numeric_limits<float>::lowest();
    _shr_theta = std::numeric_limits<float>::lowest();
    _shr_phi = std::numeric_limits<float>::lowest();
    _shr_px = 0;
    _shr_py = 0;
    _shr_pz = 0;
    _shr_pca_0 = std::numeric_limits<float>::lowest();
    _shr_pca_1 = std::numeric_limits<float>::lowest();
    _shr_pca_2 = std::numeric_limits<float>::lowest();

    _shr_openangle = std::numeric_limits<float>::lowest();
    _shr_bkt_pdg = 0;
    _shr_bkt_purity = std::numeric_limits<float>::lowest();
    _shr_bkt_completeness = std::numeric_limits<float>::lowest();
    _shr_bkt_E = std::numeric_limits<float>::lowest();
    _shr_tkfit_start_x = std::numeric_limits<float>::lowest();
    _shr_tkfit_start_y = std::numeric_limits<float>::lowest();
    _shr_tkfit_start_z = std::numeric_limits<float>::lowest();
    _shr_start_x = std::numeric_limits<float>::lowest();
    _shr_start_y = std::numeric_limits<float>::lowest();
    _shr_start_z = std::numeric_limits<float>::lowest();

    _shrsubclusters0 = std::numeric_limits<int>::lowest();
    _shrsubclusters1 = std::numeric_limits<int>::lowest();
    _shrsubclusters2 = std::numeric_limits<int>::lowest();

    _shrclusfrac0 = std::numeric_limits<float>::lowest();
    _shrclusfrac1 = std::numeric_limits<float>::lowest();
    _shrclusfrac2 = std::numeric_limits<float>::lowest();

    _shrclusdir0 = std::numeric_limits<float>::lowest();
    _shrclusdir1 = std::numeric_limits<float>::lowest();
    _shrclusdir2 = std::numeric_limits<float>::lowest();

    _shr_tkfit_phi = std::numeric_limits<float>::lowest();
    _shr_tkfit_theta = std::numeric_limits<float>::lowest();
    _shr_tkfit_dedx_Y = std::numeric_limits<float>::lowest();
    _shr_tkfit_dedx_U = std::numeric_limits<float>::lowest();
    _shr_tkfit_dedx_V = std::numeric_limits<float>::lowest();
    _shr_tkfit_nhits_Y = 0;
    _shr_tkfit_nhits_U = 0;
    _shr_tkfit_nhits_V = 0;
    _shr_tkfit_dedx_Y_alt = std::numeric_limits<float>::lowest();
    _shr_tkfit_dedx_U_alt = std::numeric_limits<float>::lowest();
    _shr_tkfit_dedx_V_alt = std::numeric_limits<float>::lowest();
    _shr_tkfit_nhits_Y_alt = 0;
    _shr_tkfit_nhits_U_alt = 0;
    _shr_tkfit_nhits_V_alt = 0;


    _shr_tkfit_npoints = std::numeric_limits<int>::lowest();
    _shr_tkfit_npointsvalid = std::numeric_limits<int>::lowest();

    _shr_trkfitmedangle = std::numeric_limits<float>::lowest();

    _shrmoliereavg = std::numeric_limits<float>::lowest();
    _shrmoliererms = std::numeric_limits<float>::lowest();

    _shr_tkfit_2cm_dedx_Y = std::numeric_limits<float>::lowest();
    _shr_tkfit_2cm_dedx_U = std::numeric_limits<float>::lowest();
    _shr_tkfit_2cm_dedx_V = std::numeric_limits<float>::lowest();
    _shr_tkfit_2cm_nhits_Y = 0;
    _shr_tkfit_2cm_nhits_U = 0;
    _shr_tkfit_2cm_nhits_V = 0;

    _shr_tkfit_gap05_dedx_Y = std::numeric_limits<float>::lowest();
    _shr_tkfit_gap05_dedx_U = std::numeric_limits<float>::lowest();
    _shr_tkfit_gap05_dedx_V = std::numeric_limits<float>::lowest();
    _shr_tkfit_gap05_nhits_Y = 0;
    _shr_tkfit_gap05_nhits_U = 0;
    _shr_tkfit_gap05_nhits_V = 0;

    _shr_tkfit_gap10_dedx_Y = std::numeric_limits<float>::lowest();
    _shr_tkfit_gap10_dedx_U = std::numeric_limits<float>::lowest();
    _shr_tkfit_gap10_dedx_V = std::numeric_limits<float>::lowest();
    _shr_tkfit_gap10_nhits_Y = 0;
    _shr_tkfit_gap10_nhits_U = 0;
    _shr_tkfit_gap10_nhits_V = 0;

    _shr_pidchipr = std::numeric_limits<float>::lowest();
    _shr_pidchimu = std::numeric_limits<float>::lowest();
    _shr_bragg_p = std::numeric_limits<float>::lowest();
    _shr_bragg_mu = std::numeric_limits<float>::lowest();
    _shr_bragg_mip = std::numeric_limits<float>::lowest();
    _shr_bragg_pion = std::numeric_limits<float>::lowest();
    _shr_bragg_kaon = std::numeric_limits<float>::lowest();

    _matched_E = 0;
    _total_hits_y = 0;
    _extra_energy_y = 0;
    _trk_energy_hits_tot = 0;
    _hits_outfv = 0;
    _contained_fraction = 0;
    _sps_contained_fraction = 0;
}

void CC0piNpSelection::setBranches(TTree *_tree)
{

    _tree->Branch("trk_id", &_trk_pfp_id, "trk_pfp_id/i");
    _tree->Branch("shr_id", &_shr_pfp_id, "shr_pfp_id/i");

    _tree->Branch("shr_energy_tot", &_shr_energy_tot, "shr_energy_tot/F");
    _tree->Branch("shr_energy", &_shr_energy, "shr_energy/F");
    _tree->Branch("shr_energy_tot_cali", &_shr_energy_tot_cali, "shr_energy_tot_cali/F");
    _tree->Branch("shr_energy_cali", &_shr_energy_cali, "shr_energy_cali/F");
    _tree->Branch("shr_theta", &_shr_theta, "shr_theta/F");
    _tree->Branch("shr_phi", &_shr_phi, "shr_phi/F");
    _tree->Branch("shr_pca_0", &_shr_pca_0, "shr_pca_0/F");
    _tree->Branch("shr_pca_1", &_shr_pca_1, "shr_pca_1/F");
    _tree->Branch("shr_pca_2", &_shr_pca_2, "shr_pca_2/F");

    _tree->Branch("shr_px", &_shr_px, "shr_px/F");
    _tree->Branch("shr_py", &_shr_py, "shr_py/F");
    _tree->Branch("shr_pz", &_shr_pz, "shr_pz/F");
    _tree->Branch("shr_openangle", &_shr_openangle, "shr_openangle/F");
    _tree->Branch("shr_tkfit_start_x", &_shr_tkfit_start_x, "shr_tkfit_start_x/F");
    _tree->Branch("shr_tkfit_start_y", &_shr_tkfit_start_y, "shr_tkfit_start_y/F");
    _tree->Branch("shr_tkfit_start_z", &_shr_tkfit_start_z, "shr_tkfit_start_z/F");
    _tree->Branch("shr_tkfit_theta", &_shr_tkfit_theta, "shr_tkfit_theta/F");
    _tree->Branch("shr_tkfit_phi", &_shr_tkfit_phi, "shr_tkfit_phi/F");
    _tree->Branch("shr_start_x", &_shr_start_x, "shr_start_x/F");
    _tree->Branch("shr_start_y", &_shr_start_y, "shr_start_y/F");
    _tree->Branch("shr_start_z", &_shr_start_z, "shr_start_z/F");
    _tree->Branch("shr_dedx_Y", &_shr_dedx_Y, "shr_dedx_Y/F");
    _tree->Branch("shr_dedx_V", &_shr_dedx_V, "shr_dedx_V/F");
    _tree->Branch("shr_dedx_U", &_shr_dedx_U, "shr_dedx_U/F");
    _tree->Branch("shr_dedx_Y_cali", &_shr_dedx_Y_cali, "shr_dedx_Y_cali/F");
    _tree->Branch("shr_dedx_V_cali", &_shr_dedx_V_cali, "shr_dedx_V_cali/F");
    _tree->Branch("shr_dedx_U_cali", &_shr_dedx_U_cali, "shr_dedx_U_cali/F");
    _tree->Branch("shr_tkfit_dedx_Y", &_shr_tkfit_dedx_Y, "shr_tkfit_dedx_Y/F");
    _tree->Branch("shr_tkfit_dedx_V", &_shr_tkfit_dedx_V, "shr_tkfit_dedx_V/F");
    _tree->Branch("shr_tkfit_dedx_U", &_shr_tkfit_dedx_U, "shr_tkfit_dedx_U/F");
    _tree->Branch("shr_tkfit_nhits_Y", &_shr_tkfit_nhits_Y, "shr_tkfit_nhits_Y/i");
    _tree->Branch("shr_tkfit_nhits_V", &_shr_tkfit_nhits_V, "shr_tkfit_nhits_V/i");
    _tree->Branch("shr_tkfit_nhits_U", &_shr_tkfit_nhits_U, "shr_tkfit_nhits_U/i");
    _tree->Branch("shr_tkfit_dedx_Y_alt", &_shr_tkfit_dedx_Y_alt, "shr_tkfit_dedx_Y_alt/F");
    _tree->Branch("shr_tkfit_dedx_V_alt", &_shr_tkfit_dedx_V_alt, "shr_tkfit_dedx_V_alt/F");
    _tree->Branch("shr_tkfit_dedx_U_alt", &_shr_tkfit_dedx_U_alt, "shr_tkfit_dedx_U_alt/F");
    _tree->Branch("shr_tkfit_nhits_Y_alt", &_shr_tkfit_nhits_Y_alt, "shr_tkfit_nhits_Y_alt/i");
    _tree->Branch("shr_tkfit_nhits_V_alt", &_shr_tkfit_nhits_V_alt, "shr_tkfit_nhits_V_alt/i");
    _tree->Branch("shr_tkfit_nhits_U_alt", &_shr_tkfit_nhits_U_alt, "shr_tkfit_nhits_U_alt/i");
    _tree->Branch("shr_tkfit_npoints", &_shr_tkfit_npoints, "shr_tkfit_npoints/i");
    _tree->Branch("shr_tkfit_npointsvalid", &_shr_tkfit_npointsvalid, "shr_tkfit_npointsvalid/i");
    _tree->Branch("shr_trkfitmedangle", &_shr_trkfitmedangle, "shr_trkfitmedangle/f");
    _tree->Branch("shrmoliereavg", &_shrmoliereavg, "shrmoliereavg/f");
    _tree->Branch("shrmoliererms", &_shrmoliererms, "shrmoliererms/f");
    if (fSaveMoreDedx)
    {
        _tree->Branch("shr_tkfit_2cm_dedx_Y", &_shr_tkfit_2cm_dedx_Y, "shr_tkfit_2cm_dedx_Y/F");
        _tree->Branch("shr_tkfit_2cm_dedx_V", &_shr_tkfit_2cm_dedx_V, "shr_tkfit_2cm_dedx_V/F");
        _tree->Branch("shr_tkfit_2cm_dedx_U", &_shr_tkfit_2cm_dedx_U, "shr_tkfit_2cm_dedx_U/F");
        _tree->Branch("shr_tkfit_2cm_nhits_Y", &_shr_tkfit_2cm_nhits_Y, "shr_tkfit_2cm_nhits_Y/i");
        _tree->Branch("shr_tkfit_2cm_nhits_V", &_shr_tkfit_2cm_nhits_V, "shr_tkfit_2cm_nhits_V/i");
        _tree->Branch("shr_tkfit_2cm_nhits_U", &_shr_tkfit_2cm_nhits_U, "shr_tkfit_2cm_nhits_U/i");
        _tree->Branch("shr_tkfit_gap05_dedx_Y", &_shr_tkfit_gap05_dedx_Y, "shr_tkfit_gap05_dedx_Y/F");
        _tree->Branch("shr_tkfit_gap05_dedx_V", &_shr_tkfit_gap05_dedx_V, "shr_tkfit_gap05_dedx_V/F");
        _tree->Branch("shr_tkfit_gap05_dedx_U", &_shr_tkfit_gap05_dedx_U, "shr_tkfit_gap05_dedx_U/F");
        _tree->Branch("shr_tkfit_gap05_nhits_Y", &_shr_tkfit_gap05_nhits_Y, "shr_tkfit_gap05_nhits_Y/i");
        _tree->Branch("shr_tkfit_gap05_nhits_V", &_shr_tkfit_gap05_nhits_V, "shr_tkfit_gap05_nhits_V/i");
        _tree->Branch("shr_tkfit_gap05_nhits_U", &_shr_tkfit_gap05_nhits_U, "shr_tkfit_gap05_nhits_U/i");
        _tree->Branch("shr_tkfit_gap10_dedx_Y", &_shr_tkfit_gap10_dedx_Y, "shr_tkfit_gap10_dedx_Y/F");
        _tree->Branch("shr_tkfit_gap10_dedx_V", &_shr_tkfit_gap10_dedx_V, "shr_tkfit_gap10_dedx_V/F");
        _tree->Branch("shr_tkfit_gap10_dedx_U", &_shr_tkfit_gap10_dedx_U, "shr_tkfit_gap10_dedx_U/F");
        _tree->Branch("shr_tkfit_gap10_nhits_Y", &_shr_tkfit_gap10_nhits_Y, "shr_tkfit_gap10_nhits_Y/i");
        _tree->Branch("shr_tkfit_gap10_nhits_V", &_shr_tkfit_gap10_nhits_V, "shr_tkfit_gap10_nhits_V/i");
        _tree->Branch("shr_tkfit_gap10_nhits_U", &_shr_tkfit_gap10_nhits_U, "shr_tkfit_gap10_nhits_U/i");
    }
    _tree->Branch("shr_chipr", &_shr_pidchipr, "shr_chipr/F");
    _tree->Branch("shr_chimu", &_shr_pidchimu, "shr_chimu/F");
    _tree->Branch("shr_bragg_p", &_shr_bragg_p, "shr_bragg_p/F");
    _tree->Branch("shr_bragg_mu", &_shr_bragg_mu, "shr_bragg_mu/F");
    _tree->Branch("shr_bragg_mip", &_shr_bragg_mip, "shr_bragg_mip/F");
    _tree->Branch("shr_bragg_kaon", &_shr_bragg_kaon, "shr_bragg_kaon/F");
    _tree->Branch("shr_bragg_pion", &_shr_bragg_pion, "shr_bragg_pion/F");

    _tree->Branch("tksh_distance", &_tksh_distance, "tksh_distance/F");
    _tree->Branch("tksh_angle", &_tksh_angle, "tksh_angle/F");
    _tree->Branch("shr_distance", &_shr_distance, "shr_distance/F");
    _tree->Branch("shr_score", &_shr_score, "shr_score/F");
    _tree->Branch("shr_bkt_pdg", &_shr_bkt_pdg, "shr_bkt_pdg/I");
    _tree->Branch("shr_bkt_purity", &_shr_bkt_purity, "shr_bkt_purity/F");
    _tree->Branch("shr_bkt_completeness", &_shr_bkt_completeness, "shr_bkt_completeness/F");
    _tree->Branch("shr_bkt_E", &_shr_bkt_E, "shr_bkt_E/F");

    _tree->Branch("trk_len", &_trk_len, "trk_len/F");
    _tree->Branch("trk_theta", &_trk_theta, "trk_theta/F");
    _tree->Branch("trk_phi", &_trk_phi, "trk_phi/F");
    _tree->Branch("trk_energy", &_trk_energy, "trk_energy/F");
    _tree->Branch("trk_energy_muon", &_trk_energy_muon, "trk_energy_muon/F");
    _tree->Branch("trk_energy_muon_mcs", &_trk_energy_muon_mcs, "trk_energy_muon_mcs/F");
    _tree->Branch("trk_energy_tot", &_trk_energy_tot, "trk_energy_tot/F");
    _tree->Branch("trk_energy_muon_tot", &_trk_energy_muon_tot, "trk_energy_muon_tot/F");
    _tree->Branch("trk_distance", &_trk_distance, "trk_distance/F");
    _tree->Branch("trk_score", &_trk_score, "trk_score/F");
    _tree->Branch("trk_bkt_pdg", &_trk_bkt_pdg, "trk_bkt_pdg/I");
    _tree->Branch("trk_bkt_purity", &_trk_bkt_purity, "trk_bkt_purity/F");
    _tree->Branch("trk_bkt_completeness", &_trk_bkt_completeness, "trk_bkt_completeness/F");
    _tree->Branch("trk_bkt_E", &_trk_bkt_E, "trk_bkt_E/F");
    _tree->Branch("trk_chipr_best", &_trk_pidchipr_best, "trk_chipr_best/F");
    _tree->Branch("trk_chipr_worst", &_trk_pidchipr_worst, "trk_chipr_worst/F");
    _tree->Branch("trk_chimu_best", &_trk_pidchimu_best, "trk_chimu_best/F");
    _tree->Branch("trk_chimu_worst", &_trk_pidchimu_worst, "trk_chimu_worst/F");

    _tree->Branch("trk_chipr", &_trk_pidchipr, "trk_chipr/F");
    _tree->Branch("trk_chimu", &_trk_pidchimu, "trk_chimu/F");
    _tree->Branch("trk_pida", &_trk_pida, "trk_pida/F");
    _tree->Branch("trk_bragg_p", &_trk_bragg_p, "trk_bragg_p/F");
    _tree->Branch("trk_bragg_mu", &_trk_bragg_mu, "trk_bragg_mu/F");
    _tree->Branch("trk_bragg_mip", &_trk_bragg_mip, "trk_bragg_mip/F");
    _tree->Branch("trk_bragg_kaon", &_trk_bragg_kaon, "trk_bragg_kaon/F");
    _tree->Branch("trk_bragg_pion", &_trk_bragg_pion, "trk_bragg_pion/F");

    _tree->Branch("trk_hits_max", &_trk_hits_max, "trk_hits_max/i");
    _tree->Branch("shr_hits_max", &_shr_hits_max, "shr_hits_max/i");

    _tree->Branch("trkshrhitdist0", &_trkshrhitdist0, "trkshrhitdist0/F");
    _tree->Branch("trkshrhitdist1", &_trkshrhitdist1, "trkshrhitdist1/F");
    _tree->Branch("trkshrhitdist2", &_trkshrhitdist2, "trkshrhitdist2/F");

    _tree->Branch("total_hits_y", &_total_hits_y, "total_hits_y/i");
    _tree->Branch("extra_energy_y", &_extra_energy_y, "extra_energy_y/F");
    _tree->Branch("trk_energy_hits_tot", &_trk_energy_hits_tot, "trk_energy_hits_tot/F");

    _tree->Branch("shrsubclusters0", &_shrsubclusters0, "shrsubclusters0/i");
    _tree->Branch("shrsubclusters1", &_shrsubclusters1, "shrsubclusters1/i");
    _tree->Branch("shrsubclusters2", &_shrsubclusters2, "shrsubclusters2/i");

    _tree->Branch("shrclusfrac0", &_shrclusfrac0, "shrclusfrac0/f");
    _tree->Branch("shrclusfrac1", &_shrclusfrac1, "shrclusfrac1/f");
    _tree->Branch("shrclusfrac2", &_shrclusfrac2, "shrclusfrac2/f");

    _tree->Branch("shrclusdir0", &_shrclusdir0, "shrclusdir0/f");
    _tree->Branch("shrclusdir1", &_shrclusdir1, "shrclusdir1/f");
    _tree->Branch("shrclusdir2", &_shrclusdir2, "shrclusdir2/f");

    _tree->Branch("shr_hits_tot", &_shr_hits_tot, "shr_hits_tot/i");
    _tree->Branch("shr_hits_y_tot", &_shr_hits_y_tot, "shr_hits_y_tot/i");
    _tree->Branch("shr_hits_u_tot", &_shr_hits_u_tot, "shr_hits_u_tot/i");
    _tree->Branch("shr_hits_v_tot", &_shr_hits_v_tot, "shr_hits_v_tot/i");

    _tree->Branch("trk_hits_tot", &_trk_hits_tot, "trk_hits_tot/i");
    _tree->Branch("trk_hits_y_tot", &_trk_hits_y_tot, "trk_hits_y_tot/i");
    _tree->Branch("trk_hits_u_tot", &_trk_hits_u_tot, "trk_hits_u_tot/i");
    _tree->Branch("trk_hits_v_tot", &_trk_hits_v_tot, "trk_hits_v_tot/i");

    _tree->Branch("n_tracks_contained", &_n_tracks_contained, "n_tracks_contained/i");
    _tree->Branch("n_showers_contained", &_n_showers_contained, "n_showers_contained/i");

    _tree->Branch("matched_E", &_matched_E, "matched_E/F");

    _tree->Branch("hits_ratio", &_hits_ratio, "hits_ratio/F");
    _tree->Branch("contained_fraction", &_contained_fraction, "contained_fraction/F");
    _tree->Branch("sps_contained_fraction", &_sps_contained_fraction, "sps_contained_fraction/F");
    _tree->Branch("pt", &_pt, "pt/F");
    _tree->Branch("p", &_p, "p/F");
    _tree->Branch("pt_assume_muon", &_pt_assume_muon, "pt_assume_muon/F");
    _tree->Branch("p_assume_muon", &_p_assume_muon, "p_assume_muon/F");
}

DEFINE_ART_CLASS_TOOL(CC0piNpSelection)
} // namespace selection

#endif
