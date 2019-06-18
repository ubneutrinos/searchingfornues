#ifndef ANALYSIS_TRACKANALYSIS_CXX
#define ANALYSIS_TRACKANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "../CommonDefs/Typedefs.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/PIDFuncs.h"

#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       TrackAnalysis
// File:        TrackAnalysis.cc
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

class TrackAnalysis : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  TrackAnalysis(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~TrackAnalysis(){};

  // provide for initialization
  void configure(fhicl::ParameterSet const &pset);

  /**
     * @brief Analysis function
     */
  void analyzeEvent(art::Event const &e, bool fData) override;

  /**
     * @brief Analyze slice
     */
  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

  /**
     * @brief Save truth info for event associated to neutrino
     */
  void SaveTruth(art::Event const &e);

  /**
     * @brief Fill Default info for event associated to neutrino
     */
  void fillDefault();

  /**
     * @brief set branches for TTree
     */
  void setBranches(TTree *_tree) override;

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree) override;


private:
  trkf::TrackMomentumCalculator _trkmom;

  TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
  TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);

  art::InputTag fCALOproducer;
  art::InputTag fPIDproducer;
  art::InputTag fTRKproducer;
  float fTrkShrScore; /**< Threshold on the Pandora track score (default 0.5) */

  unsigned int _n_tracks;

  std::vector<size_t> _trk_pfp_id;

  std::vector<float> _trk_score_v;

  std::vector<float> _trk_start_x;
  std::vector<float> _trk_start_y;
  std::vector<float> _trk_start_z;

  std::vector<float> _trk_distance;

  std::vector<float> _trk_theta;
  std::vector<float> _trk_phi;

  std::vector<float> _trk_dir_x;
  std::vector<float> _trk_dir_y;
  std::vector<float> _trk_dir_z;

  std::vector<float> _trk_end_x;
  std::vector<float> _trk_end_y;
  std::vector<float> _trk_end_z;

  std::vector<float> _trk_length;

  std::vector<float> _trk_bragg_p;
  std::vector<float> _trk_bragg_mu;
  std::vector<float> _trk_bragg_mip;

  std::vector<float> _trk_pidchi;
  std::vector<float> _trk_pidchipr;
  std::vector<float> _trk_pidchika;
  std::vector<float> _trk_pidchipi;
  std::vector<float> _trk_pidchimu;
  std::vector<float> _trk_pida;
  std::vector<float> _trk_energy_proton;
  std::vector<float> _trk_energy_muon;


};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
TrackAnalysis::TrackAnalysis(const fhicl::ParameterSet &p)
{

  fCALOproducer = p.get< art::InputTag > ("CALOproducer");
  fPIDproducer  = p.get< art::InputTag > ("PIDproducer" );
  fTRKproducer  = p.get< art::InputTag > ("TRKproducer" );
  fTrkShrScore  = p.get<float>("TrkShrScore", 0.5);

}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void TrackAnalysis::configure(fhicl::ParameterSet const &p)
{
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void TrackAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
}

void TrackAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
  searchingfornues::ProxyCaloColl_t const& calo_proxy = proxy::getCollection<std::vector<recob::Track> >(e, fTRKproducer,
													 proxy::withAssociated<anab::Calorimetry>(fCALOproducer));

  searchingfornues::ProxyPIDColl_t const& pid_proxy = proxy::getCollection<std::vector<recob::Track> >(e, fTRKproducer,
													 proxy::withAssociated<anab::ParticleID>(fPIDproducer));

	TVector3 nuvtx;

  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
  {
    auto PDG = fabs(slice_pfp_v[i_pfp]->PdgCode());

    if (PDG == 12 || PDG == 14)
    {
      // grab vertex
      float_t xyz[3] = {};

      auto vtx = slice_pfp_v[i_pfp].get<recob::Vertex>();
      if (vtx.size() != 1)
      {
        std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
      }
      else
      {
        // save vertex to array
        vtx.at(0)->XYZ(xyz);
        nuvtx.SetXYZ(xyz[0], xyz[1], xyz[2]);
      }

      break;
    }
  }

  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
  {
    auto PDG = fabs(slice_pfp_v[i_pfp]->PdgCode());

    if (PDG == 12 || PDG == 14)
      continue;

    auto trk_v = slice_pfp_v[i_pfp].get<recob::Track>();

    auto ntrk = trk_v.size();
    if (ntrk == 0)
    {
      fillDefault();
    }
    else if (ntrk == 1)
    {
      auto trk = trk_v.at(0);
      float trkscore = searchingfornues::GetTrackShowerScore(slice_pfp_v[i_pfp]);
      if (trkscore > fTrkShrScore)
      {
        _n_tracks++;
      }
      _trk_score_v.push_back(trkscore);

      // if (trkscore > fTrkShrScore) {
      // get track proxy in order to fetch calorimtry
      // auto trkpxy1 = calo_proxy[trk.key()];
      // auto calopxy_v = trkpxy1.get<anab::Calorimetry>();

      // get trk proxy in order to fetch PID
      auto trkpxy2 = pid_proxy[trk.key()];
      auto pidpxy_v = trkpxy2.get<anab::ParticleID>();

      float bragg_p = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 2),
                                searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 2));

      float bragg_mu = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 2),
                                  searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 2));

      float bragg_mip = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 2);

      float pidchipr = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 0);
      float pidchimu = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 0);
      float pidchipi = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 0);
      float pidchika = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 0);

      float pida_mean = searchingfornues::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 2);



      _trk_bragg_p.push_back(bragg_p);
      _trk_bragg_mu.push_back(bragg_mu);
      _trk_bragg_mip.push_back(bragg_mip);
      _trk_pida.push_back(pida_mean);
      _trk_pidchipr.push_back(pidchipr);
      _trk_pidchimu.push_back(pidchimu);
      _trk_pidchipi.push_back(pidchipi);
      _trk_pidchika.push_back(pidchika);

      // Kinetic energy using tabulated stopping power (GeV)
      float energy_proton = std::sqrt(std::pow(_trkmom.GetTrackMomentum(trk->Length(), 2212), 2) + std::pow(proton->Mass(), 2)) - proton->Mass();
      float energy_muon = std::sqrt(std::pow(_trkmom.GetTrackMomentum(trk->Length(), 13), 2) + std::pow(muon->Mass(), 2)) - muon->Mass();

      _trk_energy_proton.push_back(energy_proton);
      _trk_energy_muon.push_back(energy_muon);

      _trk_dir_x.push_back(trk->StartDirection().X());
      _trk_dir_y.push_back(trk->StartDirection().Y());
      _trk_dir_z.push_back(trk->StartDirection().Z());

      _trk_start_x.push_back(trk->Start().X());
      _trk_start_y.push_back(trk->Start().Y());
      _trk_start_z.push_back(trk->Start().Z());

      _trk_end_x.push_back(trk->End().X());
      _trk_end_y.push_back(trk->End().Y());
      _trk_end_z.push_back(trk->End().Z());

      _trk_theta.push_back(trk->Theta());
      _trk_phi.push_back(trk->Phi());

      _trk_length.push_back(trk->Length());

      TVector3 trk_vtx;
      trk_vtx.SetXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
      trk_vtx -= nuvtx;
      _trk_distance.push_back(trk_vtx.Mag());

      _trk_pfp_id.push_back(i_pfp);
    }
  } // for all PFParticles
}

void ShowerAnalysis::fillDefault()
{
  _trk_pfp_id.push_back(std::numeric_limits<int>::lowest());

  _trk_score_v.push_back(std::numeric_limits<float>::lowest());

  _trk_start_x.push_back(std::numeric_limits<float>::lowest());
  _trk_start_y.push_back(std::numeric_limits<float>::lowest());
  _trk_start_z.push_back(std::numeric_limits<float>::lowest());

  _trk_distance.push_back(std::numeric_limits<float>::lowest());

  _trk_theta.push_back(std::numeric_limits<float>::lowest());
  _trk_phi.push_back(std::numeric_limits<float>::lowest());

  _trk_dir_x.push_back(std::numeric_limits<float>::lowest());
  _trk_dir_y.push_back(std::numeric_limits<float>::lowest());
  _trk_dir_z.push_back(std::numeric_limits<float>::lowest());

  _trk_end_x.push_back(std::numeric_limits<float>::lowest());
  _trk_end_y.push_back(std::numeric_limits<float>::lowest());
  _trk_end_z.push_back(std::numeric_limits<float>::lowest());

  _trk_length.push_back(std::numeric_limits<float>::lowest());

  _trk_bragg_p.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_mu.push_back(std::numeric_limits<float>::lowest());
  _trk_bragg_mip.push_back(std::numeric_limits<float>::lowest());

  _trk_pidchi.push_back(std::numeric_limits<float>::lowest());
  _trk_pidchipr.push_back(std::numeric_limits<float>::lowest());
  _trk_pidchika.push_back(std::numeric_limits<float>::lowest());
  _trk_pidchipi.push_back(std::numeric_limits<float>::lowest());
  _trk_pidchimu.push_back(std::numeric_limits<float>::lowest());
  _trk_pida.push_back(std::numeric_limits<float>::lowest());
  _trk_energy_proton.push_back(std::numeric_limits<float>::lowest());
  _trk_energy_muon.push_back(std::numeric_limits<float>::lowest());
}

void TrackAnalysis::setBranches(TTree *_tree)
{
  _tree->Branch("n_tracks", &_n_tracks, "n_tracks/i");
  _tree->Branch("trk_score_v", "std::vector<float>",  &_trk_score_v);
  _tree->Branch("trk_bragg_p_v", "std::vector< float >", &_trk_bragg_p);
  _tree->Branch("trk_bragg_mu_v", "std::vector< float >", &_trk_bragg_mu);
  _tree->Branch("trk_bragg_mip_v", "std::vector< float >", &_trk_bragg_mip);
  _tree->Branch("trk_pida_v", "std::vector< float >", &_trk_pida);
  _tree->Branch("trk_pid_chipr_v", "std::vector< float >", &_trk_pidchipr);
  _tree->Branch("trk_pid_chipi_v", "std::vector< float >", &_trk_pidchipi);
  _tree->Branch("trk_pid_chika_v", "std::vector< float >", &_trk_pidchika);
  _tree->Branch("trk_pid_chimu_v", "std::vector< float >", &_trk_pidchimu);
  _tree->Branch("trk_pfp_id", "std::vector< size_t >", &_trk_pfp_id);
  _tree->Branch("trk_dir_x", "std::vector< float >", &_trk_dir_x);
  _tree->Branch("trk_dir_y", "std::vector< float >", &_trk_dir_y);
  _tree->Branch("trk_dir_z", "std::vector< float >", &_trk_dir_z);

  _tree->Branch("trk_start_x", "std::vector< float >", &_trk_start_x);
  _tree->Branch("trk_start_y", "std::vector< float >", &_trk_start_y);
  _tree->Branch("trk_start_z", "std::vector< float >", &_trk_start_z);

  _tree->Branch("trk_end_x", "std::vector< float >", &_trk_end_x);
  _tree->Branch("trk_end_y", "std::vector< float >", &_trk_end_y);
  _tree->Branch("trk_end_z", "std::vector< float >", &_trk_end_z);
  _tree->Branch("trk_dist_v", "std::vector< float >", &_trk_distance);

  _tree->Branch("trk_theta_v", "std::vector< float >", &_trk_theta);
  _tree->Branch("trk_phi_v", "std::vector< float >", &_trk_phi);

  _tree->Branch("trk_len_v", "std::vector< float >", &_trk_length);
  _tree->Branch("trk_energy_proton", "std::vector< float >", &_trk_energy_proton);
  _tree->Branch("trk_energy_muon", "std::vector< float >", &_trk_energy_muon);
}

void TrackAnalysis::resetTTree(TTree *_tree)
{
  _n_tracks = 0;

  _trk_score_v.clear();
  _trk_bragg_p.clear();
  _trk_bragg_mu.clear();
  _trk_bragg_mip.clear();
  _trk_pida.clear();
  _trk_pidchipr.clear();
  _trk_pidchika.clear();
  _trk_pidchipi.clear();
  _trk_pidchimu.clear();
  _trk_pidchi.clear();
  _trk_pfp_id.clear();

  _trk_start_x.clear();
  _trk_start_y.clear();
  _trk_start_z.clear();

  _trk_end_x.clear();
  _trk_end_y.clear();
  _trk_end_z.clear();

  _trk_dir_x.clear();
  _trk_dir_y.clear();
  _trk_dir_z.clear();
  _trk_distance.clear();

  _trk_theta.clear();
  _trk_phi.clear();

  _trk_length.clear();

  _trk_energy_muon.clear();
  _trk_energy_proton.clear();
}

DEFINE_ART_CLASS_TOOL(TrackAnalysis)
} // namespace analysis

#endif
