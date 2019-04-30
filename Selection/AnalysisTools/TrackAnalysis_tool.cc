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
     * @brief set branches for TTree
     */
  void setBranches(TTree *_tree) override;

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree) override;

  /**
     * @brief return PID information
     */
  double PID(art::Ptr<anab::ParticleID> selected_pid,
             std::string AlgName,
             anab::kVariableType VariableType,
             anab::kTrackDir TrackDirection,
             int pdgCode,
             int selectedPlane);

private:
  trkf::TrackMomentumCalculator _trkmom;

  TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
  TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);

  art::InputTag fCALOproducer;
  art::InputTag fPIDproducer;
  art::InputTag fTRKproducer;

  std::vector<size_t> _trk_pfp_id;

  std::vector<float> _trkscore_v;

  std::vector<double> _trk_start_x;
  std::vector<double> _trk_start_y;
  std::vector<double> _trk_start_z;

  std::vector<double> _trk_theta;
  std::vector<double> _trk_phi;

  std::vector<double> _trk_dir_x;
  std::vector<double> _trk_dir_y;
  std::vector<double> _trk_dir_z;

  std::vector<double> _trk_end_x;
  std::vector<double> _trk_end_y;
  std::vector<double> _trk_end_z;

  std::vector<double> _trk_length;

  std::vector<double> _trk_bragg_p;
  std::vector<double> _trk_bragg_mu;
  std::vector<double> _trk_bragg_mip;

  std::vector<double> _trk_pidchi;
  std::vector<double> _trk_pidchipr;
  std::vector<double> _trk_pidchika;
  std::vector<double> _trk_pidchipi;
  std::vector<double> _trk_pidchimu;
  std::vector<double> _trk_pida;
  std::vector<double> _trk_energy_proton;
  std::vector<double> _trk_energy_muon;


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



  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
  {
    auto PDG = fabs(slice_pfp_v[i_pfp]->PdgCode());

    if (PDG == 12 || PDG == 14)
      continue;

    auto trk_v = slice_pfp_v[i_pfp].get<recob::Track>(); 

    auto ntrk = trk_v.size();

    if (ntrk != 1) continue;

    auto trk = trk_v.at(0);
    // get track proxy in order to fetch calorimtry
    auto trkpxy1 = calo_proxy[trk.key()];
    auto calopxy_v = trkpxy1.get<anab::Calorimetry>();
    std::cout << "There are "<< calopxy_v.size() << " associated calo objects " << std::endl;

    // get trk proxy in order to fetch PID
    auto trkpxy2 = pid_proxy[trk.key()];
    auto pidpxy_v = trkpxy2.get<anab::ParticleID>();
    std::cout << "There are "<< pidpxy_v.size() << " associated PID objects " << std::endl;
    double bragg_p = std::max(PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 2),
                              PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 2));

    double bragg_mu = std::max(PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 2),
                               PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 2));

    double bragg_mip = PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 2);

    double pidchipr = PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 0);
    double pidchimu = PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 0);
    double pidchipi = PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 0);
    double pidchika = PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 0);

    double pida_mean = PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 2);

    _trk_bragg_p.push_back(bragg_p);
    _trk_bragg_mu.push_back(bragg_mu);
    _trk_bragg_mip.push_back(bragg_mip);
    _trk_pida.push_back(pida_mean);
    _trk_pidchipr.push_back(pidchipr);
    _trk_pidchimu.push_back(pidchimu);
    _trk_pidchipi.push_back(pidchipi);
    _trk_pidchika.push_back(pidchika);

    // Kinetic energy using tabulated stopping power (GeV)
    double energy_proton = std::sqrt(std::pow(_trkmom.GetTrackMomentum(trk->Length(), 2212), 2) + std::pow(proton->Mass(), 2)) - proton->Mass();
    double energy_muon = std::sqrt(std::pow(_trkmom.GetTrackMomentum(trk->Length(), 13), 2) + std::pow(muon->Mass(), 2)) - muon->Mass();

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

    _trkscore_v.push_back(searchingfornues::GetTrackShowerScore(slice_pfp_v[i_pfp]));
    _trk_pfp_id.push_back(i_pfp);

  }// for all PFParticles

}

double TrackAnalysis::PID(art::Ptr<anab::ParticleID> selected_pid,
                             std::string AlgName,
                             anab::kVariableType VariableType,
                             anab::kTrackDir TrackDirection,
                             int pdgCode,
                             int selectedPlane)
{

  std::vector<anab::sParticleIDAlgScores> AlgScoresVec = selected_pid->ParticleIDAlgScores();
  for (size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
  {
    anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
    int planeid = -1;
    // std::cout << "AlgScore index " << i_algscore << " plane " << AlgScore.fPlaneMask<< std::endl;
    if (AlgScore.fPlaneMask.none() || AlgScore.fPlaneMask.count() > 1 || (AlgScore.fPlaneMask.count() == 1 && !(AlgScore.fPlaneMask.test(0) || AlgScore.fPlaneMask.test(1) || AlgScore.fPlaneMask.test(2))))
    {
      std::cout << "[uB_PlaneIDBitsetHelper] Cannot return a single MicroBooNE plane for bitset " << AlgScore.fPlaneMask << ". Returning -1 (invalid planeID)." << std::endl;
      continue;
    }
    else if (AlgScore.fPlaneMask.test(0))
      planeid = 2;
    else if (AlgScore.fPlaneMask.test(1))
      planeid = 1;
    else if (AlgScore.fPlaneMask.test(2))
      planeid = 0;
    else
      planeid = -1;

    if (selectedPlane != planeid) {
      continue;
    }

    if (AlgScore.fAlgName == AlgName)
    {
        if (anab::kVariableType(AlgScore.fVariableType) == VariableType && anab::kTrackDir(AlgScore.fTrackDir) == TrackDirection)
        {
            if (AlgScore.fAssumedPdg == pdgCode)
            {
                double alg_value = AlgScore.fValue;
                return alg_value;
            }
        }
    }
  }
  return std::numeric_limits<double>::lowest();
}

void TrackAnalysis::setBranches(TTree *_tree)
{
  _tree->Branch("trkscore_v", "std::vector<float>",  &_trkscore_v);
  _tree->Branch("trk_bragg_p", "std::vector< double >", &_trk_bragg_p);
  _tree->Branch("trk_bragg_mu", "std::vector< double >", &_trk_bragg_mu);
  _tree->Branch("trk_bragg_mip", "std::vector< double >", &_trk_bragg_mip);
  _tree->Branch("trk_pida", "std::vector< double >", &_trk_pida);
  _tree->Branch("trk_pid_chipr", "std::vector< double >", &_trk_pidchipr);
  _tree->Branch("trk_pid_chipi", "std::vector< double >", &_trk_pidchipi);
  _tree->Branch("trk_pid_chika", "std::vector< double >", &_trk_pidchika);
  _tree->Branch("trk_pid_chimu", "std::vector< double >", &_trk_pidchimu);
  _tree->Branch("trk_pfp_id", "std::vector< size_t >", &_trk_pfp_id);
  _tree->Branch("trk_dir_x", "std::vector< double >", &_trk_dir_x);
  _tree->Branch("trk_dir_y", "std::vector< double >", &_trk_dir_y);
  _tree->Branch("trk_dir_z", "std::vector< double >", &_trk_dir_z);

  _tree->Branch("trk_start_x", "std::vector< double >", &_trk_start_x);
  _tree->Branch("trk_start_y", "std::vector< double >", &_trk_start_y);
  _tree->Branch("trk_start_z", "std::vector< double >", &_trk_start_z);

  _tree->Branch("trk_end_x", "std::vector< double >", &_trk_end_x);
  _tree->Branch("trk_end_y", "std::vector< double >", &_trk_end_y);
  _tree->Branch("trk_end_z", "std::vector< double >", &_trk_end_z);

  _tree->Branch("trk_theta", "std::vector< double >", &_trk_theta);
  _tree->Branch("trk_phi", "std::vector< double >", &_trk_phi);

  _tree->Branch("trk_length", "std::vector< double >", &_trk_length);
  _tree->Branch("trk_energy_proton", "std::vector< double >", &_trk_energy_proton);
  _tree->Branch("trk_energy_muon", "std::vector< double >", &_trk_energy_muon);
}

void TrackAnalysis::resetTTree(TTree *_tree)
{
  _trkscore_v.clear();
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

  _trk_theta.clear();
  _trk_phi.clear();

  _trk_length.clear();

  _trk_energy_muon.clear();
  _trk_energy_proton.clear();
}

DEFINE_ART_CLASS_TOOL(TrackAnalysis)
} // namespace analysis

#endif
