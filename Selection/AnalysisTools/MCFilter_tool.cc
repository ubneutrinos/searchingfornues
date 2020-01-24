#ifndef ANALYSIS_MCFILTER_CXX
#define ANALYSIS_MCFILTER_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "larcore/Geometry/Geometry.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       MCFilter
// File:        MCFilter.cc
//
//              Analysis Tool to save information about MC filters
//
// Configuration parameters:
//
// TBD
//
// Created by Giuseppe Cerati (cerati@fnal.gov) on 03/15/2019
//
////////////////////////////////////////////////////////////////////////

class MCFilter : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  MCFilter(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~MCFilter(){};

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

private:

  art::InputTag fMCTproducer;
  art::InputTag fMCRproducer;

  float _mcf_nu_e; // neutrino energy
  float _mcf_lep_e; // lepton energy
  int _mcf_actvol; // is neutrino interaction in active volume
  int _mcf_nmm; // number of mu-
  int _mcf_nmp; // number of mu+
  int _mcf_nem; // number of e-
  int _mcf_nep; // number of e+
  int _mcf_np0; // number of pi0
  int _mcf_npp; // number of pi+
  int _mcf_npm; // number of pi-
  float _mcf_mcshr_elec_etot; // electron MCShower max energy
  //filter decisions
  int _mcf_pass_ccpi0; // pass ccpi0 filter
  int _mcf_pass_ncpi0; // pass ncpi0 filter
  int _mcf_pass_ccnopi; // pass ccnopi filter
  int _mcf_pass_ncnopi; // pass ncnopi filter
  int _mcf_pass_cccpi; // pass cccpi filter
  int _mcf_pass_nccpi; // pass nccpi filter
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
MCFilter::MCFilter(const fhicl::ParameterSet &p)
{
  fMCTproducer = p.get<art::InputTag>("MCTproducer");
  fMCRproducer = p.get<art::InputTag>("MCRproducer");
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void MCFilter::configure(fhicl::ParameterSet const &p) {}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void MCFilter::analyzeEvent(art::Event const &e, bool fData)
{
  if (fData) return;

  // load MCTruth
  auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth> >(fMCTproducer);

  auto mct = mct_h->at(0);
  auto neutrino = mct.GetNeutrino();
  auto nu = neutrino.Nu();

  _mcf_nu_e = nu.Trajectory().E(0);
  _mcf_lep_e = neutrino.Lepton().E();

  art::ServiceHandle<geo::Geometry const> geo;
  auto const& tpcActiveBox = geo->TPC().ActiveBoundingBox();
  _mcf_actvol = tpcActiveBox.ContainsPosition(geo::Point_t(nu.Vx(),nu.Vy(),nu.Vz()));

  size_t npart = mct.NParticles();
  for (size_t i = 0; i < npart; i++) {
    //
    auto const &part = mct.GetParticle(i);
    //
    if (part.StatusCode() != 1) continue;
    //
    // if muon
    if (part.PdgCode() == 13) _mcf_nmm += 1;
    if (part.PdgCode() == -13) _mcf_nmp += 1;
    // if electron
    if (part.PdgCode() == 11) _mcf_nem += 1;
    if (part.PdgCode() == -11) _mcf_nep += 1;
    // if pi0
    if (part.PdgCode() == 111) _mcf_np0 += 1;
    // if pi+/-
    if (part.PdgCode() == 211) _mcf_npp += 1;
    if (part.PdgCode() == -211) _mcf_npm += 1;
  } // for all MCParticles

  float maxElecMCShwEMeV = 0.;
  const std::vector<sim::MCShower> &inputMCShower = *(e.getValidHandle<std::vector<sim::MCShower> >(fMCRproducer));
  for (auto& mcs : inputMCShower) {
    if (std::abs(mcs.PdgCode())==11) {
      if (mcs.Start().E()>maxElecMCShwEMeV) {
	_mcf_mcshr_elec_etot = mcs.Start().E();
	maxElecMCShwEMeV = _mcf_mcshr_elec_etot;
      }
    }
  }

  // now replicate the filter conditions
  _mcf_pass_ccpi0 = 0;
  if (_mcf_actvol==1 && _mcf_nmm==1 && _mcf_nem==0 && _mcf_nep==0 && _mcf_np0==1) _mcf_pass_ccpi0 = 1;
  //
  _mcf_pass_ncpi0 = 0;
  if (_mcf_actvol==1 && _mcf_nmm==0 && _mcf_nmp==0 && _mcf_nem==0 && _mcf_nep==0 && _mcf_np0==1) _mcf_pass_ncpi0 = 1;
  //
  _mcf_pass_ccnopi = 0;
  if (_mcf_actvol==1 && _mcf_nmm==1 && _mcf_nem==0 && _mcf_nep==0 &&
      _mcf_np0==0 && _mcf_npp==0 && _mcf_npm==0 &&
      (((_mcf_lep_e-0.105)>0.02 && _mcf_lep_e<0.3) || _mcf_mcshr_elec_etot>15)) _mcf_pass_ccnopi = 1;
  //
  _mcf_pass_ncnopi = 0;
  if (_mcf_actvol==1 && _mcf_nmm==0 && _mcf_nmp==0 && _mcf_nem==0 && _mcf_nep==0 &&
      _mcf_np0==0 && _mcf_npp==0 && _mcf_npm==0 && _mcf_nu_e>0.9) _mcf_pass_ncnopi = 1;
  //
  _mcf_pass_cccpi = 0;
  if (_mcf_actvol==1 && _mcf_nmm==1 && _mcf_nem==0 && _mcf_nep==0 &&
      _mcf_np0==0 && (_mcf_npp==1 || _mcf_npm==1) &&
      (((_mcf_lep_e-0.105)>0.02 && _mcf_lep_e<0.4) || _mcf_mcshr_elec_etot>35)) _mcf_pass_cccpi = 1;
  //
  _mcf_pass_nccpi = 0;
  if (_mcf_actvol==1 && _mcf_nmm==0 && _mcf_nmp==0 && _mcf_nem==0 && _mcf_nep==0 &&
      _mcf_np0==0 && (_mcf_npp==1 || _mcf_npm==1)) _mcf_pass_nccpi = 1;
  //
}

void MCFilter::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected){}

void MCFilter::setBranches(TTree *_tree)
{
  _tree->Branch("mcf_nu_e", &_mcf_nu_e, "mcf_nu_e/F");
  _tree->Branch("mcf_lep_e", &_mcf_lep_e, "mcf_lep_e/F");
  _tree->Branch("mcf_actvol", &_mcf_actvol, "mcf_actvol/I");
  _tree->Branch("mcf_nmm", &_mcf_nmm, "mcf_nmm/I");
  _tree->Branch("mcf_nmp", &_mcf_nmp, "mcf_nmp/I");
  _tree->Branch("mcf_nem", &_mcf_nem, "mcf_nem/I");
  _tree->Branch("mcf_nep", &_mcf_nep, "mcf_nep/I");
  _tree->Branch("mcf_np0", &_mcf_np0, "mcf_np0/I");
  _tree->Branch("mcf_npp", &_mcf_npp, "mcf_npp/I");
  _tree->Branch("mcf_npm", &_mcf_npm, "mcf_npm/I");
  _tree->Branch("mcf_mcshr_elec_etot", &_mcf_mcshr_elec_etot, "mcf_mcshr_elec_etot/F");
  _tree->Branch("mcf_pass_ccpi0", &_mcf_pass_ccpi0, "mcf_pass_ccpi0/I");
  _tree->Branch("mcf_pass_ncpi0", &_mcf_pass_ncpi0, "mcf_pass_ncpi0/I");
  _tree->Branch("mcf_pass_ccnopi", &_mcf_pass_ccnopi, "mcf_pass_ccnopi/I");
  _tree->Branch("mcf_pass_ncnopi", &_mcf_pass_ncnopi, "mcf_pass_ncnopi/I");
  _tree->Branch("mcf_pass_cccpi", &_mcf_pass_cccpi, "mcf_pass_cccpi/I");
  _tree->Branch("mcf_pass_nccpi", &_mcf_pass_nccpi, "mcf_pass_nccpi/I");
}

void MCFilter::resetTTree(TTree *_tree)
{
  _mcf_nu_e = -1.;
  _mcf_lep_e = -1;
  _mcf_actvol = -1;
  _mcf_nmm = 0;
  _mcf_nmp = 0;
  _mcf_nem = 0;
  _mcf_nep = 0;
  _mcf_np0 = 0;
  _mcf_npp = 0;
  _mcf_npm = 0;
  _mcf_mcshr_elec_etot = -1.;
  _mcf_pass_ccpi0 = -1;
  _mcf_pass_ncpi0 = -1;
  _mcf_pass_ccnopi = -1;
  _mcf_pass_ncnopi = -1;
  _mcf_pass_cccpi = -1;
  _mcf_pass_nccpi = -1;
}

DEFINE_ART_CLASS_TOOL(MCFilter)
} // namespace analysis

#endif
