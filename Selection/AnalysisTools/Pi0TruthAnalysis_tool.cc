#ifndef ANALYSIS_SHOWERANALYSIS_CXX
#define ANALYSIS_SHOWERANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       Pi0TruthAnalysis
// File:        Pi0TruthAnalysis.cc
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

class Pi0TruthAnalysis : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  Pi0TruthAnalysis(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~Pi0TruthAnalysis(){};

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
  art::InputTag fMCSproducer;

  float _gamma1_edep;
  float _gamma1_etot;
  float _gamma2_edep;
  float _gamma2_etot;
  
  int   _gamma_parent;
  float _gamma_edep;
  float _gamma_etot;
  float _gamma_dist;

  int   _elec_parent;
  float _elec_edep;
  float _elec_etot;
  float _elec_dist;

};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
Pi0TruthAnalysis::Pi0TruthAnalysis(const fhicl::ParameterSet &p)
{
  fMCTproducer = p.get<art::InputTag>("MCTproducer", "");
  fMCSproducer = p.get<art::InputTag>("MCSproducer", "");
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void Pi0TruthAnalysis::configure(fhicl::ParameterSet const &p)
{
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void Pi0TruthAnalysis::analyzeEvent(art::Event const &e, bool fData)
{

  if (fData == true) {
    return;
  }

  // load MCTruth
  auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);

  // load MCShowers
  auto const &mcs_h = e.getValidHandle<std::vector<sim::MCShower>>(fMCSproducer);

  auto mct = mct_h->at(0);
  auto neutrino = mct.GetNeutrino();
  auto nu = neutrino.Nu();

  float vtx_x = nu.EndX();
  float vtx_y = nu.EndY();
  float vtx_z = nu.EndZ();

  // nucleus induced pi0 gammas
  _gamma1_edep = 0;
  _gamma1_etot = 0;
  _gamma2_edep = 0;
  _gamma2_etot = 0;
  
  _gamma_parent = 0;
  _gamma_edep = 0;
  _gamma_etot = 0;
  _gamma_dist = 0;

  _elec_parent = 0;
  _elec_edep = 0;
  _elec_etot = 0;
  _elec_dist = 0;

  for (size_t i=0; i < mcs_h->size(); i++){
    auto const& mcs = mcs_h->at(i);
    // distance from vertex                                                                
    double x = mcs.Start().X();
    double y = mcs.Start().Y();
    double z = mcs.Start().Z();
    double d = sqrt( ( (vtx_x - x) * (vtx_x - x) ) +
		     ( (vtx_y - y) * (vtx_y - y) ) +
		     ( (vtx_z - z) * (vtx_z - z) ) );

    auto etot = mcs.Start().E();
    auto edep = mcs.DetProfile().E();

    // if originating from the nucleus
    if ( d < 0.01 ){

      if (etot > _gamma1_etot) {
	_gamma2_etot = _gamma1_etot;
	_gamma2_edep = _gamma1_edep;
	_gamma1_etot = etot;
	_gamma1_edep = edep;
      }
      else if (etot > _gamma2_etot) {
	_gamma2_etot = etot;
	_gamma2_edep = edep;
      }

    }// if originating from the neutrino

    // check for all gammas / electrons
    if (fabs(mcs.PdgCode()) == 11) {
      if (edep > _elec_edep) {
	_elec_edep = edep;
	_elec_etot = etot;
	_elec_parent = mcs.MotherPdgCode();
	_elec_dist = d;
      }
    }// if electron
    if (fabs(mcs.PdgCode()) == 22) {
      if (edep > _gamma_edep) {
	_gamma_edep = edep;
	_gamma_etot = etot;
	_gamma_parent = mcs.MotherPdgCode();
	_gamma_dist = d;
      }
    }// if photon
	
  }// for all MCShowers                                                               

}

void Pi0TruthAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{

  return;
}

void Pi0TruthAnalysis::setBranches(TTree *_tree)
{

  _tree->Branch("gamma_edep", &_gamma_edep, "gamma_edep/F");
  _tree->Branch("gamma_etot", &_gamma_etot, "gamma_etot/F");
  _tree->Branch("gamma_dist", &_gamma_dist, "gamma_dist/F");
  _tree->Branch("gamma_parent", &_gamma_parent, "gamma_parent/I");

  _tree->Branch("elec_edep", &_elec_edep, "elec_edep/F");
  _tree->Branch("elec_etot", &_elec_etot, "elec_etot/F");
  _tree->Branch("elec_dist", &_elec_dist, "elec_dist/F");
  _tree->Branch("elec_parent", &_elec_parent, "elec_parent/I");

  _tree->Branch("gamma1_edep", &_gamma1_edep, "gamma1_edep/F");
  _tree->Branch("gamma1_etot", &_gamma1_etot, "gamma1_etot/F");

  _tree->Branch("gamma2_edep", &_gamma2_edep, "gamma2_edep/F");
  _tree->Branch("gamma2_etot", &_gamma2_etot, "gamma2_etot/F");

}

void Pi0TruthAnalysis::resetTTree(TTree *_tree)
{

  _gamma1_edep = 0;
  _gamma1_etot = 0;
  _gamma2_edep = 0;
  _gamma2_etot = 0;
  
  _gamma_parent = 0;
  _gamma_edep = 0;
  _gamma_etot = 0;
  _gamma_dist = 0;

  _elec_parent = 0;
  _elec_edep = 0;
  _elec_etot = 0;
  _elec_dist = 0;

}

DEFINE_ART_CLASS_TOOL(Pi0TruthAnalysis)
} // namespace analysis

#endif
