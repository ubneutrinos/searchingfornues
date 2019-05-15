////////////////////////////////////////////////////////////////////////
// Class:       TruthFilter
// Plugin Type: filter (art v3_01_02)
// File:        TruthFilter_module.cc
//
// Generated at Thu May  2 21:20:54 2019 by David Caratelli using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"

#include "larcore/Geometry/Geometry.h"

#include <memory>

class TruthFilter;


class TruthFilter : public art::EDFilter {
public:
  explicit TruthFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TruthFilter(TruthFilter const&) = delete;
  TruthFilter(TruthFilter&&) = delete;
  TruthFilter& operator=(TruthFilter const&) = delete;
  TruthFilter& operator=(TruthFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  float fEnulow, fEnuhigh;
  float fProtonThreshold;
  bool fFiducialVolume;
  bool fCCNC;
  bool fNeutralCurrent;
  bool fChargedCurrent;
  int  fNpi0;

};


TruthFilter::TruthFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  fEnulow  = p.get<float>("Enulow" ,0  ); // GeV
  fEnuhigh = p.get<float>("Enuhigh",1e6); // GeV
  fProtonThreshold = p.get<float>("ProtonThreshold",0.04); // GeV
  fCCNC    = p.get<bool>("CCNC",false);
  fFiducialVolume = p.get<bool>("FiducialVolume",false);
  fNeutralCurrent = p.get<bool>("NeutralCurrent",false);
  fChargedCurrent = p.get<bool>("ChargedCurrent",false);
  fNpi0 = p.get<int>("Npi0",-1);
}

bool TruthFilter::filter(art::Event& e)
{

  art::ServiceHandle<geo::Geometry> geo;

  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");

  auto mct      = mct_h->at(0);
  auto neutrino = mct.GetNeutrino();
  auto nu       = neutrino.Nu();

  auto nu_e  = nu.Trajectory().E(0);
  auto vtx_x = nu.EndX();
  auto vtx_y = nu.EndY();
  auto vtx_z = nu.EndZ();

  // cut on energy
  if (nu_e < fEnulow) return false;
  if (nu_e > fEnuhigh) return false;
  // cut on FV
  if (fFiducialVolume == true) {
    if (vtx_x < 0) return false;
    if (vtx_x > 2 * geo->DetHalfWidth()) return false;
    if (vtx_y < -geo->DetHalfHeight()) return false;
    if (vtx_y > geo->DetHalfHeight()) return false;
    if (vtx_z < 0) return false;
    if (vtx_z > geo->DetLength()) return false;
  }
  // cut on final state
  if (fCCNC) {
    if (fNeutralCurrent && (neutrino.CCNC() == 0)) return false;
    if (fChargedCurrent && (neutrino.CCNC() == 1)) return false;
  }// if cut on CC or NC

  // loop through particles

  int nelec = 0;
  int nmuon = 0;
  int npi0 = 0;
  int nproton = 0; // with 40 MeV KE threshold
  int npion = 0;

  float protonHighE = 0;

  size_t npart = mct.NParticles();
  for (size_t i = 0; i < npart; i++)
  {
    auto const &part = mct.GetParticle(i);
    // if muon
    if ((std::abs(part.PdgCode()) == 13) and (part.StatusCode() == 1))
      {
	nmuon += 1;
      } // if muon
    // if electron
    if ((std::abs(part.PdgCode()) == 11) and (part.StatusCode() == 1))
      {
	nelec += 1;
      } // if electron
    // if pi0
    if ((part.PdgCode() == 111) and (part.StatusCode() == 1))
      {
	npi0 += 1;
      } // if pi0
    // if proton
    if ((part.PdgCode() == 2212) and (part.StatusCode() == 1))
      {
	// if highest energy, update energy
	if (part.Momentum(0).E() > protonHighE)
	  protonHighE = part.Momentum(0).E();
	if (part.Momentum(0).E() > fProtonThreshold)
	  nproton += 1;
      } // if proton
    // if pion
    if ((std::abs(part.PdgCode()) == 211) and (part.StatusCode() == 1))
      {
	npion += 1;
      } // if pion
  }// for all MCParticles
  
  if ( (fNpi0 > 0) && (npi0 == 0) ) return false;
  if ( (fNpi0 == 0) && (npi0 > 0) ) return false;

  return true;
}

void TruthFilter::beginJob()
{
  // Implementation of optional member function here.
}

void TruthFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(TruthFilter)
