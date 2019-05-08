#ifndef ANALYSIS_DEFAULTANALYSIS_CXX
#define ANALYSIS_DEFAULTANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcore/Geometry/Geometry.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       DefaultAnalysis
// File:        DefaultAnalysis.cc
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

class DefaultAnalysis : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  DefaultAnalysis(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~DefaultAnalysis(){};

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
  /**
     * @brief apply SCE corections from reconstruction maps
     */
  void ApplySCECorrection(const float &vtx_x, const float &vtx_y, const float &vtx_z,
                          float &sce_x, float &sce_y, float &sce_z);

  /**
   * @brief Determine if the specified point is in the fiducial volume
   *        Not recommended, no array size checking is done.
   *
   * @param x array of 3D location
   * @return True if the point is inside the fiducial volume
   */
  bool isFiducial(const double x[3]) const;

  art::InputTag fCRTVetoproducer; // producer for CRT veto ass tag [anab::T0 <-> recob::OpFlash]
  art::InputTag fCLSproducer;     // cluster associated to PFP
  art::InputTag fMCTproducer;     // MCTruth from neutrino generator
  art::InputTag fMCPproducer;     // MCParticle from Geant4 stage
  art::InputTag fBacktrackTag;
  art::InputTag fHproducer;
  art::InputTag fMCRproducer;
  art::InputTag fSLCproducer; // slice associated to PFP

  float fFidvolXstart;
  float fFidvolXend;
  float fFidvolYstart;
  float fFidvolYend;
  float fFidvolZstart;
  float fFidvolZend;

  const int k_nu_e_other = 1;
  const int k_nu_e_cc0pinp = 11;
  const int k_nu_mu_other = 2;
  const int k_nu_mu_pi0 = 21;
  const int k_nc = 3;
  const int k_nc_pi0 = 31;
  const int k_cosmic = 4;
  const int k_outfv = 5;
  const int k_other = 6;
  const int k_data = 0;
  // kinematic thresholds to define signal
  float fProtonThreshold;

  int _category; // event category

  // neutrino vertx (reco)
  float _nu_vtx_x, _nu_vtx_y, _nu_vtx_z;
  // neutrino vertex SCE position corrections (reco)
  float _nu_sce_x, _nu_sce_y, _nu_sce_z;

  int _run, _sub, _evt; // event info
  // neutrino information
  float _nu_e;                  // neutrino energy [GeV]
  int _nu_pdg;                  // neutrino PDG code
  int _ccnc;                    // CC or NC tag from GENIE
  float _vtx_x, _vtx_y, _vtx_z; // neutrino interaction vertex coordinates [cm]
  float _vtx_t;                 // neutrino generation time
  bool _isVtxInActive;          // true if neutrino in active volume, 0 < x < 256 -116 < y < 116;  0 < z <  1036
  bool _isVtxInFiducial;        // true if neutrino in fiducial volume

  // final state particle information
  int _nmuon;                            // is there a final-state muon from the neutrino? [1=yes 0=no]
  float _muon_e, _muon_p, _muon_c;       // energy, purity, completeness.
  int _nelec;                            // is there a final-state electron from the neutrino? [1=yes 0=no]
  float _elec_e, _elec_p, _elec_c;       // energy, purity, completeness.
  int _npi0;                             // how many pi0s are there?
  int _pi0;                              // is there a final-state pi0 from the neutrino? [1=yes 0=no]
  float _pi0_e, _pi0_p, _pi0_c;          // energy, purity, completeness.
  int _nproton;                          // how many protons are there?
  int _proton;                           // is there a final-state proton from the neutrino? [1=yes 0=no]
  float _proton_e, _proton_p, _proton_c; // energy, purity, completeness.
  int _npion;                            // how many pions are there?
  int _pion;                             // is there a final-state charged pion from the neutrino? [1=yes 0=no]
  float _pion_e, _pion_p, _pion_c;       // energy, purity, completeness.

  // end muon process
  std::string _endmuonprocess;
  float _endmuonmichel;

  // number of slices in the event
  int _nslice;
  int _crtveto;    // is the event vetoed by the CRT Veto?
  float _crthitpe; // pe associated to CRT hit

  //
  std::vector<int> _pfp_slice_idx; // index of PFP is vector of PFPs in nu slice

  // reco PFParticle backtracking. One entry for PFParticle in the slice
  // std::vector<int>   _backtracked_idx;    // index of PFP [key]
  // std::vector<int>   _backtracked_tid;    // TrackID of backtracked MCParticle
  std::vector<int> _backtracked_pdg;            // PDG code of backtracked particle
  std::vector<float> _backtracked_e;            // energy of backtracked particle
  std::vector<float> _backtracked_purity;       // purity of backtracking
  std::vector<float> _backtracked_completeness; // completeness of backtracking

  float _lep_e;                                              // lepton energy (if one exists) [GeV]
  int _pass;                                                 // does the slice pass the selection
  float _xtimeoffset, _xsceoffset, _ysceoffset, _zsceoffset; // offsets for generation time and SCE

  int evnhits;                                // number of hits in event
  int slpdg;                                  // PDG code of primary pfp in slice
  int slnhits;                                // number of hits in slice
  std::vector<int> pfpdg;                     // PDG code of pfp in slice
  std::vector<int> pfnhits;                   // number of hits in pfp
  std::vector<std::vector<int>> pfnplanehits; // number of hits in pfp
  unsigned int _hits_u;
  unsigned int _hits_v;
  unsigned int _hits_y;


  std::vector<int> _mc_pdg;
  std::vector<double> _mc_E;

  std::vector<double> _mc_px;
  std::vector<double> _mc_py;
  std::vector<double> _mc_pz;

  std::vector<double> _mc_vx;
  std::vector<double> _mc_vy;
  std::vector<double> _mc_vz;

  std::vector<double> _mc_endx;
  std::vector<double> _mc_endy;
  std::vector<double> _mc_endz;
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
DefaultAnalysis::DefaultAnalysis(const fhicl::ParameterSet &p)
{
  fCRTVetoproducer = p.get<art::InputTag>("CRTVetoproducer", ""); // default is no CRT veto
  fCLSproducer = p.get<art::InputTag>("CLSproducer");
  fMCTproducer = p.get<art::InputTag>("MCTproducer");
  fMCPproducer = p.get<art::InputTag>("MCPproducer");
  fBacktrackTag = p.get<art::InputTag>("BacktrackTag");
  fHproducer = p.get<art::InputTag>("Hproducer");
  fMCRproducer = p.get<art::InputTag>("MCRproducer");
  fSLCproducer = p.get<art::InputTag>("SLCproducer");
  // kinematic thresholds for defining signal
  fProtonThreshold = p.get<float>("ProtonThreshold");

  fFidvolXstart = p.get<double>("fidvolXstart", 10);
  fFidvolXend = p.get<double>("fidvolXend", 10);

  fFidvolYstart = p.get<double>("fidvolYstart", 15);
  fFidvolYend = p.get<double>("fidvolYend", 15);

  fFidvolZstart = p.get<double>("fidvolZstart", 10);
  fFidvolZend = p.get<double>("fidvolZend", 50);
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void DefaultAnalysis::configure(fhicl::ParameterSet const &p)
{
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void DefaultAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();

  if (!fData)
  {
    // SaveTruth
    SaveTruth(e);
  }

  // Grab CRT veto information if available - CRT should probably have its own tool?
  if (fCRTVetoproducer != "")
  {
    art::Handle<art::Assns<crt::CRTHit, recob::OpFlash, void>> crtveto_h;
    e.getByLabel(fCRTVetoproducer, crtveto_h);
    _crtveto = crtveto_h->size();
    if (_crtveto == 1)
      _crthitpe = crtveto_h->at(0).first->peshit;
  } // if the CRT veto label has been defined

  art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
  evnhits = inputHits->size();
}

bool DefaultAnalysis::isFiducial(const double x[3]) const
{

  art::ServiceHandle<geo::Geometry> geo;
  std::vector<double> bnd = {
      0., 2. * geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(),
      0., geo->DetLength()};

  bool is_x =
      x[0] > (bnd[0] + fFidvolXstart) && x[0] < (bnd[1] - fFidvolXend);
  bool is_y =
      x[1] > (bnd[2] + fFidvolYstart) && x[1] < (bnd[3] - fFidvolYend);
  bool is_z =
      x[2] > (bnd[4] + fFidvolZstart) && x[2] < (bnd[5] - fFidvolZend);

  return is_x && is_y && is_z;
}

void DefaultAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{

  ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
                                                                                        proxy::withAssociated<recob::Hit>(fCLSproducer));
  // somehow proxies don't work for the slice-hit association, so go back to old assns
  art::ValidHandle<std::vector<recob::Slice>> inputSlice = e.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
  auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSlice, e, fSLCproducer));

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

  size_t pfpidx = 0;
  for (auto pfp : slice_pfp_v)
  {

    if (pfp->IsPrimary())
    {
      slpdg = pfp->PdgCode();
      auto slice_pxy_v = pfp.get<recob::Slice>();
      if (slice_pxy_v.size() != 1)
      {
        std::cout << "WRONG!!! n slices = " << slice_pxy_v.size() << " " << __FILE__ << " " << __LINE__ << std::endl;
        return;
      }
      auto slicehits = assocSliceHit->at(slice_pxy_v[0].key());
      slnhits = slicehits.size();
    }

    if ((pfp->PdgCode() == 12) || (pfp->PdgCode() == 14))
    {

      // grab vertex
      Double_t xyz[3] = {};

      auto vtx = pfp.get<recob::Vertex>();
      if (vtx.size() != 1)
      {
        std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
      }

      else
      {
        // save vertex to array
        vtx.at(0)->XYZ(xyz);
        auto nuvtx = TVector3(xyz[0], xyz[1], xyz[2]);

        _nu_vtx_x = nuvtx.X();
        _nu_vtx_y = nuvtx.Y();
        _nu_vtx_z = nuvtx.Z();

        ApplySCECorrection(_nu_vtx_x, _nu_vtx_y, _nu_vtx_z, _nu_sce_x, _nu_sce_y, _nu_sce_z);
      }
    } // if neutrino PFParticle

    _pfp_slice_idx.push_back(pfpidx++);
    pfpdg.push_back(pfp->PdgCode());
    std::vector<int> nplanehits(3, 0);
    // get hits associated to this PFParticle through the clusters
    std::vector<art::Ptr<recob::Hit>> hit_v;
    auto clus_pxy_v = pfp.get<recob::Cluster>();
    for (auto ass_clus : clus_pxy_v)
    {
      // get cluster proxy
      const auto &clus = clus_proxy[ass_clus.key()];
      auto clus_hit_v = clus.get<recob::Hit>();
      nplanehits[clus->Plane().Plane] = clus_hit_v.size();


      if (clus->Plane().Plane == 0) {
        _hits_u += clus_hit_v.size();
      } else if (clus->Plane().Plane == 1) {
        _hits_v += clus_hit_v.size();
      } else if (clus->Plane().Plane == 2) {
        _hits_y += clus_hit_v.size();
      }

      for (const auto &hit : clus_hit_v)
      {
        hit_v.push_back(hit);
      }
    } // for all clusters associated to PFP
    pfnhits.push_back(hit_v.size());
    pfnplanehits.push_back(nplanehits);

    // backtrack PFParticles in the slice
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
          //_backtracked_idx.push_back(pfp->Self());
          //_backtracked_tid.push_back(mcp->TrackId());
          _backtracked_e.push_back(mcp.e);
          _backtracked_pdg.push_back(PDG);
          _backtracked_purity.push_back(purity);
          _backtracked_completeness.push_back(completeness);
          // if this is an interesting particle, save it to the TTree
          if (fabs(PDG) == 13)
          {
            if (fabs(mcp.e - _muon_e) < 0.01)
            {
              _muon_p = purity;
              _muon_c = completeness;
            }
          }
          if (fabs(PDG) == 11)
          {
            if (fabs(mcp.e - _elec_e) < 0.01)
            {
              _elec_p = purity;
              _elec_c = completeness;
            }
          }
          if (fabs(PDG) == 2212)
          {
            if (fabs(mcp.e - _proton_e) < 0.01)
            {
              _proton_p = purity;
              _proton_c = completeness;
            }
          }
        }
        else
        {
          // _backtracked_idx.push_back(0);
          // _backtracked_tid.push_back(0);
          _backtracked_e.push_back(0);
          _backtracked_pdg.push_back(0);
          _backtracked_purity.push_back(0.);
          _backtracked_completeness.push_back(0.);
        }
      } // if there are associated clusters
    }   // if MC
  }
  _nslice += 1;

  if (!fData)
  {

    // Check if there is a PFParticle associated to an electron
    std::vector<int>::iterator e_reco_id = std::find(_backtracked_pdg.begin(), _backtracked_pdg.end(), 11);
    bool there_is_reco_electron = e_reco_id != _backtracked_pdg.end();

    // Check if there is a PFParticle associated to an overlay cosmic
    std::vector<int>::iterator cosmic_id = std::find(_backtracked_pdg.begin(), _backtracked_pdg.end(), 0);
    bool there_is_reco_cosmic = cosmic_id != _backtracked_pdg.end();

    std::vector<int>::iterator e_true_id = std::find(_mc_pdg.begin(), _mc_pdg.end(), 11);
    bool there_is_true_electron = e_true_id != _mc_pdg.end();

    bool there_is_true_proton = false;
    bool there_is_true_pi = false;
    bool there_is_true_mu = false;
    bool there_is_true_pi0 = false;

    for (size_t i_pdg = 0; i_pdg < _mc_pdg.size(); i_pdg++)
    {
      if (abs(_mc_pdg[i_pdg]) == 2212)
      {
        double ke = _mc_E[i_pdg] - 0.938;
        if (ke > 0.06)
        {
          there_is_true_proton = true;
        }
      }

      if (abs(_mc_pdg[i_pdg]) == 211 || _mc_pdg[i_pdg] == 111)
      {
        double ke = _mc_E[i_pdg] - 0.135;
        if (ke > 0.04)
        {
          if (_mc_pdg[i_pdg] == 111)
            there_is_true_pi0 = true;
          there_is_true_pi = true;
        }
      }

      if (abs(_mc_pdg[i_pdg]) == 13)
      {
        double ke = _mc_E[i_pdg] - 0.105;
        if (ke > 0.04)
        {
          there_is_true_mu = true;
        }
      }
    }

    if (!there_is_reco_cosmic && !_isVtxInFiducial)
    {
      _category = k_outfv;
    }
    else if (_nu_pdg == 12 && !there_is_reco_cosmic)
    {
      if (there_is_true_electron)
      {
        if (there_is_reco_electron)
        {
          if (!there_is_true_pi && there_is_true_proton)
          {
            _category = k_nu_e_cc0pinp;
          }
          else
          {
            _category = k_nu_e_other;
          }
        }
        else
        {
          _category = k_other;
        }
      }
      else
      {
        if (!there_is_true_pi)
        {
          _category = k_nc;
        }
        else
        {
          _category = k_nc_pi0;
        }
      }
    }
    else if (_nu_pdg == 14 && !there_is_reco_cosmic)
    {
      if (there_is_true_mu)
      {
        if (there_is_true_pi0)
        {
          _category = k_nu_mu_pi0;
        }
        else
        {
          _category = k_nu_mu_other;
        }
      }
      else
      {
        if (!there_is_true_pi)
        {
          _category = k_nc;
        }
        else
        {
          _category = k_nc_pi0;
        }
      }
    }
    else
    {
      _category = k_cosmic;
    }
  }
  else
  {
    _category = k_data;
  }
  if (selected)
    _pass = 1;
}

void DefaultAnalysis::setBranches(TTree *_tree)
{
  // reconstructed neutrino vertex
  _tree->Branch("nu_vtx_x", &_nu_vtx_x, "nu_vtx_x/F");
  _tree->Branch("nu_vtx_y", &_nu_vtx_y, "nu_vtx_y/F");
  _tree->Branch("nu_vtx_z", &_nu_vtx_z, "nu_vtx_z/F");
  _tree->Branch("nu_sce_x", &_nu_sce_x, "nu_sce_x/F");
  _tree->Branch("nu_sce_y", &_nu_sce_y, "nu_sce_y/F");
  _tree->Branch("nu_sce_z", &_nu_sce_z, "nu_sce_z/F");

  // neutrino information
  _tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
  _tree->Branch("ccnc", &_ccnc, "ccnc/I");
  _tree->Branch("nu_e", &_nu_e, "nu_e/F");
  _tree->Branch("vtx_x", &_vtx_x, "vtx_x/F");
  _tree->Branch("vtx_y", &_vtx_y, "vtx_y/F");
  _tree->Branch("vtx_z", &_vtx_z, "vtx_z/F");
  _tree->Branch("isVtxInActive", &_isVtxInActive, "isVtxInActive/O");
  _tree->Branch("isVtxInFiducial", &_isVtxInFiducial, "isVtxInFiducial/O");
  // individual particles in the neutrino slice
  // legend:
  // _e -> energy of particle in GeV
  // _c -> completeness from back-tracking [0,1]
  // _p -> purity from back-tracking [0,1]
  // muon
  _tree->Branch("nmuon", &_nmuon, "nmuon/I");
  _tree->Branch("muon_e", &_muon_e, "muon_e/F");
  _tree->Branch("muon_c", &_muon_c, "muon_c/F");
  _tree->Branch("muon_p", &_muon_p, "muon_p/F");
  // electron
  _tree->Branch("nelec", &_nelec, "nelec/I");
  _tree->Branch("elec_e", &_elec_e, "elec_e/F");
  _tree->Branch("elec_c", &_elec_c, "elec_c/F");
  _tree->Branch("elec_p", &_elec_p, "elec_p/F");
  // pi0
  _tree->Branch("npi0", &_npi0, "npi0/I");
  _tree->Branch("pi0_e", &_pi0_e, "pi0_e/F");
  _tree->Branch("pi0_c", &_pi0_c, "pi0_c/F");
  _tree->Branch("pi0_p", &_pi0_p, "pi0_p/F");
  // first [highest momentum] proton
  _tree->Branch("nproton", &_nproton, "nproton/I");
  _tree->Branch("proton_e", &_proton_e, "proton_e/F");
  _tree->Branch("proton_c", &_proton_c, "proton_c/F");
  _tree->Branch("proton_p", &_proton_p, "proton_p/F");
  // charged pions
  _tree->Branch("npion", &_npion, "npion/I");
  _tree->Branch("pion_e", &_pion_e, "pion_e/F");
  _tree->Branch("pion_c", &_pion_c, "pion_c/F");
  _tree->Branch("pion_p", &_pion_p, "pion_p/F");

  _tree->Branch("nslice", &_nslice, "nslice/I");
  _tree->Branch("crtveto", &_crtveto, "crtveto/I");
  _tree->Branch("crthitpe", &_crthitpe, "crthitpe/F");

  _tree->Branch("pfp_slice_idx", "std::vector<int>", &_pfp_slice_idx);

  // PFParticle backtracking
  // _tree->Branch("backtracked_idx"   ,"std::vector<int>"  ,&_backtracked_idx   );
  // _tree->Branch("backtracked_tid"   ,"std::vector<int>"  ,&_backtracked_tid   );

  _tree->Branch("category", &_category, "category/I");

  _tree->Branch("backtracked_pdg", "std::vector<int>", &_backtracked_pdg);
  _tree->Branch("backtracked_e", "std::vector<float>", &_backtracked_e);
  _tree->Branch("backtracked_purity", "std::vector<float>", &_backtracked_purity);
  _tree->Branch("backtracked_completeness", "std::vector<float>", &_backtracked_completeness);

  _tree->Branch("lep_e", &_lep_e, "lep_e/F");
  _tree->Branch("pass", &_pass, "pass/I");
  _tree->Branch("run", &_run, "run/I");
  _tree->Branch("sub", &_sub, "sub/I");
  _tree->Branch("evt", &_evt, "evt/I");

  _tree->Branch("xtimeoffset", &_xtimeoffset, "xtimeoffset/F");
  _tree->Branch("xsceoffset", &_xsceoffset, "xsceoffset/F");
  _tree->Branch("ysceoffset", &_ysceoffset, "ysceoffset/F");
  _tree->Branch("zsceoffset", &_zsceoffset, "zsceoffset/F");

  _tree->Branch("evnhits", &evnhits, "evnhits/I");
  _tree->Branch("slpdg", &slpdg, "slpdg/I");
  _tree->Branch("slnhits", &slnhits, "slnhits/I");
  _tree->Branch("pfpdg", &pfpdg);
  _tree->Branch("pfnhits", &pfnhits);
  _tree->Branch("pfnplanehits", &pfnplanehits);
  _tree->Branch("hits_u", &_hits_u, "hits_u/i");
  _tree->Branch("hits_v", &_hits_v, "hits_v/i");
  _tree->Branch("hits_y", &_hits_y, "hits_y/i");

  _tree->Branch("mc_pdg", "std::vector< int >", &_mc_pdg);
  _tree->Branch("mc_E", "std::vector< double >", &_mc_E);

  _tree->Branch("mc_vx", "std::vector< double >",
                &_mc_vx);
  _tree->Branch("mc_vy", "std::vector< double >",
                &_mc_vy);
  _tree->Branch("mc_vz", "std::vector< double >",
                &_mc_vz);

  _tree->Branch("mc_endx", "std::vector< double >",
                &_mc_endx);
  _tree->Branch("mc_endy", "std::vector< double >",
                &_mc_endy);
  _tree->Branch("mc_endz", "std::vector< double >",
                &_mc_endz);

  _tree->Branch("mc_px", "std::vector< double >",
                &_mc_px);
  _tree->Branch("mc_py", "std::vector< double >",
                &_mc_py);
  _tree->Branch("mc_pz", "std::vector< double >",
                &_mc_pz);

  _tree->Branch("endmuonprocess", &_endmuonprocess);

  _tree->Branch("endmuonmichel", &_endmuonmichel, "endmuonmichel/F");
}

void DefaultAnalysis::resetTTree(TTree *_tree)
{
  _run = std::numeric_limits<int>::min();
  _sub = std::numeric_limits<int>::min();
  _evt = std::numeric_limits<int>::min();
  _nu_e = std::numeric_limits<float>::min();
  _nu_pdg = std::numeric_limits<int>::min();
  _ccnc = std::numeric_limits<int>::min();
  _pass = 0;
  _vtx_x = std::numeric_limits<float>::min();
  _vtx_y = std::numeric_limits<float>::min();
  _vtx_z = std::numeric_limits<float>::min();

  _category = 0;

  _nu_vtx_x = std::numeric_limits<float>::min();
  _nu_vtx_y = std::numeric_limits<float>::min();
  _nu_vtx_z = std::numeric_limits<float>::min();

  _nu_sce_x = std::numeric_limits<float>::min();
  _nu_sce_y = std::numeric_limits<float>::min();
  _nu_sce_z = std::numeric_limits<float>::min();

  _isVtxInActive = false;
  _isVtxInFiducial = false;

  _nslice = 0;
  _crtveto = 0;
  _crthitpe = 0;

  _nmuon = 0;
  _muon_e = 0;
  _muon_p = 0;
  _muon_c = 0;

  _nelec = 0;
  _elec_e = 0;
  _elec_p = 0;
  _elec_c = 0;

  _npi0 = 0;
  _pi0_e = 0;
  _pi0_p = 0;
  _pi0_c = 0;

  _npion = 0;
  _pion_e = 0;
  _pion_p = 0;
  _pion_c = 0;

  _nproton = 0;
  _proton_e = 0;
  _proton_p = 0;
  _proton_c = 0;

  _endmuonprocess = "";
  _endmuonmichel = 0;

  // _backtracked_idx.clear();
  // _backtracked_tid.clear();
  _backtracked_e.clear();
  _backtracked_pdg.clear();
  _backtracked_purity.clear();
  _backtracked_completeness.clear();

  evnhits = std::numeric_limits<int>::min();
  slpdg = std::numeric_limits<int>::min();
  slnhits = std::numeric_limits<int>::min();
  pfpdg.clear();
  pfnhits.clear();
  pfnplanehits.clear();

  _hits_u = 0;
  _hits_v = 0;
  _hits_y = 0;

  _mc_E.clear();
  _mc_pdg.clear();

  _mc_px.clear();
  _mc_py.clear();
  _mc_pz.clear();

  _mc_vx.clear();
  _mc_vy.clear();
  _mc_vz.clear();

  _mc_endx.clear();
  _mc_endy.clear();
  _mc_endz.clear();
}

void DefaultAnalysis::SaveTruth(art::Event const &e)
{

  // load MCTruth
  auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);

  // load MCTruth [from geant]
  auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);

  auto mct = mct_h->at(0);
  auto neutrino = mct.GetNeutrino();
  auto nu = neutrino.Nu();

  _ccnc = neutrino.CCNC();
  _nu_pdg = nu.PdgCode();
  _nu_e = nu.Trajectory().E(0);
  _vtx_x = nu.EndX();
  _vtx_y = nu.EndY();
  _vtx_z = nu.EndZ();
  _vtx_t = nu.T();

  art::ServiceHandle<geo::Geometry> geo;

  if (_vtx_x < 2. * geo->DetHalfWidth() && _vtx_x > 0. &&
      _vtx_y < geo->DetHalfHeight() && _vtx_y > -geo->DetHalfHeight() &&
      _vtx_z < geo->DetLength() && _vtx_z > 0.)
  {
    _isVtxInActive = true;
  }
  else
    _isVtxInActive = false;

  double vtx[3] = {_vtx_x, _vtx_y, _vtx_z};
  _isVtxInFiducial = isFiducial(vtx);

  _nelec = 0;
  _nmuon = 0;
  _npi0 = 0;
  _nproton = 0;
  _npion = 0;

  // save muon trackID [ needed to find Michel ]
  float muonMomentum = 0;

  size_t npart = mct.NParticles();
  for (size_t i = 0; i < npart; i++)
  {

    auto const &part = mct.GetParticle(i);
    if (part.StatusCode() != 1)
    {
      continue;
    }

    // if muon
    if ((std::abs(part.PdgCode()) == 13) and (part.StatusCode() == 1))
    {
      muonMomentum = part.Momentum(0).E();
      _nmuon += 1;
      _muon_e = part.Momentum(0).E();
    } // if muon
    // if electron
    if ((std::abs(part.PdgCode()) == 11) and (part.StatusCode() == 1))
    {
      _nelec += 1;
      _elec_e = part.Momentum(0).E();
    } // if electron
    // if pi0
    if ((part.PdgCode() == 111) and (part.StatusCode() == 1))
    {
      _npi0 += 1;
      _pi0_e = part.Momentum(0).E();
    } // if pi0
    // if proton
    if ((part.PdgCode() == 2212) and (part.StatusCode() == 1))
    {
      // if highest energy, update energy
      if (part.Momentum(0).E() > _proton_e)
        _proton_e = part.Momentum(0).E();
      if (part.Momentum(0).E() > fProtonThreshold)
        _nproton += 1;
    } // if proton
    // if pion
    if ((std::abs(part.PdgCode()) == 211) and (part.StatusCode() == 1))
    {
      if (part.Momentum(0).E() > _pion_e)
        _pion_e = part.Momentum(0).E();
      _npion += 1;
    } // if pion
  }   // for all MCParticles

  for (size_t p = 0; p < mcp_h->size(); p++)
  {
    auto mcp = mcp_h->at(p);
    if (!(mcp.Process() == "primary" &&
          mcp.StatusCode() == 1)) {
      continue;
    }

    _mc_E.push_back(mcp.E());

    _mc_pdg.push_back(mcp.PdgCode());

    _mc_px.push_back(mcp.Px());
    _mc_py.push_back(mcp.Py());
    _mc_pz.push_back(mcp.Pz());

    _mc_vx.push_back(mcp.Vx());
    _mc_vy.push_back(mcp.Vy());
    _mc_vz.push_back(mcp.Vz());

    _mc_endx.push_back(mcp.EndX());
    _mc_endy.push_back(mcp.EndY());
    _mc_endz.push_back(mcp.EndZ());
  }

  searchingfornues::ApplyDetectorOffsets(_vtx_t, _vtx_x, _vtx_y, _vtx_z, _xtimeoffset, _xsceoffset, _ysceoffset, _zsceoffset);

  // find if mu -> michel
  _endmuonprocess = "";
  _endmuonmichel = 0;
  bool containedMu = false;
  float muendpointX = 0;

  if (muonMomentum > 0)
  {

    int muonTrackId = -1;

    // loop through all MCParticles, find the muon
    for (size_t p = 0; p < mcp_h->size(); p++)
    {
      auto mcp = mcp_h->at(p);

      if (mcp.StatusCode() != 1)
      {
        continue;
      }

      if ((mcp.Momentum(0).E() - muonMomentum) < 0.0001)
      {
        muonTrackId = mcp.TrackId();
        // stops in the detector?
        if ((mcp.EndPosition().X() > 0) && (mcp.EndPosition().X() < 2 * geo->DetHalfWidth()) &&
            (mcp.EndPosition().Y() > -geo->DetHalfHeight()) && (mcp.EndPosition().Y() < geo->DetHalfHeight()) &&
            (mcp.EndPosition().Z() > 0) && (mcp.EndPosition().Z() < geo->DetLength()))
        {
          _endmuonprocess = mcp.EndProcess();
          containedMu = true;
          muendpointX = mcp.EndPosition().X();
        } // contained muons
        break;
      } // if we find the muon
    }   // for all MCParticles

    // loop again through MCParticles searching for the Michel, if the muon is contained
    if (containedMu)
    {

      // loop through all MCParticles, find the michel
      for (size_t p = 0; p < mcp_h->size(); p++)
      {

        auto mcp = mcp_h->at(p);

        if (fabs(mcp.PdgCode()) != 11)
          continue;

        // if parent of this particle is the muon
        if (mcp.Mother() == muonTrackId)
        {

          // do they start at the same position?
          if ((mcp.Vx() - muendpointX) < 0.0001)
          {
            _endmuonprocess = mcp.Process();
            _endmuonmichel = mcp.Momentum(0).E();
            //std::cout << "Found MCParticle coincident with Michel of PDG "
          } // if start point coincides with muon end point
        }   // if parent is muon
      }     // for all particles, second loop
    }       // if muon is contained
  }         // if there is a muon, we are searching for a possible Michel decay

  return;
}

void DefaultAnalysis::ApplySCECorrection(const float &vtx_x, const float &vtx_y, const float &vtx_z,
                                         float &sce_x, float &sce_y, float &sce_z)
{

  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

  if (SCE->EnableCalSpatialSCE() == true)
  {

    auto offset = SCE->GetPosOffsets(geo::Point_t(vtx_x, vtx_y, vtx_z));
    sce_x = offset.X();
    sce_y = offset.Y();
    sce_z = offset.Z();

  } // if spatial offset calibrations are enabled

  return;
}

DEFINE_ART_CLASS_TOOL(DefaultAnalysis)
} // namespace analysis

#endif
