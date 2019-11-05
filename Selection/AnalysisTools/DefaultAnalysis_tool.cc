#ifndef ANALYSIS_DEFAULTANALYSIS_CXX
#define ANALYSIS_DEFAULTANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/Containment.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"

#include "canvas/Persistency/Common/TriggerResults.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

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

  TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
  TParticlePDG *neutron = TDatabasePDG::Instance()->GetParticle(2112);
  TParticlePDG *electron = TDatabasePDG::Instance()->GetParticle(11);
  TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);
  TParticlePDG *pion = TDatabasePDG::Instance()->GetParticle(211);
  TParticlePDG *pi0 = TDatabasePDG::Instance()->GetParticle(111);
  TParticlePDG *electron_neutrino = TDatabasePDG::Instance()->GetParticle(12);
  TParticlePDG *muon_neutrino = TDatabasePDG::Instance()->GetParticle(14);

  art::InputTag fCRTVetoproducer; // producer for CRT veto ass tag [anab::T0 <-> recob::OpFlash]
  art::InputTag fCLSproducer;     // cluster associated to PFP
  art::InputTag fMCTproducer;     // MCTruth from neutrino generator
  art::InputTag fMCPproducer;     // MCParticle from Geant4 stage
  art::InputTag fBacktrackTag;
  art::InputTag fHproducer;
  art::InputTag fMCRproducer;
  art::InputTag fSLCproducer; // slice associated to PFP
  float fTrkShrScore;         /**< Threshold on the Pandora track score (default 0.5) */

  float fFidvolXstart;
  float fFidvolXend;
  float fFidvolYstart;
  float fFidvolYend;
  float fFidvolZstart;
  float fFidvolZend;

  const int k_nu_e_other = 1;
  const int k_nu_e_cc0pi0p = 10;
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
  float fPionThreshold;
  //float fPi0Threshold;
  float fElectronThreshold;
  float fMuonThreshold;

  int _category; // event category

  float _true_nu_vtx_t, _true_nu_vtx_x, _true_nu_vtx_y, _true_nu_vtx_z;
  float _true_nu_vtx_sce_x, _true_nu_vtx_sce_y, _true_nu_vtx_sce_z;
  float _reco_nu_vtx_x, _reco_nu_vtx_y, _reco_nu_vtx_z;
  float _reco_nu_vtx_sce_x, _reco_nu_vtx_sce_y, _reco_nu_vtx_sce_z;

  int _run, _sub, _evt; // event info
  // has the swtrigger fired?
  int _swtrig;
  // neutrino information
  float _nu_e;  /**< neutrino energy [GeV] */
  float _nu_pt; /**< transverse momentum of interaction [GeV/c] */
  float _theta; /**< angle between incoming and outgoing leptons, in radians */

  int _nu_pdg;                  /**< neutrino PDG code */
  int _ccnc;                    /**< CC or NC tag from GENIE */
  int _interaction;             /**< Interaction code from GENIE */
  bool _isVtxInFiducial;        /**< true if neutrino in fiducial volume */

  // final state particle information
  int _nmuon;                         /**< is there a final-state muon from the neutrino? [1=yes 0=no] */
  float _muon_e, _muon_p, _muon_c;    /**< energy, purity, completeness. */
  int _nelec;                         /**< is there a final-state electron from the neutrino? [1=yes 0=no] */
  float _elec_e, _elec_p, _elec_c;    /**< energy, purity, completeness. */
  float _elec_vx, _elec_vy, _elec_vz; /**< electron vertex. */
  int _npi0;                          /**< how many pi0s are there? */
  //int _pi0;                              /**< is there a final-state pi0 from the neutrino? [1=yes 0=no] */
  float _pi0_e, _pi0_p, _pi0_c; /**< energy, purity, completeness. */
  int _nneutron;                /**< how many neutrons are there? */
  int _nproton;                 /**< how many protons are there? */
  //int _proton;                           /**< is there a final-state proton from the neutrino? [1=yes 0=no] */
  float _proton_e, _proton_p, _proton_c; /**< energy, purity, completeness. */
  int _npion;                            /**< how many pions are there? */
  //int _pion;                             /**< is there a final-state charged pion from the neutrino? [1=yes 0=no] */
  float _pion_e, _pion_p, _pion_c; /**< energy, purity, completeness. */

  std::string _endmuonprocess; /**< End muon process name */
  float _endmuonmichel;        /**< End muon Michel electron energy */

  int _nslice;     /**< number of slice in the event */
  int _crtveto;    /**< is the event vetoed by the CRT Veto? */
  float _crthitpe; /**< pe associated to CRT hit */

  std::vector<int> _pfp_slice_idx; /**< index of PFP is vector of PFPs in nu slice */

  // reco PFParticle backtracking. One entry for PFParticle in the slice
  // std::vector<int>   _backtracked_idx;    // index of PFP [key]
  // std::vector<int>   _backtracked_tid;    // TrackID of backtracked MCParticle
  std::vector<int> _backtracked_pdg;              // PDG code of backtracked particle
  std::vector<float> _backtracked_e;              // energy of backtracked particle
  std::vector<float> _backtracked_purity;         // purity of backtracking
  std::vector<float> _backtracked_completeness;   // completeness of backtracking
  std::vector<float> _backtracked_overlay_purity; // purity of overlay

  std::vector<float> _backtracked_px;
  std::vector<float> _backtracked_py;
  std::vector<float> _backtracked_pz;

  std::vector<float> _backtracked_start_x;
  std::vector<float> _backtracked_start_y;
  std::vector<float> _backtracked_start_z;
  std::vector<float> _backtracked_start_t;
  std::vector<float> _backtracked_start_U;
  std::vector<float> _backtracked_start_V;
  std::vector<float> _backtracked_start_Y;
  std::vector<float> _backtracked_sce_start_x;
  std::vector<float> _backtracked_sce_start_y;
  std::vector<float> _backtracked_sce_start_z;
  std::vector<float> _backtracked_sce_start_U;
  std::vector<float> _backtracked_sce_start_V;
  std::vector<float> _backtracked_sce_start_Y;

  float _lep_e;                                              // lepton energy (if one exists) [GeV]
  int _pass;                                                 // does the slice pass the selection

  int evnhits;                     // number of hits in event
  int slpdg;                       // PDG code of primary pfp in slice
  int slnhits;                     // number of hits in slice
  float _topo_score;               /**< topological score of the slice */
  std::vector<int> pfpdg;          // PDG code of pfp in slice
  std::vector<int> pfnhits;        // number of hits in pfp
  std::vector<int> pfnplanehits_U; // number of hits in pfp plane U
  std::vector<int> pfnplanehits_V; // number of hits in pfp plane V
  std::vector<int> pfnplanehits_Y; // number of hits in pfp plane Y
  float slclustfrac;               //fraction of clustered hits in the slice

  std::vector<uint> _generation;    // generation, 1 is primary
  std::vector<uint> _shr_daughters; // number of shower daughters
  std::vector<uint> _trk_daughters; // number of track daughters

  unsigned int _n_pfps;
  std::vector<float> _trk_score_v;
  unsigned int _n_tracks;
  unsigned int _n_showers;

  unsigned int _hits_u;
  unsigned int _hits_v;
  unsigned int _hits_y;

  std::vector<int> _mc_pdg;
  std::vector<float> _mc_E;

  std::vector<float> _mc_px;
  std::vector<float> _mc_py;
  std::vector<float> _mc_pz;

  std::vector<float> _mc_vx;
  std::vector<float> _mc_vy;
  std::vector<float> _mc_vz;

  std::vector<float> _mc_endx;
  std::vector<float> _mc_endy;
  std::vector<float> _mc_endz;

  std::vector<float> _mc_completeness;
  std::vector<float> _mc_purity;

  float _true_pt;
  float _true_pt_visible;
  float _true_p;
  float _true_p_visible;
  float _true_e_visible;
  float _leeweight;
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
  fTrkShrScore = p.get<float>("TrkShrScore", 0.5);
  // kinematic thresholds for defining signal
  fProtonThreshold = p.get<float>("ProtonThreshold", 0.04);
  fMuonThreshold = p.get<float>("MuonThreshold", 0.02);
  fPionThreshold = p.get<float>("PionThreshold", 0.04);
  fElectronThreshold = p.get<float>("ElectronThreshold", 0.03);

  fFidvolXstart = p.get<double>("fidvolXstart", 10);
  fFidvolXend = p.get<double>("fidvolXend", 10);

  fFidvolYstart = p.get<double>("fidvolYstart", 10);
  fFidvolYend = p.get<double>("fidvolYend", 10);

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

  std::cout << "[DefaultAnalysis::analyzeEvent] Run: " << _run << ", SubRun: " << _sub << ", Event: " << _evt << std::endl;

  // storing trigger result output for software trigger
  art::InputTag swtrig_tag("TriggerResults", "", "DataOverlayOptical");
  art::Handle<art::TriggerResults> swtrig_handle;
  e.getByLabel(swtrig_tag, swtrig_handle);
  if (swtrig_handle.isValid())
  {
    if (swtrig_handle->accept() == true)
      _swtrig = 1;
    else
      _swtrig = 0;
  } // if software trigger run by this producer

  if (!fData)
  {
    art::InputTag eventweight_tag("eventweightLEE");
    art::Handle<std::vector<evwgh::MCEventWeight>> eventweights_handle;
    e.getByLabel(eventweight_tag, eventweights_handle);
    if (eventweights_handle.isValid())
    {
      std::vector<art::Ptr<evwgh::MCEventWeight>> eventweights;
      art::fill_ptr_vector(eventweights, eventweights_handle);
      std::map<std::string, std::vector<double>> evtwgt_map = eventweights.at(0)->fWeight;
      for (std::map<std::string, std::vector<double>>::iterator it = evtwgt_map.begin(); it != evtwgt_map.end(); ++it)
      {
        if (it->second.size() == 1)
        {
          _leeweight = it->second[0];
        }
      }
    }
    else
    {
      std::cout << "[DefaultAnalysis::analyzeEvent] LEE MCEventWeight not present" << std::endl;
    }

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

void DefaultAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
  ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
                                                                                        proxy::withAssociated<recob::Hit>(fCLSproducer));
  // somehow proxies don't work for the slice-hit association, so go back to old assns
  art::ValidHandle<std::vector<recob::Slice>> inputSlice = e.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
  auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSlice, e, fSLCproducer));

  // Build larpandora info:
  lar_pandora::LArPandoraHelper larpandora;
  lar_pandora::PFParticleVector pfparticles;
  lar_pandora::PFParticleMap particleMap;
  larpandora.CollectPFParticles(e, "pandora", pfparticles);
  larpandora.BuildPFParticleMap(pfparticles, particleMap);

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
  _n_pfps = 0;
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

      auto metadata_pxy_v = pfp.get<larpandoraobj::PFParticleMetadata>();

      if (metadata_pxy_v.size() != 0)
      {
        for (unsigned int j = 0; j < metadata_pxy_v.size(); ++j)
        {
          const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(metadata_pxy_v.at(j));
          auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
          if (!pfParticlePropertiesMap.empty())
          {
            _topo_score = pfParticlePropertiesMap.at("NuScore");
          } // if PFP metadata exists!
        }
      }

      // grab vertex
      double xyz[3] = {};

      auto vtx = pfp.get<recob::Vertex>();
      if (vtx.size() == 1)
      {
        // save vertex to array
        vtx.at(0)->XYZ(xyz);
        auto nuvtx = TVector3(xyz[0], xyz[1], xyz[2]);

        _reco_nu_vtx_x = nuvtx.X();
        _reco_nu_vtx_y = nuvtx.Y();
        _reco_nu_vtx_z = nuvtx.Z();

        float _reco_nu_vtx_sce[3];
        searchingfornues::ApplySCECorrectionXYZ(_reco_nu_vtx_x, _reco_nu_vtx_y, _reco_nu_vtx_z, _reco_nu_vtx_sce);
        _reco_nu_vtx_sce_x = _reco_nu_vtx_sce[0];
        _reco_nu_vtx_sce_y = _reco_nu_vtx_sce[1];
        _reco_nu_vtx_sce_z = _reco_nu_vtx_sce[2];
      }
      else
      {
        std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
      }
      continue;
    } // if neutrino PFParticle

    _n_pfps++;
    _pfp_slice_idx.push_back(pfpidx++);
    pfpdg.push_back(pfp->PdgCode());

    // Hieracrchy information:
    _generation.push_back(larpandora.GetGeneration(particleMap, particleMap.at(pfp->Self())));
    uint this_num_trk_d = 0;
    uint this_num_shr_d = 0;
    for (size_t daughter : pfp->Daughters())
    {
      if (larpandora.IsTrack(particleMap.at(daughter)))
      {
        this_num_trk_d++; // Track daughter
      }
      else
      {
        this_num_shr_d++; // Shower daughter
      }
    }
    _shr_daughters.push_back(this_num_shr_d);
    _trk_daughters.push_back(this_num_trk_d);

    // store track score
    float trkscore = searchingfornues::GetTrackShowerScore(pfp);
    if ((trkscore >= 0) && (trkscore >= fTrkShrScore))
    {
      _n_tracks++;
    }
    else if ((trkscore >= 0) && (trkscore < fTrkShrScore))
    {
      _n_showers++;
    }
    _trk_score_v.push_back(trkscore);

    // get hits associated to this PFParticle through the clusters
    std::vector<art::Ptr<recob::Hit>> hit_v;
    auto clus_pxy_v = pfp.get<recob::Cluster>();

    pfnplanehits_U.push_back(0);
    pfnplanehits_V.push_back(0);
    pfnplanehits_Y.push_back(0);

    for (auto ass_clus : clus_pxy_v)
    {
      // get cluster proxy
      const auto &clus = clus_proxy[ass_clus.key()];
      auto clus_hit_v = clus.get<recob::Hit>();
      auto nhits = clus_hit_v.size();

      if (clus->Plane().Plane == 0)
      {
        _hits_u += nhits;
        pfnplanehits_U.back() += nhits;
      }
      else if (clus->Plane().Plane == 1)
      {
        _hits_v += nhits;
        pfnplanehits_V.back() += nhits;
      }
      else if (clus->Plane().Plane == 2)
      {
        _hits_y += nhits;
        pfnplanehits_Y.back() += nhits;
      }
      for (const auto &hit : clus_hit_v)
      {
        hit_v.push_back(hit);
      }
    } // for all clusters associated to PFP
    pfnhits.push_back(hit_v.size());

    // backtrack PFParticles in the slice
    if (!fData)
    {
      if (clus_pxy_v.size() != 0)
      {
        float purity = 0., completeness = 0., overlay_purity = 0.;
        int ibt = searchingfornues::getAssocBtPart(hit_v, assocMCPart, btparts_v, purity, completeness, overlay_purity);
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
          _backtracked_overlay_purity.push_back(overlay_purity);

          _backtracked_px.push_back(mcp.px);
          _backtracked_py.push_back(mcp.py);
          _backtracked_pz.push_back(mcp.pz);
          _backtracked_start_x.push_back(mcp.start_x);
          _backtracked_start_y.push_back(mcp.start_y);
          _backtracked_start_z.push_back(mcp.start_z);
          _backtracked_start_t.push_back(mcp.start_t);

          _backtracked_start_U.push_back(searchingfornues::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 0));
          _backtracked_start_V.push_back(searchingfornues::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 1));
          _backtracked_start_Y.push_back(searchingfornues::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 2));

          for (size_t i = 0; i < _mc_E.size(); i++)
          {
            if (_mc_E[i] == mcp.e)
            {
              _mc_completeness[i] = completeness;
              _mc_purity[i] = purity;
            }
          }

          float reco_st[3] = {mcp.start_x, mcp.start_y, mcp.start_z};

          if (PDG == 11 || PDG == 22)
          {
            reco_st[0] += searchingfornues::x_offset(mcp.start_t);
          }
          else
          {
            searchingfornues::True2RecoMappingXYZ(mcp.start_t, mcp.start_x, mcp.start_y, mcp.start_z, reco_st);
          }

          _backtracked_sce_start_x.push_back(reco_st[0]);
          _backtracked_sce_start_y.push_back(reco_st[1]);
          _backtracked_sce_start_z.push_back(reco_st[2]);

          _backtracked_sce_start_U.push_back(searchingfornues::YZtoPlanecoordinate(reco_st[1], reco_st[2], 0));
          _backtracked_sce_start_V.push_back(searchingfornues::YZtoPlanecoordinate(reco_st[1], reco_st[2], 1));
          _backtracked_sce_start_Y.push_back(searchingfornues::YZtoPlanecoordinate(reco_st[1], reco_st[2], 2));

          // if this is an interesting particle, save it to the TTree
          if (fabs(PDG) == muon->PdgCode())
          {
            if (fabs(mcp.e - _muon_e) < 0.01)
            {
              _muon_p = purity;
              _muon_c = completeness;
            }
          }
          if (fabs(PDG) == electron->PdgCode())
          {
            if (fabs(mcp.e - _elec_e) < 0.01)
            {
              _elec_p = purity;
              _elec_c = completeness;
            }
          }
          if (fabs(PDG) == proton->PdgCode())
          {
            if (fabs(mcp.e - _proton_e) < 0.0001)
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
          _backtracked_e.push_back(std::numeric_limits<float>::lowest());
          _backtracked_pdg.push_back(0);
          _backtracked_purity.push_back(std::numeric_limits<float>::lowest());
          _backtracked_completeness.push_back(std::numeric_limits<float>::lowest());
          _backtracked_overlay_purity.push_back(std::numeric_limits<float>::lowest());

          _backtracked_px.push_back(std::numeric_limits<float>::lowest());
          _backtracked_py.push_back(std::numeric_limits<float>::lowest());
          _backtracked_pz.push_back(std::numeric_limits<float>::lowest());
          _backtracked_start_x.push_back(std::numeric_limits<float>::lowest());
          _backtracked_start_y.push_back(std::numeric_limits<float>::lowest());
          _backtracked_start_z.push_back(std::numeric_limits<float>::lowest());
          _backtracked_start_t.push_back(std::numeric_limits<float>::lowest());

          _backtracked_start_U.push_back(std::numeric_limits<float>::lowest());
          _backtracked_start_V.push_back(std::numeric_limits<float>::lowest());
          _backtracked_start_Y.push_back(std::numeric_limits<float>::lowest());

          _backtracked_sce_start_x.push_back(std::numeric_limits<float>::lowest());
          _backtracked_sce_start_y.push_back(std::numeric_limits<float>::lowest());
          _backtracked_sce_start_z.push_back(std::numeric_limits<float>::lowest());

          _backtracked_sce_start_U.push_back(std::numeric_limits<float>::lowest());
          _backtracked_sce_start_V.push_back(std::numeric_limits<float>::lowest());
          _backtracked_sce_start_Y.push_back(std::numeric_limits<float>::lowest());
        }
      } // if there are associated clusters
    }   // if MC
  }
  if (slnhits > 0)
  {
    slclustfrac = 0.;
    for (auto n : pfnhits)
      slclustfrac += n;
    slclustfrac /= float(slnhits);
  }
  _nslice += 1;

  if (!fData)
  {
    bool there_is_true_proton = _nproton > 0;
    bool there_is_true_pi = _npion > 0;
    bool there_is_true_mu = _nmuon > 0;
    bool there_is_true_pi0 = _npi0 > 0;
    bool there_is_true_electron = _nelec > 0;

    if (!_isVtxInFiducial)
    {
      _category = k_outfv;
    }
    else if (abs(_nu_pdg) == electron_neutrino->PdgCode())
    {
      if (there_is_true_electron)
      {
        if (!there_is_true_pi && there_is_true_proton && !there_is_true_pi0)
        {
          _category = k_nu_e_cc0pinp;
        }
        else if (!there_is_true_pi && !there_is_true_proton && !there_is_true_pi0)
        {
          _category = k_nu_e_cc0pi0p;
        }
        else
        {
          _category = k_nu_e_other;
        }
      }
      else
      {
        if (!there_is_true_pi0)
        {
          _category = k_nc;
        }
        else
        {
          _category = k_nc_pi0;
        }
      }
    }
    else if (abs(_nu_pdg) == muon_neutrino->PdgCode())
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
        if (!there_is_true_pi0)
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
  _tree->Branch("leeweight", &_leeweight, "leeweight/F");

  _tree->Branch("true_pt", &_true_pt, "true_pt/F");
  _tree->Branch("true_pt_visible", &_true_pt_visible, "true_pt_visible/F");
  _tree->Branch("true_p", &_true_p, "true_p/F");
  _tree->Branch("true_p_visible", &_true_p_visible, "true_p_visible/F");

  _tree->Branch("true_e_visible", &_true_e_visible, "true_e_visible/F");

  // neutrino information
  _tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
  _tree->Branch("ccnc", &_ccnc, "ccnc/I");
  _tree->Branch("interaction", &_interaction, "interaction/I");
  _tree->Branch("nu_e", &_nu_e, "nu_e/F");
  _tree->Branch("nu_pt", &_nu_pt, "nu_pt/F");
  _tree->Branch("theta", &_theta, "theta/F");
  _tree->Branch("isVtxInFiducial", &_isVtxInFiducial, "isVtxInFiducial/O");

  _tree->Branch("true_nu_vtx_t", &_true_nu_vtx_t, "true_nu_vtx_t/F");
  _tree->Branch("true_nu_vtx_x", &_true_nu_vtx_x, "true_nu_vtx_x/F");
  _tree->Branch("true_nu_vtx_y", &_true_nu_vtx_y, "true_nu_vtx_y/F");
  _tree->Branch("true_nu_vtx_z", &_true_nu_vtx_z, "true_nu_vtx_z/F");
  _tree->Branch("true_nu_vtx_sce_x", &_true_nu_vtx_sce_x, "true_nu_vtx_sce_x/F");
  _tree->Branch("true_nu_vtx_sce_y", &_true_nu_vtx_sce_y, "true_nu_vtx_sce_y/F");
  _tree->Branch("true_nu_vtx_sce_z", &_true_nu_vtx_sce_z, "true_nu_vtx_sce_z/F");
  _tree->Branch("reco_nu_vtx_x", &_reco_nu_vtx_x, "reco_nu_vtx_x/F");
  _tree->Branch("reco_nu_vtx_y", &_reco_nu_vtx_y, "reco_nu_vtx_y/F");
  _tree->Branch("reco_nu_vtx_z", &_reco_nu_vtx_z, "reco_nu_vtx_z/F");
  _tree->Branch("reco_nu_vtx_sce_x", &_reco_nu_vtx_sce_x, "reco_nu_vtx_sce_x/F");
  _tree->Branch("reco_nu_vtx_sce_y", &_reco_nu_vtx_sce_y, "reco_nu_vtx_sce_y/F");
  _tree->Branch("reco_nu_vtx_sce_z", &_reco_nu_vtx_sce_z, "reco_nu_vtx_sce_z/F");

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
  _tree->Branch("elec_vx", &_elec_vx, "elec_vx/F");
  _tree->Branch("elec_vy", &_elec_vy, "elec_vy/F");
  _tree->Branch("elec_vz", &_elec_vz, "elec_vz/F");
  // pi0
  _tree->Branch("npi0", &_npi0, "npi0/I");
  _tree->Branch("pi0_e", &_pi0_e, "pi0_e/F");
  _tree->Branch("pi0_c", &_pi0_c, "pi0_c/F");
  _tree->Branch("pi0_p", &_pi0_p, "pi0_p/F");

  _tree->Branch("nneutron", &_nneutron, "nneutron/I");

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
  _tree->Branch("backtracked_overlay_purity", "std::vector<float>", &_backtracked_overlay_purity);

  _tree->Branch("backtracked_px", "std::vector<float>", &_backtracked_px);
  _tree->Branch("backtracked_py", "std::vector<float>", &_backtracked_py);
  _tree->Branch("backtracked_pz", "std::vector<float>", &_backtracked_pz);

  _tree->Branch("backtracked_start_x", "std::vector<float>", &_backtracked_start_x);
  _tree->Branch("backtracked_start_y", "std::vector<float>", &_backtracked_start_y);
  _tree->Branch("backtracked_start_z", "std::vector<float>", &_backtracked_start_z);
  _tree->Branch("backtracked_start_t", "std::vector<float>", &_backtracked_start_t);
  _tree->Branch("backtracked_start_U", "std::vector<float>", &_backtracked_start_U);
  _tree->Branch("backtracked_start_V", "std::vector<float>", &_backtracked_start_V);
  _tree->Branch("backtracked_start_Y", "std::vector<float>", &_backtracked_start_Y);
  _tree->Branch("backtracked_sce_start_x", "std::vector<float>", &_backtracked_sce_start_x);
  _tree->Branch("backtracked_sce_start_y", "std::vector<float>", &_backtracked_sce_start_y);
  _tree->Branch("backtracked_sce_start_z", "std::vector<float>", &_backtracked_sce_start_z);
  _tree->Branch("backtracked_sce_start_U", "std::vector<float>", &_backtracked_sce_start_U);
  _tree->Branch("backtracked_sce_start_V", "std::vector<float>", &_backtracked_sce_start_V);
  _tree->Branch("backtracked_sce_start_Y", "std::vector<float>", &_backtracked_sce_start_Y);

  _tree->Branch("lep_e", &_lep_e, "lep_e/F");
  _tree->Branch("pass", &_pass, "pass/I");
  _tree->Branch("run", &_run, "run/I");
  _tree->Branch("sub", &_sub, "sub/I");
  _tree->Branch("evt", &_evt, "evt/I");

  _tree->Branch("swtrig", &_swtrig, "swtrig/I");

  _tree->Branch("evnhits", &evnhits, "evnhits/I");
  _tree->Branch("slpdg", &slpdg, "slpdg/I");
  _tree->Branch("slnhits", &slnhits, "slnhits/I");
  _tree->Branch("n_pfps", &_n_pfps, "n_pfps/I");
  _tree->Branch("n_tracks", &_n_tracks, "n_tracks/I");
  _tree->Branch("n_showers", &_n_showers, "n_showers/I");

  _tree->Branch("pfp_generation_v", "std::vector< uint >", &_generation);
  _tree->Branch("pfp_trk_daughters_v", "std::vector< uint >", &_trk_daughters);
  _tree->Branch("pfp_shr_daughters_v", "std::vector< uint >", &_shr_daughters);

  _tree->Branch("trk_score_v", "std::vector< float >", &_trk_score_v);

  _tree->Branch("pfpdg", "std::vector<int>", &pfpdg);
  _tree->Branch("pfnhits", "std::vector<int>", &pfnhits);
  _tree->Branch("pfnplanehits_U", "std::vector<int>", &pfnplanehits_U);
  _tree->Branch("pfnplanehits_V", "std::vector<int>", &pfnplanehits_V);
  _tree->Branch("pfnplanehits_Y", "std::vector<int>", &pfnplanehits_Y);
  _tree->Branch("hits_u", &_hits_u, "hits_u/i");
  _tree->Branch("hits_v", &_hits_v, "hits_v/i");
  _tree->Branch("hits_y", &_hits_y, "hits_y/i");
  _tree->Branch("topological_score", &_topo_score, "topological_score/F");
  _tree->Branch("slclustfrac", &slclustfrac, "slclustfrac/F");

  _tree->Branch("mc_pdg", "std::vector< int >", &_mc_pdg);
  _tree->Branch("mc_E", "std::vector< float >", &_mc_E);

  _tree->Branch("mc_vx", "std::vector< float >", &_mc_vx);
  _tree->Branch("mc_vy", "std::vector< float >", &_mc_vy);
  _tree->Branch("mc_vz", "std::vector< float >", &_mc_vz);

  _tree->Branch("mc_endx", "std::vector< float >", &_mc_endx);
  _tree->Branch("mc_endy", "std::vector< float >", &_mc_endy);
  _tree->Branch("mc_endz", "std::vector< float >", &_mc_endz);

  _tree->Branch("mc_px", "std::vector< float >", &_mc_px);
  _tree->Branch("mc_py", "std::vector< float >", &_mc_py);
  _tree->Branch("mc_pz", "std::vector< float >", &_mc_pz);

  _tree->Branch("mc_completeness", "std::vector< float >", &_mc_completeness);
  _tree->Branch("mc_purity", "std::vector< float >", &_mc_purity);

  _tree->Branch("endmuonprocess", &_endmuonprocess);

  _tree->Branch("endmuonmichel", &_endmuonmichel, "endmuonmichel/F");
}

void DefaultAnalysis::resetTTree(TTree *_tree)
{
  _leeweight = 0;
  _run = std::numeric_limits<int>::lowest();
  _sub = std::numeric_limits<int>::lowest();
  _evt = std::numeric_limits<int>::lowest();
  _nu_e = std::numeric_limits<float>::lowest();
  _theta = std::numeric_limits<float>::lowest();
  _nu_pt = std::numeric_limits<float>::lowest();

  _swtrig = std::numeric_limits<int>::lowest();

  _nu_pdg = std::numeric_limits<int>::lowest();
  _ccnc = std::numeric_limits<int>::lowest();
  _interaction = std::numeric_limits<int>::lowest();
  _pass = 0;

  _category = 0;

  _true_nu_vtx_t = std::numeric_limits<float>::lowest();
  _true_nu_vtx_x = std::numeric_limits<float>::lowest();
  _true_nu_vtx_y = std::numeric_limits<float>::lowest();
  _true_nu_vtx_z = std::numeric_limits<float>::lowest();
  _true_nu_vtx_sce_x = std::numeric_limits<float>::lowest();
  _true_nu_vtx_sce_y = std::numeric_limits<float>::lowest();
  _true_nu_vtx_sce_z = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_x = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_y = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_z = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_sce_x = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_sce_y = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_sce_z = std::numeric_limits<float>::lowest();

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
  _elec_vx = std::numeric_limits<float>::lowest();
  _elec_vy = std::numeric_limits<float>::lowest();
  _elec_vz = std::numeric_limits<float>::lowest();
  _npi0 = 0;
  _pi0_e = 0;
  _pi0_p = 0;
  _pi0_c = 0;

  _npion = 0;
  _pion_e = 0;
  _pion_p = 0;
  _pion_c = 0;

  _nneutron = 0;

  _nproton = 0;
  _proton_e = 0;
  _proton_p = 0;
  _proton_c = 0;

  _endmuonprocess = "";
  _endmuonmichel = 0;

  _n_pfps = 0;
  _n_tracks = 0;
  _n_showers = 0;
  _trk_score_v.clear();
  // _backtracked_idx.clear();
  // _backtracked_tid.clear();
  _backtracked_e.clear();
  _backtracked_pdg.clear();
  _backtracked_purity.clear();
  _backtracked_completeness.clear();
  _backtracked_overlay_purity.clear();

  _backtracked_px.clear();
  _backtracked_py.clear();
  _backtracked_pz.clear();

  _backtracked_start_x.clear();
  _backtracked_start_y.clear();
  _backtracked_start_z.clear();
  _backtracked_start_t.clear();
  _backtracked_start_U.clear();
  _backtracked_start_V.clear();
  _backtracked_start_Y.clear();
  _backtracked_sce_start_x.clear();
  _backtracked_sce_start_y.clear();
  _backtracked_sce_start_z.clear();
  _backtracked_sce_start_U.clear();
  _backtracked_sce_start_V.clear();
  _backtracked_sce_start_Y.clear();

  evnhits = std::numeric_limits<int>::lowest();
  slpdg = std::numeric_limits<int>::lowest();
  _topo_score = std::numeric_limits<float>::lowest();
  slnhits = std::numeric_limits<int>::lowest();
  pfpdg.clear();
  pfnhits.clear();
  pfnplanehits_U.clear();
  pfnplanehits_V.clear();
  pfnplanehits_Y.clear();
  _generation.clear();
  _shr_daughters.clear();
  _trk_daughters.clear();
  slclustfrac = std::numeric_limits<float>::lowest();

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

  _mc_completeness.clear();
  _mc_purity.clear();

  _true_pt = 0;
  _true_pt_visible = 0;
  _true_p = 0;
  _true_p_visible = 0;

  _true_e_visible = 0;
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
  _interaction = neutrino.Mode();
  _nu_pdg = nu.PdgCode();
  _nu_e = nu.Trajectory().E(0);

  _true_nu_vtx_t = nu.T();
  _true_nu_vtx_x = nu.Vx();
  _true_nu_vtx_y = nu.Vy();
  _true_nu_vtx_z = nu.Vz();

  float _true_nu_vtx_sce[3];
  searchingfornues::True2RecoMappingXYZ(_true_nu_vtx_t, _true_nu_vtx_x, _true_nu_vtx_y, _true_nu_vtx_z, _true_nu_vtx_sce);

  _true_nu_vtx_sce_x = _true_nu_vtx_sce[0];
  _true_nu_vtx_sce_y = _true_nu_vtx_sce[1];
  _true_nu_vtx_sce_z = _true_nu_vtx_sce[2];

  _theta = neutrino.Theta();
  _nu_pt = neutrino.Pt();

  double vtx[3] = {_true_nu_vtx_x, _true_nu_vtx_y, _true_nu_vtx_z};
  _isVtxInFiducial = searchingfornues::isFiducial(vtx,
						  fFidvolXstart, fFidvolYstart, fFidvolZstart,
						  fFidvolXend, fFidvolYend, fFidvolZend);
  
  _nelec = 0;
  _nmuon = 0;
  _npi0 = 0;
  _nproton = 0;
  _npion = 0;
  _nneutron = 0;
  // save muon trackID [ needed to find Michel ]
  float muonMomentum = 0;

  TLorentzVector total_p;
  TLorentzVector total_p_visible;

  size_t npart = mct.NParticles();
  for (size_t i = 0; i < npart; i++)
  {

    auto const &part = mct.GetParticle(i);
    if (part.StatusCode() != 1)
    {
      continue;
    }

    total_p += part.Momentum(0);
    // if muon
    if ((std::abs(part.PdgCode()) == muon->PdgCode()) and (part.StatusCode() == 1))
    {
      muonMomentum = part.Momentum(0).E();

      if (part.Momentum(0).E() - muon->Mass() > fMuonThreshold)
      {
        _true_e_visible += part.Momentum(0).E() - muon->Mass();
        total_p_visible += part.Momentum(0);
        _nmuon += 1;
      }

      if (part.Momentum(0).E() > _muon_e)
        _muon_e = part.Momentum(0).E();
    } // if muon

    // if electron
    else if ((std::abs(part.PdgCode()) == electron->PdgCode()) and (part.StatusCode() == 1))
    {
      if (part.Momentum(0).E() - electron->Mass() > fElectronThreshold)
      {
        _nelec += 1;
        total_p_visible += part.Momentum(0);
        _true_e_visible += part.Momentum(0).E() - electron->Mass();
      }
      if (part.Momentum(0).E() > _elec_e)
        _elec_e = part.Momentum(0).E();
    } // if electron

    // if pi0
    else if ((part.PdgCode() == pi0->PdgCode()) and (part.StatusCode() == 1))
    {
      _npi0 += 1;
      total_p_visible += part.Momentum(0);
      _true_e_visible += part.Momentum(0).E();
      if (part.Momentum(0).E() > _pi0_e)
        _pi0_e = part.Momentum(0).E();
    } // if pi0

    // if proton
    else if ((part.PdgCode() == proton->PdgCode()) and (part.StatusCode() == 1))
    {
      if (part.Momentum(0).E() - proton->Mass() > fProtonThreshold)
      {
        total_p_visible += part.Momentum(0);
        _true_e_visible += part.Momentum(0).E() - proton->Mass();
        _nproton += 1;
      }
      if (part.Momentum(0).E() > _proton_e)
        _proton_e = part.Momentum(0).E();

    } // if proton

    // if neutron
    else if ((part.PdgCode() == neutron->PdgCode()) and (part.StatusCode() == 1))
    {
      _nneutron += 1;
    }

    // if pion
    else if ((std::abs(part.PdgCode()) == pion->PdgCode()) and (part.StatusCode() == 1))
    {
      if (part.Momentum(0).E() - pion->Mass() > fPionThreshold)
      {
        _npion += 1;
        total_p_visible += part.Momentum(0);
        _true_e_visible += part.Momentum(0).E() - pion->Mass();
      }
      if (part.Momentum(0).E() > _pion_e)
        _pion_e = part.Momentum(0).E();
    } // if pion

    else
    {
      TParticlePDG *particle_pdg = TDatabasePDG::Instance()->GetParticle(part.PdgCode());
      if (particle_pdg != NULL) // PDG codes corresponding to ions e.g. 2000000101 are not in the database
        _true_e_visible += part.Momentum(0).E() - particle_pdg->Mass();
    }

  } // for all MCParticles

  _true_pt = total_p.Perp();
  _true_pt_visible = total_p_visible.Perp();
  _true_p = total_p.Mag();
  _true_p_visible = total_p_visible.Mag();

  for (size_t p = 0; p < mcp_h->size(); p++)
  {
    auto mcp = mcp_h->at(p);
    if (!(mcp.Process() == "primary" &&
          mcp.StatusCode() == 1))
    {
      continue;
    }

    if (mcp.PdgCode() == electron->PdgCode())
    {
      _elec_vx = mcp.Vx();
      _elec_vy = mcp.Vy();
      _elec_vz = mcp.Vz();
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

    _mc_completeness.push_back(std::numeric_limits<float>::lowest());
    _mc_purity.push_back(std::numeric_limits<float>::lowest());
  }

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
        art::ServiceHandle<geo::Geometry> geo;
        geo::TPCGeo const &thisTPC = geo->TPC();
        geo::BoxBoundedGeo theTpcGeo = thisTPC.ActiveBoundingBox();
        if ((mcp.EndPosition().X() > theTpcGeo.MinX()) && (mcp.EndPosition().X() < theTpcGeo.MaxX()) &&
            (mcp.EndPosition().Y() > theTpcGeo.MinY()) && (mcp.EndPosition().Y() < theTpcGeo.MaxY()) &&
            (mcp.EndPosition().Z() > theTpcGeo.MinZ()) && (mcp.EndPosition().Z() < theTpcGeo.MaxZ()))
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

        if (fabs(mcp.PdgCode()) != electron->PdgCode())
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


DEFINE_ART_CLASS_TOOL(DefaultAnalysis)
} // namespace analysis

#endif
