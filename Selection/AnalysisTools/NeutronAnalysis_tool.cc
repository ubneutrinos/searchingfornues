#ifndef ANALYSIS_NEUTRONANALYSIS_CXX
#define ANALYSIS_NEUTRONANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "../CommonDefs/Typedefs.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/Containment.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/ProximityClustering.h"
#include "../CommonDefs/PIDFuncs.h"
#include "../CommonDefs/LLR_PID.h"
#include "../CommonDefs/LLRPID_proton_muon_lookup.h"
#include "../CommonDefs/LLRPID_correction_lookup.h"
#include "../CommonDefs/CalibrationFuncs.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "../SelectionTools/SelectionToolBase.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
////
//// Class:       NeutronAnalysis
//// File:        NeutronAnalysis.cc
////
////              Neutron analysis tools: This tool tracks non-primary
////		  pfps and backtracks them to the mcparticle
////
//// Configuration parameters:
////
//// TBD
////
//// Created by Burke Irwin (irwin175@umn.edu) on 04/18/2023
////
//////////////////////////////////////////////////////////////////////////
//
class NeutronAnalysis : public AnalysisToolBase
{

public:
  /**
 *      *  @brief  Constructor
 *           *
 *                *  @param  pset
 *                     */
  NeutronAnalysis(const fhicl::ParameterSet &pset);

 /**
 *      *  @brief  Destructor
 *           */
  ~NeutronAnalysis(){};

  // provide for initialization
  void configure(fhicl::ParameterSet const &pset);
  
  /**
     * @brief Analysis function
     */
  void analyzeEvent(art::Event const &e, bool fData) override;
  
  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;  

  void setBranches(TTree *_tree) override;

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree) override;
  using ProxyPfpColl_t = selection::ProxyPfpColl_t;

private:
  int _run, _sub, _evt;

  bool  _is_nu_reco;

  std::vector<int>   _nonprim_backtracked_pdg;            // PDG code of nonprim_backtracked particle
  std::vector<float> _nonprim_backtracked_e;              // energy of nonprim_backtracked particle
  std::vector<int>   _nonprim_backtracked_tid;            // track-id of nonprim_backtracked particle
  std::vector<float> _nonprim_backtracked_purity;         // purity of backtracking
  std::vector<float> _nonprim_backtracked_completeness;   // completeness of backtracking
  std::vector<float> _nonprim_backtracked_overlay_purity; // purity of overlay

  std::vector<float> _nonprim_backtracked_px;
  std::vector<float> _nonprim_backtracked_py;
  std::vector<float> _nonprim_backtracked_pz;

  std::vector<float> _nonprim_backtracked_start_x;
  std::vector<float> _nonprim_backtracked_start_y;
  std::vector<float> _nonprim_backtracked_start_z;
  std::vector<float> _nonprim_backtracked_start_t;
  std::vector<float> _nonprim_backtracked_start_U;
  std::vector<float> _nonprim_backtracked_start_V;
  std::vector<float> _nonprim_backtracked_start_Y;
  std::vector<float> _nonprim_backtracked_sce_start_x;
  std::vector<float> _nonprim_backtracked_sce_start_y;
  std::vector<float> _nonprim_backtracked_sce_start_z;
  std::vector<float> _nonprim_backtracked_sce_start_U;
  std::vector<float> _nonprim_backtracked_sce_start_V;
  std::vector<float> _nonprim_backtracked_sce_start_Y;
  std::vector<std::__cxx11::string> _nonprim_backtracked_process;

  std::vector<int> nonprim_pfpdg;          // PDG code of non primary pfp
  std::vector<uint> _nonprim_generation;    // generation, 1 is primary
  std::vector<int> _nonprim_slc_id_v;
  std::vector<int> pfnhits;        // number of hits in pfp

  unsigned int _n_nonprim_pfps;		// Number of non-prime pfps in event

  std::vector<int> _all_mc_pdg;
  std::vector<float> _all_mc_E;

  std::vector<float> _all_mc_px;
  std::vector<float> _all_mc_py;
  std::vector<float> _all_mc_pz;

  std::vector<float> _all_mc_vx;
  std::vector<float> _all_mc_vy;
  std::vector<float> _all_mc_vz;

  std::vector<float> _all_mc_endx;
  std::vector<float> _all_mc_endy;
  std::vector<float> _all_mc_endz;

  std::vector<float> _all_mc_completeness;
  std::vector<float> _all_mc_purity;
  std::vector<std::__cxx11::string>  _all_mc_process;

  std::vector<int> _all_mc_mother;
  std::vector<int> _all_mc_trkid;
  std::vector<float> _all_mc_distance;

  // PID Tools
  searchingfornues::LLRPID llr_pid_calculator;
  searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;
  searchingfornues::CorrectionLookUpParameters correction_parameters;

  std::vector<float>  _nonprim_topo_score_v;
  std::vector<float>  _nonprim_trk_score_v;

  std::vector<float>  _nonprim_trk_calo_energy_u_v;
  std::vector<float>  _nonprim_trk_calo_energy_v_v;
  std::vector<float>  _nonprim_trk_calo_energy_y_v;

  std::vector<float>  _nonprim_trk_llr_pid_u_v;
  std::vector<float>  _nonprim_trk_llr_pid_v_v;
  std::vector<float>  _nonprim_trk_llr_pid_y_v;
  std::vector<float>  _nonprim_trk_llr_pid_v;
  std::vector<float>  _nonprim_trk_llr_pid_score_v;

  std::vector<int>  _nonprim_trk_nhits_u_v;
  std::vector<int>  _nonprim_trk_nhits_v_v;
  std::vector<int>  _nonprim_trk_nhits_y_v;

  std::vector<float> _nonprim_trk_start_x_v;
  std::vector<float> _nonprim_trk_start_y_v;
  std::vector<float> _nonprim_trk_start_z_v;

  std::vector<float> _nonprim_trk_sce_start_x_v;
  std::vector<float> _nonprim_trk_sce_start_y_v;
  std::vector<float> _nonprim_trk_sce_start_z_v;

  std::vector<float> _nonprim_trk_distance_v;
  std::vector<float> _nonprim_trk_theta_v;
  std::vector<float> _nonprim_trk_phi_v;

  std::vector<float> _nonprim_trk_end_x_v;
  std::vector<float> _nonprim_trk_end_y_v;
  std::vector<float> _nonprim_trk_end_z_v;

  std::vector<float> _nonprim_trk_sce_end_x_v;
  std::vector<float> _nonprim_trk_sce_end_y_v;
  std::vector<float> _nonprim_trk_sce_end_z_v;

  std::vector<float> _nonprim_trk_len_v;

  bool fRecalibrateHits;
  float fEnergyThresholdForMCHits;
  std::vector<float> fADCtoE; // vector of ADC to # of e- conversion [to be taken from production reco2 fhicl files]

  art::InputTag fPCAproducer;
  art::InputTag fPFPproducer;
  art::InputTag fTRKproducer;
  art::InputTag fVTXproducer;
  art::InputTag fSHRproducer;

  art::InputTag fCALOproducer;
  art::InputTag fPIDproducer;

  art::InputTag fCRTVetoproducer; // producer for CRT veto ass tag [anab::T0 <-> recob::OpFlash]
  art::InputTag fCLSproducer;     // cluster associated to PFP
  art::InputTag fMCTproducer;     // MCTruth from neutrino generator
  art::InputTag fMCPproducer;     // MCParticle from Geant4 stage
  art::InputTag fMCFluxproducer;  // MCFlux producer
  art::InputTag fBacktrackTag;
  art::InputTag fHproducer;
  art::InputTag fMCRproducer;
  art::InputTag fSLCproducer; // slice associated to PFP

  void BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col);

  std::map<unsigned int, unsigned int> _pfpmap;
};

NeutronAnalysis::NeutronAnalysis(const fhicl::ParameterSet &p)
{
  fPFPproducer = p.get<art::InputTag>("PFPproducer");
  fPCAproducer = p.get<art::InputTag>("PCAproducer");
  fTRKproducer = p.get<art::InputTag>("TRKproducer");
  fVTXproducer = p.get<art::InputTag>("VTXproducer");
  fSHRproducer = p.get<art::InputTag>("SHRproducer");

  fCRTVetoproducer = p.get<art::InputTag>("CRTVetoproducer", ""); // default is no CRT veto
  fCLSproducer = p.get<art::InputTag>("CLSproducer");
  fMCTproducer = p.get<art::InputTag>("MCTproducer");
  fMCPproducer = p.get<art::InputTag>("MCPproducer");
  fMCFluxproducer = p.get<art::InputTag>("MCFluxproducer");
  fBacktrackTag = p.get<art::InputTag>("BacktrackTag");
  fHproducer = p.get<art::InputTag>("Hproducer");
  fMCRproducer = p.get<art::InputTag>("MCRproducer");
  fSLCproducer = p.get<art::InputTag>("SLCproducer");
  fCALOproducer = p.get<art::InputTag>("CALOproducer");
  fPIDproducer = p.get<art::InputTag>("PIDproducer");
  fRecalibrateHits = p.get<bool>("RecalibrateHits", false);
  fEnergyThresholdForMCHits = p.get<float>("EnergyThresholdForMCHits", 0.1);
  fADCtoE = p.get<std::vector<float>>("ADCtoE");

  // set dedx pdf parameters
  llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
  llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
  llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);

  llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
  llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
  llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);

  llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
  llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
  llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);

  // set correction parameters
  if (fRecalibrateHits)
  {
    llr_pid_calculator.set_corr_par_binning(0, correction_parameters.parameter_correction_edges_pl_0);
    llr_pid_calculator.set_correction_tables(0, correction_parameters.correction_table_pl_0);

    llr_pid_calculator.set_corr_par_binning(1, correction_parameters.parameter_correction_edges_pl_1);
    llr_pid_calculator.set_correction_tables(1, correction_parameters.correction_table_pl_1);

    llr_pid_calculator.set_corr_par_binning(2, correction_parameters.parameter_correction_edges_pl_2);
    llr_pid_calculator.set_correction_tables(2, correction_parameters.correction_table_pl_2);
  }
}


void NeutronAnalysis::configure(fhicl::ParameterSet const &p)
{
}

void NeutronAnalysis::analyzeEvent(art::Event const &e, bool fData)// std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData)
{

std::cout<<"RUNNING NEUTRON ANALYSIS TOOL ANALYZE EVENT"<<std::endl;

  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();

  ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, fPFPproducer,
                                            proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
                                            proxy::withAssociated<recob::Cluster>(fCLSproducer),
                                            proxy::withAssociated<recob::Slice>(fSLCproducer),
                                            proxy::withAssociated<recob::Track>(fTRKproducer),
                                            proxy::withAssociated<recob::Vertex>(fVTXproducer),
                                            proxy::withAssociated<recob::PCAxis>(fPCAproducer),
                                            proxy::withAssociated<recob::Shower>(fSHRproducer),
                                            proxy::withAssociated<recob::SpacePoint>(fPFPproducer));

  // Build larpandora info:
  lar_pandora::LArPandoraHelper larpandora;
  lar_pandora::PFParticleVector pfparticles;
  lar_pandora::PFParticleMap particleMap;
  larpandora.CollectPFParticles(e, "pandora", pfparticles);
  larpandora.BuildPFParticleMap(pfparticles, particleMap);

  _is_nu_reco = false; //flag that shows whether a neutrino has been reconstructed in this event

  TVector3 nuvtx;
  for (auto &pfp : pfp_proxy)
  {
    if (pfp->PdgCode() == 12 || pfp->PdgCode() == 14)
    {
      std::cout << "pfp pdg: " << pfp->PdgCode() <<"\n"<< std::endl;
      double xyz[3] = {};
      auto vtx = pfp.get<recob::Vertex>();
      auto slc = pfp.get<recob::Slice>();
      std::cout<<"slice.ID(): "<<slc.at(0)->ID()<<std::endl;
      if (vtx.size() != 1)
      {
        std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
      }
      else
      {
        // save vertex to array	
        vtx.at(0)->XYZ(xyz);
        nuvtx.SetXYZ(xyz[0], xyz[1], xyz[2]);
	std::cout << "nuvtx position: " << nuvtx[0] << " " << nuvtx[1] << " " << nuvtx[2] <<"\n"<< std::endl;
      }
      _is_nu_reco = true;
      break;
    }
  }

  if (!fData)
  {
  auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
  if (_is_nu_reco)
  {
    std::cout<<"Neutrino Reconstructed"<<std::endl;
    ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
                                                                                        proxy::withAssociated<recob::Hit>(fCLSproducer));

    searchingfornues::ProxyCaloColl_t const &calo_proxy = proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer,
                                                                                                        proxy::withAssociated<anab::Calorimetry>(fCALOproducer));

    searchingfornues::ProxyPIDColl_t const &pid_proxy = proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer,
                                                                                                      proxy::withAssociated<anab::ParticleID>(fPIDproducer));

    BuildPFPMap(pfp_proxy);
    
    // Initialize Backtracker vector and associated MC Particles for Event
    std::vector<searchingfornues::BtPart> btparts_v;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;

    // Fill Backtracker vector with all primary showers and tracks plus all showers and tracks created by secondary protons with process neutronInelastic
    const std::vector<sim::MCShower> &inputMCShower = *(e.getValidHandle<std::vector<sim::MCShower>>(fMCRproducer));
    const std::vector<sim::MCTrack> &inputMCTrack = *(e.getValidHandle<std::vector<sim::MCTrack>>(fMCRproducer));
    art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
    assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
    btparts_v = searchingfornues::initBacktrackingParticleVec(inputMCShower, inputMCTrack, *inputHits, assocMCPart,true);

    // load MCTruth [from geant]

    auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);
    auto mct = mct_h->at(0);
    auto neutrino = mct.GetNeutrino();
    auto nu = neutrino.Nu();
    
    TVector3 mc_nu_vtx;
    std::cout<<"mc nuvtx: "<< mc_nu_vtx[0] << " " << mc_nu_vtx[1] << " " << mc_nu_vtx[2] <<"\n"<< std::endl;
    mc_nu_vtx.SetXYZ( nu.Vx(), nu.Vy(), nu.Vz());
    std::cout<< "True Neutrino PDG Code:"<< nu.PdgCode() << std::endl;
    std::cout<<"mc nuvtx: "<< mc_nu_vtx[0] << " " << mc_nu_vtx[1] << " " << mc_nu_vtx[2] <<"\n"<< std::endl; 
    if (21.5<=nu.Vx() && nu.Vx() <= 234.85 && nu.Vy() >= -95. && nu.Vy() <= 95. && nu.Vz() >= 21.5 && nu.Vz() <= 966.8){std::cout<<"NuVTX in fiducial volume"<<std::endl;}

    int pfp_counter = 0;
    _n_nonprim_pfps = 0;
    for (auto &pfp : pfp_proxy)
    {
      auto trk_v = pfp.get<recob::Track>();
      auto slc_v = pfp.get<recob::Slice>();
      if (slc_v.size() == 1)
      {
	_nonprim_slc_id_v.push_back(slc_v.at(0)->ID());
      }
      else
      {
	_nonprim_slc_id_v.push_back(std::numeric_limits<int>::lowest());
      }

      std::__cxx11::string backtracked_process;

      nonprim_pfpdg.push_back(pfp->PdgCode());
      _nonprim_generation.push_back(larpandora.GetGeneration(particleMap, particleMap.at(pfp->Self())));
      // get hits associated to this PFParticle through the clusters
      std::vector<art::Ptr<recob::Hit>> hit_v;
      auto clus_pxy_v = pfp.get<recob::Cluster>();
      for (auto ass_clus : clus_pxy_v)
      {
        // get cluster proxy
        const auto &clus = clus_proxy[ass_clus.key()];
        auto clus_hit_v = clus.get<recob::Hit>();
        for (const auto &hit : clus_hit_v)
        {
          hit_v.push_back(hit);
        } // for all hits in this cluster
      } // for all clusters associated to PFP

      if (clus_pxy_v.size() != 0)
      {
        float purity = 0., completeness = 0., overlay_purity = 0.;
        int ibt = searchingfornues::getAssocBtPart(hit_v, assocMCPart, btparts_v, purity, completeness, overlay_purity);
        if (ibt >= 0)
        {
	  std::cout<<"PFP number: "<<pfp_counter+1<<std::endl;
          auto &mcp = btparts_v.at(ibt);  //replace square braces with .at()
          int PDG = mcp.pdg;
	  std::cout<<"pdg: "<<PDG<<std::endl;
	  // This loop checks to see if the backtracked particle is a neutron Inelasticly produced particle
	  for (size_t p = 0; p < mcp_h->size(); p++)
	  {
	    int mcptrkid = mcp_h->at(p).TrackId();
	    int ibttrkid = mcp.tids.at(0);
	    if (ibttrkid == mcptrkid)
	    {
	      backtracked_process = mcp_h->at(p).Process();
	      break;
	    }
	  }
          _nonprim_backtracked_e.push_back(mcp.e);
          _nonprim_backtracked_tid.push_back(mcp.tids.at(0));
          _nonprim_backtracked_pdg.push_back(PDG);
          _nonprim_backtracked_purity.push_back(purity);
          _nonprim_backtracked_completeness.push_back(completeness);
          _nonprim_backtracked_overlay_purity.push_back(overlay_purity);
          _nonprim_backtracked_px.push_back(mcp.px);
          _nonprim_backtracked_py.push_back(mcp.py);
          _nonprim_backtracked_pz.push_back(mcp.pz);
          _nonprim_backtracked_start_x.push_back(mcp.start_x);
          _nonprim_backtracked_start_y.push_back(mcp.start_y);
          _nonprim_backtracked_start_z.push_back(mcp.start_z);
          _nonprim_backtracked_start_t.push_back(mcp.start_t);
          _nonprim_backtracked_start_U.push_back(searchingfornues::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 0));
          _nonprim_backtracked_start_V.push_back(searchingfornues::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 1));
          _nonprim_backtracked_start_Y.push_back(searchingfornues::YZtoPlanecoordinate(mcp.start_y, mcp.start_z, 2));
	  _nonprim_backtracked_process.push_back(backtracked_process);
	}
	else
	{
	  _nonprim_backtracked_e.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_tid.push_back(std::numeric_limits<int>::lowest());
          _nonprim_backtracked_pdg.push_back(0);
          _nonprim_backtracked_purity.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_completeness.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_overlay_purity.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_px.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_py.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_pz.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_start_x.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_start_y.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_start_z.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_start_t.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_start_U.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_start_V.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_start_Y.push_back(std::numeric_limits<float>::lowest());
          _nonprim_backtracked_process.push_back("N/A");
	}
      } //Cluster Loop
      if (trk_v.size() == 1){
        _nonprim_trk_score_v.push_back(searchingfornues::GetTrackShowerScore(pfp));
	auto trk = trk_v.at(0); 
	auto calo_v = calo_proxy[trk.key()].get<anab::Calorimetry>();
	
	_nonprim_trk_calo_energy_u_v.push_back(-1); 
	_nonprim_trk_calo_energy_v_v.push_back(-1);
	_nonprim_trk_calo_energy_y_v.push_back(-1);

        _nonprim_trk_len_v.push_back(searchingfornues::GetSCECorrTrackLength(trk));
        _nonprim_trk_start_x_v.push_back(trk->Start().X());
        _nonprim_trk_start_y_v.push_back(trk->Start().Y());
        _nonprim_trk_start_z_v.push_back(trk->Start().Z());
        _nonprim_trk_end_x_v.push_back(trk->End().X());
        _nonprim_trk_end_y_v.push_back(trk->End().Y());
        _nonprim_trk_end_z_v.push_back(trk->End().Z());
	  
        float _nonprim_trk_start_sce[3];
        searchingfornues::ApplySCECorrectionXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z(), _nonprim_trk_start_sce);
	_nonprim_trk_sce_start_x_v.push_back(_nonprim_trk_start_sce[0]);
        _nonprim_trk_sce_start_y_v.push_back(_nonprim_trk_start_sce[1]);
        _nonprim_trk_sce_start_z_v.push_back(_nonprim_trk_start_sce[2]);

        float _nonprim_trk_end_sce[3];
        searchingfornues::ApplySCECorrectionXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z(), _nonprim_trk_end_sce);
        _nonprim_trk_sce_end_x_v.push_back(_nonprim_trk_end_sce[0]);
        _nonprim_trk_sce_end_y_v.push_back(_nonprim_trk_end_sce[1]);
        _nonprim_trk_sce_end_z_v.push_back(_nonprim_trk_end_sce[2]);

        _nonprim_trk_theta_v.push_back(trk->Theta());
        _nonprim_trk_phi_v.push_back(trk->Phi());

        TVector3 nonprim_trk_vtx_v;
        nonprim_trk_vtx_v.SetXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
        nonprim_trk_vtx_v -= nuvtx;
        _nonprim_trk_distance_v.push_back(nonprim_trk_vtx_v.Mag());

	//PID LLR calculator
	_nonprim_trk_llr_pid_u_v.push_back(0);
      	_nonprim_trk_llr_pid_v_v.push_back(0);
      	_nonprim_trk_llr_pid_y_v.push_back(0);
      	_nonprim_trk_llr_pid_v.push_back(0);
      	_nonprim_trk_llr_pid_score_v.push_back(0);
	_nonprim_trk_nhits_u_v.push_back(0);
	_nonprim_trk_nhits_v_v.push_back(0);
	_nonprim_trk_nhits_y_v.push_back(0);

      	for (auto const &calo : calo_v)
        {
	  auto const &plane = calo->PlaneID().Plane;
	  if (plane == 4294967295){continue;}
		
          auto const &dedx_values = calo->dEdx();
	  auto const &dqdx_values = calo->dQdx();
          auto const &rr = calo->ResidualRange();
          auto const &pitch = calo->TrkPitchVec();
          auto const& xyz_v = calo->XYZ();
          std::vector<std::vector<float>> par_values;
          par_values.push_back(rr);
	  par_values.push_back(pitch);

          float calo_energy = 0;
	  std::vector<float> dqdx_values_corrected, dedx_values_corrected;

	  if (fData || !fRecalibrateHits)
          {
            dqdx_values_corrected = dqdx_values;
          }
          else
          {
            dqdx_values_corrected = llr_pid_calculator.correct_many_hits_one_plane(calo, trk.value(), assocMCPart, fRecalibrateHits, fEnergyThresholdForMCHits, false);
          }

          for (size_t i = 0; i < dqdx_values_corrected.size(); i++)
          {
            float aux_dedx;
            aux_dedx = searchingfornues::ModBoxCorrection(dqdx_values_corrected[i]*fADCtoE[plane], xyz_v[i].X(), xyz_v[i].Y(), xyz_v[i].Z());
            dedx_values_corrected.push_back(aux_dedx);
            calo_energy += aux_dedx * pitch[i];
          }
	  float trk_nhits = dedx_values.size();
	  float llr_pid = llr_pid_calculator.LLR_many_hits_one_plane(dedx_values_corrected, par_values, plane);
	  if (plane == 0)
          {
	    _nonprim_trk_llr_pid_u_v.back() = llr_pid;
            _nonprim_trk_calo_energy_u_v.back() = calo_energy;
	    _nonprim_trk_nhits_u_v.back() = trk_nhits;
	  }
          else if (plane == 1)
	  {
            _nonprim_trk_llr_pid_v_v.back() = llr_pid;
	    _nonprim_trk_calo_energy_v_v.back() = calo_energy;
            _nonprim_trk_nhits_v_v.back() = trk_nhits;
          }
	  else if (plane == 2)
          {
            _nonprim_trk_llr_pid_y_v.back() = llr_pid;
	    _nonprim_trk_calo_energy_y_v.back() = calo_energy;
            _nonprim_trk_nhits_y_v.back() = trk_nhits;
          }
	  _nonprim_trk_llr_pid_v.back() += llr_pid;
	  } // Calo loop
	  _nonprim_trk_llr_pid_score_v.back() = atan(_nonprim_trk_llr_pid_v.back() / 100.) * 2 / 3.14159266;
	  }//if trk_v.size() == 1
       //Cluster Loop 
      pfp_counter += 1;
    } //End pfp loop
  std::cout<<"Number of PFP Particles: "<<pfp_counter<<std::endl;
  _n_nonprim_pfps = pfp_counter;

  for (size_t p = 0; p < mcp_h->size(); p++)
  {
    auto mcp = mcp_h->at(p);
    if (mcp.PdgCode() != 22 && mcp.PdgCode() != 11 && mcp.PdgCode() != -11)
    {

      TVector3 all_mc_vtx_v;
      all_mc_vtx_v.SetXYZ(mcp.Vx(), mcp.Vy(), mcp.Vz());
      all_mc_vtx_v -= nuvtx;
      _all_mc_distance.push_back(all_mc_vtx_v.Mag());      

      _all_mc_mother.push_back(mcp.Mother());
      _all_mc_trkid.push_back(mcp.TrackId());
      _all_mc_pdg.push_back(mcp.PdgCode());
      _all_mc_E.push_back(mcp.E());
      _all_mc_px.push_back(mcp.Px());
      _all_mc_py.push_back(mcp.Py());
      _all_mc_pz.push_back(mcp.Pz());
      _all_mc_vx.push_back(mcp.Vx());
      _all_mc_vy.push_back(mcp.Vy());
      _all_mc_vz.push_back(mcp.Vz());
      _all_mc_endx.push_back(mcp.EndX());
      _all_mc_endy.push_back(mcp.EndY());
      _all_mc_endz.push_back(mcp.EndZ());
      _all_mc_process.push_back(mcp.Process());
    }
  }
  }  // Neutrino Flag

  else  // If the event does not have a reconstructed neutrino, fill all reco variables with default values
  {
  std::cout<<"Neutrino Not Reconstructed"<<std::endl;
  for (size_t i_pfp = 0; i_pfp < pfp_proxy.size(); i_pfp++)
  {
    _nonprim_backtracked_e.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_tid.push_back(std::numeric_limits<int>::lowest());
    _nonprim_backtracked_pdg.push_back(0);
    _nonprim_backtracked_purity.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_completeness.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_overlay_purity.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_px.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_py.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_pz.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_start_x.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_start_y.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_start_z.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_start_t.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_start_U.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_start_V.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_start_Y.push_back(std::numeric_limits<float>::lowest());
    _nonprim_backtracked_process.push_back("N/A");

    _nonprim_trk_score_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_calo_energy_u_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_calo_energy_v_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_calo_energy_y_v.push_back(std::numeric_limits<float>::lowest()); 
    _nonprim_trk_len_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_start_x_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_start_y_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_start_z_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_end_x_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_end_y_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_end_z_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_sce_start_x_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_sce_start_y_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_sce_start_z_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_sce_end_x_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_sce_end_y_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_sce_end_z_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_theta_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_phi_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_distance_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_llr_pid_u_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_llr_pid_v_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_llr_pid_y_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_llr_pid_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_llr_pid_score_v.push_back(std::numeric_limits<float>::lowest());
    _nonprim_trk_nhits_u_v.push_back(std::numeric_limits<int>::lowest());
    _nonprim_trk_nhits_v_v.push_back(std::numeric_limits<int>::lowest());
    _nonprim_trk_nhits_y_v.push_back(std::numeric_limits<int>::lowest());
  }

  for (size_t p = 0; p < mcp_h->size(); p++)
  {
    auto mcp = mcp_h->at(p);
    if (mcp.PdgCode() != 22 && mcp.PdgCode() != 11 && mcp.PdgCode() != -11)
    {

      TVector3 all_mc_vtx_v;
      all_mc_vtx_v.SetXYZ(mcp.Vx(), mcp.Vy(), mcp.Vz());
      all_mc_vtx_v -= nuvtx;
      _all_mc_distance.push_back(all_mc_vtx_v.Mag());  //The displacement of the mcparticle from the reconstructed neutrino vertex.

      _all_mc_mother.push_back(mcp.Mother());
      _all_mc_trkid.push_back(mcp.TrackId());
      _all_mc_pdg.push_back(mcp.PdgCode());
      _all_mc_E.push_back(mcp.E());
      _all_mc_px.push_back(mcp.Px());
      _all_mc_py.push_back(mcp.Py());
      _all_mc_pz.push_back(mcp.Pz());
      _all_mc_vx.push_back(mcp.Vx());
      _all_mc_vy.push_back(mcp.Vy());
      _all_mc_vz.push_back(mcp.Vz());
      _all_mc_endx.push_back(mcp.EndX());
      _all_mc_endy.push_back(mcp.EndY());
      _all_mc_endz.push_back(mcp.EndZ());
      _all_mc_process.push_back(mcp.Process());
    }
  }
  }  // Neutrino not reconstructed loop

  }  // End MC

  std::cout<<"\nEND NEUTRON ANALYSIS TOOL\n"<<std::endl;

}  // End AnalyzeEvent

void NeutronAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
}

void NeutronAnalysis::BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col)
{

  _pfpmap.clear();

  unsigned int p = 0;
  for (const auto &pfp_pxy : pfp_pxy_col)
  {
    _pfpmap[pfp_pxy->Self()] = p;
    p++;
  }
  return;
}

void NeutronAnalysis::setBranches(TTree *_tree)
{
  _tree->Branch("is_nu_reco", &_is_nu_reco, "_is_nu_reco/B");

  _tree->Branch("nonprim_backtracked_pdg", "std::vector<int>", &_nonprim_backtracked_pdg);
  _tree->Branch("nonprim_backtracked_e", "std::vector<float>", &_nonprim_backtracked_e);
  _tree->Branch("nonprim_backtracked_tid", "std::vector<int>", &_nonprim_backtracked_tid);
  _tree->Branch("nonprim_backtracked_purity", "std::vector<float>", &_nonprim_backtracked_purity);
  _tree->Branch("nonprim_backtracked_completeness", "std::vector<float>", &_nonprim_backtracked_completeness);
  _tree->Branch("nonprim_backtracked_overlay_purity", "std::vector<float>", &_nonprim_backtracked_overlay_purity);

  _tree->Branch("nonprim_backtracked_px", "std::vector<float>", &_nonprim_backtracked_px);
  _tree->Branch("nonprim_backtracked_py", "std::vector<float>", &_nonprim_backtracked_py);
  _tree->Branch("nonprim_backtracked_pz", "std::vector<float>", &_nonprim_backtracked_pz);

  _tree->Branch("nonprim_backtracked_start_x", "std::vector<float>", &_nonprim_backtracked_start_x);
  _tree->Branch("nonprim_backtracked_start_y", "std::vector<float>", &_nonprim_backtracked_start_y);
  _tree->Branch("nonprim_backtracked_start_z", "std::vector<float>", &_nonprim_backtracked_start_z);
  _tree->Branch("nonprim_backtracked_start_t", "std::vector<float>", &_nonprim_backtracked_start_t);
  _tree->Branch("nonprim_backtracked_start_U", "std::vector<float>", &_nonprim_backtracked_start_U);
  _tree->Branch("nonprim_backtracked_start_V", "std::vector<float>", &_nonprim_backtracked_start_V);
  _tree->Branch("nonprim_backtracked_start_Y", "std::vector<float>", &_nonprim_backtracked_start_Y);
  _tree->Branch("nonprim_backtracked_sce_start_x", "std::vector<float>", &_nonprim_backtracked_sce_start_x);
  _tree->Branch("nonprim_backtracked_sce_start_y", "std::vector<float>", &_nonprim_backtracked_sce_start_y);
  _tree->Branch("nonprim_backtracked_sce_start_z", "std::vector<float>", &_nonprim_backtracked_sce_start_z);
  _tree->Branch("nonprim_backtracked_sce_start_U", "std::vector<float>", &_nonprim_backtracked_sce_start_U);
  _tree->Branch("nonprim_backtracked_sce_start_V", "std::vector<float>", &_nonprim_backtracked_sce_start_V);
  _tree->Branch("nonprim_backtracked_sce_start_Y", "std::vector<float>", &_nonprim_backtracked_sce_start_Y);
  _tree->Branch("nonprim_backtracked_process","std::vector<std::__cxx11::string>", &_nonprim_backtracked_process);

  _tree->Branch("n_nonprim_pfps", &_n_nonprim_pfps, "n_nonprim_pfps/I");

  _tree->Branch("nonprim_pfpdg", "std::vector<int>", &nonprim_pfpdg);
  _tree->Branch("nonprim_pfp_generation_v", "std::vector< uint >", &_nonprim_generation);
  _tree->Branch("nonprim_topo_score_v","std::vector< float >", &_nonprim_topo_score_v);
  _tree->Branch("nonprim_slc_id_v","std::vector< int >",&_nonprim_slc_id_v);

  _tree->Branch("all_mc_pdg", "std::vector< int >", &_all_mc_pdg);
  _tree->Branch("all_mc_E", "std::vector< float >", &_all_mc_E);

  _tree->Branch("all_mc_vx", "std::vector< float >", &_all_mc_vx);
  _tree->Branch("all_mc_vy", "std::vector< float >", &_all_mc_vy);
  _tree->Branch("all_mc_vz", "std::vector< float >", &_all_mc_vz);

  _tree->Branch("all_mc_endx", "std::vector< float >", &_all_mc_endx);
  _tree->Branch("all_mc_endy", "std::vector< float >", &_all_mc_endy);
  _tree->Branch("all_mc_endz", "std::vector< float >", &_all_mc_endz);

  _tree->Branch("all_mc_px", "std::vector< float >", &_all_mc_px);
  _tree->Branch("all_mc_py", "std::vector< float >", &_all_mc_py);
  _tree->Branch("all_mc_pz", "std::vector< float >", &_all_mc_pz);

  _tree->Branch("all_mc_mother", "std::vector< int >", &_all_mc_mother);
  _tree->Branch("all_mc_trkid", "std::vector< int >", &_all_mc_trkid);
  _tree->Branch("all_mc_process", "std::vector< std::__cxx11::string >", &_all_mc_process);
  _tree->Branch("all_mc_distance", "std::vector< float >", &_all_mc_distance);

  _tree->Branch("nonprim_trk_score_v", "std::vector<float>", &_nonprim_trk_score_v);
  
  _tree->Branch("nonprim_trk_llr_pid_u_v", "std::vector<float>", &_nonprim_trk_llr_pid_u_v);
  _tree->Branch("nonprim_trk_llr_pid_v_v", "std::vector<float>", &_nonprim_trk_llr_pid_v_v);
  _tree->Branch("nonprim_trk_llr_pid_y_v", "std::vector<float>", &_nonprim_trk_llr_pid_y_v);
  _tree->Branch("nonprim_trk_llr_pid_v", "std::vector<float>", &_nonprim_trk_llr_pid_v);
  _tree->Branch("nonprim_trk_llr_pid_score_v", "std::vector<float>", &_nonprim_trk_llr_pid_score_v);

  _tree->Branch("nonprim_trk_calo_energy_u_v", "std::vector< float >", &_nonprim_trk_calo_energy_u_v);
  _tree->Branch("nonprim_trk_calo_energy_v_v", "std::vector< float >", &_nonprim_trk_calo_energy_v_v);
  _tree->Branch("nonprim_trk_calo_energy_y_v", "std::vector< float >", &_nonprim_trk_calo_energy_y_v);

  _tree->Branch("nonprim_trk_nhits_u_v", "std::vector<int>", &_nonprim_trk_nhits_u_v);
  _tree->Branch("nonprim_trk_nhits_v_v", "std::vector<int>", &_nonprim_trk_nhits_v_v);
  _tree->Branch("nonprim_trk_nhits_y_v", "std::vector<int>", &_nonprim_trk_nhits_y_v);

  _tree->Branch("nonprim_trk_start_x_v", "std::vector< float >", &_nonprim_trk_start_x_v);
  _tree->Branch("nonprim_trk_start_y_v", "std::vector< float >", &_nonprim_trk_start_y_v);
  _tree->Branch("nonprim_trk_start_z_v", "std::vector< float >", &_nonprim_trk_start_z_v);

  _tree->Branch("nonprim_trk_sce_start_x_v", "std::vector< float >", &_nonprim_trk_sce_start_x_v);
  _tree->Branch("nonprim_trk_sce_start_y_v", "std::vector< float >", &_nonprim_trk_sce_start_y_v);
  _tree->Branch("nonprim_trk_sce_start_z_v", "std::vector< float >", &_nonprim_trk_sce_start_z_v);

  _tree->Branch("nonprim_trk_distance_v", "std::vector< float >", &_nonprim_trk_distance_v);
  _tree->Branch("nonprim_trk_theta_v", "std::vector< float >", &_nonprim_trk_theta_v);
  _tree->Branch("nonprim_trk_phi_v", "std::vector< float >", &_nonprim_trk_phi_v);

  _tree->Branch("nonprim_trk_end_x_v", "std::vector< float >", &_nonprim_trk_end_x_v);
  _tree->Branch("nonprim_trk_end_y_v", "std::vector< float >", &_nonprim_trk_end_y_v);
  _tree->Branch("nonprim_trk_end_z_v", "std::vector< float >", &_nonprim_trk_end_z_v);

  _tree->Branch("nonprim_trk_sce_end_x_v", "std::vector< float >", &_nonprim_trk_sce_end_x_v);
  _tree->Branch("nonprim_trk_sce_end_y_v", "std::vector< float >", &_nonprim_trk_sce_end_y_v);
  _tree->Branch("nonprim_trk_sce_end_z_v", "std::vector< float >", &_nonprim_trk_sce_end_z_v);

  _tree->Branch("nonprim_trk_len_v", "std::vector< float >", &_nonprim_trk_len_v);

} //End setBranches

void NeutronAnalysis::resetTTree(TTree *_tree)
{
  _n_nonprim_pfps = 0;

  _is_nu_reco = false;

  _nonprim_backtracked_e.clear();
  _nonprim_backtracked_tid.clear();
  _nonprim_backtracked_pdg.clear();
  _nonprim_backtracked_completeness.clear();
  _nonprim_backtracked_purity.clear();
  _nonprim_backtracked_overlay_purity.clear();

  _nonprim_backtracked_px.clear();
  _nonprim_backtracked_py.clear();
  _nonprim_backtracked_pz.clear();

  _nonprim_backtracked_start_x.clear();
  _nonprim_backtracked_start_y.clear();
  _nonprim_backtracked_start_z.clear();
  _nonprim_backtracked_start_t.clear();
  _nonprim_backtracked_start_U.clear();
  _nonprim_backtracked_start_V.clear();
  _nonprim_backtracked_start_Y.clear();
  _nonprim_backtracked_sce_start_x.clear();
  _nonprim_backtracked_sce_start_y.clear();
  _nonprim_backtracked_sce_start_z.clear();
  _nonprim_backtracked_sce_start_U.clear();
  _nonprim_backtracked_sce_start_V.clear();
  _nonprim_backtracked_sce_start_Y.clear();
  _nonprim_backtracked_process.clear();

  nonprim_pfpdg.clear();
  _nonprim_generation.clear();
  _nonprim_topo_score_v.clear();
  _nonprim_slc_id_v.clear();

  _all_mc_E.clear();
  _all_mc_pdg.clear();

  _all_mc_px.clear();
  _all_mc_py.clear();
  _all_mc_pz.clear();

  _all_mc_vx.clear();
  _all_mc_vy.clear();
  _all_mc_vz.clear();

  _all_mc_endx.clear();
  _all_mc_endy.clear();
  _all_mc_endz.clear();

  _all_mc_mother.clear();
  _all_mc_trkid.clear();
  _all_mc_process.clear();
  _all_mc_distance.clear();

  _nonprim_trk_score_v.clear();

  _nonprim_trk_llr_pid_u_v.clear();
  _nonprim_trk_llr_pid_v_v.clear();
  _nonprim_trk_llr_pid_y_v.clear();
  _nonprim_trk_llr_pid_v.clear();
  _nonprim_trk_llr_pid_score_v.clear();

  _nonprim_trk_calo_energy_u_v.clear();
  _nonprim_trk_calo_energy_v_v.clear();
  _nonprim_trk_calo_energy_y_v.clear();

  _nonprim_trk_nhits_u_v.clear();
  _nonprim_trk_nhits_v_v.clear();
  _nonprim_trk_nhits_y_v.clear();

  _nonprim_trk_start_x_v.clear();
  _nonprim_trk_start_y_v.clear();
  _nonprim_trk_start_z_v.clear();

  _nonprim_trk_sce_start_x_v.clear();
  _nonprim_trk_sce_start_y_v.clear();
  _nonprim_trk_sce_start_z_v.clear();

  _nonprim_trk_distance_v.clear();
  _nonprim_trk_theta_v.clear();
  _nonprim_trk_phi_v.clear();

  _nonprim_trk_end_x_v.clear();
  _nonprim_trk_end_y_v.clear();
  _nonprim_trk_end_z_v.clear();

  _nonprim_trk_sce_end_x_v.clear();
  _nonprim_trk_sce_end_y_v.clear();
  _nonprim_trk_sce_end_z_v.clear();

  _nonprim_trk_len_v.clear();
} //End resetTree

DEFINE_ART_CLASS_TOOL(NeutronAnalysis)
} // namespace analysis

#endif

