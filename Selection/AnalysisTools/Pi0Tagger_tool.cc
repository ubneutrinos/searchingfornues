#ifndef SELECTION_SELECTIONEXAMPLE_CXX
#define SELECTION_SELECTIONEXAMPLE_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/TrackFitterFunctions.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/LLR_PID.h"
#include "../CommonDefs/LLRPID_correction_lookup.h"
#include "../CommonDefs/LLRPID_electron_photon_lookup.h"

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       Pi0Tagger
    // File:        Pi0Tagger.cc
    //
    //              A basic selection example
    //
    // Configuration parameters:
    //
    // TBD
    //
    // Created by David Caratelli (davidc@fnal.gov) on 01/30/2019
    //
    ////////////////////////////////////////////////////////////////////////
    
  class Pi0Tagger : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    Pi0Tagger(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~Pi0Tagger(){};
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset);
    
    /**
     * @brief Analysis function
     */
    void analyzeEvent(art::Event const& e, bool fData) override;

    /**
     * @brief Analyze slice
     */
  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

    /**
     * @brief set branches for TTree
     */
    void setBranches(TTree* _tree);

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree* _tree);
    
  private:

    std::pair<double,double> VtxCompatibility(const TVector3& nuvtx, const TVector3& shrvtx, const TVector3& shrdir);

    /**
     * @brief calculate PFP energy based on hits associated to clusters
     */
    template <typename T> float PFPEnergy(const T& ass_clus_v);

    /**
     * @brief read truth info associated to pi0
     */
    void ReadTruth(art::Event const& e);

    void Reset();

    // TTree variables
    
    // truth variables
    int _pi0_ispi0;
    TVector3 _pi0_mcgamma0_mom, _pi0_mcgamma1_mom;
    float _pi0_mcgamma0_e, _pi0_mcgamma0_px, _pi0_mcgamma0_py, _pi0_mcgamma0_pz;
    float _pi0_mcgamma1_e, _pi0_mcgamma1_px, _pi0_mcgamma1_py, _pi0_mcgamma1_pz;
    float _pi0_mcrcdot0, _pi0_mcrcdot1; // dot product between MC gamma and RC gamma. Each MC gamma matched to best RC gamma
    float _pi0_mcrce0, _pi0_mcrce1; // energy of reconstructed photon associated to the MC one

    // reco variables
    int _pi0_nshower;
    int _pi0_ntrack;
    int _pi0_ngamma;
    float _pi0_radlen1, _pi0_radlen2;
    float _pi0_dot1, _pi0_dot2;
    float _pi0_energy1_Y, _pi0_energy2_Y;
    float _pi0_dedx1_Y, _pi0_dedx2_Y;
    float _pi0_dedx1_fit_Y, _pi0_dedx2_fit_Y;
    float _pi0_energy1_V, _pi0_energy2_V;
    float _pi0_dedx1_V, _pi0_dedx2_V;
    float _pi0_dedx1_fit_V, _pi0_dedx2_fit_V;
    float _pi0_energy1_U, _pi0_energy2_U;
    float _pi0_dedx1_U, _pi0_dedx2_U;
    float _pi0_dedx1_fit_U, _pi0_dedx2_fit_U;
    float _pi0_shrscore1, _pi0_shrscore2;
    float _pi0_gammadot;
    float _pi0_mass_Y, _pi0_mass_V, _pi0_mass_U;
    float _pi0_rc_vtx_x, _pi0_rc_vtx_y, _pi0_rc_vtx_z; // reco neutrino vertex
    

    // module-specific settings
    bool _pi0_onlyshower;   // should we use only showers to reconstruct pi0s?
    float _pi0_dmin;        // what is the minimum distance of the trk/shr vertex to the neutrino vertex?
    float _pi0_dotmin;      // maximum dot product between shower direction and vtx->start vector
    float _pi0_trkshrscore; // score on which to cut for track/shower classification

    bool fLocaldEdx;               // use local dE/dx from calorimetry (true) or globally convert Q -> MeV (false)
    std::vector<float> fADCtoE; // vector of ADC to # of e- conversion [to be taken from production reco2 fhicl files]
  bool fRecalibrateHits;
  float fEnergyThresholdForMCHits;

    art::InputTag fTRKproducer;
    art::InputTag fCALproducer;
    art::InputTag fHitproducer;
    art::InputTag fBacktrackTag;

  // re-calibration tools
  searchingfornues::LLRPID llr_pid_calculator_shr;
  searchingfornues::ElectronPhotonLookUpParameters electronphoton_parameters;
  searchingfornues::CorrectionLookUpParameters correction_parameters;
    
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  Pi0Tagger::Pi0Tagger(const fhicl::ParameterSet& pset)
  {
    configure(pset);

    // configure shower PID tools
    llr_pid_calculator_shr.set_dedx_binning(0, electronphoton_parameters.dedx_edges_pl_0);
    llr_pid_calculator_shr.set_par_binning(0, electronphoton_parameters.parameters_edges_pl_0);
    llr_pid_calculator_shr.set_lookup_tables(0, electronphoton_parameters.dedx_pdf_pl_0);
    
    llr_pid_calculator_shr.set_dedx_binning(1, electronphoton_parameters.dedx_edges_pl_1);
    llr_pid_calculator_shr.set_par_binning(1, electronphoton_parameters.parameters_edges_pl_1);
    llr_pid_calculator_shr.set_lookup_tables(1, electronphoton_parameters.dedx_pdf_pl_1);
    
    llr_pid_calculator_shr.set_dedx_binning(2, electronphoton_parameters.dedx_edges_pl_2);
    llr_pid_calculator_shr.set_par_binning(2, electronphoton_parameters.parameters_edges_pl_2);
    llr_pid_calculator_shr.set_lookup_tables(2, electronphoton_parameters.dedx_pdf_pl_2);

    if (fRecalibrateHits){
      llr_pid_calculator_shr.set_corr_par_binning(0, correction_parameters.parameter_correction_edges_pl_0);
      llr_pid_calculator_shr.set_correction_tables(0, correction_parameters.correction_table_pl_0);
      llr_pid_calculator_shr.set_corr_par_binning(1, correction_parameters.parameter_correction_edges_pl_1);
      llr_pid_calculator_shr.set_correction_tables(1, correction_parameters.correction_table_pl_1);
      llr_pid_calculator_shr.set_corr_par_binning(2, correction_parameters.parameter_correction_edges_pl_2);
      llr_pid_calculator_shr.set_correction_tables(2, correction_parameters.correction_table_pl_2);
    }
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void Pi0Tagger::configure(fhicl::ParameterSet const & pset)
  {
    _pi0_onlyshower  = pset.get< bool >  ("onlyshower" );
    _pi0_dotmin      = pset.get< float > ("dotmin"     );
    _pi0_dmin        = pset.get< float > ("dmin"       );
    _pi0_trkshrscore = pset.get< float > ("trkshrscore");

    fADCtoE = pset.get<std::vector<float>>("ADCtoE");
    fLocaldEdx = pset.get<bool>("LocaldEdx", false);       // use dE/dx from calo?
    fRecalibrateHits = pset.get<bool>("RecalibrateHits", false);
    fEnergyThresholdForMCHits = pset.get<float>("EnergyThresholdForMCHits", 0.1);

    fTRKproducer = pset.get< art::InputTag > ("TRKproducer", "");
    fCALproducer = pset.get< art::InputTag > ("CALproducer", "");
    fHitproducer     = pset.get<art::InputTag>("Hitproducer" , "");
    fBacktrackTag = pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");

}

  
void Pi0Tagger::analyzeEvent(art::Event const &e, bool fData)
{

  if (!fData)
    ReadTruth(e);
  
}

  
  //----------------------------------------------------------------------------
  /// selectEvent
  ///
  /// Arguments:
  ///
  /// art::Event
  /// slice track pointer vector
  /// slice shower pointer vector
  ///
  void Pi0Tagger::analyzeSlice(art::Event const& e, std::vector<ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected)
  {

    TVector3 nuvtx;
    Double_t xyz[3] = {};

    Reset();

    // load backtrack information
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
    if (!fData)
      {
        art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHitproducer);
        assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
    }


    searchingfornues::ProxyCaloColl_t const* tkcalo_proxy = NULL;
    if (fTRKproducer!="") {
      tkcalo_proxy = new searchingfornues::ProxyCaloColl_t( proxy::getCollection<std::vector<recob::Track> >(e,fTRKproducer,proxy::withAssociated<anab::Calorimetry>(fCALproducer)) );
    }
    


    TVector3 gammadir1, gammadir2;

    // vector of pfp indices for showers
    std::vector<float> shr_energy_v;
    // vector of shower energies
    std::vector<short> pfp_idx_v;

    short pfp_ctr = 0;

    // loop over all PFPs and get metadata, find the neutrino and save its vertex
    for (const auto& pfp_pxy : slice_pfp_v) {

      pfp_ctr += 1;

      auto PDG = fabs(pfp_pxy->PdgCode());

      if ( (PDG == 12) || (PDG == 14) ) {

	// grab vertex
	auto vtx = pfp_pxy.get<recob::Vertex>();
	if (vtx.size() != 1) {
	  std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
	  return;
	}
	
	// save vertex to array
	vtx.at(0)->XYZ(xyz);
	nuvtx = TVector3(xyz[0],xyz[1],xyz[2]);

	_pi0_rc_vtx_x = nuvtx.X();
	_pi0_rc_vtx_y = nuvtx.Y();
	_pi0_rc_vtx_z = nuvtx.Z();

      }// if neutrino PFP

      else { // if not neutrino PFP

	// grab shower/track score
	auto trkshrscore = searchingfornues::GetTrackShowerScore(pfp_pxy);

	auto nshr = pfp_pxy.get<recob::Shower>().size();
	//auto ntrk = pfp_pxy.get<recob::Track>().size();

	// 1 -> track-like
	if (trkshrscore > _pi0_trkshrscore)  continue;

	if (nshr != 1) continue;

	auto const& shr = pfp_pxy.get<recob::Shower>().at(0);

	auto energy = shr->Energy()[2];
	
	auto vtxcompat = VtxCompatibility(nuvtx, shr->ShowerStart(), shr->Direction());

	// if blank result, continue
	if ( (vtxcompat.first == -1) && (vtxcompat.second == -1) ) continue;
	// if too close to vertex or too mis-algined, continue
	//if ( (vtxcompat.second < _pi0_dmin) || (vtxcompat.first < _pi0_dotmin) ) continue;

	shr_energy_v.push_back( energy );
	pfp_idx_v.push_back( pfp_ctr - 1 );

      }// if not the neutrino PFP

    }// for all PFP

    // if we did not find any shower, return
    if (shr_energy_v.size() == 0) return;

    // if a single shower, backtrack and quit
    if (shr_energy_v.size() == 1) {

      auto shr = slice_pfp_v.at(pfp_idx_v[0]).get<recob::Shower>().at(0);
      
      float dot = shr->Direction().Dot(_pi0_mcgamma0_mom);
      if (dot > _pi0_mcrcdot0) {
	_pi0_mcrcdot0 = dot;
	_pi0_mcrce0   = shr->Energy()[2];
      }// if new most aligned shower
      dot = shr->Direction().Dot(_pi0_mcgamma1_mom);
      if (dot > _pi0_mcrcdot1) {
	_pi0_mcrcdot1 = dot;
	_pi0_mcrce1   = shr->Energy()[2];
      }// if new most aligned shower

      return;
    }// if a single shower

    // if two or more, sort by energy
    double e1 = 0; // energy of highest energy reco shower
    short  i1 = 0; // pfp index for highest energy reco shower
    double e2 = 0; // energy of 2nd highest energy reco shower
    short  i2 = 0; // pfp index for 2nd highest energy reco shower

    // find highest energy shower
    for (size_t n0=0; n0 < shr_energy_v.size(); n0++) {
      if (shr_energy_v[n0] > e1) {
	e1 = shr_energy_v[n0];
	i1 = pfp_idx_v[n0];
      } 
    }
    // find second highest energy shower
    for (size_t n0=0; n0 < shr_energy_v.size(); n0++) {
      if (pfp_idx_v[n0] == i1) continue;
      if (shr_energy_v[n0] > e2) {
	e2 = shr_energy_v[n0];
	i2 = pfp_idx_v[n0];
      }
    }

    auto shr1 = slice_pfp_v.at(i1).get<recob::Shower>().at(0);
    auto shr2 = slice_pfp_v.at(i2).get<recob::Shower>().at(0);

    _pi0_shrscore1 = searchingfornues::GetTrackShowerScore(slice_pfp_v.at(i1));
    _pi0_shrscore2 = searchingfornues::GetTrackShowerScore(slice_pfp_v.at(i2));
    
    auto vtxcompat1 = VtxCompatibility(nuvtx, shr1->ShowerStart(), shr1->Direction());
    
    if (tkcalo_proxy!=NULL) {
      for (const searchingfornues::ProxyCaloElem_t& tk : *tkcalo_proxy) {
	// find track with ID matching the pfp index (this convention apparently works only for shower fits...)
	if (tk->ID()==int(slice_pfp_v[i1].index())) {


	  float shr_tkfit_start_x = tk->Start().X();
	  float shr_tkfit_start_y = tk->Start().Y();
	  float shr_tkfit_start_z = tk->Start().Z();
	  float shr_tkfit_start_sce[3];
	  searchingfornues::ApplySCECorrectionXYZ(shr_tkfit_start_x,shr_tkfit_start_y,shr_tkfit_start_z,shr_tkfit_start_sce);

	  auto const trkcalos = tk.get<anab::Calorimetry>();

	  for (const auto& tkcalo : trkcalos) {
	    
	    if (tkcalo->ResidualRange().size() == 0)
	      continue;
	    
	    auto const& xyz_v = tkcalo->XYZ();
	    
	    // collect XYZ coordinates of track-fitted shower
	    std::vector<float> x_v, y_v, z_v;
	    std::vector<float> dist_from_start_v;
	    for (auto xyz : xyz_v){
	      x_v.push_back(xyz.X());
	      y_v.push_back(xyz.Y());
	      z_v.push_back(xyz.Z());
	      float dist_from_start = searchingfornues::distance3d(xyz.X(), xyz.Y(), xyz.Z(),
								   shr_tkfit_start_sce[0],shr_tkfit_start_sce[1],shr_tkfit_start_sce[2]);

	      dist_from_start_v.push_back(dist_from_start);
	    }// collect XYZ coordinates of track-fitted shower
	    
	    std::vector<float> dqdx_values_corrected;
	    
	    if (fData || !fRecalibrateHits) {
	      if (!fLocaldEdx)
		dqdx_values_corrected = tkcalo->dQdx();
	      else
		dqdx_values_corrected = tkcalo->dEdx();
	    }// if re-calibration is not necessary
	    
	    else 
	      dqdx_values_corrected = llr_pid_calculator_shr.correct_many_hits_one_plane(tkcalo, tk, assocMCPart, fRecalibrateHits, fEnergyThresholdForMCHits, fLocaldEdx);
	    
	    
	    int nhits = 0;


	    
	    
	    if (tkcalo->PlaneID().Plane == 0) {
	      searchingfornues::GetTrackFitdEdx(dqdx_values_corrected,tkcalo->ResidualRange(), 0., 4., _pi0_dedx1_fit_U, nhits);
	      searchingfornues::GetTrackFitdEdx(tkcalo, 0., 4., false, _pi0_dedx1_fit_U, nhits);
	      _pi0_dedx1_fit_U = searchingfornues::GetdEdxfromdQdx(_pi0_dedx1_fit_U, shr1->ShowerStart()[0], shr1->ShowerStart()[2], shr1->ShowerStart()[2], 2.1, fADCtoE[0] );
	    }
	    if (tkcalo->PlaneID().Plane == 1) {
	      searchingfornues::GetTrackFitdEdx(dqdx_values_corrected,tkcalo->ResidualRange(), 0., 4., _pi0_dedx1_fit_V, nhits);
	      searchingfornues::GetTrackFitdEdx(tkcalo, 0., 4., false, _pi0_dedx1_fit_V, nhits);
	      _pi0_dedx1_fit_V = searchingfornues::GetdEdxfromdQdx(_pi0_dedx1_fit_V, shr1->ShowerStart()[0], shr1->ShowerStart()[2], shr1->ShowerStart()[2], 2.1, fADCtoE[1] );
	    }
	    if (tkcalo->PlaneID().Plane == 2) {
	      searchingfornues::GetTrackFitdEdx(dqdx_values_corrected,tkcalo->ResidualRange(), 0., 4., _pi0_dedx1_fit_Y, nhits);
	      searchingfornues::GetTrackFitdEdx(tkcalo, 0., 4., false, _pi0_dedx1_fit_Y, nhits);
	      _pi0_dedx1_fit_Y = searchingfornues::GetdEdxfromdQdx(_pi0_dedx1_fit_Y, shr1->ShowerStart()[0], shr1->ShowerStart()[2], shr1->ShowerStart()[2], 2.1, fADCtoE[2] );
	    }
	    
	  }// for all calorimetry objects
	}// if track matches shower index -> this is the track-fitted to the shower
      }// for all track fits to showers
    }// if track-fits to showers exist
    
    _pi0_radlen1    = vtxcompat1.second;
    _pi0_dot1       = vtxcompat1.first;
    _pi0_energy1_Y  = shr1->Energy()[2];
    _pi0_dedx1_Y    = shr1->dEdx()[2];
    _pi0_energy1_V  = shr1->Energy()[1];
    _pi0_dedx1_V    = shr1->dEdx()[1];
    _pi0_energy1_U  = shr1->Energy()[0];
    _pi0_dedx1_U    = shr1->dEdx()[0];
    
    auto vtxcompat2 = VtxCompatibility(nuvtx, shr2->ShowerStart(), shr2->Direction());
    
    if (tkcalo_proxy!=NULL) {
      for (const searchingfornues::ProxyCaloElem_t& tk : *tkcalo_proxy) {
	// find track with ID matching the pfp index (this convention apparently works only for shower fits...)
	if (tk->ID()==int(slice_pfp_v[i2].index())) {
	  
	  float shr_tkfit_start_x = tk->Start().X();
	  float shr_tkfit_start_y = tk->Start().Y();
	  float shr_tkfit_start_z = tk->Start().Z();
	  float shr_tkfit_start_sce[3];
	  searchingfornues::ApplySCECorrectionXYZ(shr_tkfit_start_x,shr_tkfit_start_y,shr_tkfit_start_z,shr_tkfit_start_sce);
	  
	  auto const trkcalos = tk.get<anab::Calorimetry>();
	  
	  for (const auto& tkcalo : trkcalos) {
	    
	    if (tkcalo->ResidualRange().size() == 0)
	      continue;
	    
	    auto const& xyz_v = tkcalo->XYZ();
	    
	    // collect XYZ coordinates of track-fitted shower
	    std::vector<float> x_v, y_v, z_v;
	    std::vector<float> dist_from_start_v;
	    for (auto xyz : xyz_v){
	      x_v.push_back(xyz.X());
	      y_v.push_back(xyz.Y());
	      z_v.push_back(xyz.Z());
	      float dist_from_start = searchingfornues::distance3d(xyz.X(), xyz.Y(), xyz.Z(),
								   shr_tkfit_start_sce[0],shr_tkfit_start_sce[1],shr_tkfit_start_sce[2]);
	      
	      dist_from_start_v.push_back(dist_from_start);
	    }// collect XYZ coordinates of track-fitted shower
	    
	    std::vector<float> dqdx_values_corrected;
	    
	    if (fData || !fRecalibrateHits) {
	      if (!fLocaldEdx)
		dqdx_values_corrected = tkcalo->dQdx();
	      else
		dqdx_values_corrected = tkcalo->dEdx();
	    }// if re-calibration is not necessary
	    
	    else 
	      dqdx_values_corrected = llr_pid_calculator_shr.correct_many_hits_one_plane(tkcalo, tk, assocMCPart, fRecalibrateHits, fEnergyThresholdForMCHits, fLocaldEdx);
	    
	    
	    int nhits = 0;

	    if (dqdx_values_corrected.size() > 0)
	      std::cout << "[DEDX] dqdx is " << dqdx_values_corrected[0] << std::endl;
	    
	    if (tkcalo->PlaneID().Plane == 0) {
	      searchingfornues::GetTrackFitdEdx(dqdx_values_corrected,tkcalo->ResidualRange(), 0., 4., _pi0_dedx2_fit_U, nhits);
	      std::cout << "[DEDX] 0 : " << _pi0_dedx2_fit_U << std::endl;
	      searchingfornues::GetTrackFitdEdx(tkcalo, 0., 4., false, _pi0_dedx2_fit_U, nhits);
	      std::cout << "[DEDX] 0 : " << _pi0_dedx2_fit_U << std::endl;
	      _pi0_dedx2_fit_U = searchingfornues::GetdEdxfromdQdx(_pi0_dedx2_fit_U, shr2->ShowerStart()[0], shr2->ShowerStart()[2], shr2->ShowerStart()[2], 2.1, fADCtoE[0] );
	    }
	    if (tkcalo->PlaneID().Plane == 1) {
	      searchingfornues::GetTrackFitdEdx(dqdx_values_corrected,tkcalo->ResidualRange(), 0., 4., _pi0_dedx2_fit_V, nhits);
	      std::cout << "[DEDX] 1 : " << _pi0_dedx2_fit_V << std::endl;
	      searchingfornues::GetTrackFitdEdx(tkcalo, 0., 4., false, _pi0_dedx2_fit_V, nhits);
	      std::cout << "[DEDX] 1 : " << _pi0_dedx2_fit_V << std::endl;
	      
	      _pi0_dedx2_fit_V = searchingfornues::GetdEdxfromdQdx(_pi0_dedx2_fit_V, shr2->ShowerStart()[0], shr2->ShowerStart()[2], shr2->ShowerStart()[2], 2.1, fADCtoE[1] );
	    }
	    if (tkcalo->PlaneID().Plane == 2) {
	      searchingfornues::GetTrackFitdEdx(dqdx_values_corrected,tkcalo->ResidualRange(), 0., 4., _pi0_dedx2_fit_Y, nhits);
	      std::cout << "[DEDX] 2 : " << _pi0_dedx2_fit_Y << std::endl;
	      searchingfornues::GetTrackFitdEdx(tkcalo, 0., 4., false, _pi0_dedx2_fit_Y, nhits);
	      std::cout << "[DEDX] 2 : " << _pi0_dedx2_fit_Y << std::endl;
	      _pi0_dedx2_fit_Y = searchingfornues::GetdEdxfromdQdx(_pi0_dedx2_fit_Y, shr2->ShowerStart()[0], shr2->ShowerStart()[2], shr2->ShowerStart()[2], 2.1, fADCtoE[2] );
	    }
	    
	  }// for all calorimetry objects
	  
	}// if track matches shower index -> this is the track-fitted to the shower
      }// for all track fits to showers
    }// if track-fits to showers exist
    
    _pi0_radlen2   = vtxcompat2.second;
    _pi0_dot2      = vtxcompat2.first;
    _pi0_energy2_Y = shr2->Energy()[2];
    _pi0_dedx2_Y   = shr2->dEdx()[2];
    _pi0_energy2_V = shr2->Energy()[1];
    _pi0_dedx2_V   = shr2->dEdx()[1];
    _pi0_energy2_U = shr2->Energy()[0];
    _pi0_dedx2_U   = shr2->dEdx()[0];
    
    _pi0_gammadot = shr1->Direction().Dot(shr2->Direction());
    _pi0_mass_Y = sqrt( 2 * _pi0_energy1_Y * _pi0_energy2_Y * (1 - _pi0_gammadot ) );
    _pi0_mass_V = sqrt( 2 * _pi0_energy1_V * _pi0_energy2_V * (1 - _pi0_gammadot ) );
    _pi0_mass_U = sqrt( 2 * _pi0_energy1_U * _pi0_energy2_U * (1 - _pi0_gammadot ) );

    // backtracking for shower 1
    float dot1 = shr1->Direction().Dot(_pi0_mcgamma0_mom);
    if (dot1 > _pi0_mcrcdot0) {
      _pi0_mcrcdot0 = dot1;
      _pi0_mcrce0   = shr1->Energy()[2];
    }// if new most aligned shower
    dot1 = shr1->Direction().Dot(_pi0_mcgamma1_mom);
    if (dot1 > _pi0_mcrcdot1) {
      _pi0_mcrcdot1 = dot1;
      _pi0_mcrce1   = shr1->Energy()[2];
    }// if new most aligned shower

    // backtracking for shower 1
    float dot2 = shr2->Direction().Dot(_pi0_mcgamma0_mom);
    if (dot2 > _pi0_mcrcdot0) {
      _pi0_mcrcdot0 = dot2;
      _pi0_mcrce0   = shr2->Energy()[2];
    }// if new most aligned shower
    dot2 = shr2->Direction().Dot(_pi0_mcgamma1_mom);
    if (dot2 > _pi0_mcrcdot1) {
      _pi0_mcrcdot1 = dot2;
      _pi0_mcrce1   = shr2->Energy()[2];
    }// if new most aligned shower
    
    return;
  }
  
  
  void Pi0Tagger::setBranches(TTree* _tree) {

    // truth
    _tree->Branch("pi0_mcgamma0_e" ,&_pi0_mcgamma0_e ,"pi0_mcgamma0_e/F" );
    _tree->Branch("pi0_mcgamma0_px",&_pi0_mcgamma0_px,"pi0_mcgamma0_px/F");
    _tree->Branch("pi0_mcgamma0_py",&_pi0_mcgamma0_py,"pi0_mcgamma0_py/F");
    _tree->Branch("pi0_mcgamma0_pz",&_pi0_mcgamma0_pz,"pi0_mcgamma0_pz/F");
    _tree->Branch("pi0_mcrcdot0",&_pi0_mcrcdot0,"pi0_mcrcdot0/F");
    _tree->Branch("pi0_mcrce0",&_pi0_mcrce0,"pi0_mcrce0/F");
    _tree->Branch("pi0_mcgamma1_e" ,&_pi0_mcgamma1_e ,"pi0_mcgamma1_e/F" );
    _tree->Branch("pi0_mcgamma1_px",&_pi0_mcgamma1_px,"pi0_mcgamma1_px/F");
    _tree->Branch("pi0_mcgamma1_py",&_pi0_mcgamma1_py,"pi0_mcgamma1_py/F");
    _tree->Branch("pi0_mcgamma1_pz",&_pi0_mcgamma1_pz,"pi0_mcgamma1_pz/F");
    _tree->Branch("pi0_mcrcdot1",&_pi0_mcrcdot1,"pi0_mcrcdot1/F");
    _tree->Branch("pi0_mcrce1",&_pi0_mcrce1,"pi0_mcrce1/F");

    // reco
    _tree->Branch("pi0_nshower",&_pi0_nshower,"pi0_nshower/I");
    _tree->Branch("pi0_ntrack" ,&_pi0_ntrack ,"pi0_ntrack/I" );
    _tree->Branch("pi0_ngamma",&_pi0_ngamma,"pi0_ngamma/I");
    _tree->Branch("pi0_radlen1",&_pi0_radlen1,"pi0_radlen1/F");
    _tree->Branch("pi0_radlen2",&_pi0_radlen2,"pi0_radlen2/F");
    _tree->Branch("pi0_dot1",&_pi0_dot1,"pi0_dot1/F");
    _tree->Branch("pi0_dot2",&_pi0_dot2,"pi0_dot2/F");
    _tree->Branch("pi0_energy1_Y",&_pi0_energy1_Y,"pi0_energy1_Y/F");
    _tree->Branch("pi0_energy2_Y",&_pi0_energy2_Y,"pi0_energy2_Y/F");
    _tree->Branch("pi0_dedx1_Y",&_pi0_dedx1_Y,"pi0_dedx1_Y/F");
    _tree->Branch("pi0_dedx2_Y",&_pi0_dedx2_Y,"pi0_dedx2_Y/F");
    _tree->Branch("pi0_dedx1_fit_Y",&_pi0_dedx1_fit_Y,"pi0_dedx1_fit_Y/F");
    _tree->Branch("pi0_dedx2_fit_Y",&_pi0_dedx2_fit_Y,"pi0_dedx2_fit_Y/F");
    _tree->Branch("pi0_energy1_V",&_pi0_energy1_V,"pi0_energy1_V/F");
    _tree->Branch("pi0_energy2_V",&_pi0_energy2_V,"pi0_energy2_V/F");
    _tree->Branch("pi0_dedx1_V",&_pi0_dedx1_V,"pi0_dedx1_V/F");
    _tree->Branch("pi0_dedx2_V",&_pi0_dedx2_V,"pi0_dedx2_V/F");
    _tree->Branch("pi0_dedx1_fit_V",&_pi0_dedx1_fit_V,"pi0_dedx1_fit_V/F");
    _tree->Branch("pi0_dedx2_fit_V",&_pi0_dedx2_fit_V,"pi0_dedx2_fit_V/F");
    _tree->Branch("pi0_energy1_U",&_pi0_energy1_U,"pi0_energy1_U/F");
    _tree->Branch("pi0_energy2_U",&_pi0_energy2_U,"pi0_energy2_U/F");
    _tree->Branch("pi0_dedx1_U",&_pi0_dedx1_U,"pi0_dedx1_U/F");
    _tree->Branch("pi0_dedx2_U",&_pi0_dedx2_U,"pi0_dedx2_U/F");
    _tree->Branch("pi0_dedx1_fit_U",&_pi0_dedx1_fit_U,"pi0_dedx1_fit_U/F");
    _tree->Branch("pi0_dedx2_fit_U",&_pi0_dedx2_fit_U,"pi0_dedx2_fit_U/F");
    _tree->Branch("pi0_shrscore1",&_pi0_shrscore1,"pi0_shrscore1/F");
    _tree->Branch("pi0_shrscore2",&_pi0_shrscore2,"pi0_shrscore2/F");
    _tree->Branch("pi0_gammadot",&_pi0_gammadot,"pi0_gammadot/F");
    _tree->Branch("pi0_mass_Y",&_pi0_mass_Y,"pi0_mass_Y/F");
    _tree->Branch("pi0_mass_V",&_pi0_mass_V,"pi0_mass_V/F");
    _tree->Branch("pi0_mass_U",&_pi0_mass_U,"pi0_mass_U/F");
    _tree->Branch("pi0_rc_vtx_x",&_pi0_rc_vtx_x,"pi0_rc_vtx_x/F");
    _tree->Branch("pi0_rc_vtx_y",&_pi0_rc_vtx_y,"pi0_rc_vtx_y/F");
    _tree->Branch("pi0_rc_vtx_z",&_pi0_rc_vtx_z,"pi0_rc_vtx_z/F");

    return;
  }

  void Pi0Tagger::resetTTree(TTree* _tree) {

    _pi0_nshower = std::numeric_limits<int>::min();
    _pi0_ntrack  = std::numeric_limits<int>::min();

    return;
  }

  std::pair<double,double> Pi0Tagger::VtxCompatibility(const TVector3& nuvtx, const TVector3& shrvtx, const TVector3& shrdir) {

    // grab shower start point and direction
    //auto shrvtx = shr->ShowerStart();
    //auto shrdir = shr->Direction();
    //gammadir = shrdir;
    
    // assess compatibility
    auto nuvtx2shrvtx = (shrvtx - nuvtx).Unit();
    auto shrdirnormed = shrdir.Unit();
    
    double dot  = nuvtx2shrvtx.Dot(shrdirnormed);
    double dist = (nuvtx-shrvtx).Mag();
    
    return std::make_pair(dot,dist);
    
    
  }// end of vertex compatibility

  template <typename T> float Pi0Tagger::PFPEnergy(const T& ass_clus_v) {

    float energy0 = 0;
    float energy1 = 0;
    float energy2 = 0;
    
    for (auto ass_clus : ass_clus_v) {

      //auto clus = clus_v[ass_clus.key()];
      
      if (ass_clus->View() == 0) { energy0 = ass_clus->Integral(); }
      if (ass_clus->View() == 1) { energy1 = ass_clus->Integral(); }
      if (ass_clus->View() == 2) { energy2 = ass_clus->Integral(); }

    }// for all clusters

    if (energy2 != 0) return energy2;
    return (energy1+energy0)/2.;
  }// calculte PFP energy based on associated hit charge

  void Pi0Tagger::ReadTruth(art::Event const& e) {

    auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");
    auto const& mcp_h = e.getValidHandle<std::vector<simb::MCParticle> >("largeant");
    
    _pi0_ispi0 = 0;
    
    auto mct = mct_h->at(0);
    size_t npart = mct.NParticles();
    
    int pi0trkid = -1;
    double pi0_vtx_x, pi0_vtx_y, pi0_vtx_z;
    
    for (size_t i=0; i < npart; i++){
      auto const& part = mct.GetParticle(i);
      if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
	pi0_vtx_x = part.Trajectory().X(0);
	pi0_vtx_y = part.Trajectory().Y(0);
	pi0_vtx_z = part.Trajectory().Z(0);
	pi0trkid = part.TrackId();
	break;
      }
    }
    
    if (pi0trkid < 0) return;
    
    // how many photons from the pi0 did we find?
    int ngamma = 0;
    
    // go through MCParticles and find photons
    for (size_t i=0; i < mcp_h->size(); i++) {
      
      auto mcp = mcp_h->at(i);
      
      if (mcp.PdgCode() != 22) continue;
      
      double x = mcp.Trajectory().X(0);
      double y = mcp.Trajectory().Y(0);
      double z = mcp.Trajectory().Z(0);
      double d = sqrt( ( (pi0_vtx_x - x) * (pi0_vtx_x - x) ) +
		       ( (pi0_vtx_y - y) * (pi0_vtx_y - y) ) +
		       ( (pi0_vtx_z - z) * (pi0_vtx_z - z) ) );
      
      if ( d < 0.01 ) {
	
	
	if (ngamma == 0) {
	  _pi0_mcgamma0_e = mcp.E(0);
	  _pi0_mcgamma0_mom = mcp.Momentum(0).Vect().Unit();
	  _pi0_mcgamma0_px = _pi0_mcgamma0_mom.X();
	  _pi0_mcgamma0_py = _pi0_mcgamma0_mom.Y();
	  _pi0_mcgamma0_pz = _pi0_mcgamma0_mom.Z();
	}
	
	if (ngamma == 1) {
	  _pi0_mcgamma1_e = mcp.E(0);
	  _pi0_mcgamma1_mom = mcp.Momentum(0).Vect().Unit();
	  _pi0_mcgamma1_px = _pi0_mcgamma1_mom.X();
	  _pi0_mcgamma1_py = _pi0_mcgamma1_mom.Y();
	  _pi0_mcgamma1_pz = _pi0_mcgamma1_mom.Z();
	}
	
	ngamma += 1; // update photon count
	
	
      }// if a pi0 induced gamma
      
    }// for all mcparticles
    

  return;
  }

  void Pi0Tagger::Reset() {

    _pi0_ispi0 = 0;
    _pi0_mcgamma0_e = 0;
    _pi0_mcgamma0_px = 0;
    _pi0_mcgamma0_py = 0;
    _pi0_mcgamma0_pz = 0;
    _pi0_mcgamma1_e = 0;
    _pi0_mcgamma1_px = 0;
    _pi0_mcgamma1_py = 0;
    _pi0_mcgamma1_pz = 0;
    _pi0_mcrce0 = 0;
    _pi0_mcrce1 = 0;
    _pi0_mcrcdot0 = -1;
    _pi0_mcrcdot1 = -1;

    _pi0_dot1 = -1;
    _pi0_dot2 = -1;
    _pi0_radlen1 = -1;
    _pi0_radlen2 = -1;
    _pi0_energy1_Y = -1;
    _pi0_energy2_Y = -1;
    _pi0_dedx1_Y = -1;
    _pi0_dedx2_Y = -1;
    _pi0_energy1_V = -1;
    _pi0_energy2_V = -1;
    _pi0_dedx1_V = -1;
    _pi0_dedx2_V = -1;
    _pi0_energy1_U = -1;
    _pi0_energy2_U = -1;
    _pi0_dedx1_U = -1;
    _pi0_dedx2_U = -1;
    _pi0_ngamma = 0;
    _pi0_ntrack = 0;
    _pi0_nshower = 0;
    _pi0_shrscore1 = -1;
    _pi0_shrscore2 = -1;
    _pi0_gammadot = -1;
    _pi0_mass_Y = -1;
    _pi0_mass_V = -1;
    _pi0_mass_U = -1;
    
    return;
  }// end of reset function

  
  DEFINE_ART_CLASS_TOOL(Pi0Tagger)
} // namespace selection

#endif
