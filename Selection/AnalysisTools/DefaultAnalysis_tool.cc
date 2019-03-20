#ifndef ANALYSIS_DEFAULTANALYSIS_CXX
#define ANALYSIS_DEFAULTANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"

#include "nusimdata/SimulationBase/MCTruth.h"

// backtracking tools
#include "ubana/ubana/searchingfornues/Selection/CommonDefs/BacktrackingFuncs.h"

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

  class DefaultAnalysis : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    DefaultAnalysis(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~DefaultAnalysis(){ };
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset);
    
    /**
     * @brief Analysis function
     */
    void analyzeEvent(art::Event const& e, bool fData) override;

    /**
     * @brief Analyze slice
     */
    void analyzeSlice(art::Event const& e, std::vector<ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) override;

    /**
     * @brief Save truth info for event associated to neutrino
     */
    void SaveTruth(art::Event const& e);
    
    /**
     * @brief set branches for TTree
     */
    void setBranches(TTree* _tree) override;

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree* _tree) override;
    
  private:

    art::InputTag fCRTVetoproducer; // producer for CRT veto ass tag [anab::T0 <-> recob::OpFlash]
    art::InputTag fCLSproducer; // cluster associated to PFP
    art::InputTag fMCTproducer;
    art::InputTag fBacktrackTag;

    // kinematic thresholds to define signal
    float fProtonThreshold;

    int   _run, _sub, _evt;       // event info
    // neutrino information
    float _nu_e;                  // neutrino energy [GeV]
    int   _nu_pdg;                // neutrino PDG code
    int   _ccnc;                  // CC or NC tag from GENIE
    float _vtx_x, _vtx_y, _vtx_z; // neutrino interaction vertex coordinates [cm]
    float _vtx_t;                 // neutrino generation time 
    bool  _isVtxInFiducial;       // true if neutrino in fiducial volume, 0 < x < 256 -116 < y < 116;  0 < z <  1036
    // final state particle information
    int   _nmuon;                    // is there a final-state muon from the neutrino? [1=yes 0=no]
    float _muon_e, _muon_p, _muon_c; // energy, purity, completeness.
    int   _nelec;                    // is there a final-state electron from the neutrino? [1=yes 0=no]
    float _elec_e, _elec_p, _elec_c; // energy, purity, completeness.
    int   _npi0;                     // how many pi0s are there?
    int   _pi0;                      // is there a final-state pi0 from the neutrino? [1=yes 0=no]
    float _pi0_e, _pi0_p, _pi0_c;    // energy, purity, completeness.
    int   _nproton;                  // how many protons are there?
    int   _proton;                       // is there a final-state proton from the neutrino? [1=yes 0=no]
    float _proton_e, _proton_p, _proton_c;       // energy, purity, completeness.
    int   _npion;                    // how many pions are there?
    int   _pion;                     // is there a final-state charged pion from the neutrino? [1=yes 0=no]
    float _pion_e, _pion_p, _pion_c; // energy, purity, completeness.

    // number of slices in the event
    int _nslice;
    int _crtveto; // is the event vetoed by the CRT Veto?
    float _crthitpe; // pe associated to CRT hit

    // reco PFParticle backtracking. One entry for PFParticle in the slice
    std::vector<int>   _backtracked_idx;    // index of PFP [key]
    std::vector<int>   _backtracked_tid;    // TrackID of backtracked MCParticle
    std::vector<int>   _backtracked_pdg;    // PDG code of backtracked particle
    std::vector<float> _backtracked_e;      // energy of backtracked particle
    std::vector<float> _backtracked_purity; // purity of backtracking

    float _lep_e;                 // lepton energy (if one exists) [GeV]
    int  _pass;                   // does the slice pass the selection
    float _xtimeoffset, _xsceoffset, _ysceoffset, _zsceoffset; // offsets for generation time and SCE
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  DefaultAnalysis::DefaultAnalysis(const fhicl::ParameterSet& p)
  {
    fCRTVetoproducer = p.get< art::InputTag > ("CRTVetoproducer",""); // default is no CRT veto
    fCLSproducer = p.get< art::InputTag > ("CLSproducer");
    fMCTproducer = p.get< art::InputTag > ("MCTproducer");
    fBacktrackTag = p.get< art::InputTag > ("BacktrackTag");
    // kinematic thresholds for defining signal
    fProtonThreshold = p.get< float > ("ProtonThreshold");
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void DefaultAnalysis::configure(fhicl::ParameterSet const & p)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void DefaultAnalysis::analyzeEvent(art::Event const& e, bool fData)
  {
    _evt = e.event();
    _sub = e.subRun();
    _run = e.run();

    if (!fData) {
      // SaveTruth
      SaveTruth(e);
    }

    // Grab CRT veto information if available - CRT should probably have its own tool?
    if (fCRTVetoproducer != "") {
      art::Handle< art::Assns<crt::CRTHit,recob::OpFlash,void> > crtveto_h;
      e.getByLabel(fCRTVetoproducer,crtveto_h);
      _crtveto = crtveto_h->size();
      if (_crtveto == 1) 
	_crthitpe = crtveto_h->at(0).first->peshit;
    }// if the CRT veto label has been defined

    return;
  }

  void DefaultAnalysis::analyzeSlice(art::Event const& e, std::vector<ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected)
  {

    ProxyClusColl_t const& clus_proxy = proxy::getCollection<std::vector<recob::Cluster> >(e, fCLSproducer,
											   proxy::withAssociated<recob::Hit>(fCLSproducer));

    // load backtrack information
    art::InputTag BacktrackTag { fBacktrackTag };
    auto const& gaushit_h = e.getValidHandle<std::vector<recob::Hit> > ("gaushit");
    art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> hittruth(gaushit_h,e,BacktrackTag);
    //const auto& hittruth = std::unique_ptr<art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> > (new art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>(gaushit_h,e,BacktrackTag));

    for (auto pfp :  slice_pfp_v) {
      // backtrack PFParticles in the slice
      if (!fData) {
	// get h  its associated to this PFParticle through the clusters
	std::vector<art::Ptr<recob::Hit> > hit_v;
	auto clus_pxy_v = pfp.get<recob::Cluster>();
	if (clus_pxy_v.size() != 0) {
	  for (auto ass_clus : clus_pxy_v) {
	    // get cluster proxy
	    const auto& clus = clus_proxy[ass_clus.key()];
	    auto clus_hit_v = clus.get<recob::Hit>();
	    for (const auto& hit : clus_hit_v)
	      hit_v.push_back(hit);
	  }// for all clusters associated to PFP
	  
	  float purity, completeness;
	  auto mcp = searchingfornues::getAssocMCParticle(hittruth, hit_v, purity, completeness);
	  if (mcp){
	    auto PDG = mcp->PdgCode();
	    _backtracked_idx.push_back(0);
	    _backtracked_tid.push_back(mcp->TrackId());
	    _backtracked_e.push_back(mcp->Momentum().E());
	    _backtracked_pdg.push_back(PDG);
	    _backtracked_purity.push_back(purity);
	    // if this is n interesting particle, save it to the TTree
	    if (fabs(PDG) == 13) {
	      if ( fabs(mcp->Momentum().E() - _muon_e) < 0.01 )
		_muon_p = purity;
	    }
	    if (fabs(PDG) == 11) {
	      if ( fabs(mcp->Momentum().E() - _elec_e) < 0.01 )
		_elec_p = purity;
	    }
	    if (fabs(PDG) == 2212) {
	      if ( fabs(mcp->Momentum().E() - _proton_e) < 0.01 )
		_proton_p = purity;
	    }
	  }
	  else{
	    _backtracked_idx.push_back(0);
	    _backtracked_tid.push_back(0);
	    _backtracked_e.push_back(0);
	    _backtracked_pdg.push_back(0);
	    _backtracked_purity.push_back(0.);
	  }
	}// if there are associated clusters
      }// if MC
    }
    _nslice += 1;
    if (selected) _pass = 1;
  }
  
  void DefaultAnalysis::setBranches(TTree* _tree) 
  {
    // neutrino information
    _tree->Branch("_nu_pdg",&_nu_pdg,"nu_pdg/I");
    _tree->Branch("_ccnc"  ,&_ccnc  ,"ccnc/I"  );
    _tree->Branch("_nu_e"  ,&_nu_e  ,"nu_e/F"  );
    _tree->Branch("_vtx_x" ,&_vtx_x ,"vtx_x/F" );
    _tree->Branch("_vtx_y" ,&_vtx_y ,"vtx_y/F" );
    _tree->Branch("_vtx_z" ,&_vtx_z ,"vtx_z/F" );
    _tree->Branch("_isVtxInFiducial" ,&_isVtxInFiducial ,"isVtxInFiducial/O" );
    // individual particles in the neutrino slice
    // legend:
    // _e -> energy of particle in GeV
    // _c -> completeness from back-tracking [0,1]
    // _p -> purity from back-tracking [0,1]
    // muon
    _tree->Branch("_nmuon",&_nmuon,"nmuon/I");
    _tree->Branch("_muon_e",&_muon_e,"muon_e/F");
    _tree->Branch("_muon_c",&_muon_c,"muon_c/F");
    _tree->Branch("_muon_p",&_muon_p,"muon_p/F");
    // electron
    _tree->Branch("_nelec",&_nelec,"nelec/I");
    _tree->Branch("_elec_e",&_elec_e,"elec_e/F");
    _tree->Branch("_elec_c",&_elec_c,"elec_c/F");
    _tree->Branch("_elec_p",&_elec_p,"elec_p/F");
    // pi0
    _tree->Branch("_pi0"   ,&_pi0   ,"pi0/I");
    _tree->Branch("_pi0_e",&_pi0_e,"pi0_e/F");
    _tree->Branch("_pi0_c",&_pi0_c,"pi0_c/F");
    _tree->Branch("_pi0_p",&_pi0_p,"pi0_p/F");
    // first [highest momentum] proton
    _tree->Branch("_proton"   ,&_proton   ,"proton/I");
    _tree->Branch("_proton_e",&_proton_e,"proton_e/F");
    _tree->Branch("_proton_c",&_proton_c,"proton_c/F");
    _tree->Branch("_proton_p",&_proton_p,"proton_p/F");

    _tree->Branch("_nslice"  ,&_nslice  ,"nslice/I"  );
    _tree->Branch("_crtveto" ,&_crtveto ,"crtveto/I" );
    _tree->Branch("_crthitpe",&_crthitpe,"crthitpe/F");

    // PFParticle backtracking
    _tree->Branch("_backtracked_idx"   ,"std::vector<int>"  ,&_backtracked_idx   );
    _tree->Branch("_backtracked_tid"   ,"std::vector<int>"  ,&_backtracked_tid   );
    _tree->Branch("_backtracked_pdg"   ,"std::vector<int>"  ,&_backtracked_pdg   );
    _tree->Branch("_backtracked_e"     ,"std::vector<float>",&_backtracked_e     );
    _tree->Branch("_backtracked_purity","std::vector<float>",&_backtracked_purity);

    _tree->Branch("_lep_e" ,&_lep_e ,"lep_e/F" );
    _tree->Branch("_pass"  ,&_pass  ,"pass/I"  );
    _tree->Branch("_run"   ,&_run   ,"run/I"   );
    _tree->Branch("_sub"   ,&_sub   ,"sub/I"   );
    _tree->Branch("_evt"   ,&_evt   ,"evt/I"   );

    _tree->Branch("_xtimeoffset",&_xtimeoffset,"xtimeoffset/F");
    _tree->Branch("_xsceoffset" ,&_xsceoffset ,"xsceoffset/F" );
    _tree->Branch("_ysceoffset" ,&_ysceoffset ,"ysceoffset/F" );
    _tree->Branch("_zsceoffset" ,&_zsceoffset ,"zsceoffset/F" );
  }

  void DefaultAnalysis::resetTTree(TTree* _tree)
  {
    _run     = std::numeric_limits<int>::min();
    _sub     = std::numeric_limits<int>::min();
    _evt     = std::numeric_limits<int>::min();
    _nu_e    = std::numeric_limits<float>::min();
    _nu_pdg  = std::numeric_limits<int>::min();
    _ccnc    = std::numeric_limits<int>::min();
    _pass    = 0;
    _vtx_x   = std::numeric_limits<float>::min();
    _vtx_y   = std::numeric_limits<float>::min();
    _vtx_z   = std::numeric_limits<float>::min();
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

    _backtracked_idx.clear();
    _backtracked_tid.clear();
    _backtracked_e.clear();
    _backtracked_pdg.clear();
    _backtracked_purity.clear();
  }

void DefaultAnalysis::SaveTruth(art::Event const& e) {

  // load MCTruth
  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >(fMCTproducer);

  auto mct      = mct_h->at(0);
  auto neutrino = mct.GetNeutrino();
  auto nu       = neutrino.Nu();

  _ccnc = neutrino.CCNC();
  _nu_pdg = nu.PdgCode();
  _nu_e  = nu.Trajectory().E(0);
  _vtx_x = nu.EndX();
  _vtx_y = nu.EndY();
  _vtx_z = nu.EndZ();
  _vtx_t = nu.T();

  if (_vtx_x <  256. && _vtx_x >    0. && 
      _vtx_y < -116. && _vtx_y >  116. && 
      _vtx_z <    0. && _vtx_z > 1036.   ) {_isVtxInFiducial = true;}
  else _isVtxInFiducial = false;


  _nelec   = 0;
  _nmuon   = 0;
  _npi0    = 0;
  _nproton = 0;
  _npion   = 0;

  size_t npart = mct.NParticles();
  for (size_t i=0; i < npart; i++){

    auto const& part = mct.GetParticle(i);

    // if muon
    if ( (fabs(part.PdgCode()) == 13) and (part.StatusCode() == 1) ){
      _nmuon += 1;
      _muon_e = part.Momentum(0).E();
    }// if muon
    // if electron
    if ( (fabs(part.PdgCode()) == 11) and (part.StatusCode() == 1) ){
      _nelec += 1;
      _elec_e = part.Momentum(0).E();
    }// if electron
    // if pi0
    if ( (part.PdgCode() == 111) and (part.StatusCode() == 1) ){
      _npi0 += 1;
      _pi0_e = part.Momentum(0).E();
    }// if pi0
    // if proton
    if ( (part.PdgCode() == 2212) and (part.StatusCode() == 1) ){
      if (_nproton == 0) {
	// if highest energy, update energy
	if (part.Momentum(0).E() > _proton_e)
	  _proton_e = part.Momentum(0).E();
	if (part.Momentum(0).E() > fProtonThreshold)
	  _nproton += 1;
      }
    }// if proton
    // if pion
    if ( (fabs(part.PdgCode()) == 211) and (part.StatusCode() == 1) ){
      if (part.Momentum(0).E() > _pion_e)
	_pion_e = part.Momentum(0).E();
      _npion += 1;
    }// if proton
  }// for all MCParticles

  searchingfornues::ApplyDetectorOffsets(_vtx_t,_vtx_x,_vtx_y,_vtx_z,_xtimeoffset,_xsceoffset,_ysceoffset,_zsceoffset);
  
  return;
}

  DEFINE_ART_CLASS_TOOL(DefaultAnalysis)
} // namespace analysis

#endif
