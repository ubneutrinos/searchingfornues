////////////////////////////////////////////////////////////////////////
// Class:       NeutrinoSelectionFilter
// Plugin Type: analyzer (art v3_00_00)
// File:        NeutrinoSelectionFilter_module.cc
//
// Generated on March 7 2019 by Elena Gramellini copy-pasting and modifiying
// NeutrinoSelection, analyzer 
// Generated at Wed Jan 30 21:51:32 2019 by David Caratelli using cetskelgen
// from cetlib version v3_04_00.
// Scope of the code is to filter neutrino events based on the NeutrinoSelection
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

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardata/Utilities/FindManyInChainP.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

// backtracking tools
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

// services for detector properties
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// selection tool
#include "SelectionTools/SelectionToolBase.h"

// saving output
#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"
#include "TVector3.h"

class NeutrinoSelectionFilter;


class NeutrinoSelectionFilter : public art::EDFilter {
public:
  explicit NeutrinoSelectionFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NeutrinoSelectionFilter(NeutrinoSelectionFilter const&) = delete;
  NeutrinoSelectionFilter(NeutrinoSelectionFilter&&) = delete;
  NeutrinoSelectionFilter& operator=(NeutrinoSelectionFilter const&) = delete;
  NeutrinoSelectionFilter& operator=(NeutrinoSelectionFilter&&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  bool endSubRun(art::SubRun &subrun);

  using ProxyPfpColl_t = selection::ProxyPfpColl_t;
  using ProxyPfpElem_t =  selection::ProxyPfpElem_t;
  //using AssocPfpMeta_t = decltype(proxy::CollectionProxyBase< ProxyPfpElem_t, ProxyPfpColl_t, larpandoraobj::PFParticleMetadata >::aux());

private:

  art::InputTag fPFPproducer;
  art::InputTag fCLSproducer; // cluster associated to PFP
  art::InputTag fHITproducer; // hit associated to cluster
  art::InputTag fSHRproducer; // shower associated to PFP
  art::InputTag fVTXproducer; // vertex associated to PFP
  art::InputTag fTRKproducer; // track associated to PFP
  art::InputTag fMCTproducer;
  art::InputTag fBacktrackTag;
  bool fVerbose;
  bool fData;

  // TTree
  TTree* _tree;
  int   _run, _sub, _evt;       // event info
  // neutrino information
  float _nu_e;                  // neutrino energy [GeV]
  int   _nu_pdg;                // neutrino PDG code
  int   _ccnc;                  // CC or NC tag from GENIE
  float _vtx_x, _vtx_y, _vtx_z; // neutrino interaction vertex coordinates [cm]
  float _vtx_t;                 // neutrino generation time 
  // final state particle information
  int   _nmuon;                    // is there a final-state muon from the neutrino? [1=yes 0=no]
  float _muon_e, _muon_p, _muon_c; // energy, purity, completeness.
  int   _nelec;                    // is there a final-state electron from the neutrino? [1=yes 0=no]
  float _elec_e, _elec_p, _elec_c; // energy, purity, completeness.
  int   _npi0;                     // how many pi0s are there?
  int   _pi0;                      // is there a final-state pi0 from the neutrino? [1=yes 0=no]
  float _pi0_e, _pi0_p, _pi0_c;    // energy, purity, completeness.
  int   _nproton;                  // how many protons are there?
  int   _p0;                       // is there a final-state proton from the neutrino? [1=yes 0=no]
  float _p0_e, _p0_p, _p0_c;       // energy, purity, completeness.
  int   _p1;                       // is there a final-state proton from the neutrino? [1=yes 0=no]
  float _p1_e, _p1_p, _p1_c;       // energy, purity, completeness.
  int   _p2;                       // is there a final-state proton from the neutrino? [1=yes 0=no]
  float _p2_e, _p2_p, _p2_c;       // energy, purity, completeness.  
  int   _npion;                    // how many pions are there?
  int   _pion;                     // is there a final-state charged pion from the neutrino? [1=yes 0=no]
  float _pion_e, _pion_p, _pion_c; // energy, purity, completeness.

  // reco PFParticle backtracking. One entry for PFParticle in the slice
  std::vector<int>   _backtracked_idx;    // index of PFP [key]
  std::vector<int>   _backtracked_tid;    // TrackID of backtracked MCParticle
  std::vector<int>   _backtracked_pdg;    // PDG code of backtracked particle
  std::vector<float> _backtracked_e;      // energy of backtracked particle
  std::vector<float> _backtracked_purity; // purity of backtracking

  float _lep_e;                 // lepton energy (if one exists) [GeV]
  int  _pass;                   // does the slice pass the selection
  float _xtimeoffset, _xsceoffset, _ysceoffset, _zsceoffset; // offsets for generation time and SCE

  TTree* _subrun_tree;
  int _run_sr;                  // The run number
  int _sub_sr;                  // The subRun number
  float _pot;                   // The total amount of POT for the current sub run
  

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;

  // selection tool
  std::unique_ptr<::selection::SelectionToolBase> _selectionTool;

  /**
   * @brief function to builf a map linking PFParticle index to Self() attribute
   *
   * @input handle to event pfparticle record
   */
  void BuildPFPMap(const ProxyPfpColl_t& pfp_pxy_col);

  
  /**
   * @brief print PFParticle metadata information
   */
  template <typename T> void printPFParticleMetadata(const ProxyPfpElem_t& pfp_pxy,
						     const T& pfParticleMetadataList);

  /**
   * @brief build PFParticle hierarchy (i.e. slice) from parent [recursive function]
   *
   * @input pfp_pxy : parent pfparticle proxy for which to add daughters
   * @input pfp_pxy_col : evnt PFP proxy collection
   * @input slice_v : passed by reference, slice containing all PFParticles in hierarchy
   *
   */
  void AddDaughters(const ProxyPfpElem_t& pfp_pxy,
		    const ProxyPfpColl_t& pfp_pxy_col,
		    std::vector<ProxyPfpElem_t>& slice_v);

  /**
   * @brief backtracking function
   * backtrack_h stores mcparticle to hit association info
   * trackid_v stores the trackid of particles to backtrack
   * hit_idx_v stores the hit to be backtracked [i.e. from the PFParticle]
   * completeness is the particle completeness
   * purity is the particle purity
   */
  /*
  size_t BackTrack(const art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> backtrack_h,
		   const std::vector<unsigned int>& trackid_v,
		   const std::vector<unsigned int>& hit_idx_v, float& completeness, float& purity);
  */

  /**
   * @brief BackTrack a single hit collection (i.e. from a PFParticle)
   */
  art::Ptr<simb::MCParticle> getAssocMCParticle(art::Event const& e,
						const std::vector<art::Ptr<recob::Hit> >& hits,
						float& purity, float& completeness);
    
  /**
   * @brief Save truth info for event associated to neutrino
   */
  void SaveTruth(art::Event const& e);

  /**
   * @brief shift coordinates for truth particles according to SCE offsets + time offsets
   */
  void ApplyDetectorOffsets();

  /**
   * @brief reset ttree
   */
  void ResetTTree();

};


NeutrinoSelectionFilter::NeutrinoSelectionFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  
  fPFPproducer = p.get< art::InputTag > ("PFPproducer");
  fSHRproducer = p.get< art::InputTag > ("SHRproducer");
  fHITproducer = p.get< art::InputTag > ("HITproducer");
  fCLSproducer = p.get< art::InputTag > ("CLSproducer");
  fVTXproducer = p.get< art::InputTag > ("VTXproducer");
  fTRKproducer = p.get< art::InputTag > ("TRKproducer");
  fMCTproducer = p.get< art::InputTag > ("MCTproducer");
  fBacktrackTag = p.get< art::InputTag > ("BacktrackTag");
  fVerbose     = p.get< bool >          ("Verbose");
  fData     = p.get< bool >          ("IsData");

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("NeutrinoSelectionFilter", "Neutrino Selection TTree");

  // neutrino information
  _tree->Branch("_nu_pdg",&_nu_pdg,"nu_pdg/I");
  _tree->Branch("_ccnc"  ,&_ccnc  ,"ccnc/I"  );
  _tree->Branch("_nu_e"  ,&_nu_e  ,"nu_e/F"  );
  _tree->Branch("_vtx_x" ,&_vtx_x ,"vtx_x/F" );
  _tree->Branch("_vtx_y" ,&_vtx_y ,"vtx_y/F" );
  _tree->Branch("_vtx_z" ,&_vtx_z ,"vtx_z/F" );
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
  _tree->Branch("_p0"   ,&_p0   ,"p0/I");
  _tree->Branch("_p0_e",&_p0_e,"p0_e/F");
  _tree->Branch("_p0_c",&_p0_c,"p0_c/F");
  _tree->Branch("_p0_p",&_p0_p,"p0_p/F");
  // second proton
  _tree->Branch("_p1"   ,&_p1   ,"p1/I");
  _tree->Branch("_p1_e",&_p1_e,"p1_e/F");
  _tree->Branch("_p1_c",&_p1_c,"p1_c/F");
  _tree->Branch("_p1_p",&_p1_p,"p1_p/F");
  // third proton
  _tree->Branch("_p2"   ,&_p2   ,"p2/I");
  _tree->Branch("_p2_e",&_p2_e,"p2_e/F");
  _tree->Branch("_p2_c",&_p2_c,"p2_c/F");
  _tree->Branch("_p2_p",&_p2_p,"p2_p/F");

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

  _subrun_tree = tfs->make<TTree>("SubRun", "SubRun TTree");
  _subrun_tree->Branch("run"   , &_run_sr   , "run/I");
  _subrun_tree->Branch("subRun", &_sub_sr   , "subRun/I");
  if (!fData)
        _subrun_tree->Branch("pot", &_pot, "pot/F");


  // configure and construct Selection Tool
  const fhicl::ParameterSet& selection_pset = p.get<fhicl::ParameterSet>("SelectionTool");  
  _selectionTool = art::make_tool<::selection::SelectionToolBase>(selection_pset);

  // pass the TTree to the selection tool so that any branch can be added to it
  _selectionTool->setBranches(_tree);
  
}

bool NeutrinoSelectionFilter::filter(art::Event & e)
{

  ResetTTree();
  
  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();

  if (!fData) {
    // SaveTruth
    SaveTruth(e);
  }

  
  if (fVerbose) { std::cout << "new event : [run,event] : [" << e.run() << ", " << e.event() << "]" << std::endl; }
  
  // grab PFParticles in event
  selection::ProxyPfpColl_t const& pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle> >(e,fPFPproducer,
												     proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
												     proxy::withAssociated<recob::Cluster>(fCLSproducer),
												     proxy::withAssociated<recob::Track>(fTRKproducer),
												     proxy::withAssociated<recob::Vertex>(fVTXproducer),
												     proxy::withAssociated<recob::Shower>(fSHRproducer));
  
  selection::ProxyClusColl_t const& clus_proxy = proxy::getCollection<std::vector<recob::Cluster> >(e, fCLSproducer,
												    proxy::withAssociated<recob::Hit>(fCLSproducer));
  
  
  // build PFParticle map  for this event
  BuildPFPMap(pfp_proxy);
  
  // loop through PFParticles
  size_t p=0;
  // Should you keep this event? Yes, but only if there's a neutrino.
  // On the other hand, you can fill the tree even if you discart the event
  bool keepEvent = false;

  for (const ProxyPfpElem_t& pfp_pxy : pfp_proxy) {

    // get metadata for this PFP
    const auto& pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();
    
    //  find neutrino candidate
    if (pfp_pxy->IsPrimary() == false) continue;
    
    auto PDG = fabs(pfp_pxy->PdgCode());

    if ( (PDG == 12) || (PDG == 14) ) {
      
      if (fVerbose) printPFParticleMetadata(pfp_pxy,pfParticleMetadataList);
      
      // collect PFParticle hierarchy originating from this neutrino candidate
      std::vector<ProxyPfpElem_t> slice_pfp_v;
      AddDaughters(pfp_pxy, pfp_proxy, slice_pfp_v);

      if (fVerbose) { std::cout << "This slice has " << slice_pfp_v.size() << " daughter PFParticles" << std::endl; }
      
      // create list of tracks and showers associated to this slice
      std::vector<art::Ptr<recob::Track>  > sliceTracks;
      std::vector<art::Ptr<recob::Shower> > sliceShowers;

      for (auto pfp :  slice_pfp_v) {

	// backtrack PFParticles in the slice
	if (!fData) {
	  // get hits associated to this PFParticle through the clusters
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
	    auto mcp = getAssocMCParticle(e, hit_v, purity, completeness);
	    if (mcp){
	      auto PDG = mcp->PdgCode();
	      _backtracked_idx.push_back(0);
	      _backtracked_tid.push_back(mcp->TrackId());
	      _backtracked_e.push_back(mcp->Momentum().E());
	      _backtracked_pdg.push_back(PDG);
	      _backtracked_purity.push_back(purity);
	      // if this is n interesting particle, save it to the TTree
	      if (fabs(PDG) == 13) {
		_muon_p = purity;
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

	// if there is a track associated to the PFParticle, add it
	auto const& ass_trk_v = pfp.get<recob::Track>();
	if (ass_trk_v.size() == 1) sliceTracks.push_back( ass_trk_v.at(0) );
	// if there is a shower associated to the PFParticle, add it
	auto const& ass_shr_v = pfp.get<recob::Shower>();
	if (ass_shr_v.size() == 1) sliceShowers.push_back( ass_shr_v.at(0) );
      }// for all PFParticles in the slice
      
      // check that # of PFP and # trk + shr matches
      if ( (slice_pfp_v.size()-1) != (sliceTracks.size() + sliceShowers.size()) )
	std::cout << "ERROR : there are "  << slice_pfp_v.size() << " PFP but " << sliceTracks.size() << " + " << sliceShowers.size() << " tracks + showers" << std::endl;
      
      // run selection on this slice
      // bool selected = _selectionTool->selectEvent(e, sliceTracks, sliceShowers);
      bool selected = _selectionTool->selectEvent(e, slice_pfp_v);
      if (fVerbose && selected) { std::cout << "SLICE was selected!" << std::endl; }
      
      if (selected) { _pass = 1; keepEvent = true; }
      
      _tree->Fill();


    }// if a neutrino PFParticle
    p++;
  }// for all PFParticles
  
  return keepEvent;
}

template <typename T> void NeutrinoSelectionFilter::printPFParticleMetadata(const ProxyPfpElem_t& pfp_pxy,
									    const T& pfParticleMetadataList) {
  
  if (pfParticleMetadataList.size()!=0) {
    
    for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j)
      {
	const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
	auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
	if (!pfParticlePropertiesMap.empty()) {
	  std::cout << " Found PFParticle " << pfp_pxy->Self() << " with: " << std::endl;
	  for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	    std::cout << "  - " << it->first << " = " << it->second << std::endl;
	  }
	}
      }
  }// if PFP metadata exists!
  return;
}

void NeutrinoSelectionFilter::BuildPFPMap(const ProxyPfpColl_t& pfp_pxy_col) {
  
  _pfpmap.clear();

  unsigned int p=0;
  for (const auto& pfp_pxy : pfp_pxy_col) {
    _pfpmap[pfp_pxy->Self()] = p;
    p++;
  }

  return;
}// BuildPFPMap

void NeutrinoSelectionFilter:: AddDaughters(const ProxyPfpElem_t& pfp_pxy,
				      const ProxyPfpColl_t& pfp_pxy_col,
				      std::vector<ProxyPfpElem_t>& slice_v) {
  
  auto daughters = pfp_pxy->Daughters();
  
  slice_v.push_back(pfp_pxy);
  
  if (fVerbose) std::cout << "\t PFP w/ PdgCode " << pfp_pxy->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;
  
  for(auto const& daughterid : daughters) {

    if (_pfpmap.find(daughterid) == _pfpmap.end()) {
      std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
      continue;
    }
    
    // const art::Ptr<recob::PFParticle> pfp_pxy(pfp_pxy_col, _pfpmap.at(daughterid) );
    auto pfp_pxy2 = pfp_pxy_col.begin();
    for (size_t j=0; j<_pfpmap.at(daughterid); ++j) ++pfp_pxy2;
    // const T& pfp_pxy2 = (pfp_pxy_col.begin()+_pfpmap.at(daughterid));
    
    AddDaughters(*pfp_pxy2, pfp_pxy_col, slice_v);
    
  }// for all daughters
  
  return;
}// AddDaughters

void NeutrinoSelectionFilter::ResetTTree() {

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
  _p0_e = 0;
  _p0_p = 0;
  _p0_c = 0;

  _p1_e = 0;
  _p1_p = 0;
  _p1_c = 0;

  _p2_e = 0;
  _p2_p = 0;
  _p2_c = 0;

  _backtracked_idx.clear();
  _backtracked_tid.clear();
  _backtracked_e.clear();
  _backtracked_pdg.clear();
  _backtracked_purity.clear();

  _selectionTool->resetTTree(_tree);

}

void NeutrinoSelectionFilter::SaveTruth(art::Event const& e) {

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
	_p0_e = part.Momentum(0).E();
      }
      if (_nproton == 1) {
	_p1_e = part.Momentum(0).E();
      }
      if (_nproton == 2) {
	_p2_e = part.Momentum(0).E();
      }
      _nproton += 1;
    }// if proton
    // if pion
    if ( (fabs(part.PdgCode()) == 211) and (part.StatusCode() == 1) ){
      _pion_e = part.Momentum(0).E();
      _npion += 1;
    }// if proton
  }// for all MCParticles

  ApplyDetectorOffsets();
  
  return;
}

// calculate offsets for truth neutrino vertex
void NeutrinoSelectionFilter::ApplyDetectorOffsets() {

  auto const& detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
  auto const& detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  double g4Ticks = detClocks->TPCG4Time2Tick(_vtx_t) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
  _xtimeoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);

  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto offset = SCE->GetPosOffsets(geo::Point_t(_vtx_x,_vtx_y,_vtx_z));
  _xsceoffset = offset.X();
  _ysceoffset = offset.Y();
  _zsceoffset = offset.Z();

}// apply SCE corrections + time-shifts
 
art::Ptr<simb::MCParticle> NeutrinoSelectionFilter::getAssocMCParticle(art::Event const& e,
								 const std::vector<art::Ptr<recob::Hit> >& hits,
								 float& purity, float& completeness) {

  // load backtrack information
  art::InputTag BacktrackTag { fBacktrackTag };
  auto const& gaushit_h = e.getValidHandle<std::vector<recob::Hit> > ("gaushit");
  art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> hittruth(gaushit_h,e,BacktrackTag);
  //const auto& hittruth = std::unique_ptr<art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> > (new art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>(gaushit_h,e,BacktrackTag));
  
  // store total charge from hits
  float pfpcharge = 0; // total hit charge from clusters
  float maxcharge = 0; // charge backtracked to best match
  
  //credit: Wes Ketchum
  std::unordered_map<int,double> trkide;
  std::unordered_map<int,float> trkq;
  double maxe=-1, tote=0;
  art::Ptr<simb::MCParticle> maxp_me; //pointer for the particle match we will calculate
  //simb::MCParticle* maxp_me; //pointer for the particle match we will calculate
  for (auto h : hits) {
    pfpcharge += h->Integral();
    //const std::vector<const simb::MCParticle*> particle_vec = hittruth.at(h.key());
    std::vector<art::Ptr<simb::MCParticle> > particle_vec = hittruth.at(h.key());
    //auto particle_vec = hittruth.at(h.key());
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec = hittruth.data(h.key());;
    //loop over particles
    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
      trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id
      trkq  [ particle_vec[i_p]->TrackId() ] += h->Integral()  * match_vec[i_p]->ideFraction; //store hit integral associated to this hit
      tote += match_vec[i_p]->energy; //calculate total energy deposited
      if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
	maxe = trkide[ particle_vec[i_p]->TrackId() ];
	maxp_me = particle_vec[i_p];
	maxcharge = trkq[ particle_vec[i_p]->TrackId() ];
      }
    }//end loop over particles per hit
  }

  purity       = maxcharge / pfpcharge;
  completeness = 0;

  return maxp_me;
}

 /*
 size_t NeutrinoSelectionFilter::BackTrack(const art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> backtrack_h,
				     const std::vector<unsigned int>& trackid_v,
				     const std::vector<unsigned int>& hit_idx_v, float& completeness, float& purity)
 {
   
   float BackTrackEnergy       = 0;
   float BackTrackShowerEnergy = 0;
   float BackTrackCharge       = 0;
   float BackTrackShowerCharge = 0;
   
   //std::cout << "\t ANCESTOR comparing with MCShower of energy " << mcs.Start().E() << std::endl;
   //std::cout << "\t ANCESTOR start is [" << mcs.Start().X() << ", " << mcs.Start().Y() << ", " << mcs.Start().Z() << "]" << std::endl;
   //std::cout << "\t ANCESTOR HAS " << shrtrackIDs.size() << " particles" << std::endl;
   
   std::vector<simb::MCParticle const*> particle_vec;
   std::vector<anab::BackTrackerHitMatchingData const*> match_vec;  
   
   for (auto const& hit_idx : hit_idx_v) {
     
     particle_vec.clear(); match_vec.clear();
     
     backtrack_h.get(hit_idx, particle_vec, match_vec);
     
     // does this hit match to the mcshower?
     for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){            
       
       auto mctrkid = particle_vec.at(i_p)->TrackId();
       auto charge  = match_vec[i_p]->numElectrons;
       auto energy  = match_vec[i_p]->energy;
       
       BackTrackCharge += charge;
       BackTrackEnergy += energy;
       // does this trackID match that of the MCShower?
       for (auto const& shrtrkid : shrtrackIDs) {
	 if ( shrtrkid == (unsigned int)mctrkid ){
	   BackTrackShowerCharge += charge;
	   BackTrackShowerEnergy += energy;
	   break;
	 }
       }
     }// for all particles associated to this hit
   }// for all hits
   
   purity       = BackTrackShowerCharge / BackTrackCharge;
   completeness = BackTrackShowerEnergy / mcs.Start().E();
   
   return;
 }// BackTrack function end
 */
 

bool NeutrinoSelectionFilter::endSubRun(art::SubRun &subrun)
{
    if (!fData)
    {
        art::Handle<sumdata::POTSummary> potSummaryHandle;
        _pot = subrun.getByLabel(fMCTproducer, potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
        std::cout << "[LArPandoraExternalEventBuilding::endSubRun] Storing POT info!" << std::endl;
    }

    _run_sr = subrun.run();
    _sub_sr = subrun.subRun();
    _subrun_tree->Fill();  
    return true;
}


DEFINE_ART_MODULE(NeutrinoSelectionFilter)
