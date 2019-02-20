////////////////////////////////////////////////////////////////////////
// Class:       NeutrinoSelection
// Plugin Type: analyzer (art v3_00_00)
// File:        NeutrinoSelection_module.cc
//
// Generated at Wed Jan 30 21:51:32 2019 by David Caratelli using cetskelgen
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
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
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardata/Utilities/FindManyInChainP.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

// selection tool
#include "SelectionTools/SelectionToolBase.h"

// saving output
#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"
#include "TVector3.h"

class NeutrinoSelection;


class NeutrinoSelection : public art::EDAnalyzer {
public:
  explicit NeutrinoSelection(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NeutrinoSelection(NeutrinoSelection const&) = delete;
  NeutrinoSelection(NeutrinoSelection&&) = delete;
  NeutrinoSelection& operator=(NeutrinoSelection const&) = delete;
  NeutrinoSelection& operator=(NeutrinoSelection&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void endSubRun(art::SubRun &subrun);

  using ProxyPfpColl_t = selection::ProxyPfpColl_t;
  using ProxyPfpElem_t =  selection::ProxyPfpElem_t;
  //using AssocPfpMeta_t = decltype(proxy::CollectionProxyBase< ProxyPfpElem_t, ProxyPfpColl_t, larpandoraobj::PFParticleMetadata >::aux());

private:

  art::InputTag fPFPproducer;
  art::InputTag fSHRproducer;
  art::InputTag fTRKproducer;
  art::InputTag fMCTproducer;
  bool fVerbose;
  bool fData;

  // TTree
  TTree* _tree;
  int   _run, _sub, _evt;       // event info
  float _nu_e;                  // neutrino energy [GeV]
  float _vtx_x, _vtx_y, _vtx_z; // neutrino interaction vertex coordinates [cm]
  int   _nu_pdg;                // neutrino PDG code
  int   _ccnc;                  // CC or NC tag from GENIE
  float _lep_e;                 // lepton energy (if one exists) [GeV]
  int  _pass;                   // does the slice pass the selection

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
   * @brief Save truth info for event associated to neutrino
   */
  void SaveTruth(art::Event const& e);

  /**
   * @brief reset ttree
   */
  void ResetTTree();

};


NeutrinoSelection::NeutrinoSelection(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  
  fPFPproducer = p.get< art::InputTag > ("PFPproducer");
  fSHRproducer = p.get< art::InputTag > ("SHRproducer");
  fTRKproducer = p.get< art::InputTag > ("TRKproducer");
  fMCTproducer = p.get< art::InputTag > ("MCTproducer");
  fVerbose     = p.get< bool >          ("Verbose");
  fData     = p.get< bool >          ("IsData");

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("NeutrinoSelection", "Neutrino Selection TTree");
  _tree->Branch("_nu_e"  ,&_nu_e  ,"nu_e/F"  );
  _tree->Branch("_vtx_x" ,&_vtx_x ,"vtx_x/F" );
  _tree->Branch("_vtx_y" ,&_vtx_y ,"vtx_y/F" );
  _tree->Branch("_vtx_z" ,&_vtx_z ,"vtx_z/F" );
  _tree->Branch("_nu_pdg",&_nu_pdg,"nu_pdg/I");
  _tree->Branch("_ccnc"  ,&_ccnc  ,"ccnc/I"  );
  _tree->Branch("_lep_e" ,&_lep_e ,"lep_e/F" );
  _tree->Branch("_pass"  ,&_pass  ,"pass/I"  );
  _tree->Branch("_run"   ,&_run   ,"run/I"   );
  _tree->Branch("_sub"   ,&_sub   ,"sub/I"   );
  _tree->Branch("_evt"   ,&_evt   ,"evt/I"   );

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

void NeutrinoSelection::analyze(art::Event const& e)
{

  ResetTTree();
  
  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();
  
  // SaveTruth
  SaveTruth(e);

  
  if (fVerbose) { std::cout << "new event : [run,event] : [" << e.run() << ", " << e.event() << "]" << std::endl; }
  
  // grab PFParticles in event
  ProxyPfpColl_t const& pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle> >(e,fPFPproducer,
								       proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
								       proxy::withAssociated<recob::Track>(fTRKproducer),
								       proxy::withAssociated<recob::Shower>(fSHRproducer));

  // build PFParticle map  for this event
  BuildPFPMap(pfp_proxy);
  
  // loop through PFParticles
  size_t p=0;
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
      
      if (selected) { _pass = 1; }
      
      _tree->Fill();

    }// if a neutrino PFParticle
    p++;
  }// for all PFParticles
  
  return;
}

template <typename T> void NeutrinoSelection::printPFParticleMetadata(const ProxyPfpElem_t& pfp_pxy,
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

void NeutrinoSelection::BuildPFPMap(const ProxyPfpColl_t& pfp_pxy_col) {
  
  _pfpmap.clear();

  unsigned int p=0;
  for (const auto& pfp_pxy : pfp_pxy_col) {
    _pfpmap[pfp_pxy->Self()] = p;
    p++;
  }

  return;
}// BuildPFPMap

void NeutrinoSelection:: AddDaughters(const ProxyPfpElem_t& pfp_pxy,
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

void NeutrinoSelection::ResetTTree() {

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

  _selectionTool->resetTTree(_tree);

}

void NeutrinoSelection::SaveTruth(art::Event const& e) {

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
  
  return;
}

void NeutrinoSelection::endSubRun(art::SubRun &subrun)
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
}


DEFINE_ART_MODULE(NeutrinoSelection)
