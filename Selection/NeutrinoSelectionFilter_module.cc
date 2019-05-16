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

#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

// selection tool
#include "SelectionTools/SelectionToolBase.h"

// analysis tool
#include "AnalysisTools/AnalysisToolBase.h"

// saving output
#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"
#include "TVector3.h"

class NeutrinoSelectionFilter;

class NeutrinoSelectionFilter : public art::EDFilter
{
public:
  explicit NeutrinoSelectionFilter(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NeutrinoSelectionFilter(NeutrinoSelectionFilter const &) = delete;
  NeutrinoSelectionFilter(NeutrinoSelectionFilter &&) = delete;
  NeutrinoSelectionFilter &operator=(NeutrinoSelectionFilter const &) = delete;
  NeutrinoSelectionFilter &operator=(NeutrinoSelectionFilter &&) = delete;

  // Required functions.
  bool filter(art::Event &e) override;

  // Selected optional functions.
  bool endSubRun(art::SubRun &subrun);

  using ProxyPfpColl_t = selection::ProxyPfpColl_t;
  using ProxyPfpElem_t = selection::ProxyPfpElem_t;
  //using AssocPfpMeta_t = decltype(proxy::CollectionProxyBase< ProxyPfpElem_t, ProxyPfpColl_t, larpandoraobj::PFParticleMetadata >::aux());

private:
  art::InputTag fPFPproducer;
  art::InputTag fCLSproducer; // cluster associated to PFP
  art::InputTag fSLCproducer; // slice associated to PFP
  art::InputTag fHITproducer; // hit associated to cluster
  art::InputTag fSHRproducer; // shower associated to PFP
  art::InputTag fVTXproducer; // vertex associated to PFP
  art::InputTag fTRKproducer; // track associated to PFP
  art::InputTag fMCTproducer;
  bool fVerbose;
  bool fData;
  bool fFilter;

  // TTree
  TTree *_tree;
  int _selected;

  TTree *_subrun_tree;
  int _run_sr; // The run number
  int _sub_sr; // The subRun number
  float _pot;  // The total amount of POT for the current sub run

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;

  // selection tool
  std::unique_ptr<::selection::SelectionToolBase> _selectionTool;

  // analysis tool
  std::vector<std::unique_ptr<::analysis::AnalysisToolBase>> _analysisToolsVec;

  /**
   * @brief function to builf a map linking PFParticle index to Self() attribute
   *
   * @input handle to event pfparticle record
   */
  void BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col);

  /**
   * @brief print PFParticle metadata information
   */
  template <typename T>
  void printPFParticleMetadata(const ProxyPfpElem_t &pfp_pxy,
                               const T &pfParticleMetadataList);

  /**
   * @brief build PFParticle hierarchy (i.e. slice) from parent [recursive function]
   *
   * @input pfp_pxy : parent pfparticle proxy for which to add daughters
   * @input pfp_pxy_col : evnt PFP proxy collection
   * @input slice_v : passed by reference, slice containing all PFParticles in hierarchy
   *
   */
  void AddDaughters(const ProxyPfpElem_t &pfp_pxy,
                    const ProxyPfpColl_t &pfp_pxy_col,
                    std::vector<ProxyPfpElem_t> &slice_v);

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
   * @brief reset ttree
   */
  void ResetTTree();
};

NeutrinoSelectionFilter::NeutrinoSelectionFilter(fhicl::ParameterSet const &p)
    : EDFilter{p} // ,
// More initializers here.
{

  fPFPproducer = p.get<art::InputTag>("PFPproducer");
  fSHRproducer = p.get<art::InputTag>("SHRproducer");
  fHITproducer = p.get<art::InputTag>("HITproducer");
  fCLSproducer = p.get<art::InputTag>("CLSproducer");
  fSLCproducer = p.get<art::InputTag>("SLCproducer");
  fVTXproducer = p.get<art::InputTag>("VTXproducer");
  fTRKproducer = p.get<art::InputTag>("TRKproducer");
  fMCTproducer = p.get<art::InputTag>("MCTproducer");
  fVerbose = p.get<bool>("Verbose");
  fData = p.get<bool>("IsData");
  fFilter = p.get<bool>("Filter", false);

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("NeutrinoSelectionFilter", "Neutrino Selection TTree");
  _tree->Branch("selected", &_selected, "selected/I");

  _subrun_tree = tfs->make<TTree>("SubRun", "SubRun TTree");
  _subrun_tree->Branch("run", &_run_sr, "run/I");
  _subrun_tree->Branch("subRun", &_sub_sr, "subRun/I");

  if (!fData)
    _subrun_tree->Branch("pot", &_pot, "pot/F");

  // configure and construct Selection Tool
  const fhicl::ParameterSet &selection_pset = p.get<fhicl::ParameterSet>("SelectionTool");
  _selectionTool = art::make_tool<::selection::SelectionToolBase>(selection_pset);

  // pass the TTree to the selection tool so that any branch can be added to it
  _selectionTool->setBranches(_tree);

  // set whether data
  _selectionTool->SetData(fData);

  // configure and construct Analysis Tool
  auto const tool_psets = p.get<fhicl::ParameterSet>("AnalysisTools");
  for (auto const &tool_pset_labels : tool_psets.get_pset_names())
  {
    auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_labels);
    _analysisToolsVec.push_back(art::make_tool<::analysis::AnalysisToolBase>(tool_pset));
  }

  // pass the TTree to the analysis tool so that any branch can be added to it
  for (size_t i = 0; i < _analysisToolsVec.size(); i++)
    _analysisToolsVec[i]->setBranches(_tree);
}

bool NeutrinoSelectionFilter::filter(art::Event &e)
{

  ResetTTree();

  if (fVerbose)
  {
    std::cout << "new event : [run,event] : [" << e.run() << ", " << e.event() << "]" << std::endl;
  }

  // grab PFParticles in event
  selection::ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, fPFPproducer,
                                                        proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
                                                        proxy::withAssociated<recob::Cluster>(fCLSproducer),
                                                        proxy::withAssociated<recob::Slice>(fSLCproducer),
                                                        proxy::withAssociated<recob::Track>(fTRKproducer),
                                                        proxy::withAssociated<recob::Vertex>(fVTXproducer),
                                                        proxy::withAssociated<recob::Shower>(fSHRproducer),
                                                        proxy::withAssociated<recob::SpacePoint>(fPFPproducer));

  BuildPFPMap(pfp_proxy);

  for (size_t i = 0; i < _analysisToolsVec.size(); i++)
  {
    _analysisToolsVec[i]->analyzeEvent(e, fData); //fixme add more arguments, make other functions down the line...
  }
  // loop through PFParticles
  size_t p = 0;
  // Should you keep this event? Yes, but only if there's a neutrino.
  // On the other hand, you can fill the tree even if you discart the event
  bool keepEvent = false;

  for (const ProxyPfpElem_t &pfp_pxy : pfp_proxy)
  {

    // get metadata for this PFP
    const auto &pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();

    //  find neutrino candidate
    if (pfp_pxy->IsPrimary() == false)
      continue;

    auto PDG = fabs(pfp_pxy->PdgCode());

    if ((PDG == 12) || (PDG == 14))
    {

      if (fVerbose)
        printPFParticleMetadata(pfp_pxy, pfParticleMetadataList);

      // collect PFParticle hierarchy originating from this neutrino candidate
      std::vector<ProxyPfpElem_t> slice_pfp_v;
      AddDaughters(pfp_pxy, pfp_proxy, slice_pfp_v);

      if (fVerbose)
      {
        std::cout << "This slice has " << slice_pfp_v.size() << " daughter PFParticles" << std::endl;
      }

      // create list of tracks and showers associated to this slice
      std::vector<art::Ptr<recob::Track>> sliceTracks;
      std::vector<art::Ptr<recob::Shower>> sliceShowers;

      for (auto pfp : slice_pfp_v)
      {
        // if there is a track associated to the PFParticle, add it
        auto const &ass_trk_v = pfp.get<recob::Track>();
        if (ass_trk_v.size() == 1)
          sliceTracks.push_back(ass_trk_v.at(0));
        // if there is a shower associated to the PFParticle, add it
        auto const &ass_shr_v = pfp.get<recob::Shower>();
        if (ass_shr_v.size() == 1)
          sliceShowers.push_back(ass_shr_v.at(0));
      } // for all PFParticles in the slice

      // check that # of PFP and # trk + shr matches
      if ((slice_pfp_v.size() - 1) != (sliceTracks.size() + sliceShowers.size()))
        std::cout << "ERROR : there are " << slice_pfp_v.size() << " PFP but " << sliceTracks.size() << " + " << sliceShowers.size() << " tracks + showers" << std::endl;

      // run selection on this slice
      bool selected = _selectionTool->selectEvent(e, slice_pfp_v);
      if (fVerbose && selected)
      {
        std::cout << "SLICE was selected!" << std::endl;
      }

      if (selected)
      {
        keepEvent = true;
        _selected = 1;
      }

      for (size_t i = 0; i < _analysisToolsVec.size(); i++) {
        _analysisToolsVec[i]->analyzeSlice(e, slice_pfp_v, fData, selected);
      }

    } // if a neutrino PFParticle
    p++;
  } // for all PFParticles

  _tree->Fill();

  if (fFilter == true)
    return keepEvent;

  return true;
}

template <typename T>
void NeutrinoSelectionFilter::printPFParticleMetadata(const ProxyPfpElem_t &pfp_pxy,
                                                      const T &pfParticleMetadataList)
{

  if (pfParticleMetadataList.size() != 0)
  {

    for (unsigned int j = 0; j < pfParticleMetadataList.size(); ++j)
    {
      const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
      auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
      if (!pfParticlePropertiesMap.empty())
      {
        if (fVerbose)
          std::cout << " Found PFParticle " << pfp_pxy->Self() << " with: " << std::endl;
        for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
        {
          if (fVerbose)
            std::cout << "  - " << it->first << " = " << it->second << std::endl;
        }
      }
    }
  } // if PFP metadata exists!
  return;
}

void NeutrinoSelectionFilter::BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col)
{

  _pfpmap.clear();

  unsigned int p = 0;
  for (const auto &pfp_pxy : pfp_pxy_col)
  {
    _pfpmap[pfp_pxy->Self()] = p;
    p++;
  }

  return;
} // BuildPFPMap

void NeutrinoSelectionFilter::AddDaughters(const ProxyPfpElem_t &pfp_pxy,
                                           const ProxyPfpColl_t &pfp_pxy_col,
                                           std::vector<ProxyPfpElem_t> &slice_v)
{

  auto daughters = pfp_pxy->Daughters();

  slice_v.push_back(pfp_pxy);

  if (fVerbose)
    std::cout << "\t PFP w/ PdgCode " << pfp_pxy->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;

  for (auto const &daughterid : daughters)
  {

    if (_pfpmap.find(daughterid) == _pfpmap.end())
    {
      std::cout << "Did not find DAUGHTERID in map! error" << std::endl;
      continue;
    }

    // const art::Ptr<recob::PFParticle> pfp_pxy(pfp_pxy_col, _pfpmap.at(daughterid) );
    auto pfp_pxy2 = pfp_pxy_col.begin();
    for (size_t j = 0; j < _pfpmap.at(daughterid); ++j)
      ++pfp_pxy2;
    // const T& pfp_pxy2 = (pfp_pxy_col.begin()+_pfpmap.at(daughterid));

    AddDaughters(*pfp_pxy2, pfp_pxy_col, slice_v);

  } // for all daughters

  return;
} // AddDaughters

void NeutrinoSelectionFilter::ResetTTree()
{

  _selected = 0;

  _selectionTool->resetTTree(_tree);
  for (size_t i = 0; i < _analysisToolsVec.size(); i++)
    _analysisToolsVec[i]->resetTTree(_tree);
}

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
