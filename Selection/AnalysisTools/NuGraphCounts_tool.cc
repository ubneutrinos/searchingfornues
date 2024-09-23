#ifndef ANALYSIS_NUGRAPHCOUNTS_CXX
#define ANALYSIS_NUGRAPHCOUNTS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardata/Utilities/FindManyInChainP.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       NuGraphCounts
// File:        NuGraphCounts.cc
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

class NuGraphCounts : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  NuGraphCounts(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~NuGraphCounts(){};

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
  //
  art::InputTag fCLSproducer; // cluster associated to PFP
  art::InputTag fSLCproducer; // slice associated to PFP
  art::InputTag fNG2producer; // nugraph2 producer

  template <typename T, typename A>
  int arg_max(std::vector<T, A> const& vec)
  {
    return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
  }

  // TTree variables

  //nu graph slice hit counts
  int slcng2mip;
  int slcng2hip;
  int slcng2shr;
  int slcng2mcl;
  int slcng2dfs;

  //nu graph clustered hit counts
  int clung2mip;
  int clung2hip;
  int clung2shr;
  int clung2mcl;
  int clung2dfs;

  //nu graph pfp counts
  std::vector<int> pfng2semlabel;
  std::vector<float> pfng2mipfrac;
  std::vector<float> pfng2hipfrac;
  std::vector<float> pfng2shrfrac;
  std::vector<float> pfng2mclfrac;
  std::vector<float> pfng2dfsfrac;

};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
NuGraphCounts::NuGraphCounts(const fhicl::ParameterSet &p)
{
  fCLSproducer = p.get<art::InputTag>("CLSproducer");
  fSLCproducer = p.get<art::InputTag>("SLCproducer");
  fNG2producer = p.get<art::InputTag>("NG2producer");
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void NuGraphCounts::configure(fhicl::ParameterSet const &p)
{
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void NuGraphCounts::analyzeEvent(art::Event const &e, bool fData) {}

void NuGraphCounts::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{

  ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer, proxy::withAssociated<recob::Hit>(fCLSproducer));
  // somehow proxies don't work for the slice-hit association, so go back to old assns
  art::ValidHandle<std::vector<recob::Slice>> inputSlice = e.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
  auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSlice, e, fSLCproducer));

  // auto GNNDescription = e.getHandle<anab::MVADescription<5>>(art::InputTag("NuGraph", "semantic"));

  auto const& hitsWithScores = proxy::getCollection<std::vector<recob::Hit>>(
    e,
    art::InputTag("gaushit"), //tag of the hit collection we ran the GNN on
    //proxy::withParallelData<anab::FeatureVector<1>>(art::InputTag("NuGraph", "filter")),
    proxy::withParallelData<anab::FeatureVector<5>>(fNG2producer));

  std::vector<int> ng2semclucounts(5,0);
  std::vector<int> ng2semslccounts(5,0);
  for (auto& h : hitsWithScores) {
    auto scores = h.get<anab::FeatureVector<5>>();
    std::vector<float> ng2semscores;
    for (size_t i=0;i<scores.size();i++) ng2semscores.push_back(scores[i]);
    unsigned int sem_label = arg_max(ng2semscores);
    ng2semslccounts[sem_label]++;
  }

  for (auto pfp : slice_pfp_v)
  {

    // get hits associated to this PFParticle through the clusters
    std::vector<art::Ptr<recob::Hit>> hit_v;
    auto clus_pxy_v = pfp.get<recob::Cluster>();
    if (clus_pxy_v.size() != 0)
    {
      for (auto ass_clus : clus_pxy_v)
      {
        // get cluster proxy
        const auto &clus = clus_proxy[ass_clus.key()];
        auto clus_hit_v = clus.get<recob::Hit>();
        for (const auto &hit : clus_hit_v)
          hit_v.push_back(hit);
      } // for all clusters associated to PFP
    }
    //
    if (hit_v.size()>0) {
      std::vector<int> ng2sempfpcounts(5,0);
      for (auto& hit : hit_v) {
	  auto scores = hitsWithScores[hit.key()].get<anab::FeatureVector<5>>();
	  std::vector<float> ng2semscores;
	  for (size_t i=0;i<scores.size();i++) ng2semscores.push_back(scores[i]);
	  unsigned int sem_label = arg_max(ng2semscores);
	  ng2sempfpcounts[sem_label]++;
	  ng2semclucounts[sem_label]++;
      }
      pfng2semlabel.push_back(arg_max(ng2sempfpcounts));
      pfng2mipfrac.push_back(float(ng2sempfpcounts[0])/hit_v.size());
      pfng2hipfrac.push_back(float(ng2sempfpcounts[1])/hit_v.size());
      pfng2shrfrac.push_back(float(ng2sempfpcounts[2])/hit_v.size());
      pfng2mclfrac.push_back(float(ng2sempfpcounts[3])/hit_v.size());
      pfng2dfsfrac.push_back(float(ng2sempfpcounts[4])/hit_v.size());
    } else {
      pfng2semlabel.push_back(-1);
      pfng2mipfrac.push_back(-1);
      pfng2hipfrac.push_back(-1);
      pfng2shrfrac.push_back(-1);
      pfng2mclfrac.push_back(-1);
      pfng2dfsfrac.push_back(-1);
    }
    //
  }
  //
  slcng2mip = ng2semslccounts[0];
  slcng2hip = ng2semslccounts[1];
  slcng2shr = ng2semslccounts[2];
  slcng2mcl = ng2semslccounts[3];
  slcng2dfs = ng2semslccounts[4];
  //
  clung2mip = ng2semclucounts[0];
  clung2hip = ng2semclucounts[1];
  clung2shr = ng2semclucounts[2];
  clung2mcl = ng2semclucounts[3];
  clung2dfs = ng2semclucounts[4];

  return;
}

void NuGraphCounts::setBranches(TTree *_tree)
{
  //
  _tree->Branch("slcng2mip", &slcng2mip, "slcng2mip/I");
  _tree->Branch("slcng2hip", &slcng2hip, "slcng2hip/I");
  _tree->Branch("slcng2shr", &slcng2shr, "slcng2shr/I");
  _tree->Branch("slcng2mcl", &slcng2mcl, "slcng2mcl/I");
  _tree->Branch("slcng2dfs", &slcng2dfs, "slcng2dfs/I");
  //
  _tree->Branch("clung2mip", &clung2mip, "clung2mip/I");
  _tree->Branch("clung2hip", &clung2hip, "clung2hip/I");
  _tree->Branch("clung2shr", &clung2shr, "clung2shr/I");
  _tree->Branch("clung2mcl", &clung2mcl, "clung2mcl/I");
  _tree->Branch("clung2dfs", &clung2dfs, "clung2dfs/I");
  //
  _tree->Branch("pfng2semlabel", &pfng2semlabel);
  _tree->Branch("pfng2mipfrac", &pfng2mipfrac);
  _tree->Branch("pfng2hipfrac", &pfng2hipfrac);
  _tree->Branch("pfng2shrfrac", &pfng2shrfrac);
  _tree->Branch("pfng2mclfrac", &pfng2mclfrac);
  _tree->Branch("pfng2dfsfrac", &pfng2dfsfrac);
}

void NuGraphCounts::resetTTree(TTree *_tree)
{
  //nu graph slice hit counts
  slcng2mip = std::numeric_limits<int>::min();
  slcng2hip = std::numeric_limits<int>::min();
  slcng2shr = std::numeric_limits<int>::min();
  slcng2mcl = std::numeric_limits<int>::min();
  slcng2dfs = std::numeric_limits<int>::min();
  //nu graph clustered hit counts
  clung2mip = std::numeric_limits<int>::min();
  clung2hip = std::numeric_limits<int>::min();
  clung2shr = std::numeric_limits<int>::min();
  clung2mcl = std::numeric_limits<int>::min();
  clung2dfs = std::numeric_limits<int>::min();
  //nu graph pfp counts
  pfng2semlabel.clear();
  pfng2mipfrac.clear();
  pfng2hipfrac.clear();
  pfng2shrfrac.clear();
  pfng2mclfrac.clear();
  pfng2dfsfrac.clear();
}

DEFINE_ART_CLASS_TOOL(NuGraphCounts)
} // namespace analysis

#endif
