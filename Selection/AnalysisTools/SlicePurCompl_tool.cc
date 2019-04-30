#ifndef ANALYSIS_SLICEPURCOMPL_CXX
#define ANALYSIS_SLICEPURCOMPL_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       SlicePurCompl
    // File:        SlicePurCompl.cc
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

  class SlicePurCompl : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    SlicePurCompl(const fhicl::ParameterSet& pset);

    /**
     *  @brief  Destructor
     */
    ~SlicePurCompl(){ };

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
    //
    art::InputTag fCLSproducer; // cluster associated to PFP
    art::InputTag fSLCproducer; // slice associated to PFP
    art::InputTag fMCTproducer;
    art::InputTag fMCPproducer;
    art::InputTag fHproducer;
    art::InputTag fHTproducer;
    //
    std::vector<size_t> lepid; //lepton
    std::vector<size_t> proid; //protons
    std::vector<size_t> pi1id; //charged pions
    std::vector<size_t> pi0id; //neutral pions
    std::vector<size_t> neuid; //neutrons
    std::vector<size_t> gamid; //nphotons
    std::vector<size_t> othid; //other
    //
    std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> > assocMCPart;
    //
    void fillId(const std::vector<double>& p, const simb::MCParticle& mcp, std::vector<size_t>& id);
    //
    void incrementCounts(const std::vector<art::Ptr<recob::Hit> >& hits,
                         const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& assocMCPart,
                         const std::vector<size_t>& lepid, const std::vector<size_t>& proid, const std::vector<size_t>& pi1id, const std::vector<size_t>& pi0id,
                         const std::vector<size_t>& neuid, const std::vector<size_t>& gamid, const std::vector<size_t>& othid,
                         int& nlephits, int& nprohits, int& npi1hits, int& npi0hits, int& nneuhits, int& ngamhits, int& nothhits) const;

    // TTree variables
    // per event
    int evnunhits;
    int evlepnhits;
    int evpronhits;
    int evpi1nhits;
    int evpi0nhits;
    int evneunhits;
    int evgamnhits;
    int evothnhits;
    //
    // per nu slice
    int slnunhits;
    int sllepnhits;
    int slpronhits;
    int slpi1nhits;
    int slpi0nhits;
    int slneunhits;
    int slgamnhits;
    int slothnhits;
    //
    // per pfp in nu slice
    std::vector<int> pfnunhits;
    std::vector<int> pflepnhits;
    std::vector<int> pfpronhits;
    std::vector<int> pfpi1nhits;
    std::vector<int> pfpi0nhits;
    std::vector<int> pfneunhits;
    std::vector<int> pfgamnhits;
    std::vector<int> pfothnhits;
  };

  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  SlicePurCompl::SlicePurCompl(const fhicl::ParameterSet& p)
  {
    fCLSproducer = p.get< art::InputTag > ("CLSproducer");
    fSLCproducer = p.get< art::InputTag > ("SLCproducer");
    fMCTproducer = p.get< art::InputTag > ("MCTproducer");
    fMCPproducer = p.get< art::InputTag > ("MCPproducer");
    fHproducer = p.get< art::InputTag > ("Hproducer");
    fHTproducer = p.get< art::InputTag > ("HTproducer");
  }

  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void SlicePurCompl::configure(fhicl::ParameterSet const & p)
  {
  }

  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void SlicePurCompl::analyzeEvent(art::Event const& e, bool fData)
  {

    if (fData) return;

    // momentum of primary products of the neutrino interaction
    std::vector<double> plep; //lepton
    std::vector<double> ppro; //protons
    std::vector<double> ppi1; //charged pions
    std::vector<double> ppi0; //neutral pions
    std::vector<double> pneu; //neutrons
    std::vector<double> pgam; //nphotons
    std::vector<double> poth; //other

    // generator, simb::MCTruth
    art::ValidHandle<std::vector<simb::MCTruth> > inputMCTruth = e.getValidHandle<std::vector<simb::MCTruth> >(fMCTproducer);
    if (inputMCTruth->size()<1) return;
    const auto& mct = inputMCTruth->at(0);
    plep.push_back( mct.GetNeutrino().Lepton().P() );
    for ( int im=0; im<mct.NParticles(); ++im ) {
      const auto& mcp = mct.GetParticle(im);
      if (mcp.StatusCode()!=1) continue;
      if (std::abs(mcp.PdgCode())==2212)      ppro.push_back( mcp.P() );
      else if (std::abs(mcp.PdgCode())==211)  ppi1.push_back( mcp.P() );
      else if (std::abs(mcp.PdgCode())==111)  ppi0.push_back( mcp.P() );
      else if (std::abs(mcp.PdgCode())==2112) pneu.push_back( mcp.P() );
      else if (std::abs(mcp.PdgCode())==22)   pgam.push_back( mcp.P() );
      else                                    poth.push_back( mcp.P() );
    }

    lepid.clear();
    proid.clear();
    pi1id.clear();
    pi0id.clear();
    neuid.clear();
    gamid.clear();
    othid.clear();

    // largeant, simb::MCParticle
    art::ValidHandle<std::vector<simb::MCParticle> > inputMCParticle = e.getValidHandle<std::vector<simb::MCParticle> >(fMCPproducer);
    for ( unsigned int im=0; im<inputMCParticle->size(); ++im ) {
      const auto& mcp = inputMCParticle->at(im);
      if (mcp.StatusCode()==1) {
        //
        if (mcp.Mother()==0) {
          // try to match primaries
          if (abs(mcp.PdgCode())==11 || abs(mcp.PdgCode())==13) fillId(plep, mcp, lepid);
          else if (std::abs(mcp.PdgCode())==2212)               fillId(ppro, mcp, proid);
          else if (std::abs(mcp.PdgCode())==211)                fillId(ppi1, mcp, pi1id);
          else if (std::abs(mcp.PdgCode())==111)                fillId(ppi0, mcp, pi0id);
          else if (std::abs(mcp.PdgCode())==2112)               fillId(pneu, mcp, neuid);
          else if (std::abs(mcp.PdgCode())==22)                 fillId(pgam, mcp, gamid);
          else                                                  fillId(poth, mcp, othid);
        } else {
          // now deal with secondaries navigating history
          if (std::find(lepid.begin(), lepid.end(), mcp.Mother())!=lepid.end()) lepid.push_back(mcp.TrackId());
          else if (std::find(proid.begin(), proid.end(), mcp.Mother())!=proid.end()) proid.push_back(mcp.TrackId());
          else if (std::find(pi1id.begin(), pi1id.end(), mcp.Mother())!=pi1id.end()) pi1id.push_back(mcp.TrackId());
          else if (std::find(pi0id.begin(), pi0id.end(), mcp.Mother())!=pi0id.end()) pi0id.push_back(mcp.TrackId());
          else if (std::find(neuid.begin(), neuid.end(), mcp.Mother())!=neuid.end()) neuid.push_back(mcp.TrackId());
          else if (std::find(gamid.begin(), gamid.end(), mcp.Mother())!=gamid.end()) gamid.push_back(mcp.TrackId());
          else if (std::find(othid.begin(), othid.end(), mcp.Mother())!=othid.end()) othid.push_back(mcp.TrackId());
        }
      }
    }

    art::ValidHandle<std::vector<recob::Hit> > inputHits = e.getValidHandle<std::vector<recob::Hit> >(fHproducer);
    assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>(inputHits, e, fHTproducer));
    int nlephits = 0, nprohits = 0, npi1hits = 0, npi0hits = 0, nneuhits = 0, ngamhits = 0, nothhits = 0;
    std::vector<art::Ptr<recob::Hit> > inHitsPtrV;
    for (unsigned int ih=0; ih<inputHits->size(); ih++) inHitsPtrV.push_back({inputHits,ih});
    incrementCounts(inHitsPtrV, assocMCPart,
		    lepid, proid, pi1id, pi0id, neuid, gamid, othid,
		    nlephits, nprohits, npi1hits, npi0hits, nneuhits, ngamhits, nothhits);
    evnunhits = nlephits+nprohits+npi1hits+npi0hits+nneuhits+ngamhits+nothhits;
    evlepnhits = nlephits;
    evpronhits = nprohits;
    evpi1nhits = npi1hits;
    evpi0nhits = npi0hits;
    evneunhits = nneuhits;
    evgamnhits = ngamhits;
    evothnhits = nothhits;

    return;
  }

  void SlicePurCompl::analyzeSlice(art::Event const& e, std::vector<ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected)
  {

    if (fData) return;

    ProxyClusColl_t const& clus_proxy = proxy::getCollection<std::vector<recob::Cluster> >(e, fCLSproducer, proxy::withAssociated<recob::Hit>(fCLSproducer));
    // somehow proxies don't work for the slice-hit association, so go back to old assns
    art::ValidHandle<std::vector<recob::Slice> > inputSlice = e.getValidHandle<std::vector<recob::Slice> >(fSLCproducer);
    auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit> >(new art::FindManyP<recob::Hit>(inputSlice, e, fSLCproducer));

    for (auto pfp :  slice_pfp_v) {
      if (pfp->IsPrimary()) {
	// check associated slice only once per hierarchy
	auto slice_pxy_v = pfp.get<recob::Slice>();
	if (slice_pxy_v.size()!=1) {
	  std::cout << "WRONG!!! n slices = " << slice_pxy_v.size() << " " << __FILE__ << " " << __LINE__ << std::endl;
	  return;
	}
	auto slicehits = assocSliceHit->at(slice_pxy_v[0].key());

	int nsllephits=0, nslprohits=0, nslpi1hits=0, nslpi0hits=0, nslneuhits=0, nslgamhits=0, nslothhits=0;
	incrementCounts(slicehits, assocMCPart,
			lepid, proid, pi1id, pi0id, neuid, gamid, othid,
			nsllephits, nslprohits, nslpi1hits, nslpi0hits, nslneuhits, nslgamhits, nslothhits);
	slnunhits = nsllephits+nslprohits+nslpi1hits+nslpi0hits+nslneuhits+nslgamhits+nslothhits;
	sllepnhits = nsllephits;
	slpronhits = nslprohits;
	slpi1nhits = nslpi1hits;
	slpi0nhits = nslpi0hits;
	slneunhits = nslneuhits;
	slgamnhits = nslgamhits;
	slothnhits = nslothhits;
      }
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
      }
      int npflephits=0, npfprohits=0, npfpi1hits=0, npfpi0hits=0, npfneuhits=0, npfgamhits=0, npfothhits=0;
      incrementCounts(hit_v, assocMCPart,
		      lepid, proid, pi1id, pi0id, neuid, gamid, othid,
		      npflephits, npfprohits, npfpi1hits, npfpi0hits, npfneuhits, npfgamhits, npfothhits);
      //
      pfnunhits.push_back(npflephits+npfprohits+npfpi1hits+npfpi0hits+npfneuhits+npfgamhits+npfothhits);
      pflepnhits.push_back(npflephits);
      pfpronhits.push_back(npfprohits);
      pfpi1nhits.push_back(npfpi1hits);
      pfpi0nhits.push_back(npfpi0hits);
      pfneunhits.push_back(npfneuhits);
      pfgamnhits.push_back(npfgamhits);
      pfothnhits.push_back(npfothhits);
    }

    return;
  }

  void SlicePurCompl::setBranches(TTree* _tree)
  {
    _tree->Branch("evnunhits",&evnunhits,"evnunhits/I");
    _tree->Branch("evlepnhits",&evlepnhits,"evlepnhits/I");
    _tree->Branch("evpronhits",&evpronhits,"evpronhits/I");
    _tree->Branch("evpi1nhits",&evpi1nhits,"evpi1nhits/I");
    _tree->Branch("evpi0nhits",&evpi0nhits,"evpi0nhits/I");
    _tree->Branch("evneunhits",&evneunhits,"evneunhits/I");
    _tree->Branch("evgamnhits",&evgamnhits,"evgamnhits/I");
    _tree->Branch("evothnhits",&evothnhits,"evothnhits/I");
    //
    _tree->Branch("slnunhits",&slnunhits,"slnunhits/I");
    _tree->Branch("sllepnhits",&sllepnhits,"sllepnhits/I");
    _tree->Branch("slpronhits",&slpronhits,"slpronhits/I");
    _tree->Branch("slpi1nhits",&slpi1nhits,"slpi1nhits/I");
    _tree->Branch("slpi0nhits",&slpi0nhits,"slpi0nhits/I");
    _tree->Branch("slneunhits",&slneunhits,"slneunhits/I");
    _tree->Branch("slgamnhits",&slgamnhits,"slgamnhits/I");
    _tree->Branch("slothnhits",&slothnhits,"slothnhits/I");
    //
    _tree->Branch("pfnunhits",&pfnunhits);
    _tree->Branch("pflepnhits",&pflepnhits);
    _tree->Branch("pfpronhits",&pfpronhits);
    _tree->Branch("pfpi1nhits",&pfpi1nhits);
    _tree->Branch("pfpi0nhits",&pfpi0nhits);
    _tree->Branch("pfneunhits",&pfneunhits);
    _tree->Branch("pfgamnhits",&pfgamnhits);
    _tree->Branch("pfothnhits",&pfothnhits);
  }

  void SlicePurCompl::resetTTree(TTree* _tree)
  {
    evnunhits = std::numeric_limits<int>::min();
    evlepnhits = std::numeric_limits<int>::min();
    evpronhits = std::numeric_limits<int>::min();
    evpi1nhits = std::numeric_limits<int>::min();
    evpi0nhits = std::numeric_limits<int>::min();
    evneunhits = std::numeric_limits<int>::min();
    evgamnhits = std::numeric_limits<int>::min();
    evothnhits = std::numeric_limits<int>::min();
    //
    slnunhits = std::numeric_limits<int>::min();
    sllepnhits = std::numeric_limits<int>::min();
    slpronhits = std::numeric_limits<int>::min();
    slpi1nhits = std::numeric_limits<int>::min();
    slpi0nhits = std::numeric_limits<int>::min();
    slneunhits = std::numeric_limits<int>::min();
    slgamnhits = std::numeric_limits<int>::min();
    slothnhits = std::numeric_limits<int>::min();
    //
    pfnunhits.clear();
    pflepnhits.clear();
    pfpronhits.clear();
    pfpi1nhits.clear();
    pfpi0nhits.clear();
    pfneunhits.clear();
    pfgamnhits.clear();
    pfothnhits.clear();
  }

  void SlicePurCompl::fillId(const std::vector<double>& p, const simb::MCParticle& mcp, std::vector<size_t>& id) {
    for (unsigned int ip=0; ip<p.size(); ip++) {
      if (std::abs(p[ip]-mcp.P())<0.00000001) {
	id.push_back( mcp.TrackId() );
      }
    }
    return;
  }

  void SlicePurCompl::incrementCounts(const std::vector<art::Ptr<recob::Hit> >& hits,
				      const std::unique_ptr<art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> >& assocMCPart,
				      const std::vector<size_t>& lepid, const std::vector<size_t>& proid, const std::vector<size_t>& pi1id, const std::vector<size_t>& pi0id,
				      const std::vector<size_t>& neuid, const std::vector<size_t>& gamid, const std::vector<size_t>& othid,
				      int& nlephits, int& nprohits, int& npi1hits, int& npi0hits, int& nneuhits, int& ngamhits, int& nothhits) const
  {
    for (unsigned int ih=0; ih<hits.size(); ih++) {
      art::Ptr<recob::Hit> hitp = hits[ih];
      auto assmcp = assocMCPart->at(hitp.key());
      auto assmdt = assocMCPart->data(hitp.key());
      for (unsigned int ia=0; ia<assmcp.size(); ++ia) {
	auto mcp = assmcp[ia];
	auto amd = assmdt[ia];
	if (amd->isMaxIDE!=1) continue;
	if (std::find(lepid.begin(),lepid.end(),mcp->TrackId())!=lepid.end()) nlephits++;
	else if (std::find(proid.begin(),proid.end(),mcp->TrackId())!=proid.end()) nprohits++;
	else if (std::find(pi1id.begin(),pi1id.end(),mcp->TrackId())!=pi1id.end()) npi1hits++;
	else if (std::find(pi0id.begin(),pi0id.end(),mcp->TrackId())!=pi0id.end()) npi0hits++;
	else if (std::find(neuid.begin(),neuid.end(),mcp->TrackId())!=neuid.end()) nneuhits++;
	else if (std::find(gamid.begin(),gamid.end(),mcp->TrackId())!=gamid.end()) ngamhits++;
	else if (std::find(othid.begin(),othid.end(),mcp->TrackId())!=othid.end()) nothhits++;
      }
    }
  }

  DEFINE_ART_CLASS_TOOL(SlicePurCompl)
} // namespace analysis

#endif
