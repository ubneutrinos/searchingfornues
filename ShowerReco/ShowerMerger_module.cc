////////////////////////////////////////////////////////////////////////
// Class:       ShowerMerger
// Plugin Type: producer (art v3_01_02)
// File:        ShowerMerger_module.cc
//
// Generated at Sat Apr 13 17:16:11 2019 by David Caratelli using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "art/Persistency/Common/PtrMaker.h"

#include "ubana/ubana/searchingfornues/Selection/CommonDefs/Typedefs.h"
#include "ubana/ubana/searchingfornues/Selection/CommonDefs/TrackShowerScoreFuncs.h"
#include "ubana/ubana/searchingfornues/Selection/CommonDefs/BacktrackingFuncs.h"

#include <memory>

#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"

class ShowerMerger;


class ShowerMerger : public art::EDProducer {
public:
  explicit ShowerMerger(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShowerMerger(ShowerMerger const&) = delete;
  ShowerMerger(ShowerMerger&&) = delete;
  ShowerMerger& operator=(ShowerMerger const&) = delete;
  ShowerMerger& operator=(ShowerMerger&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  art::InputTag fMCTproducer;
  art::InputTag fPFPproducer;
  art::InputTag fCLSproducer; // cluster associated to PFP
  art::InputTag fSLCproducer; // slice associated to PFP
  art::InputTag fHITproducer; // hit associated to cluster
  art::InputTag fSHRproducer; // shower associated to PFP
  art::InputTag fVTXproducer; // vertex associated to PFP
  art::InputTag fTRKproducer; // track associated to PFP

  // trunk-merger 
  float fShrVtxTrkDistMax; // maximum distance between shower start and track start/end
  float fShrTrkDotMin;     // minimum alignment between track and shower direction
  // branch-merger
  float fMinBranchConeDot; // minimum dot product between shr vtx -> branch vtx vector and shr dir vector
  float fMinBranchVtxDist; // minimum distance of branch vertex to shower start point

  // TTree
  TTree* _tree;
  float _nu_e, _elec_e;
  int   _nelec, _nreco;
  float _vtx_x, _vtx_y, _vtx_z, _vtx_t;
  float _Ashr_x, _Ashr_y, _Ashr_z, _Ashr_q;
  float _Pshr_x, _Pshr_y, _Pshr_z, _Pshr_q;
  float _xtimeoffset, _xsceoffset, _ysceoffset, _zsceoffset; // offsets for generation time and SCE

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;

  bool FindShowerBranches(const size_t shr_pfp_idx,
			  const std::vector<searchingfornues::ProxyPfpElem_t>& slice_pfp_v,
			  std::vector<size_t>& branch_pfp_idx_v);
    
  bool FindShowerTrunk(const size_t shr_pfp_idx,
		       const std::vector<searchingfornues::ProxyPfpElem_t>& slice_pfp_v,
		       TVector3& vtxcandidate,
		       size_t& trktrunk_pfp_idx);
  
  /**
   * @brief function to builf a map linking PFParticle index to Self() attribute
   *
   * @input handle to event pfparticle record
   */
    void BuildPFPMap(const searchingfornues::ProxyPfpColl_t& pfp_pxy_col);

  /**
   * @brief build PFParticle hierarchy (i.e. slice) from parent [recursive function]
   *
   * @input pfp_pxy : parent pfparticle proxy for which to add daughters
   * @input pfp_pxy_col : evnt PFP proxy collection
   * @input slice_v : passed by reference, slice containing all PFParticles in hierarchy
   *
   */
  void AddDaughters(const searchingfornues::ProxyPfpElem_t& pfp_pxy,
		    const searchingfornues::ProxyPfpColl_t& pfp_pxy_col,
		    std::vector<searchingfornues::ProxyPfpElem_t>& slice_v);

  void SaveTruth(art::Event const& e);

};


ShowerMerger::ShowerMerger(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{

  fMCTproducer = p.get< art::InputTag > ("MCTproducer");
  fPFPproducer = p.get< art::InputTag > ("PFPproducer");
  fSHRproducer = p.get< art::InputTag > ("SHRproducer");
  fHITproducer = p.get< art::InputTag > ("HITproducer");
  fCLSproducer = p.get< art::InputTag > ("CLSproducer");
  fSLCproducer = p.get< art::InputTag > ("SLCproducer");
  fVTXproducer = p.get< art::InputTag > ("VTXproducer");
  fTRKproducer = p.get< art::InputTag > ("TRKproducer");

  fShrVtxTrkDistMax = p.get< float > ("ShrVtxTrkDistMax");
  fShrTrkDotMin     = p.get< float > ("ShrTrkDotMin"    );
  fMinBranchConeDot = p.get< float > ("MinBranchConeDot");
  fMinBranchVtxDist = p.get< float > ("MinBranchVtxDist");

  produces<std::vector<recob::PFParticle> >();
  produces<std::vector<recob::Vertex>     >();
  produces<std::vector<recob::Cluster>    >();

  produces<art::Assns <recob::PFParticle, recob::Vertex>  >();
  produces<art::Assns <recob::PFParticle, recob::Cluster> >();
  produces<art::Assns <recob::PFParticle, larpandoraobj::PFParticleMetadata> >();

  produces<art::Assns <recob::Cluster, recob::Hit>  >();

  // TTree
  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("ShowerMerging", "Shower Merging");
  // truth
  _tree->Branch("nu_e"  ,&_nu_e  ,"nu_e/F"  );
  _tree->Branch("vtx_x" ,&_vtx_x ,"vtx_x/F" );
  _tree->Branch("vtx_y" ,&_vtx_y ,"vtx_y/F" );
  _tree->Branch("vtx_z" ,&_vtx_z ,"vtx_z/F" );
  _tree->Branch("vtx_t" ,&_vtx_t ,"vtx_t/F" );
  _tree->Branch("nelec",&_nelec,"nelec/I");
  _tree->Branch("elec_e",&_elec_e,"elec_e/F");
  // reco
  _tree->Branch("nreco",&_nreco,"nreco/I");
  _tree->Branch("Ashr_x",&_Ashr_x,"Ashr_x/F");
  _tree->Branch("Ashr_y",&_Ashr_y,"Ashr_y/F");
  _tree->Branch("Ashr_z",&_Ashr_z,"Ashr_z/F");
  _tree->Branch("Ashr_q",&_Ashr_q,"Ashr_q/F");
  _tree->Branch("Pshr_x",&_Pshr_x,"Pshr_x/F");
  _tree->Branch("Pshr_y",&_Pshr_y,"Pshr_y/F");
  _tree->Branch("Pshr_z",&_Pshr_z,"Pshr_z/F");
  _tree->Branch("Pshr_q",&_Pshr_q,"Pshr_q/F");
  // SCE
  _tree->Branch("xtimeoffset",&_xtimeoffset,"xtimeoffset/F");
  _tree->Branch("xsceoffset" ,&_xsceoffset ,"xsceoffset/F" );
  _tree->Branch("ysceoffset" ,&_ysceoffset ,"ysceoffset/F" );
  _tree->Branch("zsceoffset" ,&_zsceoffset ,"zsceoffset/F" );
  
  return;
}

void ShowerMerger::produce(art::Event& e)
{
  
  std::cout << "DAVIDC PRODUCE" << std::endl;

  SaveTruth(e);

  // PFP pointer maker
  art::PtrMaker<recob::PFParticle> PFParticlePtrMaker(e);
  // Vtx  pointer maker
  art::PtrMaker<recob::Vertex>     VtxPtrMaker(e);
  // Cluster  pointer maker
  art::PtrMaker<recob::Cluster>    ClsPtrMaker(e);

  // produce output PFParticles
  std::unique_ptr< std::vector<recob::PFParticle> > PFParticle_v(new std::vector<recob::PFParticle> );  
  std::unique_ptr< std::vector<recob::Vertex>     > Vertex_v(new std::vector<recob::Vertex>         );  
  std::unique_ptr< std::vector<recob::Cluster>    > Cluster_v(new std::vector<recob::Cluster>       );  

  std::unique_ptr< art::Assns <recob::PFParticle, recob::Vertex>  > PFP_Vtx_assn_v    ( new art::Assns<recob::PFParticle, recob::Vertex> );
  std::unique_ptr< art::Assns <recob::PFParticle, recob::Cluster> > PFP_Cls_assn_v    ( new art::Assns<recob::PFParticle, recob::Cluster>);
  std::unique_ptr< art::Assns <recob::PFParticle, larpandoraobj::PFParticleMetadata> > PFP_Meta_assn_v    ( new art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata>);

  std::unique_ptr< art::Assns <recob::Cluster, recob::Hit>  > Cls_Hit_assn_v    ( new art::Assns<recob::Cluster, recob::Hit> );


  // grab PFParticles in event
  searchingfornues::ProxyPfpColl_t const& pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle> >(e,fPFPproducer,
													    proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
													    proxy::withAssociated<recob::Cluster>(fCLSproducer),
													    proxy::withAssociated<recob::Slice>(fSLCproducer),
													    proxy::withAssociated<recob::Track>(fTRKproducer),
													    proxy::withAssociated<recob::Vertex>(fVTXproducer),
													    proxy::withAssociated<recob::Shower>(fSHRproducer));  

  // grab cluster -> hit association
  searchingfornues::ProxyClusColl_t const& clus_proxy = proxy::getCollection<std::vector<recob::Cluster> >(e, fCLSproducer,
													   proxy::withAssociated<recob::Hit>(fCLSproducer));

  // get pfparticle vector
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);
  /*
  // get metadata vector
  auto const& meta_h = e.getValidHandle<std::vector<larpandoraobj::PFParticleMetadata> >(fPFPproducer);
  // get cluster vector
  auto const& cls_h = e.getValidHandle<std::vector<recob::Cluster> >(fCLSproducer);
  */

  // keep track of largest charge shower
  float QMax = 0;
  bool largestShower = false;
  _nreco = 0;
  
  // build PFParticle map  for this event
  BuildPFPMap(pfp_proxy);

  // collect PFParticle hierarchy originating from this neutrino candidate
  std::vector<searchingfornues::ProxyPfpElem_t> slice_pfp_v;

  for (const searchingfornues::ProxyPfpElem_t& pfp_pxy : pfp_proxy) {

    // get metadata for this PFP
    //const auto& pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();
    
    //  find neutrino candidate
    if (pfp_pxy->IsPrimary() == false) continue;
    
    auto PDG = fabs(pfp_pxy->PdgCode());

    if ( (PDG == 12) || (PDG == 14) ) {
      
      AddDaughters(pfp_pxy, pfp_proxy, slice_pfp_v);

    }// if neutrino like PFP

  }// for all PFParticles

  std::cout << "New slice w/ " << slice_pfp_v.size() << " PFParticles " << std::endl;

  // go through slice and find shower-like particles
  for (size_t p=0; p < slice_pfp_v.size(); p++) {

    auto pfp_pxy = slice_pfp_v[p];

    auto PDG = fabs(pfp_pxy->PdgCode());
    
    if ( (PDG == 12) || (PDG == 14) )
      continue;

    auto ass_shr_v = pfp_pxy.get<recob::Shower>();
    //auto ass_trk_v = pfp_pxy.get<recob::Track>();
    //auto trkscore = searchingfornues::GetTrackShowerScore(slice_pfp_v[i_pfp]);
    
    if (ass_shr_v.size() == 1) {

      // create new PFP for output
      recob::PFParticle outshrpfp( pfp_h->at( pfp_pxy.index() ) );

      // prepare new vertex location and track trunk
      TVector3 vtxcandidate;
      size_t   trktrunk_pfp_idx;
      Double_t xyz[3] = {};

      largestShower = false;

      // load shower clusters
      auto pfp_clus_v = pfp_pxy.get<recob::Cluster>();
      // create one cluster per plane (will save non-empty ones)
      std::vector<recob::Cluster> out_clus_v(3, recob::Cluster());
      // vector of hit art::Ptr for cluster -> hit associations
      std::vector< std::vector< art::Ptr<recob::Hit> > > out_clus_hit_assn_v(3, std::vector< art::Ptr<recob::Hit> >());
      for (size_t c=0; c < pfp_clus_v.size(); c++) {
	auto pl = pfp_clus_v[c]->Plane().Plane;
	recob::Cluster newclus( *(pfp_clus_v[c]) );
	out_clus_v[ pl ] = newclus;

	_nreco += 1;

	if (pl == 2) {
	  //std::cout << "cluster has " << newclus.Integral() << " charge" << std::endl;
	  if ( newclus.Integral() > QMax ) {
	    largestShower = true;
	    QMax = newclus.Integral();
	    _Ashr_q = QMax;
	    auto shr = pfp_pxy.get<recob::Shower>()[0];
	    _Ashr_x = shr->ShowerStart().X();
	    _Ashr_y = shr->ShowerStart().Y();
	    _Ashr_z = shr->ShowerStart().Z();
	  }
	}

	// clus -> hit associations
	const auto& clus_pxy = clus_proxy[pfp_clus_v[c].key()];
	auto clus_hit_v = clus_pxy.get<recob::Hit>();
	for (size_t i=0; i < clus_hit_v.size(); i++) {
	  out_clus_hit_assn_v[ pl ].push_back( clus_hit_v[i] );
	}// for all hits
	std::cout << "there are " << out_clus_hit_assn_v[pl].size() << " hits associated with this cluster" << std::endl;
      }// for all clusters
      std::cout << "There are " << out_clus_v.size() << " output clusters" << std::endl;

      if ( FindShowerTrunk(p, slice_pfp_v, vtxcandidate, trktrunk_pfp_idx) == true ) {
	
	// load new clusters to be added
	auto new_clus_v = slice_pfp_v[trktrunk_pfp_idx].get<recob::Cluster>();
	std::cout << "Will be adding " << new_clus_v.size() << " clusters!" << std::endl;
	// beacuse this is the trunk of the shower, we need to modify the cluster
	// start/end points to match these new ones
	for (size_t c=0; c < new_clus_v.size(); c++) {
	  auto newclus = new_clus_v.at(c);
	  auto newPl = newclus->Plane().Plane;
	  if ( (newPl < out_clus_v.size()) && (newPl >= 0) ) {

	    // store extra hits associated

	    // clus -> hit associations
	    const auto& new_clus_pxy = clus_proxy[newclus.key()];
	    auto new_clus_hit_v = new_clus_pxy.get<recob::Hit>();
	    for (size_t i=0; i < new_clus_hit_v.size(); i++) {
	      out_clus_hit_assn_v[ newPl ].push_back( new_clus_hit_v[i] );
	    }// for all hits

	    // if the old cluster is garbage, simply load the new one
	    if (out_clus_v[ newPl ].NHits() == 0) {
	      out_clus_v[ newPl ] = (*newclus);
	    }
	    else {
	      auto start_wire   = newclus->StartWire();
	      auto start_wireS  = newclus->SigmaStartWire();
	      auto start_tick   = newclus->StartTick();
	      auto start_tickS  = newclus->SigmaStartTick();
	      auto start_charge = newclus->StartCharge();
	      auto start_angle  = newclus->StartAngle();
	      auto start_open   = newclus->StartOpeningAngle();
	      auto end_wire     = out_clus_v[ newPl ].EndWire();
	      auto end_wireS    = out_clus_v[ newPl ].SigmaEndWire();
	      auto end_tick     = out_clus_v[ newPl ].EndTick();
	      auto end_tickS    = out_clus_v[ newPl ].SigmaEndTick();
	      auto end_charge   = out_clus_v[ newPl ].EndCharge();
	      auto end_angle    = out_clus_v[ newPl ].EndAngle();
	      auto end_open     = out_clus_v[ newPl ].EndOpeningAngle();
	      auto integral     = newclus->Integral()  + out_clus_v[ newPl ].Integral();
	      auto summedADC    = newclus->SummedADC() + out_clus_v[ newPl ].SummedADC();
	      auto nhits        = newclus->NHits() + out_clus_v[ newPl ].NHits();
	      out_clus_v[ newPl ] = recob::Cluster(start_wire, start_wireS, start_tick, start_tickS, start_charge, start_angle, start_open,
						   end_wire  , end_wireS  , end_tick  , end_tickS  , end_charge  , end_angle  , end_open,
						   integral, 0., summedADC, 0., nhits, 0., 0., 
						   out_clus_v[ newPl ].ID(), out_clus_v[ newPl ].View(), out_clus_v[ newPl ].Plane() );

	      std::cout << "New cluster has NHits " << nhits << " and " << out_clus_hit_assn_v[ newPl ].size() << " hits associated" << std::endl;

	    }// if the old cluster exists and it needs to be updated
	  }// if not out of bounds

	}// for all new clusters to be merged
      }// if we need to merge the shower-trunk
      else {
	vtxcandidate = ass_shr_v[0]->ShowerStart();
      } //if this shower is not to be merged
      
      // now find shower branches
      
      std::vector<size_t> branches_pfp_idx_v;
      FindShowerBranches(p, slice_pfp_v, branches_pfp_idx_v);

      std::cout << "Shower is being merged with " << branches_pfp_idx_v.size() << " branches " << std::endl;
      
      if (branches_pfp_idx_v.size() >= 1) {
	
	for (auto const& branch_idx : branches_pfp_idx_v) {
	  
	  // load new clusters to be added
	  auto new_clus_v = slice_pfp_v[branch_idx].get<recob::Cluster>();
	  std::cout << "Will be adding " << new_clus_v.size() << " clusters!" << std::endl;
	  // beacuse this is the trunk of the shower, we need to modify the cluster
	  // start/end points to match these new ones
	  for (size_t c=0; c < new_clus_v.size(); c++) {
	    auto newclus = new_clus_v.at(c);
	    auto newPl = newclus->Plane().Plane;
	    if ( (newPl < out_clus_v.size()) && (newPl >= 0) ) {
	      
	      // store extra hits associated
	      

	      std::cout << "PRE @ plane " << newPl << " hits are " << out_clus_hit_assn_v[newPl].size() << std::endl;
	      // clus -> hit associations
	      const auto& new_clus_pxy = clus_proxy[newclus.key()];
	      auto new_clus_hit_v = new_clus_pxy.get<recob::Hit>();
	      for (size_t i=0; i < new_clus_hit_v.size(); i++) {
		out_clus_hit_assn_v[ newPl ].push_back( new_clus_hit_v[i] );
	      }// for all hits
	      std::cout << "POST @ plane " << newPl << " hits are " << out_clus_hit_assn_v[newPl].size() << std::endl;
	      
	      // if the old cluster is garbage, simply load the new one
	      if (out_clus_v[ newPl ].NHits() == 0) {
		out_clus_v[ newPl ] = (*newclus);
	      }
	      else {
		auto start_wire   = out_clus_v[ newPl ].StartWire();
		auto start_wireS  = out_clus_v[ newPl ].SigmaStartWire();
		auto start_tick   = out_clus_v[ newPl ].StartTick();
		auto start_tickS  = out_clus_v[ newPl ].SigmaStartTick();
		auto start_charge = out_clus_v[ newPl ].StartCharge();
		auto start_angle  = out_clus_v[ newPl ].StartAngle();
		auto start_open   = out_clus_v[ newPl ].StartOpeningAngle();
		auto end_wire     = newclus->EndWire();
		auto end_wireS    = newclus->SigmaEndWire();
		auto end_tick     = newclus->EndTick();
		auto end_tickS    = newclus->SigmaEndTick();
		auto end_charge   = newclus->EndCharge();
		auto end_angle    = newclus->EndAngle();
		auto end_open     = newclus->EndOpeningAngle();
		auto integral     = newclus->Integral()  + out_clus_v[ newPl ].Integral();
		auto summedADC    = newclus->SummedADC() + out_clus_v[ newPl ].SummedADC();
		auto nhits        = newclus->NHits() + out_clus_v[ newPl ].NHits();
		out_clus_v[ newPl ] = recob::Cluster(start_wire, start_wireS, start_tick, start_tickS, start_charge, start_angle, start_open,
						     end_wire  , end_wireS  , end_tick  , end_tickS  , end_charge  , end_angle  , end_open,
						     integral, 0., summedADC, 0., nhits, 0., 0., 
						     out_clus_v[ newPl ].ID(), out_clus_v[ newPl ].View(), out_clus_v[ newPl ].Plane() );
		
		std::cout << "New cluster has NHits " << nhits << " and " << out_clus_hit_assn_v[ newPl ].size() << " hits associated" << std::endl;
		
	      }// if the old cluster exists and it needs to be updated
	    }// if not out of bounds	      
	  }// for all clusters
	  
	}// for all branches to be merged
	
      }// if there are branches to be merged
      
      
      PFParticle_v->emplace_back( outshrpfp );
      
      // save vertex to be associated to this shower
      xyz[0] = vtxcandidate.X();
      xyz[1] = vtxcandidate.Y();
      xyz[2] = vtxcandidate.Z();
      recob::Vertex vtx(xyz);
      Vertex_v->emplace_back(vtx);

      if (largestShower == true) {
	_Pshr_q = out_clus_v[2].Integral();
	_Pshr_x = xyz[0];
	_Pshr_y = xyz[1];
	_Pshr_z = xyz[2];
      }
      
      art::Ptr<recob::PFParticle> const PFParticlePtr = PFParticlePtrMaker(PFParticle_v->size()-1);
      art::Ptr<recob::Vertex>     const VtxPtr        = VtxPtrMaker(Vertex_v->size()-1);

      // step 1 save metadata association
      const art::Ptr<recob::PFParticle> PFPPtr(pfp_h, pfp_pxy.index() );
      PFP_Meta_assn_v->addSingle( PFParticlePtr, pfp_pxy.get<larpandoraobj::PFParticleMetadata>()[0] );

      // step 2 save vertex association
      PFP_Vtx_assn_v->addSingle( PFParticlePtr, VtxPtr);

      // step 3 save cluster association
      for (size_t c=0; c < out_clus_v.size(); c++) {
	// if the cluster is meaningful
	if (out_clus_v[c].NHits() == 0) continue;
	Cluster_v->emplace_back( out_clus_v[c] );
	art::Ptr<recob::Cluster> const ClsPtr = ClsPtrMaker(Cluster_v->size()-1);
	PFP_Cls_assn_v->addSingle( PFParticlePtr, ClsPtr);
	// grab associated hit
	for (size_t h=0; h < out_clus_hit_assn_v[ c ].size(); h++) 
	  Cls_Hit_assn_v->addSingle( ClsPtr, out_clus_hit_assn_v[ c ].at( h ) );
      }// for all clusters
      
    }// if shower-like PFP

  }// for all slice PFPs

  if (slice_pfp_v.size())
    _tree->Fill();

  e.put(std::move(PFParticle_v));
  e.put(std::move(Vertex_v));
  e.put(std::move(Cluster_v));

  e.put(std::move(PFP_Meta_assn_v));
  e.put(std::move(PFP_Vtx_assn_v));
  e.put(std::move(PFP_Cls_assn_v));
  
  e.put(std::move(Cls_Hit_assn_v));

}

void ShowerMerger::beginJob()
{
  // Implementation of optional member function here.
}

void ShowerMerger::endJob()
{
  // Implementation of optional member function here.
}

bool ShowerMerger::FindShowerTrunk(const size_t shr_pfp_idx,
				   const std::vector<searchingfornues::ProxyPfpElem_t>& slice_pfp_v,
				   TVector3& vtxcandidate,
				   size_t& trktrunk_pfp_idx) {

  // do any of the slices have a track-like segment aligned with the shower?

  if (slice_pfp_v.size() <= shr_pfp_idx) return false;

  auto ass_shr_v = slice_pfp_v[shr_pfp_idx].get<recob::Shower>();

  if (ass_shr_v.size() != 1) return false;

  std::cout << "\t shower @ idx " << shr_pfp_idx << std::endl;
  
  auto shr = ass_shr_v[0];
  // get shower 3D direction and starting point
  auto shrVtx = shr->ShowerStart(); // xyz coordinate [TVector3]
  auto shrDir = shr->Direction(); // unit vector [TVector3]

  // keep track of possible matching track
  // track with best alignment in dot product is the candidate
  float bestdot = 0.;
  //float bestdist = 1e6;

  for (size_t p=0; p < slice_pfp_v.size(); p++) {

    // skip the pfparticle associated to the shower itself
    if (p == shr_pfp_idx) continue;

    // shower trunk needs to be track-like to be merged by this algorithm
    // filter out non-track like PFPs

    auto ass_trk_v = slice_pfp_v[p].get<recob::Track>();
    if (ass_trk_v.size() != 1) continue;

    auto trk = ass_trk_v[0];

    // track start direction / start point / end point
    TVector3 trkDir(trk->StartDirection().X(), trk->StartDirection().Y(), trk->StartDirection().Z());
    TVector3 trkVtx(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
    TVector3 trkEnd(trk->End().X(),   trk->End().Y(),   trk->End().Z());

    // compatibility?
    double TrkShrDot = trkDir.Dot(shrDir);
    // distance between shower start and track start/end
    // whicever is smallest
    double TrkShrDist = (shrVtx - trkVtx).Mag();
    if ( ((shrVtx - trkEnd).Mag() ) < TrkShrDot )
      TrkShrDist = (shrVtx - trkEnd).Mag();

    std::cout << "Compare new track..." << std::endl;
    std::cout << "TrkShr Dist is " << TrkShrDist << std::endl;
    std::cout << "TrkShr Dot  is " << TrkShrDot << std::endl;
    std::cout << std::endl;

    if ( (TrkShrDist < fShrVtxTrkDistMax) && (fabs(TrkShrDot) > fShrTrkDotMin) ) {

      // if this is the most aligned track
      if ( fabs(TrkShrDot) > fabs(bestdot) ) {
	trktrunk_pfp_idx = p;
	bestdot  = TrkShrDot;
	//bestdist = TrkShrDist;
	// vertex candidate is track start/end depending on sign
	// of dot-product.
	if (bestdot > 0) { vtxcandidate = trkVtx; }
	if (bestdot < 0) { vtxcandidate = trkEnd; }
	// also save index of track to be mergd
      }

    }// if agrees within user-defined specs
      
  }// for all PFParticles 
  

  // if bestdot != 0 -> means we found a compatible match
  if (bestdot == 0)
    return false;

  return true;
}// end FindShowerTrunk


bool ShowerMerger::FindShowerBranches(const size_t shr_pfp_idx,
				      const std::vector<searchingfornues::ProxyPfpElem_t>& slice_pfp_v,
				      std::vector<size_t>& branch_pfp_idx_v) {

  // do any slices lie in the cone of an upstream shower?

  if (slice_pfp_v.size() <= shr_pfp_idx) return false;

  auto ass_shr_v = slice_pfp_v[shr_pfp_idx].get<recob::Shower>();

  if (ass_shr_v.size() != 1) return false;

  std::cout << "\t shower @ idx " << shr_pfp_idx << std::endl;
  
  auto shr = ass_shr_v[0];
  // get shower 3D direction and starting point
  auto shrVtx = shr->ShowerStart(); // xyz coordinate [TVector3]
  auto shrDir = shr->Direction(); // unit vector [TVector3]

  for (size_t p=0; p < slice_pfp_v.size(); p++) {

    // skip the pfparticle associated to the shower itself
    if (p == shr_pfp_idx) continue;

    // find start point and direction of downstream PFParticles
    // regardless of whether shower/track like
    TVector3 BranchVtx;
    TVector3 BranchDir;

    if (slice_pfp_v[p].get<recob::Track>().size() == 1) {

      auto trk = slice_pfp_v[p].get<recob::Track>()[0];

      BranchDir = TVector3(trk->StartDirection().X(), trk->StartDirection().Y(), trk->StartDirection().Z());
      BranchVtx = TVector3(trk->Start().X(), trk->Start().Y(), trk->Start().Z());

    }// if associated to a track
    else if (slice_pfp_v[p].get<recob::Shower>().size() == 1) {
      
      auto shrB = slice_pfp_v[p].get<recob::Shower>()[0];

      BranchVtx = shrB->ShowerStart();
      BranchDir = shrB->Direction();

    }// if associated to a shower

    // is the branch candidate vertex compatible with the shower cone?
    
    // (1) angle between shower direction and shr start -> branch start
    TVector3 Branch2ShrDir = (BranchVtx-shrVtx).Unit();
    double branchDot = Branch2ShrDir.Dot( shrDir );

    std::cout << " Branch finding dot  : " << branchDot << std::endl;


    if (branchDot < fMinBranchConeDot) continue;

    // (2) make sure branch starts enough downstream
    double BranchDist = (BranchVtx-shrVtx).Mag();

    std::cout << " Branch finding dist : " << BranchDist << std::endl;

    if ( BranchDist < fMinBranchVtxDist ) continue;

    // made it this far -> add candidate branch
    branch_pfp_idx_v.push_back( p );
      
  }// for all PFParticles 
  

  // if bestdot != 0 -> means we found a compatible match
  if (branch_pfp_idx_v.size() == 0)
    return false;
  
  return true;
}// end FindShowerTrunk

void ShowerMerger::BuildPFPMap(const searchingfornues::ProxyPfpColl_t& pfp_pxy_col) {
  
  _pfpmap.clear();

  unsigned int p=0;
  for (const auto& pfp_pxy : pfp_pxy_col) {
    _pfpmap[pfp_pxy->Self()] = p;
    p++;
  }

  return;
}// BuildPFPMap

void ShowerMerger:: AddDaughters(const searchingfornues::ProxyPfpElem_t& pfp_pxy,
				 const searchingfornues::ProxyPfpColl_t& pfp_pxy_col,
				 std::vector<searchingfornues::ProxyPfpElem_t>& slice_v) {
  
  auto daughters = pfp_pxy->Daughters();
  
  slice_v.push_back(pfp_pxy);
  
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


void ShowerMerger::SaveTruth(art::Event const& e) {

  // load MCTruth
  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >(fMCTproducer);

  auto mct      = mct_h->at(0);
  auto neutrino = mct.GetNeutrino();
  auto nu       = neutrino.Nu();

  _nu_e  = nu.Trajectory().E(0);
  _vtx_x = nu.EndX();
  _vtx_y = nu.EndY();
  _vtx_z = nu.EndZ();
  _vtx_t = nu.T();

  _nelec   = 0;

  size_t npart = mct.NParticles();
  for (size_t i=0; i < npart; i++){

    auto const& part = mct.GetParticle(i);

    // if electron
    if ( (std::abs(part.PdgCode()) == 11) and (part.StatusCode() == 1) ){
      _nelec += 1;
      _elec_e = part.Momentum(0).E();
    }// if electron
  }// for all MCParticles
  
  searchingfornues::ApplyDetectorOffsets(_vtx_t,_vtx_x,_vtx_y,_vtx_z,_xtimeoffset,_xsceoffset,_ysceoffset,_zsceoffset);

  return;
}

DEFINE_ART_MODULE(ShowerMerger)
