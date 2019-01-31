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
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardata/Utilities/FindManyInChainP.h"

// selection tool
#include "SelectionTools/SelectionToolBase.h"

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
  void beginJob() override;
  void endJob() override;

private:

  art::InputTag fPFPproducer;
  art::InputTag fSHRproducer;
  art::InputTag fTRKproducer;
  bool fVerbose;

  // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
  std::map<unsigned int, unsigned int> _pfpmap;

  // selection tool
  std::unique_ptr<::selection::SelectionToolBase> _selectionTool;

  /**
   * @brief function to builf a map linking PFParticle index to Self() attribute
   *
   * @input handle to event pfparticle record
   */
  void BuildPFPMap(const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h);

  
  /**
   * @brief print PFParticle metadata information
   */
  void printPFParticleMetadata(const art::Ptr<recob::PFParticle>& pfp_ptr,
			       const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList);

  /**
   * @brief build PFParticle hierarchy (i.e. slice) from parent [recursive function]
   *
   * @input pfp_ptr : parent pfparticle for which to add daughters
   * @input pfp_h   : evnt PFP handle
   * @input slice_v : passed by reference, slice containing all PFParticles in hierarchy
   *
   */
  void AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  
		    const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h,
		    std::vector<art::Ptr<recob::PFParticle> > &slice_v);

};


NeutrinoSelection::NeutrinoSelection(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  
  fPFPproducer = p.get< art::InputTag > ("PFPproducer");
  fSHRproducer = p.get< art::InputTag > ("SHRproducer");
  fTRKproducer = p.get< art::InputTag > ("TRKproducer");
  fVerbose     = p.get< bool >          ("Verbose");

  // configure and construct Selection Tool
  const fhicl::ParameterSet& selection_pset = p.get<fhicl::ParameterSet>("SelectionTool");  
  _selectionTool = art::make_tool<::selection::SelectionToolBase>(selection_pset);

}

void NeutrinoSelection::analyze(art::Event const& e)
{

  if (fVerbose) { std::cout << "new event : [run,event] : [" << e.run() << ", " << e.event() << "]" << std::endl; }
    
    // grab PFParticles in event
    auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);

    // build PFParticle map  for this event
    BuildPFPMap(pfp_h);
  
    // grab tracks associated with PFParticles
    art::FindManyP<recob::Track>  pfp_track_assn_v (pfp_h, e, fTRKproducer);
    // grab showers associated with PFParticles
    art::FindManyP<recob::Shower> pfp_shower_assn_v(pfp_h, e, fSHRproducer);
    // grab associated metadata
    art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfp_h, e, fPFPproducer);

    // loop through PFParticles
    for (size_t p=0; p < pfp_h->size(); p++) {
      
      const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p );
      
      // get metadata for this PFP
      const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(p));
      
      //  find neutrino candidate
      if (pfp_ptr->IsPrimary() == false) continue;
      
      auto PDG = fabs(pfp_ptr->PdgCode());
      if ( (PDG == 12) || (PDG == 14) ) {
	
	if (fVerbose) printPFParticleMetadata(pfp_ptr,pfParticleMetadataList);

	// collect PFParticle hierarchy originating from this neutrino candidate
	std::vector<art::Ptr<recob::PFParticle> > slice_pfp_v;
	AddDaughters(pfp_ptr, pfp_h, slice_pfp_v);
	if (fVerbose) { std::cout << "This slice has " << slice_pfp_v.size() << " daughter PFParticles" << std::endl; }

	// create list of tracks and showers associated to this slice
	std::vector<art::Ptr<recob::Track>  > sliceTracks;
	std::vector<art::Ptr<recob::Shower> > sliceShowers;
	
	for (size_t pfpidx = 0; pfpidx < slice_pfp_v.size(); pfpidx++) {
	  auto pfpkey = slice_pfp_v.at(pfpidx).key();
	  // if there is a track associated to the PFParticle, add it
	  auto const& ass_trk_v = pfp_track_assn_v.at( pfpkey );
	  if (ass_trk_v.size() == 1) sliceTracks.push_back( ass_trk_v.at(0) );
	  // if there is a shower associated to the PFParticle, add it
	  auto const& ass_shr_v = pfp_shower_assn_v.at( pfpkey );
	  if (ass_shr_v.size() == 1) sliceShowers.push_back( ass_shr_v.at(0) );
	}// for all PFParticles in the slice

	// check that # of PFP and # trk + shr matches
	if ( (slice_pfp_v.size()-1) != (sliceTracks.size() + sliceShowers.size()) )
	  std::cout << "ERROR : there are "  << slice_pfp_v.size() << " PFP but " << sliceTracks.size() << " + " << sliceShowers.size() << " tracks + showers" << std::endl;
	
	// run selection on this slice
	bool selected = _selectionTool->selectEvent(e, sliceTracks, sliceShowers);
	if (fVerbose && selected) { std::cout << "SLICE was selected!" << std::endl; }
	
      }// if a neutrino PFParticle
      
    }// for all PFParticles
    
    return;
}

void NeutrinoSelection::beginJob()
{
  // Implementation of optional member function here.
}

void NeutrinoSelection::endJob()
{
  // Implementation of optional member function here.
}

void NeutrinoSelection::printPFParticleMetadata(const art::Ptr<recob::PFParticle>& pfp_ptr,
						const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList) {
  
  if (!pfParticleMetadataList.empty()) {
    
    for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j)
      {
	const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
	auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
	if (!pfParticlePropertiesMap.empty()) {
	  std::cout << " Found PFParticle " << pfp_ptr->Self() << " with: " << std::endl;
	  for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	    std::cout << "  - " << it->first << " = " << it->second << std::endl;
	  }
	}
      }
  }// if PFP metadata exists!
  return;
}

void NeutrinoSelection::BuildPFPMap(const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h) {
  
  _pfpmap.clear();

  for (unsigned int p=0; p < pfp_h->size(); p++)
    _pfpmap[pfp_h->at(p).Self()] = p;

  return;
}// BuildPFPMap

void NeutrinoSelection::AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  
				     const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h,
				     std::vector<art::Ptr<recob::PFParticle> > &slice_v) {
  
  auto daughters = pfp_ptr->Daughters();
  
  slice_v.push_back(pfp_ptr);
  
  if (fVerbose) std::cout << "\t PFP w/ PdgCode " << pfp_ptr->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;
  
  for(auto const& daughterid : daughters) {

    if (_pfpmap.find(daughterid) == _pfpmap.end()) {
      std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
      continue;
    }
    
    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, _pfpmap.at(daughterid) );
    
    AddDaughters(pfp_ptr, pfp_h, slice_v);
    
  }// for all daughters
  
  return;
}// AddDaughters


DEFINE_ART_MODULE(NeutrinoSelection)
