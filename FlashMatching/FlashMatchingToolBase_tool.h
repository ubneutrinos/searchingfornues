/**
 *  @file   larpandora/LArPandoraEventBuilding/SliceIdTools/FlashNeutrinoId_tool.h
 *
 *  @brief  header of the flash based neutrino id tool
 */

#ifndef FLASHMATCHINGTOOLBASE_H
#define FLASHMATCHINGTOOLBASE_H

#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/Geometry.h"

#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"

#include "ubreco/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "ubreco/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"

#include "Objects/CartesianVector.h"

#include "TFile.h"
#include "TTree.h"
#include <iostream>

#include <numeric>

namespace flashmatch {
  
  /**
   *  @brief  Neutrino ID tool that selects the most likely neutrino slice using PMT information
   */
  class FlashMatchingToolBase
  {
  public:

    virtual ~FlashMatchingToolBase() noexcept = default;
    
    float ClassifySlice() { return 0; }

    /**
     *  @brief  Classify slices as neutrino or cosmic
     *
     *  @param  slices the input vector of slices to classify
     *  @param  evt the art event
     */
    virtual float ClassifySlice(const art::Event &evt,
				const std::vector< art::Ptr<recob::PFParticle> > &pfp_v,
				const std::vector< std::vector<art::Ptr< recob::SpacePoint> > > &spacepoint_v_v,
				const std::vector< std::vector<art::Ptr<recob::Hit> > > &hit_v_v,
				std::vector<float>& recospectrum,
				std::vector<float>& hypospectrum) = 0;


    virtual float ClassifyTrack(const art::Event &evt,
				const std::vector<art::Ptr< recob::SpacePoint> > &spacepoint_v,
				const std::vector<art::Ptr<recob::Hit> > &hit_v,
				std::vector<float>& recospectrum,
				std::vector<float>& hypospectrum) = 0;
    
  protected:

  };
  
}// end namespace


#endif
