#ifndef SELECTIONTOOLBASE_H
#define SELECTIONTOOLBASE_H
////////////////////////////////////////////////////////////////////////
//
// Class:       IHitEfficiencyHistogramTool
// Module Type: tool
// File:        IHitEfficiencyHistogramTool.h
//
//              This provides an interface for tools which do histogramming
//              of various quantities associated to recob::Hit objects
//
// Created by David Caratelli (davidc@fnal.gov) on January 30 2019
//
////////////////////////////////////////////////////////////////////////

// art TOOLS
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Event.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"

namespace selection {

class SelectionToolBase {

public:

    /**
     *  @brief  Virtual Destructor
     */
    virtual ~SelectionToolBase() noexcept = default;
    
    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    void configure(const fhicl::ParameterSet&){};

    /**
     * @brief Selection function
     *
     * @param art::Event event record for selection
     */
    virtual bool selectEvent(art::Event const& e,
			     const std::vector<art::Ptr<recob::Track>  >& trkptr_v,
			     const std::vector<art::Ptr<recob::Shower> >& shrptr_v) = 0;

};

} // selection namespace

#endif
