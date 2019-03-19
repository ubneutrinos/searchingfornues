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
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

#include "TTree.h"
#include <limits>

namespace selection {

using ProxyPfpColl_t = decltype(proxy::getCollection<std::vector<recob::PFParticle> >(
					      std::declval<art::Event>(),std::declval<art::InputTag>(),
                                              proxy::withAssociated<larpandoraobj::PFParticleMetadata>(std::declval<art::InputTag>()),
                                              proxy::withAssociated<recob::Cluster>(std::declval<art::InputTag>()),
                                              proxy::withAssociated<recob::Track>(std::declval<art::InputTag>()),
                                              proxy::withAssociated<recob::Vertex>(std::declval<art::InputTag>()),
                                              proxy::withAssociated<recob::Shower>(std::declval<art::InputTag>())) );
using ProxyPfpElem_t = ProxyPfpColl_t::element_proxy_t;

 using ProxyClusColl_t = decltype(proxy::getCollection<std::vector<recob::Cluster> >(
										     std::declval<art::Event>(),std::declval<art::InputTag>(),
										     proxy::withAssociated<recob::Hit>(std::declval<art::InputTag>()) ) );
 
 using ProxyClusElem_t = ProxyClusColl_t::element_proxy_t;

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
			     const std::vector<ProxyPfpElem_t>& pfp_pxy_v) = 0;

    /**
     * @brief set branches for TTree
     */
    virtual void setBranches(TTree* _tree) = 0;

    
    /**
     * @brief resetset TTree branches
     */
    virtual void resetTTree(TTree* _tree) = 0;


    /**
     * @brief set if data
     */
    void SetData(bool isdata) { fData = isdata; }


 protected:


    bool fData;


};

} // selection namespace

#endif
