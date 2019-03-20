#ifndef COMMONDEFS_TYPEDEFS_H
#define COMMONDEFS_TYPEDEFS_H

#include "art/Framework/Principal/Event.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

namespace searchingfornues {

  using ProxyPfpColl_t = decltype(proxy::getCollection<std::vector<recob::PFParticle> >(
					      std::declval<art::Event>(),std::declval<art::InputTag>(),
                                              proxy::withAssociated<larpandoraobj::PFParticleMetadata>(std::declval<art::InputTag>()),
                                              proxy::withAssociated<recob::Cluster>(std::declval<art::InputTag>()),
                                              proxy::withAssociated<recob::Track>(std::declval<art::InputTag>()),
                                              proxy::withAssociated<recob::Vertex>(std::declval<art::InputTag>()),
                                              proxy::withAssociated<recob::Shower>(std::declval<art::InputTag>()),
					      proxy::withAssociated<recob::Slice>(std::declval<art::InputTag>())) );
  using ProxyPfpElem_t = ProxyPfpColl_t::element_proxy_t;

  using ProxyClusColl_t = decltype(proxy::getCollection<std::vector<recob::Cluster> >(
										      std::declval<art::Event>(),std::declval<art::InputTag>(),
										      proxy::withAssociated<recob::Hit>(std::declval<art::InputTag>()) ) );
  using ProxyClusElem_t = ProxyClusColl_t::element_proxy_t;

}

#endif
