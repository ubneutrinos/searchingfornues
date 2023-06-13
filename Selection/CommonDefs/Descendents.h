#ifndef DESCENDENTSFUNCS_H
#define DESCENDENTSFUNCS_H

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

// Recursive functions for finding descendents in the PFParticle hierarchy 
//////////////////////////////////////////////////////////////////////////

namespace searchingfornues
{
    lar_pandora::PFParticleVector GetDaughters(const art::Ptr<recob::PFParticle> &particle, const lar_pandora::PFParticleMap &pfParticleMap)
    {
        lar_pandora::PFParticleVector daughters;
        for (int i = 0; i < particle->NumDaughters(); ++i)
        {
            const auto daughterIter = pfParticleMap.find(particle->Daughter(i));
            if (daughterIter != pfParticleMap.end()) daughters.push_back(daughterIter->second);
        }

        return daughters;
    }

    void GetDownstreamParticles(const art::Ptr<recob::PFParticle> &particle, const lar_pandora::PFParticleMap &pfParticleMap, lar_pandora::PFParticleVector &downstreamParticles)
    {
        downstreamParticles.push_back(particle);
        for (const auto &daughter : GetDaughters(particle, pfParticleMap))
            GetDownstreamParticles(daughter, pfParticleMap, downstreamParticles);
    }

    unsigned int GetNDescendents(const art::Ptr<recob::PFParticle> &particle, const lar_pandora::PFParticleMap &pfParticleMap)
    {
        unsigned int nDescendents = 0u;
        lar_pandora::PFParticleVector downstreamParticles;
        GetDownstreamParticles(particle, pfParticleMap, downstreamParticles);
        for (const auto &downstreamParticle : downstreamParticles)
            if (downstreamParticle != particle) nDescendents++;
        return nDescendents;
    }
} // namespace searchingfornues

#endif