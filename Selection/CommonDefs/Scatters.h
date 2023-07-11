#ifndef SCATTERSFUNCS_H
#define SCATTERSFUNCS_H

// Recursive functions for finding scatters in the MCParticle hierarchy 
///////////////////////////////////////////////////////////////////////

namespace searchingfornues
{
    std::vector< art::Ptr<simb::MCParticle>> GetDaughters(const art::Ptr<simb::MCParticle> &particle, const std::map<int, art::Ptr<simb::MCParticle> > &mcParticleMap)
    {
        std::vector< art::Ptr<simb::MCParticle>> daughters;
        for (int i = 0; i < particle->NumberDaughters(); ++i)
        {
            const auto daughterIter = mcParticleMap.find(particle->Daughter(i));
            if (daughterIter != mcParticleMap.end()) daughters.push_back(daughterIter->second);
        }
        return daughters;
    }

void GetNScatters(const art::ValidHandle<std::vector<simb::MCParticle>> &mcp_h, const art::Ptr<simb::MCParticle> &mcParticle, art::Ptr<simb::MCParticle> &mcScatteredParticle, unsigned int &nElastic, unsigned int &nInelastic)
    {
        mcScatteredParticle = mcParticle;

        // Create the MCParticleMap
        std::map<int, art::Ptr<simb::MCParticle> > mcParticleMap;
        for (size_t d = 0; d < mcp_h->size(); d++)
        {
            const art::Ptr<simb::MCParticle> mcParticle(mcp_h, d);
            if (!mcParticleMap.emplace(mcParticle->TrackId(), mcParticle).second)
                throw cet::exception("searchingfornues::GetNScatters") << " - Found repeated MCParticle with TrackId = " << mcParticle->TrackId() << "." << std::endl;
        }

        // Loop over the daughters and count the number of elastic and inelastic scatters
        art::Ptr<simb::MCParticle> finalStateParticle;
        bool foundInelasticScatter = false;
        for (const auto &daughter : GetDaughters(mcParticle, mcParticleMap))
        {
            const auto& process = daughter->Process();

            if (process == "hadElastic") nElastic++;
            else if (process.find("Inelastic") != std::string::npos)
            {
                // If the daughter has the same PDG as the incident particle, then assume it's the same particle before & after the scatter
                if (daughter->PdgCode() != mcParticle->PdgCode()) continue;

                // If there are multiple particles from the inelastic collision, then don't treat it as a scatter
                if (foundInelasticScatter) 
                {
                    // nInelastic = 0;
                    foundInelasticScatter = false;
                    break;
                }
                finalStateParticle = daughter;
                foundInelasticScatter = true;
            }
        }
        if (foundInelasticScatter)
        {
            nInelastic++;
            GetNScatters(mcp_h, finalStateParticle, mcScatteredParticle, nElastic, nInelastic);
        }
    }

} // namespace searchingfornues

#endif