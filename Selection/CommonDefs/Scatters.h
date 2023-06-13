#ifndef SCATTERSFUNCS_H
#define SCATTERSFUNCS_H

namespace searchingfornues
{
    std::vector< art::Ptr<simb::MCParticle>> GetDaughters(const art::Ptr<simb::MCParticle> &particle, const std::map<int, art::Ptr<simb::MCParticle> > &mcParticleMap)
    {
        std::vector< art::Ptr<simb::MCParticle>> daughters;
        for (int i = 0; i < particle->NumberDaughters(); ++i)
        {
            const auto daughterIter = mcParticleMap.find(particle->Daughter(i));
            if (daughterIter != mcParticleMap.end()) daughters.push_back(daughterIter->second); // Check needed because MCParticles outside of the cryostat are dropped and the hierarchy gets truncated
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

            if (process == "hadElastic")
            {
                nElastic++;
            }
            else if (process.find("Inelastic") != std::string::npos)
            {
                // If the daughter has the same PDG as the incident particle, then assume it's the same particle before & after the scatter
                if (daughter->PdgCode() != mcParticle->PdgCode()) 
                {
                    // nInelastic++;
                    continue;
                }

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

//     TruthHelper::EndState TruthHelper::GetEndState(const art::Ptr<simb::MCParticle> &particle, const MCParticleMap &mcParticleMap)
// {
//     MCParticleVector products;
//     const auto finalParticleMomentum = particle->Momentum(std::max(static_cast<unsigned int>(0), particle->NumberTrajectoryPoints() - 2)).Vect();

//     // ATTN for now, just return OTHER for all non-pions, we only use the end-state of the pions
//     if (particle->PdgCode() != 211)
//         return finalParticleMomentum;

//     // Follow the scatters to make sure we are really looking at the correct MCParticle
//     // ScatterVector scatters;
//     // art::Ptr<simb::MCParticle> scatteredParticle;
//     // TruthHelper::FollowScatters(particle, mcParticleMap, scatters, scatteredParticle);
//     // if (particle != scatteredParticle)
//     //     throw cet::exception("TruthHelper::GetEndState") << " - Input MCParticle undergos a scatter before reaching it's end-state" << std::endl;

//     bool hasPi0 = false;
//     bool hasDecayMuon = false;
//     bool hasDecayMuonNeutrino = false;

//     for (const auto &daughter : GetDaughters(particle, mcParticleMap))
//     {
//         // Ignore ionisation electrons
//         if (daughter->PdgCode() == 11 && daughter->Process() == "hIoni")
//             continue;

//         // Ignore nuclei from elastic scatters
//         if (daughter->Process() == "hadElastic")
//             continue;

//         // Treat everything else as an end-state interaction product
//         products.push_back(daughter);

//         if (daughter->PdgCode() == 111 && daughter->Process() == "pi+Inelastic")
//             hasPi0 = true;

//         if (daughter->PdgCode() == -13 && daughter->Process() == "Decay")
//             hasDecayMuon = true;

//         if (daughter->PdgCode() == 14 && daughter->Process() == "Decay")
//             hasDecayMuonNeutrino = true;
//     }

//     // Work out the end-state type
//     if (products.empty())
//     {
//         type = EndState::Type::None;
//     }
//     else if (hasDecayMuon && hasDecayMuonNeutrino && products.size() == 2)
//     {
//         type = EndState::Type::DecayToMuon;
//     }
//     else if (particle->EndProcess() == "pi+Inelastic")
//     {
//         type = hasPi0 ? EndState::Type::Pi0ChargeExchange : EndState::Type::InelasticAbsorption;
//     }

//     return TruthHelper::EndState(type, finalParticleMomentum, products);
// }


} // namespace searchingfornues

#endif