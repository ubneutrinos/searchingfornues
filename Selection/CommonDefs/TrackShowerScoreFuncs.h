#ifndef TRACKSHOWERSCOREFUNCS_H
#define TRACKSHOWERSCOREFUNCS_H

namespace searchingfornues
{

float GetTrackShowerScore(const ProxyPfpElem_t &pfp_pxy)
{

  const auto &pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();

  if (pfParticleMetadataList.size() == 0)
    return 1;

  for (unsigned int j = 0; j < pfParticleMetadataList.size(); ++j)
  {

    const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
    auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
    if (!pfParticlePropertiesMap.empty())
    {
      for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
      {
        if (it->first == "TrackScore")
          return it->second;
      } // for map elements
    }   // if pfp metadata map not empty
  }     // for list

  return 1;
}

} // namespace searchingfornues

#endif
