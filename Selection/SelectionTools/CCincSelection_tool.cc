#include "CCincSelection_tool.h"

namespace selection
{

//----------------------------------------------------------------------------
/// Event selection method.
///
/// Arguments:
///
/// e - Art Event.
/// pfp_pxy_v - Proxy for the pfpparticles in slice.
///
bool CCincSelection::selectEvent(art::Event const &e,
                                 const std::vector<ProxyPfpElem_t> &pfp_pxy_v)
{

  std::cout << "[CCincSelection::selectEvent] Number of Pfp ins slice: " << pfp_pxy_v.size() << std::endl;
  int electron_candidate_index = -1;
  float electron_candidate_E = 0;

  // Loop over Pfp in slice
  for (size_t i_pfp = 0; i_pfp < pfp_pxy_v.size(); i_pfp++)
  {
    auto const &pfp_pxy = pfp_pxy_v.at(i_pfp);
    auto PDG = fabs(pfp_pxy->PdgCode());

    // Pfp is the neutrino
    if (PDG == 12 || PDG == 14)
    {
      // Fill the field of the fiducial reconstructed vtx
    }
    // Pfp is a daughter
    else
    {
      // Decide if the event is fiducial by looking at all daughters as tracks
      // To Do

      // Decide which index corresponds to the electron candidate
      float trkscore = searchingfornues::GetTrackShowerScore(pfp_pxy);
      if (trkscore < m_trkScore)
      {
        size_t n_shw = pfp_pxy.get<recob::Shower>().size();
        if (n_shw == 1) // Our electron candidate needs to be reconstructed as a shower
        {
          auto shr = pfp_pxy.get<recob::Shower>()[0];
          if (electron_candidate_E < shr->Energy()[2])
          {
            electron_candidate_E = shr->Energy()[2];
            electron_candidate_index = i_pfp;
          }
        }
      }
    }
  }
  if (electron_candidate_index != -1)
  {
    std::cout << "[CCincSelection::selectEvent] Electron candidate found! ";
    std::cout << "pfp_id: " << electron_candidate_index;
    std::cout << ", collection plane energy: " << electron_candidate_E / 1000 << std::endl;
    // Fill the information of the electron candidate
    m_shrPfpId = electron_candidate_index;
    m_electron_candidate = FillElectronCandidate(pfp_pxy_v.at(m_shrPfpId));
  }
  else
  {
    std::cout << "[CCincSelection::selectEvent] No electron candidate found!" << std::endl;
  }
  return true;
}

bool CCincSelection::isFiducial(const double x[3]) const
{
  art::ServiceHandle<geo::Geometry> geo;
  std::vector<double> bnd = {
      0., 2. * geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(),
      0., geo->DetLength()};

  bool is_x = x[0] > (bnd[0] + m_fidvolXstart) && x[0] < (bnd[1] - m_fidvolXend);
  bool is_y = x[1] > (bnd[2] + m_fidvolYstart) && x[1] < (bnd[3] - m_fidvolYend);
  bool is_z = x[2] > (bnd[4] + m_fidvolZstart) && x[2] < (bnd[5] - m_fidvolZend);

  return is_x && is_y && is_z;
}

bool CCincSelection::FillElectronCandidate(const ProxyPfpElem_t &pfp_pxy)
{
  return true;
}

} // namespace selection
