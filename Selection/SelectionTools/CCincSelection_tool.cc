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
      // Decide which index corresponds to the electron candidate
      // Decide if the event is fiducial by looking at all daughters as tracks
    }
  }
  // Fill the information of the electron candidate

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

} // namespace selection
