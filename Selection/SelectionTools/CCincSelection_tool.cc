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

  std::cout << "[CCincSelection::selectEvent] Number of Pfp in slice: " << pfp_pxy_v.size() << std::endl;
  int electron_candidate_index = -1;
  float electron_candidate_E = 0;
  TVector3 nuvtx;

  // Loop over Pfp in slice
  for (size_t i_pfp = 0; i_pfp < pfp_pxy_v.size(); i_pfp++)
  {
    auto const &pfp_pxy = pfp_pxy_v.at(i_pfp);
    auto PDG = fabs(pfp_pxy->PdgCode());

    // Pfp is the neutrino
    if (PDG == 12 || PDG == 14)
    {
      double xyz[3] = {};
      auto vtx = pfp_pxy.get<recob::Vertex>();
      if (vtx.size() != 1)
      {
        std::cout << "[CCincSelection::selectEvent] Found neutrino PFP without associated vertex" << std::endl;
      }
      else
      {
        vtx.at(0)->XYZ(xyz);
        nuvtx.SetXYZ(xyz[0], xyz[1], xyz[2]);
      }
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
          auto shr = pfp_pxy.get<recob::Shower>().at(0);
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
    m_electron_candidate = FillElectronCandidate(e, pfp_pxy_v.at(m_shrPfpId));
    auto shr = pfp_pxy_v.at(m_shrPfpId).get<recob::Shower>().at(0);
    m_shrDistance = (shr->ShowerStart() - nuvtx).Mag();
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

bool CCincSelection::isContained(const double x[3], float tolerance) const
{
  art::ServiceHandle<geo::Geometry> geo;
  std::vector<double> bnd = {
      0., 2. * geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(),
      0., geo->DetLength()};

  bool is_x = x[0] > (bnd[0] + tolerance) && x[0] < (bnd[1] - tolerance);
  bool is_y = x[1] > (bnd[2] + tolerance) && x[1] < (bnd[3] - tolerance);
  bool is_z = x[2] > (bnd[4] + tolerance) && x[2] < (bnd[5] - tolerance);

  return is_x && is_y && is_z;
}

bool CCincSelection::FillElectronCandidate(art::Event const &e,
                                           const ProxyPfpElem_t &pfp_pxy)
{
  // Get the trackscore
  m_shrScore = searchingfornues::GetTrackShowerScore(pfp_pxy);

  // Get number of hits
  ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e,
                                                                                        fCLSproducer,
                                                                                        proxy::withAssociated<recob::Hit>(fCLSproducer));
  auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();
  for (auto ass_clus : clus_pxy_v)
  {
    // get cluster proxy
    const auto &clus = clus_proxy[ass_clus.key()];
    auto clus_hit_v = clus.get<recob::Hit>();
    auto nhits = clus_hit_v.size();

    m_shrHits += nhits;
    if (clus->Plane().Plane == 0)
    {
      m_shrHitsU += nhits;
    }
    else if (clus->Plane().Plane == 1)
    {
      m_shrHitsV += nhits;
    }
    else if (clus->Plane().Plane == 2)
    {
      m_shrHitsY += nhits;
    }
  }
  // Get number of spacepoints
  auto sps_pxy_v = pfp_pxy.get<recob::SpacePoint>();
  m_shrSps = sps_pxy_v.size();

  // Shower Principal Components
  auto &pca_pxy_v = pfp_pxy.get<recob::PCAxis>();
  if (pca_pxy_v.size() > 0)
  {
    m_shrPca0 = pca_pxy_v[0]->getEigenValues()[0];
    m_shrPca1 = pca_pxy_v[1]->getEigenValues()[1];
    m_shrPca2 = pca_pxy_v[2]->getEigenValues()[2];
  }

  // Get shower fields
  auto shr = pfp_pxy.get<recob::Shower>().at(0);
  m_shrEnergy = shr->Energy()[2] / 1000;
  m_shrStartX = shr->ShowerStart().X();
  m_shrStartY = shr->ShowerStart().Y();
  m_shrStartZ = shr->ShowerStart().Z();
  m_shrDedxU = shr->dEdx()[0];
  m_shrDedxV = shr->dEdx()[1];
  m_shrDedxY = shr->dEdx()[2];
  m_shrTheta = shr->Direction().Theta();
  m_shrPhi = shr->Direction().Phi();
  m_shrOpenangle = shr->OpenAngle();
  m_shrPx = shr->Direction().X();
  m_shrPy = shr->Direction().Y();
  m_shrPz = shr->Direction().Z();

  // Get calibration factor for energy
  art::ValidHandle<std::vector<recob::PFParticle>> inputPfParticle = e.getValidHandle<std::vector<recob::PFParticle>>(fCLSproducer);
  art::ValidHandle<std::vector<recob::SpacePoint>> inputSpacePoint = e.getValidHandle<std::vector<recob::SpacePoint>>(fCLSproducer);
  auto assocSpacePoint = std::unique_ptr<art::FindManyP<recob::SpacePoint>>(new art::FindManyP<recob::SpacePoint>(inputPfParticle, e, fCLSproducer));
  auto assocSpacePointHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSpacePoint, e, fCLSproducer));
  std::vector<art::Ptr<recob::SpacePoint>> spcpnts = assocSpacePoint->at(m_shrPfpId);
  std::vector<float> cali_corr(3);
  searchingfornues::getCali(spcpnts, *assocSpacePointHit, cali_corr);
  m_shrEnergyCali = m_shrEnergy * cali_corr[2];
  // Get calibration factor for dedx
  std::vector<float> dqdx_cali(3);
  searchingfornues::getDQdxCali(shr, dqdx_cali);
  m_shrDedxYCali = m_shrDedxY * dqdx_cali[2];
  m_shrDedxVCali = m_shrDedxV * dqdx_cali[1];
  m_shrDedxUCali = m_shrDedxU * dqdx_cali[0];

  // Get track and PID information
  auto trk_v = pfp_pxy.get<recob::Track>();
  auto ntrk = trk_v.size();
  if (ntrk == 1)
  {
    searchingfornues::ProxyPIDColl_t const &pid_proxy = proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer,
                                                                                                        proxy::withAssociated<anab::ParticleID>(fPIDproducer));

    auto trk = trk_v.at(0);
    auto trkpxy2 = pid_proxy[trk.key()];
    auto pidpxy_v = trkpxy2.get<anab::ParticleID>();

    m_shrBragg_p = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 2),
                            searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 2));

    m_shrBragg_mu = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 2),
                             searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 2));

    m_shrBragg_mip = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 2);
    m_shrPidchipr = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 0);
    m_shrPidchimu = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 0);

    m_shrTrkLength = trk->Length();
    m_shrRangeMomMuon = m_trkmom.GetTrackMomentum(m_shrTrkLength, 13);
    const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle = e.getValidHandle<std::vector<recob::MCSFitResult>>("pandoraMCSMu");
    const recob::MCSFitResult &mcsMu = MCSMu_handle->at(trk.key());
    m_shrMCSMomMuon = mcsMu.fwdMomentum();
    double trkEnd[3] = {trk->End().X(), trk->End().Y(), trk->End().Z()};
    m_shrTrkContained = isContained(trkEnd, m_fidVolContained);
  }
  std::cout << "[CCincSelection::selectEvent] Shower energy: " << m_shrEnergy << std::endl;
  std::cout << "[CCincSelection::selectEvent] Calibration factor: " << cali_corr[2] << std::endl;

  return true;
}

} // namespace selection
