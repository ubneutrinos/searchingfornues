#ifndef SELECTION_CC0PINPSELECTION_CXX
#define SELECTION_CC0PINPSELECTION_CXX

#include <iostream>
#include "SelectionToolBase.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "../CommonDefs/Typedefs.h"
#include "../CommonDefs/PIDFuncs.h"
#include "larcore/Geometry/Geometry.h"
// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "canvas/Persistency/Common/FindManyP.h"

namespace selection
{
////////////////////////////////////////////////////////////////////////
///
/// Class:       CC0piNpSelection
/// File:        CC0piNpSelection_tool.cc
///
///              Selection of electron neutrinos with 0 pions and at least on proton in the final state
///
/// Created by Stefano Roberto Soleti (srsoleti@fnal.gov) on 04/11/2019
///
////////////////////////////////////////////////////////////////////////

class CC0piNpSelection : public SelectionToolBase
{

public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset List of parameters in the FCL file
     */
    CC0piNpSelection(const fhicl::ParameterSet &pset);

    /**
     *  @brief  Destructor
     */
    ~CC0piNpSelection(){};

    // provide for initialization
    void configure(fhicl::ParameterSet const &pset);

    /**
     * @brief Selection function
     */
    bool selectEvent(art::Event const &e,
                     const std::vector<ProxyPfpElem_t> &pfp_pxy_v);

    /**
     * @brief Set branches for TTree
     */
    void setBranches(TTree *_tree);

    /**
     * @brief Reset TTree branches
     */
    void resetTTree(TTree *_tree);

    /**
     * @brief Check if inside fiducial volume
     */
    bool isFiducial(const double x[3]) const;

private:
    trkf::TrackMomentumCalculator _trkmom;
    TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
    TParticlePDG *electron = TDatabasePDG::Instance()->GetParticle(11);
    TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);

    float fTrkShrscore;
    float fFidvolZstart;
    float fFidvolZend;
    float fFidvolYstart;
    float fFidvolYend;
    float fFidvolXstart;
    float fFidvolXend;

    art::InputTag fCLSproducer;
    art::InputTag fPIDproducer;
    art::InputTag fTRKproducer;
    art::InputTag fBacktrackTag;
    art::InputTag fHproducer;
    art::InputTag fMCRproducer;
    art::InputTag fMCPproducer;
    art::InputTag fCALproducer;

    unsigned int _n_showers_contained; /**< Number of showers with a starting point within the fiducial volume */
    unsigned int _n_tracks_contained; /**< Number of tracks fully contained in the fiducial volume */
    unsigned int _shr_hits_max; /**< Number of hits of the leading shower */
    unsigned int _trk_hits_max; /**< Number of hits of the leading track */
    unsigned int _shr_hits_tot; /**< Total number of shower hits */
    unsigned int _trk_hits_tot; /**< Total number of track hits */
    unsigned int _trk_hits_y_tot; /**< Total number of track hits on the Y plane */
    unsigned int _trk_hits_v_tot; /**< Total number of track hits on the V plane */
    unsigned int _trk_hits_u_tot; /**< Total number of track hits on the U plane */
    unsigned int _shr_hits_y_tot; /**< Total number of shower hits on the Y plane */
    unsigned int _shr_hits_v_tot; /**< Total number of shower hits on the V plane */
    unsigned int _shr_hits_u_tot; /**< Total number of shower hits on the U plane */
    float _shr_energy;
    float _shr_energy_tot;
    float _shr_dedx_Y;
    float _shr_dedx_V;
    float _shr_dedx_U;
    float _shr_distance;
    float _tksh_distance;
    float _tksh_angle;
    float _shr_score;
    float _shr_theta;
    float _shr_phi;
    float _shr_px;
    float _shr_py;
    float _shr_pz;
    float _shr_openangle;

    size_t _shr_pfp_id;

    float _trk_len;
    float _trk_energy;
    float _trk_energy_tot;
    float _trk_distance;
    float _trk_theta;
    float _trk_phi;
    size_t _trk_pfp_id;

    float _hits_ratio;
    float _trk_bragg_p;
    float _trk_bragg_mu;
    float _trk_bragg_mip;
    float _trk_pidchipr;
    float _trk_pidchipr_best;
    float _trk_pidchipr_worst;

    float _trk_pidchimu;
    float _trk_pida;
    float _trk_score;

    float _pt;
    float _p;

    float _shr_bkt_purity;
    float _shr_bkt_completeness;
    float _shr_bkt_E;
    int _shr_bkt_pdg;
    float _trk_bkt_purity;
    float _trk_bkt_completeness;
    float _trk_bkt_E;
    int _trk_bkt_pdg;

    float _dep_E;
    float _matched_E;

    float _shr_tkfit_start_x;
    float _shr_tkfit_start_y;
    float _shr_tkfit_start_z;
    float _shr_start_x;
    float _shr_start_y;
    float _shr_start_z;
    float _shr_tkfit_phi;
    float _shr_tkfit_theta;
    float _shr_tkfit_dedx_Y;
    float _shr_tkfit_dedx_U;
    float _shr_tkfit_dedx_V;
    unsigned int _shr_tkfit_nhits_Y;
    unsigned int _shr_tkfit_nhits_V;
    unsigned int _shr_tkfit_nhits_U;
    unsigned int _hits_outfv;
    float _contained_fraction;

    float _trk_energy_hits_tot;

    unsigned int _total_hits_y;
    float _extra_energy_y;
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
CC0piNpSelection::CC0piNpSelection(const fhicl::ParameterSet &pset)
{
    configure(pset);
}

bool CC0piNpSelection::isFiducial(const double x[3]) const
{

    art::ServiceHandle<geo::Geometry> geo;
    std::vector<double> bnd = {
        0., 2. * geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(),
        0., geo->DetLength()};

    bool is_x =
        x[0] > (bnd[0] + fFidvolXstart) && x[0] < (bnd[1] - fFidvolXend);
    bool is_y =
        x[1] > (bnd[2] + fFidvolYstart) && x[1] < (bnd[3] - fFidvolYend);
    bool is_z =
        x[2] > (bnd[4] + fFidvolZstart) && x[2] < (bnd[5] - fFidvolZend);

    return is_x && is_y && is_z;
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void CC0piNpSelection::configure(fhicl::ParameterSet const &pset)
{
    fTrkShrscore = pset.get<float>("TrkShrscore", 0.5);
    fCLSproducer = pset.get<art::InputTag>("CLSproducer", "pandora");
    fTRKproducer = pset.get<art::InputTag>("TRKproducer", "shrreco3dKalmanShower");
    fPIDproducer = pset.get<art::InputTag>("PIDproducer", "pandoracalipidSCE");
    fCALproducer = pset.get<art::InputTag>("CALproducer", "shrreco3dKalmanShowercali");
    fHproducer = pset.get<art::InputTag>("Hproducer", "gaushit");
    fMCRproducer = pset.get<art::InputTag>("MCRproducer", "mcreco");
    fBacktrackTag = pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
    fMCPproducer = pset.get<art::InputTag>("MCPproducer", "largeant");

    fFidvolZstart = pset.get<float>("FidvolZstart", 10);
    fFidvolZend = pset.get<float>("FidvolZend", 50);
    fFidvolYstart = pset.get<float>("FidvolYstart", 15);
    fFidvolYend = pset.get<float>("FidvolYend", 15);
    fFidvolXstart = pset.get<float>("FidvolXstart", 10);
    fFidvolXend = pset.get<float>("FidvolXend", 10);
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
bool CC0piNpSelection::selectEvent(art::Event const &e,
                                   const std::vector<ProxyPfpElem_t> &pfp_pxy_v)
{
    ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
                                                                                          proxy::withAssociated<recob::Hit>(fCLSproducer));
    searchingfornues::ProxyPIDColl_t const &pid_proxy = proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer,
                                                                                                        proxy::withAssociated<anab::ParticleID>(fPIDproducer));

    art::ValidHandle<std::vector<recob::PFParticle>> inputPfParticle = e.getValidHandle<std::vector<recob::PFParticle>>(fCLSproducer);
    art::ValidHandle<std::vector<recob::Slice>> inputSlice = e.getValidHandle<std::vector<recob::Slice>>(fCLSproducer);
    auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSlice, e, fCLSproducer));
    auto assocSlice = std::unique_ptr<art::FindManyP<recob::Slice>>(new art::FindManyP<recob::Slice>(inputPfParticle, e, fCLSproducer));
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfPartToMetadataAssoc(inputPfParticle, e, fCLSproducer);
    searchingfornues::ProxyCaloColl_t const *tkcalo_proxy = NULL;
    if (fTRKproducer != "")
    {
        tkcalo_proxy = new searchingfornues::ProxyCaloColl_t(proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer, proxy::withAssociated<anab::Calorimetry>(fCALproducer)));
    }
    // load backtrack information
    std::vector<searchingfornues::BtPart> btparts_v;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
    if (!fData)
    {
        const std::vector<sim::MCShower> &inputMCShower = *(e.getValidHandle<std::vector<sim::MCShower>>(fMCRproducer));
        const std::vector<sim::MCTrack> &inputMCTrack = *(e.getValidHandle<std::vector<sim::MCTrack>>(fMCRproducer));
        art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
        assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
        btparts_v = searchingfornues::initBacktrackingParticleVec(inputMCShower, inputMCTrack, *inputHits, assocMCPart);

        auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
        for (size_t p = 0; p < mcp_h->size(); p++)
        {
            auto mcp = mcp_h->at(p);
            if (!(mcp.Process() == "primary" &&
                  mcp.StatusCode() == 1))
            {
                continue;
            }

            auto PDG = mcp.PdgCode();
            if (fabs(PDG) == proton->PdgCode())
            {
                double ke = mcp.E() - proton->Mass();
                if (ke > 0.040)
                    _dep_E += ke;
            }
            if (fabs(PDG) == electron->PdgCode())
            {
                double ke = mcp.E() - electron->Mass();
                if (ke > 0.030)
                    _dep_E += ke;
            }
        }
    }

    for (unsigned int inpf = 0; inpf < inputPfParticle->size(); ++inpf)
    {
        art::Ptr<recob::PFParticle> npfp(inputPfParticle, inpf);
        bool isTheNeutrino = false;

        const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> &pfParticleMetadataList(pfPartToMetadataAssoc.at(inpf));
        if (!pfParticleMetadataList.empty())
        {
            for (unsigned int j = 0; j < pfParticleMetadataList.size(); ++j)
            {
                const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));

                const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
                for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
                {
                    if (it->first == "IsNeutrino" && it->second == 1)
                        isTheNeutrino = true;
                }
            }
        }
        if (npfp->IsPrimary() == false)
            continue;
        //  recob::Slices are unique for the neutrino slice and for clear cosmics. Ambiguous cosmics may have >1 primary PFParticle per slice
        //(i.e. the slice is not unique in terms of primary PFParticles)
        auto slices = assocSlice->at(npfp.key());
        if (slices.size() != 1)
        {
            std::cout << "WRONG!!! n slices = " << slices.size() << std::endl;
        }
        if (isTheNeutrino == true)
        {
            auto slicehit = assocSliceHit->at(slices[0].key());
            for (unsigned int isl = 0; isl < slicehit.size(); isl++)
            {
                const auto &slhit = slicehit.at(isl);
                int slhit_pl = slhit->WireID().Plane;
                if (slhit_pl == 2)
                {
                    _total_hits_y++;
                    _extra_energy_y += 1.01 * slhit->Integral() * 238 * 23.6e-6 / 0.6 / 1000;
                }
            }
        }
    }

    double nu_vtx[3] = {};
    TVector3 nuvtx;
    for (size_t i_pfp = 0; i_pfp < pfp_pxy_v.size(); i_pfp++)
    {
        auto PDG = fabs(pfp_pxy_v[i_pfp]->PdgCode());

        if (PDG == 12 || PDG == 14)
        {

            auto vtx = pfp_pxy_v[i_pfp].get<recob::Vertex>();
            if (vtx.size() != 1)
            {
                std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
            }
            else
            {
                vtx.at(0)->XYZ(nu_vtx);
                if (!isFiducial(nu_vtx))
                {
                    return false;
                }
                nuvtx.SetXYZ(nu_vtx[0], nu_vtx[1], nu_vtx[2]);
            }

            break;
        }
    }

    TVector3 total_p;
    TVector3 trk_vtx;
    TVector3 shr_vtx;
    TVector3 trk_p;
    TVector3 shr_p;

    for (size_t i_pfp = 0; i_pfp < pfp_pxy_v.size(); i_pfp++)
    {

        auto const &pfp_pxy = pfp_pxy_v.at(i_pfp);

        auto PDG = fabs(pfp_pxy->PdgCode());

        if (PDG == 12 || PDG == 14)
            continue;

        auto trkshrscore = searchingfornues::GetTrackShowerScore(pfp_pxy);

        for (const auto &shr : pfp_pxy.get<recob::Shower>())
        {
            if (trkshrscore < fTrkShrscore)
            {
                double shr_vertex[3] = {shr->ShowerStart().X(), shr->ShowerStart().Y(), shr->ShowerStart().Z()};
                auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();
                std::vector<art::Ptr<recob::Hit>> hit_v;

                if (!isFiducial(shr_vertex))
                {
                    for (auto ass_clus : clus_pxy_v)
                    {
                        const auto &clus = clus_proxy[ass_clus.key()];
                        auto clus_hit_v = clus.get<recob::Hit>();
                        _hits_outfv += clus_hit_v.size();
                    }
                    continue;
                }

                _n_showers_contained++;

                unsigned int shr_hits = 0;

                for (auto ass_clus : clus_pxy_v)
                {
                    // get cluster proxy
                    const auto &clus = clus_proxy[ass_clus.key()];
                    auto clus_hit_v = clus.get<recob::Hit>();
                    shr_hits += clus_hit_v.size();
                    for (const auto &hit : clus_hit_v)
                    {
                        hit_v.push_back(hit);
                    }
                    if (clus->Plane().Plane == 0)
                    {
                        _shr_hits_u_tot += clus_hit_v.size();
                    }
                    else if (clus->Plane().Plane == 1)
                    {
                        _shr_hits_v_tot += clus_hit_v.size();
                    }
                    else if (clus->Plane().Plane == 2)
                    {
                        _shr_hits_y_tot += clus_hit_v.size();
                    }
                }

                _shr_energy_tot += shr->Energy()[2] / 1000;
                _shr_hits_tot += shr_hits;
                if (shr_hits > _shr_hits_max)
                {
                    if (!fData)
                    {
                        if (clus_pxy_v.size() != 0)
                        {
                            float purity = 0., completeness = 0.;
                            int ibt = searchingfornues::getAssocBtPart(hit_v, assocMCPart, btparts_v, purity, completeness);
                            if (ibt >= 0)
                            {
                                auto &mcp = btparts_v[ibt];
                                auto PDG = mcp.pdg;
                                auto E = mcp.e;
                                _shr_bkt_pdg = PDG;
                                _shr_bkt_purity = purity;
                                _shr_bkt_completeness = completeness;
                                _shr_bkt_E = E;
                                _matched_E += E;
                            }
                        }
                    }
                    _shr_distance = (shr->ShowerStart() - nuvtx).Mag();
                    shr_vtx = shr->ShowerStart();
                    shr_p.SetXYZ(shr->Direction().X(), shr->Direction().Y(), shr->Direction().Z());
                    shr_p.SetMag(sqrt(pow(shr->Energy()[2] / 1000. + electron->Mass(), 2) - pow(electron->Mass(), 2)));
                    total_p += shr_p;
                    _shr_dedx_Y = shr->dEdx()[2];
                    _shr_dedx_V = shr->dEdx()[1];
                    _shr_dedx_U = shr->dEdx()[0];
                    _shr_start_x = shr->ShowerStart().X();
                    _shr_start_y = shr->ShowerStart().Y();
                    _shr_start_z = shr->ShowerStart().Z();
                    _shr_energy = shr->Energy()[2] / 1000; // GeV
                    _shr_pfp_id = i_pfp;
                    _shr_hits_max = shr_hits;
                    _shr_score = trkshrscore;
                    _shr_theta = shr->Direction().Theta();
                    _shr_phi = shr->Direction().Phi();
                    _shr_px = shr_p.X();
                    _shr_py = shr_p.Y();
                    _shr_pz = shr_p.Z();
                    _shr_openangle = shr->OpenAngle();

                    if (tkcalo_proxy == NULL)
                    {
                        _shr_tkfit_start_x = std::numeric_limits<float>::lowest();
                        _shr_tkfit_start_y = std::numeric_limits<float>::lowest();
                        _shr_tkfit_start_z = std::numeric_limits<float>::lowest();

                        _shr_tkfit_phi = std::numeric_limits<float>::lowest();
                        _shr_tkfit_theta = std::numeric_limits<float>::lowest();
                        _shr_tkfit_dedx_Y = std::numeric_limits<float>::lowest();
                        _shr_tkfit_dedx_V = std::numeric_limits<float>::lowest();
                        _shr_tkfit_dedx_U = std::numeric_limits<float>::lowest();

                        continue;
                    }

                    for (const searchingfornues::ProxyCaloElem_t &tk : *tkcalo_proxy)
                    {
                        // find track with ID matching the pfp index (this convention apparently works only for shower fits...)
                        if (tk->ID() != int(pfp_pxy_v[i_pfp].index()))
                            continue;
                        _shr_tkfit_start_x = tk->Start().X();
                        _shr_tkfit_start_y = tk->Start().Y();
                        _shr_tkfit_start_z = tk->Start().Z();

                        _shr_tkfit_phi = tk->StartDirection().Phi();
                        _shr_tkfit_theta = tk->StartDirection().Theta();

                        auto const tkcalos = tk.get<anab::Calorimetry>();
                        for (const auto &tkcalo : tkcalos)
                        {
                            if (tkcalo->ResidualRange().size() == 0)
                                continue;
                            std::vector<float> dedx4cm;
                            for (size_t ic = 0; ic < tkcalo->ResidualRange().size(); ++ic)
                            {
                                if ((tkcalo->ResidualRange().back() - tkcalo->ResidualRange()[ic]) < 4.)
                                {
                                    dedx4cm.push_back(tkcalo->dEdx()[ic]);
                                }
                            }
                            float dedx4cm_med = -1.;
                            if (dedx4cm.size() > 0)
                            {
                                std::sort(dedx4cm.begin(), dedx4cm.end());
                                if (dedx4cm.size() % 2 == 1)
                                    dedx4cm_med = dedx4cm[dedx4cm.size() / 2];
                                else
                                    dedx4cm_med = 0.5 * (dedx4cm[dedx4cm.size() / 2] + dedx4cm[dedx4cm.size() / 2 - 1]);
                            }
                            if (tkcalo->PlaneID().Plane == 2)
                            {
                                _shr_tkfit_dedx_Y = dedx4cm_med;
                                _shr_tkfit_nhits_Y = dedx4cm.size();
                            }
                            else if (tkcalo->PlaneID().Plane == 1)
                            {
                                _shr_tkfit_dedx_V = dedx4cm_med;
                                _shr_tkfit_nhits_V = dedx4cm.size();
                            }
                            else if (tkcalo->PlaneID().Plane == 0)
                            {
                                _shr_tkfit_dedx_U = dedx4cm_med;
                                _shr_tkfit_nhits_U = dedx4cm.size();
                            }
                        }
                    }
                }
            }
        }

        for (const auto &trk : pfp_pxy.get<recob::Track>())
        {
            if (trkshrscore > fTrkShrscore)
            {
                double trk_start[3] = {trk->Start().X(), trk->Start().Y(), trk->Start().Z()};
                double trk_end[3] = {trk->End().X(), trk->End().Y(), trk->End().Z()};

                auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();

                if (!isFiducial(trk_start) || !isFiducial(trk_end))
                {
                    for (auto ass_clus : clus_pxy_v)
                    {
                        const auto &clus = clus_proxy[ass_clus.key()];
                        auto clus_hit_v = clus.get<recob::Hit>();
                        _hits_outfv += clus_hit_v.size();
                    }
                    continue;
                }
                _n_tracks_contained++;
                unsigned int trk_hits = 0;

                std::vector<art::Ptr<recob::Hit>> hit_v;

                for (auto ass_clus : clus_pxy_v)
                {
                    // get cluster proxy
                    const auto &clus = clus_proxy[ass_clus.key()];
                    auto clus_hit_v = clus.get<recob::Hit>();
                    trk_hits += clus_hit_v.size();
                    for (const auto &hit : clus_hit_v)
                    {
                        hit_v.push_back(hit);
                        if (clus->Plane().Plane == 2)
                            _trk_energy_hits_tot += 1.01 * hit->Integral() * 238 * 23.6e-6 / 0.6 / 1000;
                    }
                    if (clus->Plane().Plane == 0)
                    {
                        _trk_hits_u_tot += clus_hit_v.size();
                    }
                    else if (clus->Plane().Plane == 1)
                    {
                        _trk_hits_v_tot += clus_hit_v.size();
                    }
                    else if (clus->Plane().Plane == 2)
                    {
                        _trk_hits_y_tot += clus_hit_v.size();
                    }
                }

                // Kinetic energy from stopping power of proton in LAr
                float energy_proton = std::sqrt(std::pow(_trkmom.GetTrackMomentum(trk->Length(), proton->PdgCode()), 2) + std::pow(proton->Mass(), 2)) - proton->Mass();
                _trk_energy_tot += energy_proton;
                _trk_hits_tot += trk_hits;

                auto trkpxy2 = pid_proxy[trk.key()];
                auto pidpxy_v = trkpxy2.get<anab::ParticleID>();
                float chipr = std::numeric_limits<float>::lowest();;
                if (pidpxy_v.size() != 0)
                {
                    chipr = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, proton->PdgCode(), 2);
                    if (chipr > 0 && chipr < _trk_pidchipr_best)
                        _trk_pidchipr_best = chipr;
                    if (chipr > 0 && chipr > _trk_pidchipr_worst)
                        _trk_pidchipr_worst = chipr;
                }

                if (trk_hits > _trk_hits_max)
                {
                    if (!fData)
                    {
                        if (clus_pxy_v.size() != 0)
                        {
                            float purity = 0., completeness = 0.;
                            int ibt = searchingfornues::getAssocBtPart(hit_v, assocMCPart, btparts_v, purity, completeness);
                            if (ibt >= 0)
                            {
                                auto &mcp = btparts_v[ibt];
                                auto PDG = mcp.pdg;
                                auto E = mcp.e;
                                _trk_bkt_pdg = PDG;
                                _trk_bkt_purity = purity;
                                _trk_bkt_completeness = completeness;
                                _trk_bkt_E = E;
                                _matched_E += E;
                            }
                        }
                    }

                    trk_vtx.SetXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
                    _trk_distance = (trk_vtx - nu_vtx).Mag();
                    _trk_len = trk->Length();
                    _trk_energy = energy_proton;
                    _trk_pfp_id = i_pfp;
                    _trk_hits_max = trk_hits;
                    _trk_theta = trk->Theta();
                    _trk_phi = trk->Phi();
                    trk_p.SetXYZ(trk->StartDirection().X(), trk->StartDirection().Y(), trk->StartDirection().Z());
                    trk_p.SetMag(sqrt(pow(energy_proton + proton->Mass(), 2) - pow(proton->Mass(), 2)));
                    total_p += trk_p;

                    _trk_score = trkshrscore;

                    _trk_bragg_p = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, proton->PdgCode(), 2),
                                            searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, proton->PdgCode(), 2));
                    _trk_bragg_mu = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, muon->PdgCode(), 2),
                                             searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, muon->PdgCode(), 2));
                    _trk_bragg_mip = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 2);
                    _trk_pidchipr = chipr;
                    _trk_pidchimu = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, muon->PdgCode(), 2);
                    _trk_pida = searchingfornues::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 2);
                }
            }
        }
    }

    _extra_energy_y -= (_trk_energy_hits_tot + _shr_energy_tot);
    _pt = total_p.Perp();
    _p = total_p.Mag();
    _hits_ratio = (float)_shr_hits_tot / (_trk_hits_tot + _shr_hits_tot);
    _tksh_distance = (trk_vtx - shr_vtx).Mag();
    _tksh_angle = trk_p.Dot(shr_p) / (trk_p.Mag() * shr_p.Mag());

    _contained_fraction = ((float)(_trk_hits_tot + _shr_hits_tot)) / (_trk_hits_tot + _shr_hits_tot + _hits_outfv);

    if (_contained_fraction < 0.9)
        return false;

    if (!(_n_tracks_contained > 0 && _n_showers_contained > 0))
        return false;
    return true;
}

void CC0piNpSelection::resetTTree(TTree *_tree)
{

    _shr_pfp_id = 0;
    _trk_pfp_id = 0;

    _trk_hits_max = 0;
    _shr_hits_max = 0;

    _shr_dedx_Y = std::numeric_limits<float>::lowest();
    _shr_dedx_V = std::numeric_limits<float>::lowest();
    _shr_dedx_U = std::numeric_limits<float>::lowest();
    _shr_score = std::numeric_limits<float>::lowest();
    _shr_energy = 0;
    _shr_energy_tot = 0;
    _shr_distance = std::numeric_limits<float>::lowest();
    _tksh_distance = std::numeric_limits<float>::lowest();
    _tksh_angle = std::numeric_limits<float>::lowest();

    _n_showers_contained = 0;
    _n_tracks_contained = 0;

    _trk_distance = std::numeric_limits<float>::lowest();
    _trk_len = std::numeric_limits<float>::lowest();
    _trk_energy = 0;
    _trk_energy_tot = 0;
    _hits_ratio = 0;
    _trk_hits_tot = 0;
    _trk_hits_y_tot = 0;
    _trk_hits_u_tot = 0;
    _trk_hits_v_tot = 0;

    _shr_hits_tot = 0;
    _shr_hits_y_tot = 0;
    _shr_hits_v_tot = 0;
    _shr_hits_u_tot = 0;

    _trk_bragg_p = std::numeric_limits<float>::lowest();
    _trk_bragg_mu = std::numeric_limits<float>::lowest();
    _trk_bragg_mip = std::numeric_limits<float>::lowest();
    _trk_pidchipr = std::numeric_limits<float>::lowest();
    _trk_pidchipr_best = std::numeric_limits<float>::max();
    _trk_pidchipr_worst = std::numeric_limits<float>::lowest();

    _trk_pidchimu = std::numeric_limits<float>::lowest();
    _trk_pida = std::numeric_limits<float>::lowest();
    _trk_score = std::numeric_limits<float>::lowest();
    _trk_bkt_pdg = 0;
    _trk_bkt_purity = std::numeric_limits<float>::lowest();
    _trk_bkt_completeness = std::numeric_limits<float>::lowest();
    _trk_bkt_E = std::numeric_limits<float>::lowest();

    _pt = 0;
    _p = 0;
    _trk_theta = std::numeric_limits<float>::lowest();
    _trk_phi = std::numeric_limits<float>::lowest();
    _shr_theta = std::numeric_limits<float>::lowest();
    _shr_phi = std::numeric_limits<float>::lowest();
    _shr_px = 0;
    _shr_py = 0;
    _shr_pz = 0;
    _shr_openangle = std::numeric_limits<float>::lowest();
    _shr_bkt_pdg = 0;
    _shr_bkt_purity = std::numeric_limits<float>::lowest();
    _shr_bkt_completeness = std::numeric_limits<float>::lowest();
    _shr_bkt_E = std::numeric_limits<float>::lowest();
    _shr_tkfit_start_x = std::numeric_limits<float>::lowest();
    _shr_tkfit_start_y = std::numeric_limits<float>::lowest();
    _shr_tkfit_start_z = std::numeric_limits<float>::lowest();
    _shr_start_x = std::numeric_limits<float>::lowest();
    _shr_start_y = std::numeric_limits<float>::lowest();
    _shr_start_z = std::numeric_limits<float>::lowest();

    _shr_tkfit_phi = std::numeric_limits<float>::lowest();
    _shr_tkfit_theta = std::numeric_limits<float>::lowest();
    _shr_tkfit_dedx_Y = std::numeric_limits<float>::lowest();
    _shr_tkfit_dedx_U = std::numeric_limits<float>::lowest();
    _shr_tkfit_dedx_V = std::numeric_limits<float>::lowest();
    _shr_tkfit_nhits_Y = 0;
    _shr_tkfit_nhits_U = 0;
    _shr_tkfit_nhits_V = 0;
    _dep_E = 0;
    _matched_E = 0;
    _total_hits_y = 0;
    _extra_energy_y = 0;
    _trk_energy_hits_tot = 0;
    _hits_outfv = 0;
    _contained_fraction = 0;
}

void CC0piNpSelection::setBranches(TTree *_tree)
{

    _tree->Branch("trk_id", &_trk_pfp_id, "trk_pfp_id/i");
    _tree->Branch("shr_id", &_shr_pfp_id, "shr_pfp_id/i");

    _tree->Branch("shr_energy_tot", &_shr_energy_tot, "shr_energy_tot/F");
    _tree->Branch("shr_energy", &_shr_energy, "shr_energy/F");
    _tree->Branch("shr_theta", &_shr_theta, "shr_theta/F");
    _tree->Branch("shr_phi", &_shr_phi, "shr_phi/F");
    _tree->Branch("shr_px", &_shr_px, "shr_px/F");
    _tree->Branch("shr_py", &_shr_py, "shr_py/F");
    _tree->Branch("shr_pz", &_shr_pz, "shr_pz/F");
    _tree->Branch("shr_openangle", &_shr_openangle, "shr_openangle/F");
    _tree->Branch("shr_tkfit_start_x", &_shr_tkfit_start_x, "shr_tkfit_start_x/F");
    _tree->Branch("shr_tkfit_start_y", &_shr_tkfit_start_y, "shr_tkfit_start_y/F");
    _tree->Branch("shr_tkfit_start_z", &_shr_tkfit_start_z, "shr_tkfit_start_z/F");
    _tree->Branch("shr_tkfit_theta", &_shr_tkfit_theta, "shr_tkfit_theta/F");
    _tree->Branch("shr_tkfit_phi", &_shr_tkfit_phi, "shr_tkfit_phi/F");
    _tree->Branch("shr_start_x", &_shr_start_x, "shr_start_x/F");
    _tree->Branch("shr_start_y", &_shr_start_y, "shr_start_y/F");
    _tree->Branch("shr_start_z", &_shr_start_z, "shr_start_z/F");
    _tree->Branch("shr_dedx_Y", &_shr_dedx_Y, "shr_dedx_Y/F");
    _tree->Branch("shr_dedx_V", &_shr_dedx_V, "shr_dedx_V/F");
    _tree->Branch("shr_dedx_U", &_shr_dedx_U, "shr_dedx_U/F");
    _tree->Branch("shr_tkfit_dedx_Y", &_shr_tkfit_dedx_Y, "shr_tkfit_dedx_Y/F");
    _tree->Branch("shr_tkfit_dedx_V", &_shr_tkfit_dedx_V, "shr_tkfit_dedx_V/F");
    _tree->Branch("shr_tkfit_dedx_U", &_shr_tkfit_dedx_U, "shr_tkfit_dedx_U/F");
    _tree->Branch("shr_tkfit_nhits_Y", &_shr_tkfit_nhits_Y, "shr_tkfit_nhits_Y/i");
    _tree->Branch("shr_tkfit_nhits_V", &_shr_tkfit_nhits_V, "shr_tkfit_nhits_V/i");
    _tree->Branch("shr_tkfit_nhits_U", &_shr_tkfit_nhits_U, "shr_tkfit_nhits_U/i");
    _tree->Branch("tksh_distance", &_tksh_distance, "tksh_distance/F");
    _tree->Branch("tksh_angle", &_tksh_angle, "tksh_angle/F");
    _tree->Branch("shr_distance", &_shr_distance, "shr_distance/F");
    _tree->Branch("shr_score", &_shr_score, "shr_score/F");
    _tree->Branch("shr_bkt_pdg", &_shr_bkt_pdg, "shr_bkt_pdg/I");
    _tree->Branch("shr_bkt_purity", &_shr_bkt_purity, "shr_bkt_purity/F");
    _tree->Branch("shr_bkt_completeness", &_shr_bkt_completeness, "shr_bkt_completeness/F");
    _tree->Branch("shr_bkt_E", &_shr_bkt_E, "shr_bkt_E/F");

    _tree->Branch("trk_len", &_trk_len, "trk_len/F");
    _tree->Branch("trk_theta", &_trk_theta, "trk_theta/F");
    _tree->Branch("trk_phi", &_trk_phi, "trk_phi/F");
    _tree->Branch("trk_energy", &_trk_energy, "trk_energy/F");
    _tree->Branch("trk_energy_tot", &_trk_energy_tot, "trk_energy_tot/F");
    _tree->Branch("trk_distance", &_trk_distance, "trk_distance/F");
    _tree->Branch("trk_score", &_trk_score, "trk_score/F");
    _tree->Branch("trk_bkt_pdg", &_trk_bkt_pdg, "trk_bkt_pdg/I");
    _tree->Branch("trk_bkt_purity", &_trk_bkt_purity, "trk_bkt_purity/F");
    _tree->Branch("trk_bkt_completeness", &_trk_bkt_completeness, "trk_bkt_completeness/F");
    _tree->Branch("trk_bkt_E", &_trk_bkt_E, "trk_bkt_E/F");
    _tree->Branch("trk_chipr_best", &_trk_pidchipr_best, "trk_chipr_best/F");
    _tree->Branch("trk_chipr_worst", &_trk_pidchipr_worst, "trk_chipr_worst/F");

    _tree->Branch("trk_chipr", &_trk_pidchipr, "trk_chipr/F");
    _tree->Branch("trk_chimu", &_trk_pidchimu, "trk_chimu/F");
    _tree->Branch("trk_pida", &_trk_pida, "trk_pida/F");
    _tree->Branch("trk_bragg_p", &_trk_bragg_p, "trk_bragg_p/F");
    _tree->Branch("trk_bragg_mu", &_trk_bragg_mu, "trk_bragg_mu/F");
    _tree->Branch("trk_bragg_mip", &_trk_bragg_mip, "trk_bragg_mip/F");

    _tree->Branch("trk_hits_max", &_trk_hits_max, "trk_hits_max/i");
    _tree->Branch("shr_hits_max", &_shr_hits_max, "shr_hits_max/i");

    _tree->Branch("total_hits_y", &_total_hits_y, "total_hits_y/i");
    _tree->Branch("extra_energy_y", &_extra_energy_y, "extra_energy_y/F");
    _tree->Branch("trk_energy_hits_tot", &_trk_energy_hits_tot, "trk_energy_hits_tot/F");

    _tree->Branch("shr_hits_tot", &_shr_hits_tot, "shr_hits_tot/i");
    _tree->Branch("shr_hits_y_tot", &_shr_hits_y_tot, "shr_hits_y_tot/i");
    _tree->Branch("shr_hits_u_tot", &_shr_hits_u_tot, "shr_hits_u_tot/i");
    _tree->Branch("shr_hits_v_tot", &_shr_hits_v_tot, "shr_hits_v_tot/i");

    _tree->Branch("trk_hits_tot", &_trk_hits_tot, "trk_hits_tot/i");
    _tree->Branch("trk_hits_y_tot", &_trk_hits_y_tot, "trk_hits_y_tot/i");
    _tree->Branch("trk_hits_u_tot", &_trk_hits_u_tot, "trk_hits_u_tot/i");
    _tree->Branch("trk_hits_v_tot", &_trk_hits_v_tot, "trk_hits_v_tot/i");

    _tree->Branch("n_tracks_contained", &_n_tracks_contained, "n_tracks_contained/i");
    _tree->Branch("n_showers_contained", &_n_showers_contained, "n_showers_contained/i");

    _tree->Branch("dep_E", &_dep_E, "dep_E/F");
    _tree->Branch("matched_E", &_matched_E, "matched_E/F");

    _tree->Branch("hits_ratio", &_hits_ratio, "hits_ratio/F");
    _tree->Branch("contained_fraction", &_contained_fraction, "contained_fraction/F");
    _tree->Branch("pt", &_pt, "pt/F");
    _tree->Branch("p", &_p, "p/F");
}

DEFINE_ART_CLASS_TOOL(CC0piNpSelection)
} // namespace selection

#endif
