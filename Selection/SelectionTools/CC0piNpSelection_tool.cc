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

namespace selection
{
////////////////////////////////////////////////////////////////////////
//
// Class:       CC0piNpSelection
// File:        CC0piNpSelection_tool.cc
//
//              nu_e CC0pi-Np selection tool
//
// Configuration parameters:
//
// TBD
//
// Created by Stefano Roberto Soleti (srsoleti@fnal.gov) on 04/11/2019
//
////////////////////////////////////////////////////////////////////////

class CC0piNpSelection : public SelectionToolBase
{

public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
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
     * @brief set branches for TTree
     */
    void setBranches(TTree *_tree);

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree *_tree);

    /**
     * @brief check if inside fiducial volume
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

    unsigned int _n_showers;
    unsigned int _n_tracks;
    unsigned int _max_hits_shower;
    unsigned int _max_hits_track;
    unsigned int _shr_hits_tot;
    unsigned int _trk_hits_tot;

    float _shr_energy;
    float _shr_energy_tot;
    float _shr_dedx_Y;
    float _shr_dedx_V;
    float _shr_dedx_U;
    float _shr_distance;
    float _shr_score;
    float _shr_theta;
    float _shr_phi;
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
    float _trk_pidchimu;
    float _trk_pida;
    float _trk_score;

    float _pt;
    float _p;

    float _shr_bkt_purity;
    float _shr_bkt_completeness;
    int _shr_bkt_pdg;
    float _trk_bkt_purity;
    float _trk_bkt_completeness;
    int _trk_bkt_pdg;

    float _dep_E;
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
    fTRKproducer = pset.get<art::InputTag>("TRKproducer", "pandora");
    fPIDproducer = pset.get<art::InputTag>("PIDproducer", "pandoracalipidSCE");
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
            if (fabs(PDG) == proton->PdgCode()) {
                double ke = mcp.E() - proton->Mass();
                if (ke > 0.040)
                    _dep_E += ke;
            }
            if (fabs(PDG) == electron->PdgCode()) {
                double ke = mcp.E() - electron->Mass();
                if (ke > 0.030)
                    _dep_E += ke;
            }
        }
    }

    double nu_vtx[3] = {};
    TVector3 nuvtx;

    _max_hits_track = 0;
    _max_hits_shower = 0;

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
                if (!isFiducial(nu_vtx)) {
                    return false;
                }
                nuvtx.SetXYZ(nu_vtx[0], nu_vtx[1], nu_vtx[2]);
            }

            break;
        }
    }

    TVector3 total_p;

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
                if (!isFiducial(shr_vertex)) {
                    return false;
                }
                _n_showers++;

                unsigned int shr_hits = 0;

                auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();
                std::vector<art::Ptr<recob::Hit>> hit_v;

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
                }

                _shr_energy_tot += shr->Energy()[2] / 1000;
                _shr_hits_tot += shr_hits;
                if (shr_hits > _max_hits_shower)
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
                                _shr_bkt_pdg = PDG;
                                _shr_bkt_purity = purity;
                                _shr_bkt_completeness = completeness;
                            }
                        }
                    }
                    _shr_distance = (shr->ShowerStart() - nuvtx).Mag();
                    TVector3 shr_p = shr->ShowerStart();
                    shr_p.SetMag(sqrt(pow(shr->Energy()[2] / 1000. + electron->Mass(), 2) - pow(electron->Mass(), 2)));
                    total_p += shr_p;
                    _shr_dedx_Y = shr->dEdx()[2];
                    _shr_dedx_V = shr->dEdx()[1];
                    _shr_dedx_U = shr->dEdx()[0];
                    _shr_energy = shr->Energy()[2] / 1000; // GeV
                    _shr_pfp_id = i_pfp;
                    _max_hits_shower = shr_hits;
                    _shr_score = trkshrscore;
                    _shr_theta = shr->Direction().Theta();
                    _shr_phi = shr->Direction().Phi();
                }
            }
        }

        for (const auto &trk : pfp_pxy.get<recob::Track>())
        {
            if (trkshrscore > fTrkShrscore)
            {
                double trk_start[3] = {trk->Start().X(), trk->Start().Y(), trk->Start().Z()};
                double trk_end[3] = {trk->End().X(), trk->End().Y(), trk->End().Z()};
                if (!isFiducial(trk_start) || !isFiducial(trk_end))
                {
                    return false;
                }
                _n_tracks++;
                unsigned int trk_hits = 0;

                auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();
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
                    }
                }

                float energy_proton = std::sqrt(std::pow(_trkmom.GetTrackMomentum(trk->Length(), proton->PdgCode()), 2) + std::pow(proton->Mass(), 2)) - proton->Mass();
                _trk_energy_tot += energy_proton;
                _trk_hits_tot += trk_hits;

                if (trk_hits > _max_hits_track)
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
                                _trk_bkt_pdg = PDG;
                                _trk_bkt_purity = purity;
                                _trk_bkt_completeness = completeness;
                            }
                        }
                    }
                    TVector3 trk_vtx;
                    trk_vtx.SetXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
                    trk_vtx -= nuvtx;
                    _trk_distance = trk_vtx.Mag();
                    _trk_len = trk->Length();
                    _trk_energy = energy_proton;
                    _trk_pfp_id = i_pfp;
                    _max_hits_track = trk_hits;
                    _trk_theta = trk->Theta();
                    _trk_phi = trk->Phi();
                    TVector3 trk_p = trk_vtx;
                    trk_p.SetMag(sqrt(pow(energy_proton + proton->Mass(), 2) - pow(proton->Mass(), 2)));
                    total_p += trk_p;

                    _trk_score = trkshrscore;

                    auto trkpxy2 = pid_proxy[trk.key()];
                    auto pidpxy_v = trkpxy2.get<anab::ParticleID>();
                    if (pidpxy_v.size() == 0)
                    {
                        continue;
                    }

                    _trk_bragg_p = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, proton->PdgCode(), 2),
                                            searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, proton->PdgCode(), 2));
                    _trk_bragg_mu = std::max(searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, muon->PdgCode(), 2),
                                             searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, muon->PdgCode(), 2));
                    _trk_bragg_mip = searchingfornues::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 2);
                    _trk_pidchipr = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, proton->PdgCode(), 2);
                    _trk_pidchimu = searchingfornues::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, muon->PdgCode(), 2);
                    _trk_pida = searchingfornues::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 2);
                }
            }
        }
    }

    if (!(_n_tracks > 0 && _n_showers > 0))
        return false;

    _pt = total_p.Perp();
    _p = total_p.Mag();
    _hits_ratio = (float)_shr_hits_tot/(_trk_hits_tot+_shr_hits_tot);

    return true;
}

void CC0piNpSelection::resetTTree(TTree *_tree)
{

    _shr_pfp_id = 0;
    _trk_pfp_id = 0;

    _max_hits_track = 0;
    _max_hits_shower = 0;

    _shr_dedx_Y = std::numeric_limits<float>::min();
    _shr_dedx_V = std::numeric_limits<float>::min();
    _shr_dedx_U = std::numeric_limits<float>::min();
    _shr_score = std::numeric_limits<float>::min();
    _shr_energy = 0;
    _shr_energy_tot = 0;
    _shr_distance = 0;
    _n_showers = 0;
    _n_tracks = 0;

    _trk_distance = 0;
    _trk_len = 0;
    _trk_energy = 0;
    _trk_energy_tot = 0;
    _hits_ratio = 0;
    _trk_hits_tot = 0;
    _shr_hits_tot = 0;

    _trk_bragg_p = std::numeric_limits<float>::min();
    _trk_bragg_mu = std::numeric_limits<float>::min();
    _trk_bragg_mip = std::numeric_limits<float>::min();
    _trk_pidchipr = std::numeric_limits<float>::min();
    _trk_pidchimu = std::numeric_limits<float>::min();
    _trk_pida = std::numeric_limits<float>::min();
    _trk_score = std::numeric_limits<float>::min();
    _trk_bkt_pdg = 0;
    _trk_bkt_purity = std::numeric_limits<float>::min();
    _trk_bkt_completeness = std::numeric_limits<float>::min();

    _pt = 0;
    _p = 0;
    _trk_theta = std::numeric_limits<float>::min();
    _trk_phi = std::numeric_limits<float>::min();
    _shr_theta = std::numeric_limits<float>::min();
    _shr_phi = std::numeric_limits<float>::min();
    _shr_bkt_pdg = 0;
    _shr_bkt_purity = std::numeric_limits<float>::min();
    _shr_bkt_completeness = std::numeric_limits<float>::min();

    _dep_E = 0;
}

void CC0piNpSelection::setBranches(TTree *_tree)
{

    _tree->Branch("trk_id", &_trk_pfp_id, "trk_pfp_id/i");
    _tree->Branch("shr_id", &_shr_pfp_id, "shr_pfp_id/i");

    _tree->Branch("shr_energy_tot", &_shr_energy_tot, "shr_energy_tot/F");
    _tree->Branch("shr_energy", &_shr_energy, "shr_energy/F");
    _tree->Branch("shr_theta", &_shr_theta, "shr_theta/F");
    _tree->Branch("shr_phi", &_shr_phi, "shr_phi/F");
    _tree->Branch("shr_dedx_Y", &_shr_dedx_Y, "shr_dedx_Y/F");
    _tree->Branch("shr_dedx_V", &_shr_dedx_V, "shr_dedx_V/F");
    _tree->Branch("shr_dedx_U", &_shr_dedx_U, "shr_dedx_U/F");
    _tree->Branch("shr_distance", &_shr_distance, "shr_distance/F");
    _tree->Branch("shr_score", &_shr_score, "shr_score/F");
    _tree->Branch("shr_bkt_pdg", &_shr_bkt_pdg, "shr_bkt_pdg/I");
    _tree->Branch("shr_bkt_purity", &_shr_bkt_purity, "shr_bkt_purity/F");
    _tree->Branch("shr_bkt_completeness", &_shr_bkt_completeness, "shr_bkt_completeness/F");

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
    _tree->Branch("trk_chipr", &_trk_pidchipr, "trk_chipr/F");
    _tree->Branch("trk_chimu", &_trk_pidchimu, "trk_chimu/F");
    _tree->Branch("trk_pida", &_trk_pida, "trk_pida/F");
    _tree->Branch("trk_bragg_p", &_trk_bragg_p, "trk_bragg_p/F");
    _tree->Branch("trk_bragg_mu", &_trk_bragg_mu, "trk_bragg_mu/F");
    _tree->Branch("trk_bragg_mip", &_trk_bragg_mip, "trk_bragg_mip/F");

    _tree->Branch("trk_hits_max", &_max_hits_track, "trk_hits_max/i");
    _tree->Branch("shr_hits_max", &_max_hits_shower, "shr_hits_max/i");

    _tree->Branch("shr_hits_tot", &_shr_hits_tot, "shr_hits_tot/i");
    _tree->Branch("trk_hits_tot", &_trk_hits_tot, "trk_hits_tot/i");

    _tree->Branch("n_tracks", &_n_tracks, "n_tracks/i");
    _tree->Branch("n_showers", &_n_showers, "n_showers/i");

    _tree->Branch("dep_E", &_dep_E, "dep_E/F");

    _tree->Branch("hits_ratio", &_hits_ratio, "hits_ratio/F");
    _tree->Branch("pt", &_pt, "pt/F");
    _tree->Branch("p", &_p, "p/F");

}

DEFINE_ART_CLASS_TOOL(CC0piNpSelection)
} // namespace selection

#endif