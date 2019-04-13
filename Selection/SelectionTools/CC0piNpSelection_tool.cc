#ifndef SELECTION_CC0PINPSELECTION_CXX
#define SELECTION_CC0PINPSELECTION_CXX

#include <iostream>
#include "SelectionToolBase.h"
#include "nusimdata/SimulationBase/MCTruth.h"

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

private:
    /**
     * @brief read truth info associated to pi0
     */
    void readTruth(art::Event const &e);

    std::vector<double> _elec_vtx;
    std::vector<double> _p_e;
    std::vector<std::vector<double>> _p_vtx;

    unsigned int _n_showers;
    unsigned int _n_tracks;

    std::vector<std::vector<double>> _shr_energy;
    std::vector<std::vector<double>> _shr_dedx;
    std::vector<size_t> _shr_pfp_id;

    std::vector<double> _trk_len;
    std::vector<size_t> _trk_pfp_id;
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

    if (!fData)
        readTruth(e);

    double nu_vtx[3] = {};

    for (size_t i_pfp = 0; i_pfp < pfp_pxy_v.size(); i_pfp++)
    {
        auto PDG = fabs(pfp_pxy_v[i_pfp]->PdgCode());


        if (PDG == 12 || PDG == 14) {

            auto vtx = pfp_pxy_v[i_pfp].get<recob::Vertex>();
            if (vtx.size() != 1)
            {
                std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
            }
            else
            {
                vtx.at(0)->XYZ(nu_vtx);
            }

            break;
        }
    }

    for (size_t i_pfp = 0; i_pfp < pfp_pxy_v.size(); i_pfp++) {
        auto PDG = fabs(pfp_pxy_v[i_pfp]->PdgCode());

        if (PDG == 12 || PDG == 14)
            continue;

        for (const auto &shr : pfp_pxy_v[i_pfp].get<recob::Shower>())
        {
            _n_showers++;
            _shr_dedx.push_back(shr->dEdx());
            _shr_energy.push_back(shr->Energy());
            _shr_pfp_id.push_back(i_pfp);
        }

        for (const auto &trk : pfp_pxy_v[i_pfp].get<recob::Track>())
        {
            _n_tracks++;
            _trk_len.push_back(trk->Length());
            _trk_pfp_id.push_back(i_pfp);
        }
    }

    if (_n_tracks + _n_showers < 2)
        return false;

    return true;
}

void CC0piNpSelection::readTruth(art::Event const &e)
{

    std::cout << "CC0piNpSelection reading MCTruth" << std::endl;

    auto const &mcparticles_handle = e.getValidHandle<std::vector<simb::MCParticle>>("largeant");
    auto const &mcparticles(*mcparticles_handle);

    for (auto &mcparticle : mcparticles)
    {
        if (mcparticle.StatusCode() != 1)
        {
            continue;
        }
        if (mcparticle.PdgCode() == 13)
        {
            _elec_vtx.push_back(mcparticle.Trajectory().X(0));
            _elec_vtx.push_back(mcparticle.Trajectory().Y(0));
            _elec_vtx.push_back(mcparticle.Trajectory().Z(0));
        }

        if (mcparticle.PdgCode() == 2212)
        {
            std::vector<double> p_vertex;

            p_vertex.push_back(mcparticle.Trajectory().X(0));
            p_vertex.push_back(mcparticle.Trajectory().Y(0));
            p_vertex.push_back(mcparticle.Trajectory().Z(0));

            _p_vtx.push_back(p_vertex);
            _p_e.push_back(mcparticle.E());
        }
    }
}

void CC0piNpSelection::resetTTree(TTree *_tree)
{
    std::cout << "CC0piNpSelection reset TTree" << std::endl;
    _elec_vtx.clear();
    _p_vtx.clear();
    _p_e.clear();

    _shr_energy.clear();
    _shr_dedx.clear();
    _shr_pfp_id.clear();

    _trk_len.clear();
    _trk_pfp_id.clear();

    _n_showers = 0;
    _n_tracks = 0;
}

void CC0piNpSelection::setBranches(TTree *_tree)
{
    _tree->Branch("elec_vtx", "std::vector< double >", &_elec_vtx);
    _tree->Branch("p_vtx", "std::vector< std::vector< double > >", &_p_vtx);
    _tree->Branch("p_e", "std::vector< double >", &_p_e);

    _tree->Branch("trk_len", "std::vector< double >", &_trk_len);
    _tree->Branch("trk_pfp_id", "std::vector< size_t >", &_trk_pfp_id);

    _tree->Branch("n_tracks", &_n_tracks, "n_tracks/i");
}

DEFINE_ART_CLASS_TOOL(CC0piNpSelection)
} // namespace selection

#endif