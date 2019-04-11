#ifndef SELECTION_CC0PINPSELECTION_CXX
#define SELECTION_CC0PINPSELECTION_CXX

#include <iostream>
#include "SelectionToolBase.h"

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
    void setBranches(TTree *_tree){};

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree *_tree){};

private:
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

    return true;
}

DEFINE_ART_CLASS_TOOL(CC0piNpSelection)
} // namespace selection

#endif