#ifndef ANALYSIS_COSMICIP_CXX
#define ANALYSIS_COSMICIP_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"

#include "lardataobj/Simulation/SimPhotons.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "larcore/Geometry/Geometry.h" 
#include "larcorealg/Geometry/GeometryCore.h" 
#include "lardata/Utilities/GeometryUtilities.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"

#include <string>

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       SimPhoton
    // File:        SimPhoton.cc
    //
    //              A basic analysis example
    //
    // Configuration parameters:
    //
    // TBD
    //
    // Created by David Caratelli (davidc@fnal.gov) on 08/09/2021
    //
    ////////////////////////////////////////////////////////////////////////

  class SimPhoton : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    SimPhoton(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~SimPhoton(){ };
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset);
    
    /**
     * @brief Analysis function
     */
    void analyzeEvent(art::Event const& e, bool fData) override;

    /**
     * @brief Analyze slice
     */
    void analyzeSlice(art::Event const& e, std::vector<ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) override;

    /**
     * @brief set branches for TTree
     */
    void setBranches(TTree* _tree) override;

    /**
     * @brief reset ttree branches
     */
    void resetTTree(TTree* _tree) override;
    
  private:

    art::InputTag fSimPhotonProducer;

    std::vector<float> _simphoton_number_v;
    std::vector<float> _simphoton_tmin_v;
    std::vector<float> _simphoton_tmax_v;

  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  SimPhoton::SimPhoton(const fhicl::ParameterSet& p)
  {

    fSimPhotonProducer  = p.get< art::InputTag >("SimPhotonProducer");

    // load PMT coordinates
    art::ServiceHandle<geo::Geometry> geom;
    double xyz[3];
    for (size_t pmt=0; pmt < 32; pmt++) {
      geom->OpDetGeoFromOpDet(pmt).GetCenter(xyz);
      std::cout << "PMT OpDet " << pmt << " has coordinates [" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "]" << std::endl;
    }

  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void SimPhoton::configure(fhicl::ParameterSet const & p)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void SimPhoton::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
  {


    return;
  }
  
  void SimPhoton::analyzeEvent(art::Event const &e, bool fData)
  {


    art::Handle<std::vector<sim::SimPhotons> > simphotons_h;
    e.getByLabel( fSimPhotonProducer , simphotons_h );   
    
    for (size_t p=0; p < simphotons_h->size(); p++) {

      auto simphoton = simphotons_h->at(p);

      auto ch = simphoton.OpChannel();

      std::cout << "Channel " << ch << " has " << simphoton.size() << " SimPhotons" << std::endl;

      if (ch >= 32) continue;

      _simphoton_number_v.at(ch) = simphoton.size();
      // calculate min and max time for photon arrival
      float tmin = std::numeric_limits<float>::max();
      float tmax = std::numeric_limits<float>::lowest();
      for (size_t i=0; i < simphoton.size(); i++) {
	auto thistime = simphoton[i].Time;
	std::cout << "\t Photon with time " << thistime << " has MotherTrackID " << simphoton[i].MotherTrackID << std::endl;
	if (thistime > tmax) {tmax = thistime; }
	if (thistime < tmin) {tmin = thistime; }
      }// for all sim-photons for this channel
      _simphoton_tmin_v.at(ch) = tmin;
      _simphoton_tmax_v.at(ch) = tmax;
      
    }// for all SimPhotons

    return;
  }

  void SimPhoton::setBranches(TTree* _tree)
  {

    std::cout << "[SimPhotonTool] setBranches begin" << std::endl;

    _simphoton_number_v = std::vector<float>(32,0);
    _simphoton_tmin_v   = std::vector<float>(32,0);
    _simphoton_tmax_v   = std::vector<float>(32,0);

    _tree->Branch("simphoton_number_v","std::vector<float>",&_simphoton_number_v);
    _tree->Branch("simphoton_tmin_v"  ,"std::vector<float>",&_simphoton_tmin_v  );
    _tree->Branch("simphoton_tmax_v"  ,"std::vector<float>",&_simphoton_tmax_v  );

  }
  
  void SimPhoton::resetTTree(TTree* _tree)
  {

    std::cout << "[SimPhotonTool] resetTTree begin" << std::endl;

    _simphoton_number_v = std::vector<float>(32,0);
    _simphoton_tmin_v   = std::vector<float>(32,0);
    _simphoton_tmax_v   = std::vector<float>(32,0);

  }

  
  DEFINE_ART_CLASS_TOOL(SimPhoton)
} // namespace analysis

#endif
