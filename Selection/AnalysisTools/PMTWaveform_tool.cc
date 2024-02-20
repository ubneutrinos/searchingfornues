#ifndef ANALYSIS_COSMICIP_CXX
#define ANALYSIS_COSMICIP_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "larcore/Geometry/Geometry.h" 
#include "larcorealg/Geometry/GeometryCore.h" 
#include "lardata/Utilities/GeometryUtilities.h"

// pmt gains & LY services
#include "larevt/CalibrationDBI/Interface/PmtGainService.h" 
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"
#include "ubevt/Utilities/PMTRemapProvider.h"
#include "ubevt/Utilities/PMTRemapService.h"
#include "ubevt/Database/LightYieldService.h"
#include "ubevt/Database/LightYieldProvider.h"
#include "ubevt/Database/UbooneLightYieldProvider.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"

#include <string>

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       PMTWaveform
    // File:        PMTWaveform.cc
    //
    //              A basic analysis example
    //
    // Configuration parameters:
    //
    // TBD
    //
    // Created by Giuseppe Cerati (cerati@fnal.gov) on 03/15/2019
    //
    ////////////////////////////////////////////////////////////////////////

  class PMTWaveform : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    PMTWaveform(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~PMTWaveform(){ };
    
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

    art::InputTag fFLASHproducer;
    art::InputTag fFLASHCALIBproducer;
    art::InputTag fPMTWFproducer;

    float _flash_pe_uncalib, _flash_pe_uncalib_calib;
    float _flash_zcenter, _flash_ycenter, _flash_zwidth, _flash_ywidth;
    std::vector<float> _flash_pe_uncalib_v;
    std::vector<float> _flash_pe_uncalib_calib_v;
    std::vector<float> _gain_ampl_v, _gain_area_v;
    float _flash_time;
    std::vector< std::vector<short> > _opwf_v;
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  PMTWaveform::PMTWaveform(const fhicl::ParameterSet& p)
  {

    fFLASHproducer  = p.get< art::InputTag >("FLASHproducer");
    fFLASHCALIBproducer  = p.get< art::InputTag >("FLASHCALIBproducer");
    fPMTWFproducer  = p.get< art::InputTag >("PMTWFproducer");

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
  void PMTWaveform::configure(fhicl::ParameterSet const & p)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void PMTWaveform::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
  {


    return;
  }
  
  void PMTWaveform::analyzeEvent(art::Event const &e, bool fData)
  {

    // PMT gains
    const lariov::PmtGainProvider& gain_provider = art::ServiceHandle<lariov::PmtGainService>()->GetProvider();
    const ::util::PMTRemapProvider &pmtremap_provider = art::ServiceHandle<util::PMTRemapService>()->GetProvider();

    std::cout << "[PMTWaveformTool] get gains" << std::endl;
    
    _gain_ampl_v = std::vector<float>(32,0.);
    _gain_area_v = std::vector<float>(32,0.);
    for (size_t i=0; i < 32; i++) {
      auto oldch = pmtremap_provider.OriginalOpChannel(i);
      float amplgain = gain_provider.ExtraInfo(oldch%100).GetFloatData("amplitude_gain");
      float areagain = gain_provider.Gain(oldch%100);
      _gain_area_v[i] = areagain;
      _gain_ampl_v[i] = amplgain;
    }

    std::cout << "[PMTWaveformTool] get flash" << std::endl;

    // FLASH
    art::Handle<std::vector<recob::OpFlash> > flash_h;
    e.getByLabel( fFLASHproducer , flash_h );   
    for (size_t f=0; f < flash_h->size(); f++) {

      auto flash = flash_h->at(f);
      
      if (flash.TotalPE() > _flash_pe_uncalib) {
        _flash_pe_uncalib    = flash.TotalPE();
        _flash_time  = flash.Time();
        _flash_zcenter = flash.ZCenter();
        _flash_zwidth  = flash.ZWidth();
        _flash_ycenter = flash.YCenter();
        _flash_ywidth  = flash.YWidth();

        for (int pmt=0; pmt < 32; pmt++){
          _flash_pe_uncalib_v[pmt] = flash.PE(pmt);
        }
      }// if larger then other flashes

    }// for all flashes

    // FLASH (CALIBRATED)
    art::Handle<std::vector<recob::OpFlash> > flash_calib_h;
    e.getByLabel( fFLASHCALIBproducer , flash_calib_h );
  
    for (size_t f=0; f < flash_calib_h->size(); f++) {

      auto flash = flash_calib_h->at(f);

      if (flash.TotalPE() > _flash_pe_uncalib_calib) {
	      _flash_pe_uncalib_calib = flash.TotalPE();
        for (int pmt=0; pmt < 32; pmt++){
          _flash_pe_uncalib_calib_v[pmt] = flash.PE(pmt);
        }
      }// if larger then other flashes                                                                                                                                       
    }// for all calibrated flashes
    
    std::cout << "[PMTWaveformTool] get waveforms" << std::endl;
    
    // WAVEFORM
    //_opwf_v = std::vector< std::vector<short> >(32,std::vector<short>(700,0));
    
    art::Handle< std::vector< raw::OpDetWaveform > > wf_handle;
    e.getByLabel( fPMTWFproducer, wf_handle );
    
    for (size_t i=0; i < wf_handle->size(); i++) {

      auto const wf = wf_handle->at(i);

      //auto ch = wf.ChannelNumber();

      if (i >= _opwf_v.size()) continue;

      
      for (size_t tick=0; tick < 700; tick++){
	  _opwf_v[i][tick] = wf[tick];
      }// for all channels
    }// for all waveforms      

    std::cout << "[PMTWaveformTool] end" << std::endl;

    return;
  }

  void PMTWaveform::setBranches(TTree* _tree)
  {

    std::cout << "[PMTWaveformTool] setBranches begin" << std::endl;

    _opwf_v = std::vector< std::vector<short> >(32,std::vector<short>(700,0));
    _flash_pe_uncalib_v = std::vector<float>(32,0);
    _flash_pe_uncalib_calib_v = std::vector<float>(32,0);

    _tree->Branch("flash_pe_uncalib",&_flash_pe_uncalib,"flash_pe_uncalib/F");
    _tree->Branch("flash_pe_uncalib_calib",&_flash_pe_uncalib_calib,"flash_pe_uncalib_calib/F");
    _tree->Branch("flash_time",&_flash_time,"flash_time/F");
    _tree->Branch("flash_zcenter",&_flash_zcenter,"flash_zcenter/F"); 
    _tree->Branch("flash_ycenter",&_flash_ycenter,"flash_ycenter/F");
    _tree->Branch("flash_zwidth",&_flash_zwidth,"flash_zwidth/F");
    _tree->Branch("flash_ywidth",&_flash_ywidth,"flash_ywidth/F");
    _tree->Branch("flash_pe_uncalib_v","std::vector<float>",&_flash_pe_uncalib_v);
    _tree->Branch("flash_pe_uncalib_calib_v","std::vector<float>",&_flash_pe_uncalib_calib_v);
    _tree->Branch("gain_area_v","std::vector<float>",&_gain_area_v);
    _tree->Branch("gain_ampl_v","std::vector<float>",&_gain_ampl_v);

    for (unsigned int i=0; i < 32; i++) 
      _tree->Branch( TString::Format("waveform_%02u",i), "std::vector<short>", &(_opwf_v[i]) );
    
    std::cout << "[PMTWaveformTool] setBranches end" << std::endl;

  }
  
  void PMTWaveform::resetTTree(TTree* _tree)
  {

    std::cout << "[PMTWaveformTool] resetTTree begin" << std::endl;

    _flash_pe_uncalib   = std::numeric_limits<float>::lowest();
    _flash_pe_uncalib_calib   = std::numeric_limits<float>::lowest();
    _flash_time = std::numeric_limits<float>::lowest();
    _gain_area_v = std::vector<float>(32,0);
    _gain_ampl_v = std::vector<float>(32,0);
    _flash_pe_uncalib_v = std::vector<float>(32,0);
    _flash_pe_uncalib_calib_v = std::vector<float>(32,0);

    std::cout << "[PMTWaveformTool] resetTTree end" << std::endl;

  }

  
  DEFINE_ART_CLASS_TOOL(PMTWaveform)
} // namespace analysis

#endif
