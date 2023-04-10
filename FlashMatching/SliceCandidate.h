#ifndef _SLICECANDIDATE_
#define _SLICECANDIDATE_

#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larcore/Geometry/Geometry.h"

#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"

#include "ubevt/Database/UbooneElectronLifetimeProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeService.h"

#include "ubreco/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "ubreco/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"
#include "Objects/CartesianVector.h"

#include "ubevt/Utilities/PMTRemapService.h"
#include "ubevt/Utilities/PMTRemapProvider.h"

// -------------------------------------------------------------------------------------------------------------------------------------

namespace flashmatch {
    
  /**
   *  @brief  A candidate for the beam flash
   */
  class FlashCandidate
  {
  public:
    /**
     *  @brief  Default constructor
     */
  FlashCandidate()
    : m_run(-std::numeric_limits<int>::max()),
      m_subRun(-std::numeric_limits<int>::max()),
      m_event(-std::numeric_limits<int>::max()),
      m_time(-std::numeric_limits<float>::max()),
      m_totalPE(-std::numeric_limits<float>::max()),
      m_centerY(-std::numeric_limits<float>::max()),
      m_centerZ(-std::numeric_limits<float>::max()),
      m_widthY(-std::numeric_limits<float>::max()),
      m_widthZ(-std::numeric_limits<float>::max()),
      m_inBeamWindow(false),
      m_isBrightestInWindow(false),
      m_isBeamFlash(false)
	{
	}
      
    /**
     *  @brief  Parametrized constructor
     *
     *  @param  event the art event
     *  @param  flash the flash
     */
  FlashCandidate(const art::Event &event, const recob::OpFlash &flash)
    : m_run(event.run()),
      m_subRun(event.subRun()),
      m_event(event.event()),
      m_time(flash.Time()),
      m_totalPE(flash.TotalPE()),
      m_centerY(flash.YCenter()),
      m_centerZ(flash.ZCenter()),
      m_widthY(flash.YWidth()),
      m_widthZ(flash.ZWidth()),
      m_inBeamWindow(false),
      m_isBrightestInWindow(false),
      m_isBeamFlash(false)
	{
	
	  art::ServiceHandle<geo::Geometry> geo;
	  // gain service
	  /* const ::lariov::PmtGainProvider& gain_provider = art::ServiceHandle<lariov::PmtGainService>()->GetProvider(); */
	  /* // pmt remapping service */
	  /* const ::util::PMTRemapProvider& pmtremap_provider = art::ServiceHandle<util::PMTRemapService>()->GetProvider(); */
	
	  const auto nOpDets(geo->NOpDets());
	
	  /* // which FEM to get? */
	  /* int FEM = (flash.PEs().size() / 100) * 100; */
	  /* std::cout << "FEM=" << FEM << std::endl; */
	
	  m_peSpectrum.resize(nOpDets);
	
	  for (size_t OpChannel = 0; OpChannel < nOpDets; ++OpChannel) {
	    if (!geo->IsValidOpChannel(OpChannel)) continue;
	    size_t OpDet     = geo->OpDetFromOpChannel(OpChannel);
	    m_peSpectrum.at(OpDet)  = flash.PEs()[OpChannel];
	  }
	  /* for (size_t pmt=FEM; pmt < FEM+nOpDets; pmt++) { */
	  /*   // pmt%100 here represents the channel number */
	  /*   size_t OpChannel = pmt%100; */
	    /* if (!geo->IsValidOpChannel(OpChannel)) continue; */
	    /* size_t OpDet     = geo->OpDetFromOpChannel(OpChannel); */
	    /* auto PE   = flash.PEs()[pmt]; */

	    /* // SPE calibration performed on incorrect OpDetWaveform */
	    /* // channel nomenclature */
	    /* // OpFlash channel number is after it was remapped */
	    /* // to apply the gain for the appropriate channel */
	    /* // one needs to find the OpChannel number */
	    /* // associated to the original waveform  */
	    /* // hence apply the inverse mapping */
	    /* auto oldch = pmtremap_provider.OriginalOpChannel( OpChannel ); */
	    /* std::cout << "Remapping channel " << OpChannel << " to channel " << oldch << std::endl; */
	    /* auto gain = gain_provider.Gain( oldch ); */
	    /* //auto gain = gain_provider.Gain(OpChannel); */

	    /* std::cout << "Applying gain " << gain << " for OpChannel " << OpChannel << " corresponding to flash index " << pmt << std::endl; */
	    /* if (gain != 0) */
	    /*   m_peSpectrum.at(OpDet)  = PE * 120. / gain; */
	  /* }// or all OpChannels */

	}
      
    /**
     *  @brief  Check if the time of the flash is in the beam window, and save the information for later
     *
     *  @param  beamWindowStart the starting time of the beam window
     *  @param  beamWindowWend the end time of the beam window
     */
    bool IsInBeamWindow(const float beamWindowStart, const float beamWindowEnd)
    {
      m_inBeamWindow = (m_time > beamWindowStart && m_time < beamWindowEnd);
      return m_inBeamWindow;
    }
      
    /**
     *  @brief  Check if the flash passes the minimum PE threshold
     *
     *  @param  minBeamFlashPE the minimum number of photo electrons to pass
     */
    bool PassesPEThreshold(const float minBeamFlashPE) const
    {
      return (m_totalPE > minBeamFlashPE);
    }
      
    /**
     *  @breif  Convert to a flashana::Flash_t
     *
     *  @return the flashana::Flash_t
     */
    flashana::Flash_t ConvertFlashFormat() const
      {
	// Set the flash properties
	flashana::Flash_t flash;
	flash.x = 0;
	flash.x_err = 0;
	flash.y = m_centerY;
	flash.y_err = m_widthY;
	flash.z = m_centerZ;
	flash.z_err = m_widthZ;
	flash.time = m_time;
	flash.pe_v.clear();
	flash.pe_err_v.clear();
       
	for (auto const& pe : m_peSpectrum) {
	  flash.pe_v.push_back( pe );
	  flash.pe_err_v.push_back( std::sqrt(pe) );
	}
	return flash;
      }
      
    // Features of the flash are used when writing to file is enabled
    int                 m_run;                  ///< The run number
    int                 m_subRun;               ///< The subRun number
    int                 m_event;                ///< The event number
    float               m_time;                 ///< Time of the flash
    std::vector<float>  m_peSpectrum;           ///< The number of PEs on each PMT
    std::vector<float>  m_peHypSpectrum;           ///< The number of PEs on each PMT
    float               m_totalPE;              ///< The total number of photoelectrons over all PMTs in the flash
    float               m_centerY;              ///< The PE weighted center Y position of the flash
    float               m_centerZ;              ///< The PE weighted center Z position of the flash
    float               m_widthY;               ///< The PE weighted width of the flash in Y
    float               m_widthZ;               ///< The PE weighted width of the flash in Z
    bool                m_inBeamWindow;         ///< If the flash is in time with the beam window
    bool                m_isBrightestInWindow;  ///< If the flash is the brightest in the event
    bool                m_isBeamFlash;          ///< If the flash has been selected as the beam flash

  }; // end of FlashCandidate class

  /**
   *  @brief  A candidate for the target slice
   */
  class SliceCandidate
  {
  public:
    /**
     *  @brief  Data to describe an amount of charge deposited in a given 3D position
     */
    class Deposition
    {
    public:
      /**
       *  @brief  Default constructor
       *
       *  @param  x the x-component of the charge position
       *  @param  y the z-component of the charge position
       *  @param  z the z-component of the charge position
       *  @param  charge the charge deposited
       *  @param  nPhotons the estimated numer of photons produced
       */
    Deposition(const float x, const float y, const float z, const float charge, const float nPhotons)
      : m_x(x),
	m_y(y),
	m_z(z),
	m_charge(charge),
	m_nPhotons(nPhotons)
	{
	}
	
      float m_x;         ///< The x-component of the charge position
      float m_y;         ///< The z-component of the charge position
      float m_z;         ///< The z-component of the charge position
      float m_charge;    ///< The charge deposited
      float m_nPhotons;  ///< The estimated numer of photons produced
    };
      
    typedef std::vector<Deposition> DepositionVector;
      
    // ---------------------------------------------------------------------------------------------------------------------------------
      
    /**
     *  @brief  Default constructor
     */
    SliceCandidate(){}
      
    /**
     *  @brief  Parametrized constructor
     *
     *  @param  event the art event
     *  @param  slice the slice
     *  @param  pfParticleMap the input mapping from PFParticle ID to PFParticle
     *  @param  pfParticleToSpacePointMap the input mapping from PFParticles to SpacePoints
     *  @param  spacePointToHitMap the input mapping from SpacePoints to Hits
     */
  SliceCandidate(const std::vector< art::Ptr<recob::PFParticle>> &pfp_v,
		 const std::vector< std::vector<art::Ptr<recob::SpacePoint> > > &spacepoint_v_v,
		 const std::vector< std::vector<art::Ptr<recob::Hit> > > &hit_v_v,
		 const float chargeToNPhotonsTrack, 
		 const float chargeToNPhotonsShower,
		 bool applyLifetimeCorr = true)
    : m_sliceId(-std::numeric_limits<int>::max()),
      m_hasDeposition(false),
      m_totalCharge(-std::numeric_limits<float>::max()),
      m_centerX(-std::numeric_limits<float>::max()),
      m_centerY(-std::numeric_limits<float>::max()),
      m_centerZ(-std::numeric_limits<float>::max()),
      m_minX(-std::numeric_limits<float>::max()),
      m_deltaY(-std::numeric_limits<float>::max()),
      m_deltaZ(-std::numeric_limits<float>::max()),
      m_deltaYSigma(-std::numeric_limits<float>::max()),
      m_deltaZSigma(-std::numeric_limits<float>::max()),
      m_chargeToLightRatio(-std::numeric_limits<float>::max()),
      m_xChargeLightVariable(-std::numeric_limits<float>::max()),
      m_passesPrecuts(false),
      m_flashMatchScore(-std::numeric_limits<float>::max()),
      m_totalPEHypothesis(-std::numeric_limits<float>::max()),
      m_isTaggedAsTarget(false),
      m_targetMethod(-std::numeric_limits<int>::max()),
      m_isConsideredByFlashId(false),
      m_hasBestTopologicalScore(false),
      m_hasBestFlashMatchScore(false),
      m_chargeToNPhotonsTrack(chargeToNPhotonsTrack),
      m_chargeToNPhotonsShower(chargeToNPhotonsShower),
      m_xclCoef(-std::numeric_limits<float>::max()),
      m_applyLifetimeCorr(applyLifetimeCorr)
	{
	
	  const auto chargeDeposition(this->GetDepositionVector(pfp_v, spacepoint_v_v, hit_v_v));
	  m_lightCluster = this->GetLightCluster(chargeDeposition);
	  m_totalCharge = this->GetTotalCharge(chargeDeposition);
	  m_hasDeposition = (m_totalCharge > std::numeric_limits<float>::epsilon());
	
	  if (!m_hasDeposition)
	    return;
	
	  const auto chargeCenter(this->GetChargeWeightedCenter(chargeDeposition));
	  m_centerX = chargeCenter.GetX();
	  m_centerY = chargeCenter.GetY();
	  m_centerZ = chargeCenter.GetZ();
	
	  m_minX = this->GetMinimumXPosition(chargeDeposition);
	}
      
  SliceCandidate(const std::vector<art::Ptr<recob::SpacePoint> > &spacepoint_v,
		 const std::vector<art::Ptr<recob::Hit> > &hit_v,
		 const float chargeToNPhotonsTrack, 
		 const float chargeToNPhotonsShower,
		 bool applyLifetimeCorr = true)
    : m_sliceId(-std::numeric_limits<int>::max()),
      m_hasDeposition(false),
      m_totalCharge(-std::numeric_limits<float>::max()),
      m_centerX(-std::numeric_limits<float>::max()),
      m_centerY(-std::numeric_limits<float>::max()),
      m_centerZ(-std::numeric_limits<float>::max()),
      m_minX(-std::numeric_limits<float>::max()),
      m_deltaY(-std::numeric_limits<float>::max()),
      m_deltaZ(-std::numeric_limits<float>::max()),
      m_deltaYSigma(-std::numeric_limits<float>::max()),
      m_deltaZSigma(-std::numeric_limits<float>::max()),
      m_chargeToLightRatio(-std::numeric_limits<float>::max()),
      m_xChargeLightVariable(-std::numeric_limits<float>::max()),
      m_passesPrecuts(false),
      m_flashMatchScore(-std::numeric_limits<float>::max()),
      m_totalPEHypothesis(-std::numeric_limits<float>::max()),
      m_isTaggedAsTarget(false),
      m_targetMethod(-std::numeric_limits<int>::max()),
      m_isConsideredByFlashId(false),
      m_hasBestTopologicalScore(false),
      m_hasBestFlashMatchScore(false),
      m_chargeToNPhotonsTrack(chargeToNPhotonsTrack),
      m_chargeToNPhotonsShower(chargeToNPhotonsShower),
      m_xclCoef(-std::numeric_limits<float>::max()),
      m_applyLifetimeCorr(applyLifetimeCorr)
	{
	
	  const auto chargeDeposition(this->GetDepositionVector(spacepoint_v, hit_v));
	  m_lightCluster = this->GetLightCluster(chargeDeposition);
	  m_totalCharge = this->GetTotalCharge(chargeDeposition);
	  m_hasDeposition = (m_totalCharge > std::numeric_limits<float>::epsilon());
	
	  if (!m_hasDeposition)
	    return;
	
	  const auto chargeCenter(this->GetChargeWeightedCenter(chargeDeposition));
	  m_centerX = chargeCenter.GetX();
	  m_centerY = chargeCenter.GetY();
	  m_centerZ = chargeCenter.GetZ();
	
	  m_minX = this->GetMinimumXPosition(chargeDeposition);
	}

    /**
     *  @breif  Determine if a given slice is compatible with the beam flash by applying pre-selection cuts
     *
     *  @param  beamFlash the beam flash
     *  @param  maxDeltaY the maximum difference in Y between the beam flash center and the weighted charge center
     *  @param  maxDeltaZ the maximum difference in Z between the beam flash center and the weighted charge center
     *  @param  maxDeltaYSigma as for maxDeltaY, but measured in units of the flash width in Y
     *  @param  maxDeltaZSigma as for maxDeltaZ, but measured in units of the flash width in Z
     *  @param  minChargeToLightRatio the minimum ratio between the total charge and the total PE
     *  @param  maxChargeToLightRatio the maximum ratio between the total charge and the total PE
     *
     *  @return if the slice is compatible with the beamFlash
     */
    bool IsCompatibleWithBeamFlash(const FlashCandidate &beamFlash, const float maxDeltaY, const float maxDeltaZ,
				   const float maxDeltaYSigma, const float maxDeltaZSigma, const float minChargeToLightRatio, const float maxChargeToLightRatio)
    {
      // Check the flash is usable
      if (beamFlash.m_totalPE <= std::numeric_limits<float>::epsilon())
	{
	  return false;
	}

      if (beamFlash.m_widthY <= std::numeric_limits<float>::epsilon())
	{
	  return false;
	}
	
      if (beamFlash.m_widthZ <= std::numeric_limits<float>::epsilon())
	{
	  return false;
	}
	
      if (m_totalCharge <= std::numeric_limits<float>::epsilon())
	{
	  return false;
	}
	
      // Calculate the pre-selection variables
      m_deltaY = (m_centerY - beamFlash.m_centerY);
      m_deltaZ = (m_centerZ - beamFlash.m_centerZ);
      m_deltaYSigma = m_deltaY / beamFlash.m_widthY;
      m_deltaZSigma = m_deltaZ / beamFlash.m_widthZ;
      m_chargeToLightRatio = m_totalCharge / beamFlash.m_totalPE;
      m_xChargeLightVariable = m_xclCoef*log10(m_chargeToLightRatio)-m_centerX;
	
      // Check if the slice passes the pre-selection cuts
      m_passesPrecuts = (std::abs(m_deltaY) < maxDeltaY &&
			 std::abs(m_deltaZ) < maxDeltaZ &&
			 std::abs(m_deltaYSigma) < maxDeltaYSigma &&
			 std::abs(m_deltaZSigma) < maxDeltaZSigma &&
			 m_xChargeLightVariable > minChargeToLightRatio &&
			 m_xChargeLightVariable < maxChargeToLightRatio);
	
      return m_passesPrecuts;
    }
      
    /** 
     *  @brief  Get the flash matching score between this slice and the beam flash
     *
     *  @param  beamFlash the beam flash
     *  @param  flashMatchManager the flash matching manager
     *
     *  @return the flash matching score
     */
    float GetFlashMatchScore(FlashCandidate &beamFlash, flashana::FlashMatchManager &flashMatchManager)
    {
      flashMatchManager.Reset();

      // Convert the flash and the charge cluster into the required format for flash matching
      auto flash(beamFlash.ConvertFlashFormat());
	
      // Perform the match
      flashMatchManager.Emplace(std::move(flash));
      flashMatchManager.Emplace(std::move(m_lightCluster));
      const auto matches(flashMatchManager.Match());
	
      // Unable to match
      if (matches.empty())
	return -1.f;
	
      if (matches.size() != 1)
	throw cet::exception("FlashMatchingTool") << "Flash matching returned multiple matches!" << std::endl;
	
      // Fill the slice candidate with the details of the matching
      const auto match(matches.front());
	
      m_flashMatchScore = match.score;
      m_flashMatchX = match.tpc_point.x;
      m_totalPEHypothesis = std::accumulate(match.hypothesis.begin(), match.hypothesis.end(), 0.f);
	
      // Fill the slice with the hypothesized PE spectrum
      if (!m_peHypothesisSpectrum.empty())
	throw cet::exception("FlashMatchingTool") << "Hypothesized PE spectrum already set for this flash" << std::endl;
	
      art::ServiceHandle<geo::Geometry> geo;
      const auto nOpDets(geo->NOpDets());

      if (match.hypothesis.size() != nOpDets)
	throw cet::exception("FlashMatchingTool") << "Hypothesized PE spectrum has the wrong size" << std::endl;
	
      beamFlash.m_peHypSpectrum.clear();
      for (unsigned int i = 0; i < nOpDets; ++i){
	m_peHypothesisSpectrum.push_back(static_cast<float>(match.hypothesis.at(i)));
	beamFlash.m_peHypSpectrum.push_back(static_cast<float>(match.hypothesis.at(i)));
      }
	
      return m_flashMatchScore;
    }
      
      
  private:
    /**
     *  @breif  Get the 3D spacepoints (with charge) associated with the PFParticles in the slice that are produced from hits in the W view
     *
     *  @param  pfParticleMap the input mapping from PFParticle ID to PFParticle
     *  @param  pfParticleToSpacePointMap the input mapping from PFParticles to SpacePoints
     *  @param  spacePointToHitMap the input mapping from SpacePoints to Hits
     *  @param  slice the input slice
     *
     *  @return the output depositionVector
     */
    DepositionVector GetDepositionVector(const std::vector< art::Ptr<recob::PFParticle> > &pfp_v,
					 const std::vector< std::vector<art::Ptr<recob::SpacePoint> > > &spacepoint_v_v,
					 const std::vector< std::vector<art::Ptr<recob::Hit> > > &hit_v_v ) const
    {
      // Collect all PFParticles in the slice, including those downstream of the primaries
      // ATTN here we only use the neutrino hypothesis, in theory this should work with either (or indeed both with some thought)
      //PFParticleVector allParticlesInSlice;
      //this->CollectDownstreamPFParticles(pfParticleMap, slice.GetTargetHypothesis(), allParticlesInSlice);

      //--------------------------------------------------------------------
      // implementing electron lifetime correction [D. Caratelli 08/12/2022]
      const detinfo::DetectorProperties* detprop;
      detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();

      //handle to electron lifetime calibration provider
      const lariov::UBElectronLifetimeProvider& elifetimeCalibProvider
	= art::ServiceHandle<lariov::UBElectronLifetimeService>()->GetProvider();

      float elifetime  = elifetimeCalibProvider.Lifetime(); // [ms]
      float driftvelocity = detprop->DriftVelocity(); // [cm/us]

      //std::cout << "LIFETIMECORRECTION [FlashNeutrinoId][GetDepositionVector] lifetime is : "
      //      << elifetime << " [ms] and drift velocity is " << driftvelocity << " [cm/us]" << std::endl;
      // implementing electron lifetime correction [D. Caratelli 08/12/2022]
      //--------------------------------------------------------------------
	
      DepositionVector depositionVector;
	
      for (size_t p=0; p < pfp_v.size(); p++){
	  
	//auto pfp_ptr = pfp_v.at(p);
	  
	// Get the associated spacepoints
	auto spacepoint_ptr_v = spacepoint_v_v.at(p);
	auto hit_ptr_v = hit_v_v.at(p);
	  
	if (spacepoint_ptr_v.size() != hit_ptr_v.size()) 
	  std::cout << "HELP! different number of points and hits!" << std::endl;
	  
	for (size_t h=0; h < hit_ptr_v.size(); h++) {
	    
	  auto hit = hit_ptr_v.at(h);
	  auto spacepoint = spacepoint_ptr_v.at(h);
	    
	  if (hit->View() != geo::kZ)
	    continue;
	    
	  // Add the charged point to the vector
	  const auto &position(spacepoint->XYZ());
	  const auto charge(hit->Integral());
	  //------------------------------------------------------
	  // implement lifetime correction [D. Caratelli 08/12/22]
	  float lifetimecorrection = 1.;
	  if (m_applyLifetimeCorr) lifetimecorrection = exp( (position[0]) / (elifetime * driftvelocity * 1000.0));
	  //std::cout << "LIFETIMECORRECTION [FlashNeutrinoId][GetDepositionVector]: lifetime correction @ lifetime of " << elifetime << " [ms] "
	  //      << "@ position of " << position[0] << " [cm] is " << lifetimecorrection << std::endl;
	  // implement lifetime correction [D. Caratelli 08/12/22]
	  //------------------------------------------------------
	
	  // ADC -> MeV
	  //float ADCtoMeV = 240. * (23.6/1e6) / 0.5;

	  //depositionVector.emplace_back(position[0], position[1], position[2], charge, this->GetNPhotons(charge, particle));
	  depositionVector.emplace_back(position[0], position[1], position[2], charge * lifetimecorrection, charge * lifetimecorrection * m_chargeToNPhotonsTrack);
	  //std::cout << "adding to deposition vector @ location [" << position[0] << ", " << position[1] << ", " << position[2] 
	  //	      << " with charge " << charge * m_chargeToNPhotonsTrack << std::endl;
	}// for all hits/spacepoints
      }
      return depositionVector;
    }

    DepositionVector GetDepositionVector(const std::vector<art::Ptr<recob::SpacePoint> > &spacepoint_v,
					 const std::vector<art::Ptr<recob::Hit> > &hit_v ) const
    {
      // Collect all PFParticles in the slice, including those downstream of the primaries
      // ATTN here we only use the neutrino hypothesis, in theory this should work with either (or indeed both with some thought)
      //PFParticleVector allParticlesInSlice;
      //this->CollectDownstreamPFParticles(pfParticleMap, slice.GetTargetHypothesis(), allParticlesInSlice);

      //--------------------------------------------------------------------
      // implementing electron lifetime correction [D. Caratelli 08/12/2022]
      const detinfo::DetectorProperties* detprop;
      detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();

      //handle to electron lifetime calibration provider
      const lariov::UBElectronLifetimeProvider& elifetimeCalibProvider
	= art::ServiceHandle<lariov::UBElectronLifetimeService>()->GetProvider();

      float elifetime  = elifetimeCalibProvider.Lifetime(); // [ms]
      float driftvelocity = detprop->DriftVelocity(); // [cm/us]

      //std::cout << "LIFETIMECORRECTION [FlashNeutrinoId][GetDepositionVector] lifetime is : "
      //      << elifetime << " [ms] and drift velocity is " << driftvelocity << " [cm/us]" << std::endl;
      // implementing electron lifetime correction [D. Caratelli 08/12/2022]
      //--------------------------------------------------------------------

      DepositionVector depositionVector;
	
      if (spacepoint_v.size() != hit_v.size()) 
	std::cout << "HELP! different number of points and hits!" << std::endl;
	
      for (size_t h=0; h < hit_v.size(); h++) {
	  
	auto hit = hit_v.at(h);
	auto spacepoint = spacepoint_v.at(h);
	  
	if (hit->View() != geo::kZ)
	  continue;
	  
	// Add the charged point to the vector
	const auto &position(spacepoint->XYZ());
	const auto charge(hit->Integral());
	//------------------------------------------------------
	// implement lifetime correction [D. Caratelli 08/12/22]
	float lifetimecorrection = 1.;
	if (m_applyLifetimeCorr) lifetimecorrection = exp( (position[0]) / (elifetime * driftvelocity * 1000.0));
	//std::cout << "LIFETIMECORRECTION [FlashNeutrinoId][GetDepositionVector]: lifetime correction @ lifetime of " << elifetime << " [ms] "
	//      << "@ position of " << position[0] << " [cm] is " << lifetimecorrection << std::endl;
	// implement lifetime correction [D. Caratelli 08/12/22]
	//------------------------------------------------------
	  
	//depositionVector.emplace_back(position[0], position[1], position[2], charge, this->GetNPhotons(charge, particle));
	depositionVector.emplace_back(position[0], position[1], position[2], charge * lifetimecorrection, charge * lifetimecorrection * m_chargeToNPhotonsTrack);
      }// for all hits/spacepoints

      return depositionVector;
    }
    
    /**
     *  @brief  Convert from deposited charge to number of photons for a given particle
     *
     *  @param  charge the input charge
     *  @param  particle the input particle
     *
     *  @return the number of photons
     */
    float GetNPhotons(const float charge, const art::Ptr<recob::PFParticle> &particle) const
    {
      return charge * m_chargeToNPhotonsTrack; //* (LArPandoraHelper::IsTrack(particle) ? m_chargeToNPhotonsTrack : m_chargeToNPhotonsShower);
    }
      
    /**
     *  @brief  Get the centroid of the input charge cluster, weighted by charge
     *
     *  @param  depositionVector the input charge cluster
     *
     *  @return the charge weighted centroid
     */
    pandora::CartesianVector GetChargeWeightedCenter(const DepositionVector &depositionVector) const
      {
	pandora::CartesianVector center(0.f, 0.f, 0.f);
	float totalCharge(0.f);
	
	for (const auto &chargePoint : depositionVector)
	  {
	    center += pandora::CartesianVector(chargePoint.m_x, chargePoint.m_y, chargePoint.m_z) * chargePoint.m_charge;
	    totalCharge += chargePoint.m_charge;
	  }
	
	if (totalCharge <= std::numeric_limits<float>::epsilon())
	  throw cet::exception("FlashMatchingTool") << "Can't find charge weighted center of slice with zero total charge" << std::endl;

	center *= (1.f / totalCharge);
	
	return center;
      }
      
    /**
     *  @brief  Get the total charge from an input charge cluster
     *
     *  @param  depositionVector the input charge cluster
     *
     *  @return the total charge
     */
    float GetTotalCharge(const DepositionVector &depositionVector) const
    {
      float totalCharge(0.f);
	
      for (const auto &chargePoint : depositionVector)
	totalCharge += chargePoint.m_charge;
	
      return totalCharge;
    }
      
    /**
     *  @brief  Get the minimum X-position of all deposition points given
     *
     *  @param  depositionVector the input charge cluster
     *
     *  @return the minimum X-position
     */
    float GetMinimumXPosition(const DepositionVector &depositionVector) const
    {
      float minX(std::numeric_limits<float>::max());
	
      for (const auto &chargePoint : depositionVector)
	minX = std::min(chargePoint.m_x, minX);
	
      return minX;
    }
      
    /**
     *  @brief  Convert a charge deposition into a light cluster by applying the chargeToPhotonFactor to every point
     *
     *  @param  depositionVector the input charge cluster
     *
     *  @return the output light cluster
     */
    flashana::QCluster_t GetLightCluster(const DepositionVector &depositionVector) const
      {
	flashana::QCluster_t lightCluster;
	
	for (const auto &point : depositionVector)
	  lightCluster.emplace_back(point.m_x, point.m_y, point.m_z, point.m_nPhotons);
	
	return lightCluster;
      }
      
  public:
    // Features of the slice are used when writing to file is enabled
    int                  m_sliceId;                  ///< The sliceId   
    int                  m_run;                      ///< The run number
    int                  m_subRun;                   ///< The subRun number
    int                  m_event;                    ///< The event number
    bool                 m_hasDeposition;            ///< If the slice has any charge deposited on the collection plane which produced a spacepoint
    float                m_totalCharge;              ///< The total charge deposited on the collection plane by hits that produced spacepoints
    float                m_centerX;                  ///< The charge weighted center of the slice in X
    float                m_centerY;                  ///< The charge weighted center of the slice in Y
    float                m_centerZ;                  ///< The charge weighted center of the slice in Z
    float                m_minX;                     ///< The minimum X-coordinate of all spacepoints in the slice
    float                m_deltaY;                   ///< The distance of the slice centroid from the flash centroid in Y
    float                m_deltaZ;                   ///< The distance of the slice centroid from the flash centroid in Z
    float                m_deltaYSigma;              ///< deltaY but in units of the flash width in Y
    float                m_deltaZSigma;              ///< deltaZ but in units of the flash width in Z
    float                m_chargeToLightRatio;       ///< The ratio between the total charge and the total PE of the beam flash
    float                m_xChargeLightVariable;     ///< m_xclCoef*log10(chargeToLightRatio)- centerX
    bool                 m_passesPrecuts;            ///< If the slice passes the preselection cuts
    float                m_flashMatchScore;          ///< The flash matching score between the slice and the beam flash
    float                m_flashMatchX;              ///< The etimated X coordinate of the flashmatching
    float                m_totalPEHypothesis;        ///< The total PE of the hypothesized flash for this slice
    std::vector<float>   m_peHypothesisSpectrum;     ///< The PE of the hypothesized flash of this slice 
    bool                 m_isTaggedAsTarget;         ///< If the slice has been tagged as the target (neutrino)
    int                  m_targetMethod;             ///< 0: only one slice passed precuts, 1: has best toposcore, 2: had best flashmatchscore
    bool                 m_isConsideredByFlashId;    ///< If the slice was considered by the flash ID tool - this will be false if there wasn't a beam flash found in the event
    float                m_topologicalNeutrinoScore; ///< The topological-information-only neutrino ID score from Pandora
    bool                 m_hasBestTopologicalScore;  ///< If this slice has the highest topological score in the event
    bool                 m_hasBestFlashMatchScore;   ///< From the slices passing the precuts, if this one has the highest flashmatch score
    float                m_chargeToNPhotonsTrack;    ///< The conversion factor between charge and number of photons for tracks
    float                m_chargeToNPhotonsShower;   ///< The conversion factor between charge and number of photons for showers
    float                m_xclCoef;                  ///< m_xclCoef*log10(chargeToLightRatio)- centerX
    flashana::QCluster_t m_lightCluster;             ///< The hypothesised light produced - used by flashmatching
    bool                 m_applyLifetimeCorr;        ///< Whether we want the lifetime correction to be applied when computing the deposition vector

  }; // end of SliceCandidate class

}//end of namespace flashmatch

#endif
