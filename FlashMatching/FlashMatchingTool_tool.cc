#include "FlashMatchingToolBase_tool.h"
#include "SliceCandidate.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

namespace flashmatch {
  
  class FlashMatchingTool : public FlashMatchingToolBase {
    
  public:
    
    /**
     *  @brief  Default constructor
     *
     *  @param  pset FHiCL parameter set
     */
    explicit FlashMatchingTool(const fhicl::ParameterSet& pset)
      : m_pFlashTree(nullptr)
    {
      std::cout << "Called constructor" << std::endl;
      configure(pset);
    }
    
    // default destructor
    ~FlashMatchingTool(){};
    
    void configure(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Classify slices as neutrino or cosmic
     *
     *  @param  slices the input vector of slices to classify
     *  @param  evt the art event
     */
    float ClassifySlice(const art::Event &evt,
			const std::vector< art::Ptr<recob::PFParticle> > &pfp_v,
			const std::vector< std::vector<art::Ptr< recob::SpacePoint> > > &spacepoint_v_v,
			const std::vector< std::vector<art::Ptr<recob::Hit> > > &hit_v_v,
			std::vector<float>& recospectrum,
			std::vector<float>& hypospectrum);


    float ClassifyTrack(const art::Event &evt,
			const std::vector<art::Ptr< recob::SpacePoint> > &spacepoint_v,
			const std::vector<art::Ptr<recob::Hit> > &hit_v,
			std::vector<float>& recospectrum,
			std::vector<float>& hypospectrum);


  private:

    //float Classify(const art::Event &evt, const flashmatch::FlashMatchingTool::SliceCandidate& slice);
    
    // -------------------------------------------------------------------------------------------------------------------------------------
    /**
     *  @brief  A description of the reason the tool couldn't find a neutrino candidate
     */
    class FailureMode
    {
    public:
      /**
       *  @brief  Default constructor
       *
       *  @param  reason the reason for the failure
       */
      FailureMode(const std::string &reason)
	: m_reason(reason)
      {}
      
      /**
       *  @brief  Default destructor - explains the failure
       */
      ~FailureMode()
      {
	std::cout << "[Flash neutrino ID] Failed to find neutrino slice: ";
	std::cout << m_reason << std::endl
		  << std::endl;
      }

      void Print()
      {
	std::cout << "[Flash neutrino ID] Failed to find neutrino slice: ";
	std::cout << m_reason << std::endl
		  << std::endl;
      }
      
    private:
      std::string m_reason;  ///< The reason for the failure
    };
    
    // -------------------------------------------------------------------------------------------------------------------------------------
    
    /**
     *  @brief  Class to hold information about the event for monitoring
     */
    class OutputEvent
    {
    public:
      /**
       *  @brief  Reset the variables to default dummy values
       *
       *  @param  event the art event
       */
      void Reset(const art::Event &event)
      {
	m_run = event.run();
	m_subRun = event.subRun();
	m_event = event.event();
	m_nFlashes = -std::numeric_limits<int>::max();
	m_nFlashesInBeamWindow = -std::numeric_limits<int>::max();
	m_hasBeamFlash = false;
	m_nSlices = -std::numeric_limits<int>::max();
	m_nSlicesAfterPrecuts = -std::numeric_limits<int>::max();
	m_foundATargetSlice = false;
	m_targetSliceMethod = -1;
      }
      
      
      int  m_run;                   ///< The run number
      int  m_subRun;                ///< The subRun number
      int  m_event;                 ///< The event number
      int  m_nFlashes;              ///< The number of flashes
      int  m_nFlashesInBeamWindow;  ///< The number of flashes in the beam window
      bool m_hasBeamFlash;          ///< If a beam flash was found
      int  m_nSlices;               ///< The number of slices
      int  m_nSlicesAfterPrecuts;   ///< The number of slices remaining after the preselection cuts
      bool m_foundATargetSlice;     ///< If a slice was identified as the target (neutrino)
      int  m_targetSliceMethod;     ///< 0: only one slice passed precuts, 1: has best toposcore, 2: has best flashmatchscore
    };    
        
    // -------------------------------------------------------------------------------------------------------------------------------------
    
    typedef std::vector<SliceCandidate> SliceCandidateVector;
    typedef std::vector<FlashCandidate> FlashCandidateVector;
    
    /**
     *  @brief  Get the candidate flashes in the event
     *
     *  @param  event the art event
     *  @param  flashCandidates the output vector of flash candidates
     */
    void GetFlashCandidates(const art::Event &event, FlashCandidateVector &flashCandidates)
    {

      // Collect all flashes from the event
      art::InputTag flashTag(m_flashLabel);
      const auto flashes(*event.getValidHandle<std::vector<recob::OpFlash>>(flashTag));

      for (const auto &flash : flashes)
	flashCandidates.emplace_back(event, flash);
      
      m_outputEvent.m_nFlashes = flashCandidates.size();
    }
    
    /**
     *  @breif  Try to find the brightest flash with sufficent photoelectons that is in time with the beam
     *
     *  @param  flashCandidates the input vector of slice candidates
     *
     *  @return the beam flash
     */
    FlashCandidate& GetBeamFlash(FlashCandidateVector &flashCandidates)
    {

      bool foundFlashInBeamWindow(false);
      unsigned int brightestFlashIndex(std::numeric_limits<unsigned int>::max());
      float maxTotalPE(-std::numeric_limits<float>::max());
      m_outputEvent.m_nFlashesInBeamWindow = 0;

      // Find the brightest flash in the beam window
      for (unsigned int flashIndex = 0; flashIndex < flashCandidates.size(); ++flashIndex)
	{
	  // ATTN non const reference is required since monitoring variables are stored in the slice candidate
	  auto &flashCandidate(flashCandidates.at(flashIndex));
	  
	  if (!flashCandidate.IsInBeamWindow(m_beamWindowStart, m_beamWindowEnd))
	    continue;
	  
	  m_outputEvent.m_nFlashesInBeamWindow++;
	  
	  const auto totalPE(flashCandidate.m_totalPE);
	  if (totalPE < maxTotalPE)
	    continue;

	  foundFlashInBeamWindow = true;
	  maxTotalPE = totalPE;
	  brightestFlashIndex = flashIndex;
	}
      
      if (!foundFlashInBeamWindow)
	throw FailureMode("There were no flashes in the beam window");
      
      // Ensure it is sufficiently bright
      auto &brightestFlash(flashCandidates.at(brightestFlashIndex));
      brightestFlash.m_isBrightestInWindow = true;
      
      if (!brightestFlash.PassesPEThreshold(m_minBeamFlashPE))
	throw FailureMode("No flashes in the beam window passed the PE threshold");
      
      // Save the monitoring information
      brightestFlash.m_isBeamFlash = true;
      m_outputEvent.m_hasBeamFlash = true;
      
      return brightestFlash;
    }
    
    /**
     *  @brief  get the score
     *
     *  @param  beamFlash the beam flash
     *  @param  sliceCandidates the neutrino slice candidates
     */
    float GetSliceScore(FlashCandidate &beamFlash, SliceCandidate &slice)
    {
      //bool foundViableSlice(false);
      //bool foundHighestTopoligicalScore(false);
      //unsigned int bestFlashMatchSliceIndex(std::numeric_limits<unsigned int>::max());
      //unsigned int bestCombinedSliceIndex(std::numeric_limits<unsigned int>::max());
      m_outputEvent.m_nSlicesAfterPrecuts = 0;
      
      /*
      // Apply the pre-selection cuts to ensure that the slice is compatible with the beam flash
      if (!slice.IsCompatibleWithBeamFlash(beamFlash, m_maxDeltaY, m_maxDeltaZ, m_maxDeltaYSigma, m_maxDeltaZSigma,
      m_minChargeToLightRatio, m_maxChargeToLightRatio)) {
      
      std::cout << "Not compatible with beam flash!" << std::endl;
      return 0;
      }
      */
      
      m_outputEvent.m_nSlicesAfterPrecuts++;
      
      // ATTN if there is only one slice that passes the pre-selection cuts, then the score won't be used
      const auto &score(slice.GetFlashMatchScore(beamFlash, m_flashMatchManager));
      
      //foundViableSlice = true;
      
      return score;
    }
    
  private:
    
    std::string  m_flashLabel;    ///< The label of the flash producer
    std::string  m_pandoraLabel;  ///< The label of the allOutcomes pandora producer
    
    // Cuts for selecting the beam flash
    float        m_beamWindowStart;  ///< The start time of the beam window
    float        m_beamWindowEnd;    ///< The end time of the beam window
    float        m_minBeamFlashPE;   ///< The minimum number of photoelectrons required to consider a flash as the beam flash
    
    // Coefficient to account for the x-dependency in the charge light-ratio
    //float        m_xclCoef;          ///< m_xclCoef*log10(chargeToLightRatio)- centerX

    // Pre-selection cuts to determine if a slice is compatible with the beam flash
    float        m_maxDeltaY;              ///< The maximum difference in Y between the beam flash center and the weighted charge center
    float        m_maxDeltaZ;              ///< The maximum difference in Z between the beam flash center and the weighted charge center
    float        m_maxDeltaYSigma;         ///< As for maxDeltaY, but measured in units of the flash width in Y
    float        m_maxDeltaZSigma;         ///< As for maxDeltaZ, but measured in units of the flash width in Z
    float        m_minChargeToLightRatio;  ///< The minimum ratio between the total charge and the total PE
    float        m_maxChargeToLightRatio;  ///< The maximum ratio between the total charge and the total PE

    // Variables required for flash matching
    float                          m_chargeToNPhotonsTrack;   ///< The conversion factor between charge and number of photons for tracks
    float                          m_chargeToNPhotonsShower;  ///< The conversion factor between charge and number of photons for showers
    flashana::FlashMatchManager    m_flashMatchManager;       ///< The flash match manager

    // Debugging / testing
    bool                                    m_shouldWriteToFile;   ///< If we should write interesting information to a root file
    bool                                    m_hasMCNeutrino;       ///< If there is an MC neutrino we can use to get truth information
    //int                                     m_nuInteractionType;   ///< The interaction type code from MCTruth
    //int                                     m_nuCCNC;                ///< Charged current or neutral current?
    //float                                   m_nuEnergy;            ///< The true neutrino energy
    //float                                   m_leptonEnergy;        ///< The true energy of the lepton coming from the CC interaction
    //float                                   m_nuVertexX;           ///< The true neutrino vertex X position
    //float                                   m_nuVertexY;           ///< The true neutrino vertex Y position
    //float                                   m_nuVertexZ;           ///< The true neutrino vertex Z position
    //float                                   m_nuTime;              ///< The time of the true neutrino interaction
    //int                                     m_nuPdgCode;           ///< The true neutrino pdg code
    std::string                             m_truthLabel;          ///< The MCTruth producer label
    std::string                             m_mcParticleLabel;     ///< The MCParticle producer label
    std::string                             m_hitLabel;            ///< The Hit producer label
    std::string                             m_backtrackLabel;      ///< The Hit -> MCParticle producer label
    OutputEvent                             m_outputEvent;         ///< The output event whose address is used by the output branch
    FlashCandidate                          m_outputFlash;         ///< The output flash whose address is used by the output branch
    SliceCandidate                          m_outputSlice;         ///< The output slice whose address is used by the output branch

    TTree                                  *m_pFlashTree;          ///< The flash tree

  };
  

void FlashMatchingTool::configure(const fhicl::ParameterSet& pset) {
    m_flashLabel = pset.get<std::string>("FlashLabel");
    m_pandoraLabel = pset.get<std::string>("PandoraAllOutcomesLabel");
    m_beamWindowStart = pset.get<float>("BeamWindowStartTime");
    m_beamWindowEnd = pset.get<float>("BeamWindowEndTime");
    m_minBeamFlashPE = pset.get<float>("BeamFlashPEThreshold");
    m_maxDeltaY = pset.get<float>("MaxDeltaY");
    m_maxDeltaZ = pset.get<float>("MaxDeltaZ");
    m_maxDeltaYSigma = pset.get<float>("MaxDeltaYSigma");
    m_maxDeltaZSigma = pset.get<float>("MaxDeltaZSigma");
    m_minChargeToLightRatio = pset.get<float>("MinChargeToLightRatio");
    m_maxChargeToLightRatio = pset.get<float>("MaxChargeToLightRatio");
    m_chargeToNPhotonsTrack = pset.get<float>("ChargeToNPhotonsTrack");
    m_chargeToNPhotonsShower = pset.get<float>("ChargeToNPhotonsShower");
    m_shouldWriteToFile = pset.get<bool>("ShouldWriteToFile", false);
    m_hasMCNeutrino = m_shouldWriteToFile ? pset.get<bool>("HasMCNeutrino") : false;
    m_truthLabel = m_hasMCNeutrino ? pset.get<std::string>("MCTruthLabel") : "";
    m_mcParticleLabel = m_hasMCNeutrino ? pset.get<std::string>("MCParticleLabel") : "";
    m_hitLabel = m_hasMCNeutrino ? pset.get<std::string>("HitLabel") : "";
    m_backtrackLabel = m_hasMCNeutrino ? pset.get<std::string>("BacktrackerLabel") : "";
    m_flashMatchManager.Configure(pset.get<flashana::Config_t>("FlashMatchConfig")); 


    // Set up the output branches
    art::ServiceHandle<art::TFileService> fileService;

    m_pFlashTree = fileService->make<TTree>("flashes","");
    m_pFlashTree->Branch("run"                , &m_outputFlash.m_run                , "run/I");
    m_pFlashTree->Branch("subRun"             , &m_outputFlash.m_subRun             , "subRun/I");
    m_pFlashTree->Branch("event"              , &m_outputFlash.m_event              , "event/I");
    m_pFlashTree->Branch("time"               , &m_outputFlash.m_time               , "time/F");
    m_pFlashTree->Branch("centerY"            , &m_outputFlash.m_centerY            , "centerY/F");
    m_pFlashTree->Branch("centerZ"            , &m_outputFlash.m_centerZ            , "centerZ/F");
    m_pFlashTree->Branch("widthY"             , &m_outputFlash.m_widthY             , "widthY/F");
    m_pFlashTree->Branch("widthZ"             , &m_outputFlash.m_widthZ             , "widthZ/F");
    m_pFlashTree->Branch("totalPE"            , &m_outputFlash.m_totalPE            , "totalPE/F");
    m_pFlashTree->Branch("peSpectrum"         , "std::vector< float >"              , &m_outputFlash.m_peSpectrum);
    m_pFlashTree->Branch("peHypothesisSpectrum"   , "std::vector< float >"                   , &m_outputFlash.m_peHypSpectrum);
    m_pFlashTree->Branch("inBeamWindow"       , &m_outputFlash.m_inBeamWindow       , "inBeamWindow/O");
    m_pFlashTree->Branch("isBrightestInWindow", &m_outputFlash.m_isBrightestInWindow, "isBrightestInWindow/O");
    m_pFlashTree->Branch("isBeamFlash"        , &m_outputFlash.m_isBeamFlash        , "isBeamFlash/O");

    return;
  }

  float FlashMatchingTool::ClassifySlice(const art::Event &evt,
					 const std::vector< art::Ptr<recob::PFParticle> > &pfp_v,
					 const std::vector< std::vector<art::Ptr< recob::SpacePoint> > > &spacepoint_v_v,
					 const std::vector< std::vector<art::Ptr<recob::Hit> > > &hit_v_v,
					 std::vector<float>& recospectrum,
					 std::vector<float>& hypospectrum)

  {
    // Reset the output addresses in case we are writing monitoring details to an output file
    //m_outputEvent.Reset(evt);
    
    // create slice candidate
    SliceCandidate slice(pfp_v, spacepoint_v_v, hit_v_v, m_chargeToNPhotonsTrack, m_chargeToNPhotonsShower);

    //return Classify(evt, slice);

    //return Classify(evt, slice);

      FlashCandidateVector flashCandidates;
      FlashCandidate beamFlash;
      
      float score = 0;
      
      try
	{
	  // Find the flash, if any, in time with the beam with the largest number of photoelectrons that is sufficiently bright
	  this->GetFlashCandidates(evt, flashCandidates);
	  beamFlash = this->GetBeamFlash(flashCandidates);
	}
      catch (...)//const FailureMode &)
	{
	  return score;
	  //std::cout << "Failure!" << std::endl;
	}
      try
	{
	  if (m_outputEvent.m_hasBeamFlash)
	    {
	      // get score for match
	      score = this->GetSliceScore(beamFlash, slice);

	      recospectrum = beamFlash.m_peSpectrum;
	      hypospectrum = beamFlash.m_peHypSpectrum;

	      if (m_pFlashTree){
		m_outputFlash = beamFlash;
		m_pFlashTree->Fill();
	      }

	    }
	  else
	    {
	      std::cout << "No beam flash!" << std::endl;
	      slice.m_isConsideredByFlashId = 0;
	    }
	}
      catch (const FailureMode &)
	{
	  //std::cout << "Failure 2!" << std::endl;
	}
      
      
      return score;

  }


  float FlashMatchingTool::ClassifyTrack(const art::Event &evt,
					 const std::vector<art::Ptr< recob::SpacePoint> > &spacepoint_v,
					 const std::vector<art::Ptr<recob::Hit> > &hit_v,
					 std::vector<float>& recospectrum,
					 std::vector<float>& hypospectrum)


  {
    // Reset the output addresses in case we are writing monitoring details to an output file
    //m_outputEvent.Reset(evt);
    
    // create slice candidate
    SliceCandidate slice(spacepoint_v, hit_v, m_chargeToNPhotonsTrack, m_chargeToNPhotonsShower);

      FlashCandidateVector flashCandidates;
      FlashCandidate beamFlash;
      
      float score = 0;
      
      try
	{
	  // Find the flash, if any, in time with the beam with the largest number of photoelectrons that is sufficiently bright
	  this->GetFlashCandidates(evt, flashCandidates);
	  beamFlash = this->GetBeamFlash(flashCandidates);
	}
      catch ( FailureMode &exc)
	{
	  exc.Print();
	  return score;
	  //std::cout << "Failure!" << std::endl;
	}
      try
	{
	  if (m_outputEvent.m_hasBeamFlash)
	    {
	      // get score for match
	      score = this->GetSliceScore(beamFlash, slice);

	      recospectrum = beamFlash.m_peSpectrum;
	      hypospectrum = beamFlash.m_peHypSpectrum;

	      if (m_pFlashTree){
		m_outputFlash = beamFlash;
		m_pFlashTree->Fill();
	      }
	    }
	  else
	    {
	      std::cout << "No beam flash!" << std::endl;
	      slice.m_isConsideredByFlashId = 0;
	    }
	}
      catch (const FailureMode &)
	{
	  //std::cout << "Failure 2!" << std::endl;
	}
      


      
      return score;

    //return Classify(evt, slice);
  }

    
    
  DEFINE_ART_CLASS_TOOL(FlashMatchingTool)  
} // end namespace


