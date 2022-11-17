#ifndef ANALYSIS_EVENTWEIGHTTREE_CXX
#define ANALYSIS_EVENTWEIGHTTREE_CXX

#include "AnalysisToolBase.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"

namespace analysis
{
  ////////////////////////////////////////////////////////////////////////
  //
  // Class:       EventWeightTree
  // File:        EventWeightTree_tool.cc
  //
  //              Provide systematic uncertainty and spline event weights
  //
  // Configuration parameters:
  //
  // TBD
  //
  // Created by Sebastien Prince (sprince@fas.harvard.edu) on 05/13/2019
  //
  // 
  // 
  ////////////////////////////////////////////////////////////////////////
  class EventWeightTree : public AnalysisToolBase{

    public:
      /**
       *  @brief  Constructor
       *
       *  @param  pset
       */
      EventWeightTree(const fhicl::ParameterSet &pset);

      /**
       *  @brief  Destructor
       */
      ~EventWeightTree(){};

      // provide for initialization
      void configure(fhicl::ParameterSet const & pset);

      /**
       * @brief Analysis function
       */
      void analyzeEvent(art::Event const& e, bool fData) override;

      /**
       * @brief Analyze slice
       */
      void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

      /**
       * @brief set branches for TTree
       */
      void setBranches(TTree* _tree) override;

      /**
       * @brief reset ttree branches
       */
      void resetTTree(TTree* _tree) override;

    private:
      TTree *_weightstree;
      std::map<std::string, std::vector<double>> _mapWeight;
      std::vector<unsigned short> _vecWeightFlux;
      std::vector<double> _vecWeightFluxD;
      std::vector<unsigned short> _vecWeightsGenie;
      std::vector<double> _vecWeightsGenieD;
      std::vector<double> _vecWeightsGenie_vec;
      std::vector<unsigned short> _vecWeightsPPFX;
      std::vector<double> _vecWeightsPPFXD;
      std::vector<int> _vecWeightsGenie_nam;
      std::vector<unsigned short> _vecWeightsReint;
      std::vector<double> _vecWeightsReintD;
      std::vector<unsigned short> _vecWeightsGenieUp;//new check please FIXME
      std::vector<unsigned short> _vecWeightsGenieDn;
      double _knobRPAup;
      double _knobCCMECup;
      double _knobAxFFCCQEup;
      double _knobVecFFCCQEup;
      double _knobDecayAngMECup;
      double _knobThetaDelta2Npiup;
      double _knobThetaDelta2NRadup;
      //double _knobRPA_CCQE_Reducedup;
      double _knobNormCCCOHup;
      double _knobNormNCCOHup;
      double _knobxsr_scc_Fv3up;
      double _knobxsr_scc_Fa3up;
      double _knobRPAdn;
      double _knobCCMECdn;
      double _knobAxFFCCQEdn;
      double _knobVecFFCCQEdn;
      double _knobDecayAngMECdn;
      double _knobThetaDelta2Npidn;
      double _knobThetaDelta2NRaddn;
      //double _knobRPA_CCQE_Reduceddn;
      double _knobNormCCCOHdn;
      double _knobNormNCCOHdn;
      double _knobxsr_scc_Fv3dn;
      double _knobxsr_scc_Fa3dn;
      float _weightSpline;
      float _weightTune;
      float _weightSplineTimesTune;
      float _ppfx_cv;
      double _RootinoFix;
      bool _createDedicatedTree;
      bool _createMapBranch;
      bool _createFluxBranch;
      bool _createGenieBranch;
      bool _createReintBranch;
      bool _createSplineBranch;
      bool _createTuneBranch;
      bool _createSplineTimesTuneBranch;
      bool _createPPFXBranch;
      bool _SaveAllFlux;
      int _GenieAllUniverses;
      int _run;
      int _subRun;
      int _evt;
      bool fMakeNuMINtuple;
  };

  EventWeightTree::EventWeightTree(const fhicl::ParameterSet &p){
    _createDedicatedTree = p.get<bool>("createDedicatedTree");
    _createMapBranch = p.get<bool>("createMapBranch");
    _createFluxBranch = p.get<bool>("createFluxBranch");
    _createGenieBranch = p.get<bool>("createGenieBranch");
    _createReintBranch = p.get<bool>("createReintBranch");
    _createSplineBranch = p.get<bool>("createSplineBranch");
    _createTuneBranch = p.get<bool>("createTuneBranch");
    _createSplineTimesTuneBranch = p.get<bool>("createSplineTimesTuneBranch");
    _createPPFXBranch = p.get<bool>("createPPFXBranch",false);
    _SaveAllFlux = p.get<bool>("SaveAllFlux",false);
    _createGenieUpDnVecs = p.get<bool>("createGenieUpDnVecs", false);
    _GenieAllUniverses = p.get<int>("GenieAllUniverses",500);
    fMakeNuMINtuple = p.get<bool>("makeNuMINtuple", false);
    
    if(_createDedicatedTree){
      art::ServiceHandle<art::TFileService> tfs;
      _weightstree = tfs->make<TTree>("EventWeights", "EventWeights TTree");
    }
  }

  void EventWeightTree::configure(fhicl::ParameterSet const & p){
  }

  void EventWeightTree::analyzeEvent(art::Event const& evt, bool fData){

    _run = evt.run();
    _subRun = evt.subRun();
    _evt = evt.event();

    _vecWeightsGenie  = std::vector<unsigned short>(_GenieAllUniverses,1);
    _vecWeightsGenieD = std::vector<double>(_GenieAllUniverses,1.0);
    if (fMakeNuMINtuple)
      {
	_vecWeightsPPFX  = std::vector<unsigned short>(600,1);
	_vecWeightsPPFXD = std::vector<double>(600,1.0);
      }
    std::cout << " [ EventWeightTree ]" << " begin " << std::endl;

    std::vector<art::InputTag> vecTag;
    art::InputTag eventweight_tag_00("eventweight","","EventWeightSept24");
    //art::InputTag eventweight_tag_00("eventweight","","EventWeightMar18");
    art::InputTag eventweight_tag_01("eventweightSep24","","EventWeightSep24ExtraGENIE1");
    art::InputTag eventweight_tag_02("eventweightSep24","","EventWeightSep24ExtraGENIE2");
    art::InputTag eventweight_tag_03("eventweightSep24","","EventWeightSep24ExtraGENIE3");
    art::InputTag eventweight_tag_04("eventweightSep24","","EventWeightSep24ExtraGENIE4");
    art::InputTag eventweight_tag_knobs("eventweightGenieKnobs");
    //art::InputTag eventweight_spline_tag("eventweightSplines");
    // art::InputTag eventweight_tag("eventweight");
    vecTag.push_back(eventweight_tag_00);
    vecTag.push_back(eventweight_tag_01);
    vecTag.push_back(eventweight_tag_02);
    vecTag.push_back(eventweight_tag_03);
    vecTag.push_back(eventweight_tag_04);
    vecTag.push_back(eventweight_tag_knobs);
    //vecTag.push_back(eventweight_spline_tag);
    // vecTag.push_back(eventweight_tag);

    int ctr = 0;
    int GenieCounter = 0;
    int PPFXCounter = 0;

    for(auto& thisTag : vecTag){

      //std::cout << " [ EventWeightTree ]" << " newTag " << thisTag.label() << std::endl;
      
      art::Handle<std::vector<evwgh::MCEventWeight>> eventweights_handle;
      evt.getByLabel(thisTag, eventweights_handle);

      if(eventweights_handle.isValid()){

        std::cout << " [ EventWeightTree ]" << " isValid! " << std::endl;

	std::vector<art::Ptr<evwgh::MCEventWeight>> eventweights;
	art::fill_ptr_vector(eventweights, eventweights_handle);

	if (_createGenieUpDnVecs && thisTag.label()=="eventweightGenieKnobs") {

	  // we enforce that the index in vecWeightsGenieUp/Dn matches this vector
    std::vector<std::string> knobList = {"AGKYpT1pi_UBGenie","AGKYxF1pi_UBGenie","AhtBY_UBGenie","AxFFCCQEshape_UBGenie","BhtBY_UBGenie",
                 "CV1uBY_UBGenie","CV2uBY_UBGenie","DecayAngMEC_UBGenie","EtaNCEL_UBGenie","FrAbs_N_UBGenie",
                 "FrAbs_pi_UBGenie","FrCEx_N_UBGenie","FrCEx_pi_UBGenie","FrInel_N_UBGenie","FrInel_pi_UBGenie",
                 "FrPiProd_N_UBGenie","FrPiProd_pi_UBGenie","FracDelta_CCMEC_UBGenie","FracPN_CCMEC_UBGenie","MFP_N_UBGenie",
                 "MFP_pi_UBGenie","MaCCQE_UBGenie","MaCCRES_UBGenie","MaNCEL_UBGenie","MaNCRES_UBGenie",
                 "MvCCRES_UBGenie","MvNCRES_UBGenie","NonRESBGvbarnCC1pi_UBGenie","NonRESBGvbarnCC2pi_UBGenie","NonRESBGvbarnNC1pi_UBGenie",
                 "NonRESBGvbarnNC2pi_UBGenie","NonRESBGvbarpCC1pi_UBGenie","NonRESBGvbarpCC2pi_UBGenie","NonRESBGvbarpNC1pi_UBGenie","NonRESBGvbarpNC2pi_UBGenie",
                 "NonRESBGvnCC1pi_UBGenie","NonRESBGvnCC2pi_UBGenie","NonRESBGvnNC1pi_UBGenie","NonRESBGvnNC2pi_UBGenie","NonRESBGvpCC1pi_UBGenie",
                 "NonRESBGvpCC2pi_UBGenie","NonRESBGvpNC1pi_UBGenie","NonRESBGvpNC2pi_UBGenie","NormCCMEC_UBGenie","NormNCMEC_UBGenie",
                 "RDecBR1eta_UBGenie","RDecBR1gamma_UBGenie","RPA_CCQE_UBGenie","Theta_Delta2Npi_UBGenie","TunedCentralValue_UBGenie",
                 "VecFFCCQEshape_UBGenie","XSecShape_CCMEC_UBGenie","splines_general_Spline"};

	  std::map<std::string, std::vector<double>> evtwgt_map = eventweights.at(0)->fWeight;
	  // let's iterate over the list above, so that the order is guaranteed
	  for (size_t count=0; count<knobList.size(); count++) {
	    bool knobFound = false;
	    for(std::map<std::string, std::vector<double>>::iterator it=evtwgt_map.begin(); it!=evtwgt_map.end(); ++it){
	      if (it->first != knobList[count]) continue;
	      knobFound = true;
	      //std::cout << " [ EventWeightTree ]" << "Entering variation " << it->first << " with " << (it->second).size() << " weights " << std::endl;
	      std::vector<double>& vals = it->second;

	      if (vals.size()==1) {
		float w0 = vals[0];
		unsigned short w0short = (unsigned short)(w0*1000.);
		if (w0 > 65535)  { w0short = std::numeric_limits<unsigned short>::max(); }
		if (w0 < 0) { w0short = 1; }
		_vecWeightsGenieUp.push_back(w0short);
		_vecWeightsGenieDn.push_back((unsigned short)(0));
	      }
	      else if (vals.size()==2) {
		float w0 = vals[0];
		unsigned short w0short = (unsigned short)(w0*1000.);
		if (w0 > 65535)  { w0short = std::numeric_limits<unsigned short>::max(); }
		if (w0 < 0) { w0short = 1; }
		_vecWeightsGenieUp.push_back(w0short);
		float w1 = vals[1];
		unsigned short w1short = (unsigned short)(w1*1000.);
		if (w1 > 65535)  { w1short = std::numeric_limits<unsigned short>::max(); }
		if (w1 < 0) { w1short = 1; }
		_vecWeightsGenieDn.push_back(w1short);
	      }
	      else std::cout << "Argh, this is unexpected! Size should be either 0 or 1..." << std::endl;
	    }
	    // if the index is not in the order we expect, terminate the program
	    if (knobFound==false) {
	      std::cout << "knob " << knobList[count] << " NOT found!" << std::endl;
	      assert(0);
	    }
	  }
	}

	std::map<std::string, std::vector<double>> evtwgt_map = eventweights.at(0)->fWeight;
          
        if (evtwgt_map.find("RPA_CCQE_UBGenie") != evtwgt_map.end()) { 
          _knobRPAup = evtwgt_map.find("RPA_CCQE_UBGenie")->second[0]; 
          _knobRPAdn = evtwgt_map.find("RPA_CCQE_UBGenie")->second[1]; 
        }
        if (evtwgt_map.find("XSecShape_CCMEC_UBGenie") != evtwgt_map.end()) { 
          _knobCCMECup = evtwgt_map.find("XSecShape_CCMEC_UBGenie")->second[0]; 
          _knobCCMECdn = evtwgt_map.find("XSecShape_CCMEC_UBGenie")->second[1]; 
        }
        if (evtwgt_map.find("AxFFCCQEshape_UBGenie") != evtwgt_map.end()) { 
          _knobAxFFCCQEup = evtwgt_map.find("AxFFCCQEshape_UBGenie")->second[0]; 
          _knobAxFFCCQEdn = evtwgt_map.find("AxFFCCQEshape_UBGenie")->second[1]; 
        }
        if (evtwgt_map.find("VecFFCCQEshape_UBGenie") != evtwgt_map.end()) { 
          _knobVecFFCCQEup = evtwgt_map.find("VecFFCCQEshape_UBGenie")->second[0]; 
          _knobVecFFCCQEdn = evtwgt_map.find("VecFFCCQEshape_UBGenie")->second[1]; 
        }
        if (evtwgt_map.find("DecayAngMEC_UBGenie") != evtwgt_map.end()) { 
          _knobDecayAngMECup = evtwgt_map.find("DecayAngMEC_UBGenie")->second[0]; 
          _knobDecayAngMECdn = evtwgt_map.find("DecayAngMEC_UBGenie")->second[1]; 
        }
        if (evtwgt_map.find("Theta_Delta2Npi_UBGenie") != evtwgt_map.end()) { 
          _knobThetaDelta2Npiup = evtwgt_map.find("Theta_Delta2Npi_UBGenie")->second[0]; 
          _knobThetaDelta2Npidn = evtwgt_map.find("Theta_Delta2Npi_UBGenie")->second[1]; 
        }
        if (evtwgt_map.find("ThetaDelta2NRad_UBGenie") != evtwgt_map.end()) { 
          _knobThetaDelta2NRadup = evtwgt_map.find("ThetaDelta2NRad_UBGenie")->second[0]; 
          _knobThetaDelta2NRaddn = evtwgt_map.find("ThetaDelta2NRad_UBGenie")->second[1]; 
        }
        //if (evtwgt_map.find("RPA_CCQE_Reduced_UBGenie") != evtwgt_map.end()) { 
        //  _knobRPA_CCQE_Reducedup = evtwgt_map.find("RPA_CCQE_Reduced_UBGenie")->second[0]; 
        //  _knobRPA_CCQE_Reduceddn = evtwgt_map.find("RPA_CCQE_Reduced_UBGenie")->second[1]; 
        //}
        if (evtwgt_map.find("NormCCCOH_UBGenie") != evtwgt_map.end()) { 
          _knobNormCCCOHup = evtwgt_map.find("NormCCCOH_UBGenie")->second[0]; 
          _knobNormCCCOHdn = evtwgt_map.find("NormCCCOH_UBGenie")->second[1]; 
        }
        if (evtwgt_map.find("NormNCCOH_UBGenie") != evtwgt_map.end()) { 
          _knobNormNCCOHup = evtwgt_map.find("NormNCCOH_UBGenie")->second[0]; 
          _knobNormNCCOHdn = evtwgt_map.find("NormNCCOH_UBGenie")->second[1]; 
        }
        if (evtwgt_map.find("xsr_scc_Fv3_SCC") != evtwgt_map.end()) { 
          _knobxsr_scc_Fv3up = evtwgt_map.find("xsr_scc_Fv3_SCC")->second[0]; 
          _knobxsr_scc_Fv3dn = evtwgt_map.find("xsr_scc_Fv3_SCC")->second[1]; 
        }
        if (evtwgt_map.find("xsr_scc_Fa3_SCC") != evtwgt_map.end()) { 
          _knobxsr_scc_Fa3up = evtwgt_map.find("xsr_scc_Fa3_SCC")->second[0]; 
          _knobxsr_scc_Fa3dn = evtwgt_map.find("xsr_scc_Fa3_SCC")->second[1]; 
        }
        if (evtwgt_map.find("RootinoFix_UBGenie") != evtwgt_map.end()) { 
          _RootinoFix = evtwgt_map.find("RootinoFix_UBGenie")->second[0]; 
        }

        if(evtwgt_map.find("splines_general_Spline") != evtwgt_map.end()) _weightSpline = evtwgt_map.find("splines_general_Spline")->second[0];
        //evtwgt_map.erase("splines_general_Spline");
        if(evtwgt_map.find("TunedCentralValue_UBGenie") != evtwgt_map.end()) _weightTune = evtwgt_map.find("TunedCentralValue_UBGenie")->second[0];
        //evtwgt_map.erase("TunedCentralValue_Genie");

        //std::cout << " [ EventWeightTree ]" << " continue... " << std::endl;

        if(_weightSpline != -1 && _weightTune != -1) _weightSplineTimesTune = _weightSpline * _weightTune;

        // Get the PPFX Central Value -- this will be available only for NuMI
	if (fMakeNuMINtuple){
	  if(evtwgt_map.find("ppfx_cv_UBPPFXCV") != evtwgt_map.end()) _ppfx_cv = evtwgt_map.find("ppfx_cv_UBPPFXCV")->second[0];
	  std::cout << "ppfx cv weight: "<< _ppfx_cv<<  std::endl;
	  // evtwgt_map.erase("ppfx_cv_PPFXCV");
	}

        // old implementation -> just add all weights as they come from the input
        //_mapWeight.insert(evtwgt_map.begin(), evtwgt_map.end());


        bool isFirstVectorFlux   = true; // BNB Specific
        bool isFirstVectorReint  = true;

        // loop through all EventWeight variations
        for(std::map<std::string, std::vector<double>>::iterator it=evtwgt_map.begin(); it!=evtwgt_map.end(); ++it){

          ctr += 1;

          // variation name
          std::string keyname = it->first;

          // std::cout << " [ EventWeightTree ]" << "Entering variation " << keyname << " with " << (it->second).size() << " weights " << std::endl;

          // if this is a genie-all variation
          if (keyname.find("All") != std::string::npos) {
            std::cout << " [ EventWeightTree ]" << " Entering Genie ALL variation number " << GenieCounter << " " <<keyname << std::endl;
            for(unsigned int i = 0; i < it->second.size(); ++i) {
              if ( (i + (100 * GenieCounter) ) < _vecWeightsGenie.size()) 
                _vecWeightsGenieD[i + (100 * GenieCounter) ] *= it->second[i];
            }
            //else
            //std::cout << " [ EventWeightTree ]" << " ERROR FILLING GENIE WEIGHTS " << std::endl;
            GenieCounter += 1;
          }// if a genie-all variation
          
	  // We need to treat the flux weights differently between NuMI and BNB, so I'll try a nested condition
	  // we start with the NuMI option
          // if this is a ppfx multisim variation
          else if ((keyname.find("ppfx_ms_UBPPFX") != std::string::npos)&&fMakeNuMINtuple) {
            std::cout << " [ EventWeightTree ]" << " Entering PPFX variation number " << PPFXCounter << " " <<keyname << std::endl;
            for(unsigned int i = 0; i < it->second.size(); ++i) {
              if ( (i + (100 * PPFXCounter) ) < _vecWeightsPPFXD.size()) 
                _vecWeightsPPFXD[i + (100 * PPFXCounter) ] *= it->second[i];
            }
            //else
            //std::cout << " [ EventWeightTree ]" << " ERROR FILLING PPFX WEIGHTS " << std::endl;
            PPFXCounter += 1;
          }// if a ppfx-all variation
	  // then BNB
	  else if(!fMakeNuMINtuple &&
		  (keyname.find("horncurrent") != std::string::npos ||
		   keyname.find("expskin") != std::string::npos ||
		   keyname.find("piplus") != std::string::npos ||
		   keyname.find("piminus") != std::string::npos ||
		   keyname.find("kplus") != std::string::npos ||
		   keyname.find("kzero") != std::string::npos ||
		   keyname.find("kminus") != std::string::npos ||
		   keyname.find("pioninexsec") != std::string::npos ||
		   keyname.find("pionqexsec") != std::string::npos ||
		   keyname.find("piontotxsec") != std::string::npos ||
		   keyname.find("nucleontotxsec") != std::string::npos ||
		   keyname.find("nucleonqexsec") != std::string::npos ||
		   keyname.find("nucleoninexsec") != std::string::npos))
	    {
	      // are we storing flux variations one-by-one?
	      if (_SaveAllFlux) 
		_mapWeight.insert(*it);
	      else {
		// is this the first flux variation in the event?
		if(isFirstVectorFlux){
		  _vecWeightFluxD = it->second;
		  isFirstVectorFlux = false;
		}
		else{
		  // make sure same-size flux-weight vector
		  if ( (it->second).size() == _vecWeightFluxD.size() ) {
		    for(unsigned int i = 0; i < it->second.size(); ++i) 
		      _vecWeightFluxD[i] *= it->second[i];
		  }// if same-size flux-weight vector
		}// if not the first flux weight
	      }// if not storing flux variations one-by-one
	    }// if a flux-variation BNB
          // if a GEANT4 variation proton/piplus/piminus
          else if ( keyname == "reinteractions_piplus_Geant4" || keyname == "reinteractions_piminus_Geant4" || keyname == "reinteractions_proton_Geant4" ) {
            std::cout << "KrishReint: " << it->first << std::endl;
            // is this the first G4 variation in the event?
            if(isFirstVectorReint){
              _vecWeightsReintD = it->second;
              isFirstVectorReint = false;
            }
            else{
              // make sure same-size flux-weight vector
              if ( (it->second).size() == _vecWeightsReintD.size() ) {
                for(unsigned int i = 0; i < it->second.size(); ++i)
                  _vecWeightsReintD[i] *= it->second[i];
              }// if same-size G4-weight vector
            }// if not the first G4 weight
          }// if a G4 variation
          else{ // if not a genie-all variation
            _mapWeight.insert(*it);
          }// end if not a genie-all variation

        }// for all wegith vectors in map
      }//if the event-weight handle is valid

      //std::cout << " [ EventWeightTree ] " << "genie-all handle now has " << _vecWeightsGenie.size() << " elements" << std::endl;

    }// for all event-weight handles
    // if we are filling the flux-weights as one map:
    if (!_SaveAllFlux) {
      //std::cout << " [ EventWeightTree ] " << "flux-all weight vector now has " << _vecWeightFlux.size() << " elements" << std::endl;
      _mapWeight.insert( std::pair<std::string,std::vector<double> >("flux_all",_vecWeightFluxD) );
    }
    // store geinie-all map
    //std::cout << " [ EventWeightTree ] " << "flux-all weight vector now has " << _vecWeightsGenie.size() << " elements" << std::endl;
    _mapWeight.insert( std::pair<std::string,std::vector<double> >("All_UBGenie",_vecWeightsGenieD) );

    // store reinteraction map
    //std::cout << " [ EventWeightTree ] " << "reinteraction (G4) weight vector now has " << _vecWeightsReint.size() << " elements" << std::endl;
    _mapWeight.insert( std::pair<std::string,std::vector<double> >("reint_all",_vecWeightsReintD) );

    // store PPFX map
    //std::cout << " [ EventWeightTree ] " << "PPFX weight vector now has " << _vecWeightsPPFX.size() << " elements" << std::endl;
    _mapWeight.insert( std::pair<std::string,std::vector<double> >("ppfx_all",_vecWeightsPPFXD) );


    _vecWeightFlux = std::vector<unsigned short>(_vecWeightFluxD.size(), 1);
    for (size_t i=0; i < _vecWeightFluxD.size(); i++) {
      auto w = _vecWeightFluxD[i];
      unsigned short wshort = (unsigned short)(w*1000.);
      if (w > 65535)  { wshort = std::numeric_limits<unsigned short>::max(); }
      if (w < 0) { wshort = 1; }
      _vecWeightFlux[i] = wshort;
    }

    _vecWeightsGenie = std::vector<unsigned short>(_vecWeightsGenieD.size(), 1);
    for (size_t i=0; i < _vecWeightsGenieD.size(); i++) {
      auto w = _vecWeightsGenieD[i];
      unsigned short wshort = (unsigned short)(w*1000.);
      if (w > 65535)  { wshort = std::numeric_limits<unsigned short>::max(); }
      if (w < 0) { wshort = 1; }
      _vecWeightsGenie[i] = wshort;
    }

    _vecWeightsReint = std::vector<unsigned short>(_vecWeightsReintD.size(), 1);
    for (size_t i=0; i < _vecWeightsReintD.size(); i++) {
      auto w = _vecWeightsReintD[i];
      unsigned short wshort = (unsigned short)(w*1000.);
      if (w > 65535)  { wshort = std::numeric_limits<unsigned short>::max(); }
      if (w < 0) { wshort = 1; }
      _vecWeightsReint[i] = wshort;
    }

    if (fMakeNuMINtuple)
      {	_vecWeightsPPFX = std::vector<unsigned short>(_vecWeightsPPFXD.size(), 1);
	for (size_t i=0; i < _vecWeightsPPFXD.size(); i++) {
	  auto w = _vecWeightsPPFXD[i];
	  unsigned short wshort = (unsigned short)(w*1000.);
	  if (w > 65535)  { wshort = std::numeric_limits<unsigned short>::max(); }
	  if (w < 0) { wshort = 1; }
	  _vecWeightsPPFX[i] = wshort;
	}
      }

    // store genie-all variation

    if(_createDedicatedTree) _weightstree->Fill();
  }

  void EventWeightTree::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected){
  }

  void EventWeightTree::setBranches(TTree *_tree){
    if(_createDedicatedTree){
      _weightstree->Branch("weights", "std::map<std::string, std::vector<double>>", &_mapWeight);
      _weightstree->Branch("run",&_run,"run/I");
      _weightstree->Branch("subRun",&_subRun,"subRun/I");
      _weightstree->Branch("evt",&_evt,"evt/I");
    }
    if(_createMapBranch) _tree->Branch("weights", "std::map<std::string, std::vector<double>>", &_mapWeight);
    if(_createFluxBranch) _tree->Branch("weightsFlux", "std::vector<unsigned short>", &_vecWeightFlux);
    if(_createGenieBranch) _tree->Branch("weightsGenie", "std::vector<unsigned short>", &_vecWeightsGenie);
    //if(_createGenieBranch) _tree->Branch("weightsGenie_vec", "std::vector<double>", &_vecWeightsGenie_vec);
    //if(_createGenieBranch) _tree->Branch("weightsGenie_nam", "std::vector<int>", &_vecWeightsGenie_nam);
    if(_createReintBranch) _tree->Branch("weightsReint", "std::vector<unsigned short>", &_vecWeightsReint);
    if(_createSplineBranch) _tree->Branch("weightSpline",&_weightSpline,"weightSpline/F");
    if(_createTuneBranch) _tree->Branch("weightTune",&_weightTune,"weightTune/F");
    if(_createSplineTimesTuneBranch) _tree->Branch("weightSplineTimesTune",&_weightSplineTimesTune,"weightSplineTimesTune/F");
    if(_createPPFXBranch) _tree->Branch("ppfx_cv",&_ppfx_cv,"ppfx_cv/F");
    if(_createPPFXBranch) _tree->Branch("weightsPPFX", "std::vector<unsigned short>", &_vecWeightsPPFX);

    if (_createGenieBranch) {
      _tree->Branch("knobRPAup",&_knobRPAup,"knobRPAup/D");
      _tree->Branch("knobRPAdn",&_knobRPAdn,"knobRPAdn/D");
      _tree->Branch("knobCCMECup",&_knobCCMECup,"knobCCMECup/D");
      _tree->Branch("knobCCMECdn",&_knobCCMECdn,"knobCCMECdn/D");
      _tree->Branch("knobAxFFCCQEup",&_knobAxFFCCQEup,"knobAxFFCCQEup/D");
      _tree->Branch("knobAxFFCCQEdn",&_knobAxFFCCQEdn,"knobAxFFCCQEdn/D");
      _tree->Branch("knobVecFFCCQEup",&_knobVecFFCCQEup,"knobVecFFCCQEup/D");
      _tree->Branch("knobVecFFCCQEdn",&_knobVecFFCCQEdn,"knobVecFFCCQEdn/D");
      _tree->Branch("knobDecayAngMECup",&_knobDecayAngMECup,"knobDecayAngMECup/D");
      _tree->Branch("knobDecayAngMECdn",&_knobDecayAngMECdn,"knobDecayAngMECdn/D");
      _tree->Branch("knobThetaDelta2Npiup",&_knobThetaDelta2Npiup,"knobThetaDelta2Npiup/D");
      _tree->Branch("knobThetaDelta2Npidn",&_knobThetaDelta2Npidn,"knobThetaDelta2Npidn/D");
      _tree->Branch("knobThetaDelta2NRadup",&_knobThetaDelta2NRadup,"knobThetaDelta2NRadup/D");
      _tree->Branch("knobThetaDelta2NRaddn",&_knobThetaDelta2NRaddn,"knobThetaDelta2NRaddn/D");
      //_tree->Branch("knobRPA_CCQE_Reducedup",&_knobRPA_CCQE_Reducedup,"knobRPA_CCQE_Reducedup/D");
      //_tree->Branch("knobRPA_CCQE_Reduceddn",&_knobRPA_CCQE_Reduceddn,"knobRPA_CCQE_Reduceddn/D");
      _tree->Branch("knobNormCCCOHup",&_knobNormCCCOHup,"knobNormCCCOHup/D");
      _tree->Branch("knobNormCCCOHdn",&_knobNormCCCOHdn,"knobNormCCCOHdn/D");
      _tree->Branch("knobNormNCCOHup",&_knobNormNCCOHup,"knobNormNCCOHup/D");
      _tree->Branch("knobNormNCCOHdn",&_knobNormNCCOHdn,"knobNormNCCOHdn/D");
      _tree->Branch("knobxsr_scc_Fv3up",&_knobxsr_scc_Fv3up,"knobxsr_scc_Fv3up/D");
      _tree->Branch("knobxsr_scc_Fv3dn",&_knobxsr_scc_Fv3dn,"knobxsr_scc_Fv3dn/D");
      _tree->Branch("knobxsr_scc_Fa3up",&_knobxsr_scc_Fa3up,"knobxsr_scc_Fa3up/D");
      _tree->Branch("knobxsr_scc_Fa3dn",&_knobxsr_scc_Fa3dn,"knobxsr_scc_Fa3dn/D");
      _tree->Branch("RootinoFix",&_RootinoFix,"RootinoFix/D");
    }
    if(_createGenieUpDnVecs) _tree->Branch("weightsGenieUp", "std::vector<unsigned short>", &_vecWeightsGenieUp);
    if(_createGenieUpDnVecs) _tree->Branch("weightsGenieDn", "std::vector<unsigned short>", &_vecWeightsGenieDn);

  }

  void EventWeightTree::resetTTree(TTree *_tree){
    _mapWeight.clear();
    _vecWeightFlux.clear();
    _vecWeightFluxD.clear();
    _vecWeightsGenie.clear();
    _vecWeightsGenieD.clear();
    _vecWeightsGenie_vec.clear();
    _vecWeightsGenie_nam.clear();
    _vecWeightsReint.clear();
    _vecWeightsReintD.clear();
    _vecWeightsPPFX.clear();
    _vecWeightsPPFXD.clear();
    _weightSpline = -1;
    _weightTune = -1;
    _weightSplineTimesTune = -1;
    _ppfx_cv = -1;

    _knobRPAup = 1;
    _knobCCMECup = 1;
    _knobAxFFCCQEup = 1;
    _knobVecFFCCQEup = 1;
    _knobDecayAngMECup = 1;
    _knobThetaDelta2Npiup = 1;
    _knobThetaDelta2NRadup = 1;
    //_knobRPA_CCQE_Reducedup = 1;
    _knobNormCCCOHup = 1;
    _knobNormNCCOHup = 1;
    _knobxsr_scc_Fv3up = 1;
    _knobxsr_scc_Fa3up = 1;
    _knobRPAdn = 1;
    _knobCCMECdn = 1;
    _knobAxFFCCQEdn = 1;
    _knobVecFFCCQEdn = 1;
    _knobDecayAngMECdn = 1;
    _knobThetaDelta2Npidn = 1;
    _knobThetaDelta2NRaddn = 1;
    //_knobRPA_CCQE_Reduceddn = 1;
    _knobNormCCCOHdn = 1;
    _knobNormNCCOHdn = 1;
    _knobxsr_scc_Fv3dn = 1;
    _knobxsr_scc_Fa3dn = 1;
    _RootinoFix = 1;

    _vecWeightsGenieUp.clear();
    _vecWeightsGenieDn.clear();

    _run = -1;
    _subRun = -1;
    _evt = -1;
  }

  DEFINE_ART_CLASS_TOOL(EventWeightTree)

} // namespace analysis

#endif
