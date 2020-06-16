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
    std::vector<double> _vecWeightFlux;
    std::vector<double> _vecWeightsGenie;
    std::vector<double> _vecWeightsGenie_vec;
    std::vector<int> _vecWeightsGenie_nam;
    std::vector<double> _vecWeightsReint;
    double _knobRPAup;
    double _knobCCMECup;
    double _knobAxFFCCQEup;
    double _knobVecFFCCQEup;
    double _knobDecayAngMECup;
    double _knobThetaDelta2Npiup;
    double _knobRPAdn;
    double _knobCCMECdn;
    double _knobAxFFCCQEdn;
    double _knobVecFFCCQEdn;
    double _knobDecayAngMECdn;
    double _knobThetaDelta2Npidn;
    float _weightSpline;
    float _weightTune;
    float _weightSplineTimesTune;
    bool _createDedicatedTree;
    bool _createMapBranch;
    bool _createFluxBranch;
    bool _createGenieBranch;
    bool _createReintBranch;
    bool _createSplineBranch;
    bool _createTuneBranch;
    bool _createSplineTimesTuneBranch;
    bool _SaveAllFlux;
    int _GenieAllUniverses;
    int _run;
    int _subRun;
    int _evt;
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
    _SaveAllFlux = p.get<bool>("SaveAllFlux",false);
    _GenieAllUniverses = p.get<int>("GenieAllUniverses",500);
    
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

    _vecWeightsGenie = std::vector<double>(_GenieAllUniverses,1.0);

    //std::cout << " [ EventWeightTree ]" << " begin " << std::endl;
    
    std::vector<art::InputTag> vecTag;
    art::InputTag eventweight_tag_00("eventweight","","EventWeightMar18");
    art::InputTag eventweight_tag_01("eventweight","","EventWeightMar18ExtraGENIE1");
    art::InputTag eventweight_tag_02("eventweight","","EventWeightMar18ExtraGENIE2");
    art::InputTag eventweight_tag_03("eventweight","","EventWeightMar18ExtraGENIE3");
    art::InputTag eventweight_tag_04("eventweight","","EventWeightMar18ExtraGENIE4");
    //art::InputTag eventweight_spline_tag("eventweightSplines");
    vecTag.push_back(eventweight_tag_00);
    vecTag.push_back(eventweight_tag_01);
    vecTag.push_back(eventweight_tag_02);
    vecTag.push_back(eventweight_tag_03);
    vecTag.push_back(eventweight_tag_04);
    //vecTag.push_back(eventweight_spline_tag);

    int ctr = 0;
    
    for(auto& thisTag : vecTag){

      //std::cout << " [ EventWeightTree ]" << " newTag " << std::endl;
      
      art::Handle<std::vector<evwgh::MCEventWeight>> eventweights_handle;
      evt.getByLabel(thisTag, eventweights_handle);
      
      if(eventweights_handle.isValid()){

	//std::cout << " [ EventWeightTree ]" << " isValid! " << std::endl;

	std::vector<art::Ptr<evwgh::MCEventWeight>> eventweights;
	art::fill_ptr_vector(eventweights, eventweights_handle);
	std::map<std::string, std::vector<double>> evtwgt_map = eventweights.at(0)->fWeight;

	if (evtwgt_map.find("RPA_CCQE_UBGenie") != evtwgt_map.end()) { 
	  _knobRPAup = evtwgt_map.find("RPA_CCQE_UBGenie")->second[0]; 
	  _knobRPAup = evtwgt_map.find("RPA_CCQE_UBGenie")->second[1]; 
	}
	if (evtwgt_map.find("XSecShape_CCMEC_UBGenie") != evtwgt_map.end()) { 
	  _knobCCMECup = evtwgt_map.find("XSecShape_CCMEC_UBGenie")->second[0]; 
	  _knobCCMECup = evtwgt_map.find("XSecShape_CCMEC_UBGenie")->second[1]; 
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
        
	if(evtwgt_map.find("splines_general_Spline") != evtwgt_map.end()) _weightSpline = evtwgt_map.find("splines_general_Spline")->second[0];
	//evtwgt_map.erase("splines_general_Spline");
	if(evtwgt_map.find("TunedCentralValue_UBGenie") != evtwgt_map.end()) _weightTune = evtwgt_map.find("TunedCentralValue_UBGenie")->second[0];
	//evtwgt_map.erase("TunedCentralValue_Genie");

	//std::cout << " [ EventWeightTree ]" << " continue... " << std::endl;
        
	if(_weightSpline != -1 && _weightTune != -1) _weightSplineTimesTune = _weightSpline * _weightTune;
	
	// old implementation -> just add all weights as they come from the input
	//_mapWeight.insert(evtwgt_map.begin(), evtwgt_map.end());
	
	int GenieCounter = 0;
	bool isFirstVectorFlux   = true;
	bool isFirstVectorReint  = true;
	
	// loop through all EventWeight variations
	for(std::map<std::string, std::vector<double>>::iterator it=evtwgt_map.begin(); it!=evtwgt_map.end(); ++it){
	  
	  ctr += 1;
	  
	  // variation name
	  std::string keyname = it->first;
	  
	  //std::cout << " [ EventWeightTree ]" << "Entering variation " << keyname << " with " << (it->second).size() << " weights " << std::endl;
	  
	  // is this a flux variation?
	  if(keyname.find("horncurrent") != std::string::npos ||
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
	     keyname.find("nucleoninexsec") != std::string::npos){
	    
	    // are we storing flux variations one-by-one?
	    if (_SaveAllFlux) 
	      _mapWeight.insert(*it);
	    else {
	      // is this the first flux variation in the event?
	      if(isFirstVectorFlux){
		_vecWeightFlux = it->second;
		isFirstVectorFlux = false;
	      }
	      else{
		// make sure same-size flux-weight vector
		if ( (it->second).size() == _vecWeightFlux.size() ) {
		  for(unsigned int i = 0; i < it->second.size(); ++i)
		    _vecWeightFlux[i] *= it->second[i];
		}// if same-size flux-weight vector
	      }// if not the first flux weight
	    }// if not storing flux variations one-by-one
	  }// if a flux-variation
	  // if this is a genie-all variation
	  else if (keyname.find("All") != std::string::npos) {
	    for(unsigned int i = 0; i < it->second.size(); ++i)
	      if ( (i + (100 * GenieCounter) ) < _vecWeightsGenie.size())
		_vecWeightsGenie[i + (100 * GenieCounter) ] *= it->second[i];
	      else
		//std::cout << " [ EventWeightTree ]" << " ERROR FILLING GENIE WEIGHTS " << std::endl;
	    GenieCounter += 1;
	  }// if a genie-all variation
	  // if a GEANT4 variation
	  else if ( (keyname.find("reinteractions") != std::string::npos) ) {
	      // is this the first G4 variation in the event?
	      if(isFirstVectorReint){
		_vecWeightsReint = it->second;
		isFirstVectorReint = false;
	      }
	      else{
		// make sure same-size flux-weight vector
		if ( (it->second).size() == _vecWeightsReint.size() ) {
		  for(unsigned int i = 0; i < it->second.size(); ++i)
		    _vecWeightsReint[i] *= it->second[i];
		}// if same-size G4-weight vector
	      }// if not the first G4 weight
	  }// if a G4 variation
	  else{ // if not a genie-all variation
	    _mapWeight.insert(*it);
	  }// end if not a genie-all variation
	}// for all weith vectors in map
      }//if the event-weight handle is valid
      
      //std::cout << " [ EventWeightTree ] " << "genie-all handle now has " << _vecWeightsGenie.size() << " elements" << std::endl;
      
    }// for all event-weight handles
    
    // if we are filling the flux-weights as one map:
    if (!_SaveAllFlux) {
      //std::cout << " [ EventWeightTree ] " << "flux-all weight vector now has " << _vecWeightFlux.size() << " elements" << std::endl;
      _mapWeight.insert( std::pair<std::string,std::vector<double> >("flux_all",_vecWeightFlux) );
    }
    // store geinie-all map
    //std::cout << " [ EventWeightTree ] " << "flux-all weight vector now has " << _vecWeightsGenie.size() << " elements" << std::endl;
    _mapWeight.insert( std::pair<std::string,std::vector<double> >("All_UBGenie",_vecWeightsGenie) );

    // store reinteraction map
    //std::cout << " [ EventWeightTree ] " << "reinteraction (G4) weight vector now has " << _vecWeightsReint.size() << " elements" << std::endl;
    _mapWeight.insert( std::pair<std::string,std::vector<double> >("reint_all",_vecWeightsReint) );

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
    if(_createFluxBranch) _tree->Branch("weightsFlux", "std::vector<double>", &_vecWeightFlux);
    if(_createGenieBranch) _tree->Branch("weightsGenie", "std::vector<double>", &_vecWeightsGenie);
    //if(_createGenieBranch) _tree->Branch("weightsGenie_vec", "std::vector<double>", &_vecWeightsGenie_vec);
    //if(_createGenieBranch) _tree->Branch("weightsGenie_nam", "std::vector<int>", &_vecWeightsGenie_nam);
    if(_createReintBranch) _tree->Branch("weightsReint", "std::vector<double>", &_vecWeightsReint);
    if(_createSplineBranch) _tree->Branch("weightSpline",&_weightSpline,"weightSpline/F");
    if(_createTuneBranch) _tree->Branch("weightTune",&_weightTune,"weightTune/F");
    if(_createSplineTimesTuneBranch) _tree->Branch("weightSplineTimesTune",&_weightSplineTimesTune,"weightSplineTimesTune/F");

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
    }

  }
  
  void EventWeightTree::resetTTree(TTree *_tree){
    _mapWeight.clear();
    _vecWeightFlux.clear();
    _vecWeightsGenie.clear();
    _vecWeightsGenie_vec.clear();
    _vecWeightsGenie_nam.clear();
    _vecWeightsReint.clear();
    _weightSpline = -1;
    _weightTune = -1;
    _weightSplineTimesTune = -1;

    _knobRPAup = 1;
    _knobCCMECup = 1;
    _knobAxFFCCQEup = 1;
    _knobVecFFCCQEup = 1;
    _knobDecayAngMECup = 1;
    _knobThetaDelta2Npiup = 1;
    _knobRPAdn = 1;
    _knobCCMECdn = 1;
    _knobAxFFCCQEdn = 1;
    _knobVecFFCCQEdn = 1;
    _knobDecayAngMECdn = 1;
    _knobThetaDelta2Npidn = 1;


    _run = -1;
    _subRun = -1;
    _evt = -1;
  }
  
  DEFINE_ART_CLASS_TOOL(EventWeightTree)
  
} // namespace analysis

#endif
