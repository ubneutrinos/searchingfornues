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
            std::vector<double> _vecWeightsReint;
            float _weightSpline;
            bool _createDedicatedTree;
            bool _createMapBranch;
            bool _createFluxBranch;
            bool _createGenieBranch;
            bool _createReintBranch;
            bool _createSplineBranch;
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
        
        std::vector<art::InputTag> vecTag;
        art::InputTag eventweight_tag("eventweight");
        art::InputTag eventweight_spline_tag("eventweightSplines");
        vecTag.push_back(eventweight_tag);
        vecTag.push_back(eventweight_spline_tag);
        
        for(auto& thisTag : vecTag){
            art::Handle<std::vector<evwgh::MCEventWeight>> eventweights_handle;
            evt.getByLabel(thisTag, eventweights_handle);
            if(eventweights_handle.isValid()){
                std::vector<art::Ptr<evwgh::MCEventWeight>> eventweights;
                art::fill_ptr_vector(eventweights, eventweights_handle);
                std::map<std::string, std::vector<double>> evtwgt_map = eventweights.at(0)->fWeight;
                
                if(evtwgt_map.find("splines_general_Spline") != evtwgt_map.end()) _weightSpline = evtwgt_map.find("splines_general_Spline")->second[0];
                evtwgt_map.erase("splines_general_Spline");
                if(evtwgt_map.find("genie_all_Genie") != evtwgt_map.end()) _vecWeightsGenie = evtwgt_map.find("genie_all_Genie")->second;
                evtwgt_map.erase("genie_all_Genie");
                if(evtwgt_map.find("reinteractions_all_Reinteraction") != evtwgt_map.end()) _vecWeightsReint = evtwgt_map.find("reinteractions_all_Reinteraction")->second;
                evtwgt_map.erase("reinteractions_all_Reinteraction");

                _mapWeight.insert(evtwgt_map.begin(), evtwgt_map.end());
                
                if(_createFluxBranch || _createGenieBranch || _createReintBranch){
                    bool isFirstVector = true;

                    for(std::map<std::string, std::vector<double>>::iterator it=evtwgt_map.begin(); it!=evtwgt_map.end(); ++it){
                        std::string keyname = it->first;
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
                            if(isFirstVector){
                                _vecWeightFlux = it->second;
                                isFirstVector = false;
                            }
                            else{
                                for(unsigned int i = 0; i < it->second.size(); ++i){
                                    _vecWeightFlux[i] *= it->second[i];
                                }
                            }
                        }
                    }
                }
            }
        }

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
        if(_createReintBranch) _tree->Branch("weightsReint", "std::vector<double>", &_vecWeightsReint);
        if(_createSplineBranch) _tree->Branch("weightSpline",&_weightSpline,"weightSpline/F");
    }
    
    void EventWeightTree::resetTTree(TTree *_tree){
        _mapWeight.clear();
        _vecWeightFlux.clear();
        _vecWeightsGenie.clear();
        _vecWeightsReint.clear();
        _weightSpline = -1;
        _run = -1;
        _subRun = -1;
        _evt = -1;
    }
    
    DEFINE_ART_CLASS_TOOL(EventWeightTree)

} // namespace analysis

#endif
