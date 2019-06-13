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
    //              Provide uncertainty event weights in a dedicated tree
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
            std::map<std::string, std::vector<double>> _weights;
    };
    
    EventWeightTree::EventWeightTree(const fhicl::ParameterSet &p){
        art::ServiceHandle<art::TFileService> tfs;
        _weightstree = tfs->make<TTree>("EventWeights", "EventWeights TTree");
    }

    void EventWeightTree::configure(fhicl::ParameterSet const & p){
    }
    
    void EventWeightTree::analyzeEvent(art::Event const& evt, bool fData){
        art::InputTag eventweight_tag("eventweight");
        art::Handle<std::vector<evwgh::MCEventWeight>> eventweights_handle;
        evt.getByLabel(eventweight_tag, eventweights_handle);
        
        std::vector<art::Ptr<evwgh::MCEventWeight>> eventweights;
        art::fill_ptr_vector(eventweights, eventweights_handle);
        std::map<std::string, std::vector<double>> evtwgt_map = eventweights.at(0)->fWeight;
        _weights.insert(evtwgt_map.begin(), evtwgt_map.end());
        _weightstree->Fill();
    }
    
    void EventWeightTree::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected){
    }

    void EventWeightTree::setBranches(TTree *_tree){
       _weightstree->Branch("weights", "std::map<std::string, std::vector<double>>", &_weights);
    }
    
    void EventWeightTree::resetTTree(TTree *_tree){
       _weights.clear();
    }
    
    DEFINE_ART_CLASS_TOOL(EventWeightTree)

} // namespace analysis

#endif
