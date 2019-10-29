#ifndef ANALYSIS_BDT_CXX
#define ANALYSIS_BDT_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "art/Framework/Services/System/TriggerNamesService.h"
#include "canvas/Persistency/Common/TriggerResults.h"

#include "ubana/XGBoost/xgboost/c_api.h"

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       BDT
    // File:        BDT.cc
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

  class BDT : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    BDT(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~BDT();
    
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

    BoosterHandle booster_nuNCpi0;
    BoosterHandle booster_numuCCpi0;
    BoosterHandle booster_numuCC;
    BoosterHandle booster_ext;
    BoosterHandle booster_cosmic;
    BoosterHandle booster_global;
    float _bdt_nuNCpi0;
    float _bdt_numuCCpi0;
    float _bdt_numuCC;
    float _bdt_ext;
    float _bdt_cosmic;
    float _bdt_global;
    int _pass_antibdt_filter;
    TTree* _mytree;
    bool fVerbose;
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  BDT::BDT(const fhicl::ParameterSet& p)
  {
    fVerbose = p.get<bool>("Verbose", false);
    //documentation in xgboost/include/xgboost/c_api.h
    int xgtest = -1;
    xgtest = XGBoosterCreate(NULL, 0, &booster_nuNCpi0);
    xgtest = XGBoosterCreate(NULL, 0, &booster_numuCCpi0);
    xgtest = XGBoosterCreate(NULL, 0, &booster_numuCC);
    xgtest = XGBoosterCreate(NULL, 0, &booster_ext);
    xgtest = XGBoosterCreate(NULL, 0, &booster_cosmic);
    xgtest = XGBoosterCreate(NULL, 0, &booster_global);

    std::string _filename;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file("searchingfornues/booster_nopid_ncpi0.model",_filename);
    xgtest = XGBoosterLoadModel(booster_nuNCpi0, _filename.c_str());
    sp.find_file("searchingfornues/booster_nopid_ccpi0.model",_filename);
    xgtest = XGBoosterLoadModel(booster_numuCCpi0, _filename.c_str());
    sp.find_file("searchingfornues/booster_nopid_cc.model",_filename);
    xgtest = XGBoosterLoadModel(booster_numuCC, _filename.c_str());
    sp.find_file("searchingfornues/booster_nopid_ext.model",_filename);
    xgtest = XGBoosterLoadModel(booster_ext, _filename.c_str());
    sp.find_file("searchingfornues/booster_nopid_cosmic.model",_filename);
    xgtest = XGBoosterLoadModel(booster_cosmic, _filename.c_str());
    sp.find_file("searchingfornues/booster_nopid.model",_filename);
    xgtest = XGBoosterLoadModel(booster_global, _filename.c_str());
    assert(xgtest==0);
  }

  BDT::~BDT()
  {
    int xgtest = -1;
    xgtest = XGBoosterFree(booster_nuNCpi0);
    xgtest = XGBoosterFree(booster_numuCCpi0);
    xgtest = XGBoosterFree(booster_numuCC);
    xgtest = XGBoosterFree(booster_ext);
    xgtest = XGBoosterFree(booster_cosmic);
    xgtest = XGBoosterFree(booster_global);
    assert(xgtest==0);
  }

  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void BDT::configure(fhicl::ParameterSet const & p)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void BDT::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
  {

    art::Handle<art::TriggerResults> filter;
    e.getByLabel(art::InputTag("TriggerResults","","OverlayFiltersPostStage2"),filter);
    if (filter.isValid() && filter->size()==4) {
      _pass_antibdt_filter = filter->at(1).accept();//fixme accessing position at(1) is not robust
    }
    // std::cout << "filter size=" << filter->size() << std::endl;
    // std::cout << "filter pset=" << filter->parameterSetID().to_string() << std::endl;
    // auto tns = art::ServiceHandle<art::TriggerNamesService>();
    // size_t ntp =  tns->size();
    // std::cout << "ntp=" << ntp << std::endl;
    // // size_t ftp = ntp;
    // for (size_t itp=0;itp<filter->size();itp++) {
    //   std::cout << itp << " " << filter->at(itp).accept()  << std::endl;
    //   // std::cout << art::ServiceHandle<art::TriggerNamesService>()->getTrigPath(itp) << " " << filter->at(itp).accept()  << std::endl;
    //   // if (art::ServiceHandle<art::TriggerNamesService>()->getTrigPath(itp)=="sel2") ftp = itp;
    // }

    std::vector<float> data;
    std::vector<std::string> variables{"shr_dedx_Y", "shr_distance", "trk_distance", "pt", "hits_y",
	"shr_tkfit_dedx_Y", "shr_tkfit_dedx_U", "shr_tkfit_dedx_V", "p",
	"hits_ratio", "shr_dedx_U", "shr_dedx_V", "n_tracks_contained", "n_showers_contained",
	"shr_theta", "trk_len", "trk_score", "shr_score", "shr_energy_tot_cali", "trk_energy_tot",
	"shr_phi", "trk_theta", "trk_phi", "tksh_angle", "tksh_distance", "CosmicIP",
	"shr_pca_2", "shr_pca_1", "shr_pca_0",
	"topological_score", "slpdg"};
    std::vector<std::string> type{"F", "F", "F", "F", "i",
	"F", "F", "F", "F",
	"F", "F", "F", "I", "I",
	"F", "F", "F", "F", "F", "F",
	"F", "F", "F", "F", "F", "F",
	"F", "F", "F",
	"F", "I"};
    //retrieve variables from the tree abd fill data vector accordingly
    for (size_t iv=0;iv<variables.size();++iv) {
      if (type[iv]=="F") {
	float* tmp = (float*) _mytree->GetBranch(variables[iv].c_str())->GetAddress();
	data.push_back(*tmp);
	if (fVerbose) std::cout << type[iv] << " " << variables[iv] << "=" << data.back() << " - " << *tmp << std::endl;
      } else  if (type[iv]=="I") {
	int* tmp = (int*) _mytree->GetBranch(variables[iv].c_str())->GetAddress();
	data.push_back(*tmp);
	if (fVerbose) std::cout << type[iv] << " " << variables[iv] << "=" << data.back() << " - " << *tmp << std::endl;
      } else  if (type[iv]=="i") {
	unsigned int* tmp = (unsigned int*) _mytree->GetBranch(variables[iv].c_str())->GetAddress();
	data.push_back(*tmp);
	if (fVerbose) std::cout << type[iv] << " " << variables[iv] << "=" << data.back() << " - " << *tmp << std::endl;
      } else {exit(1);}
    }
    //get predictions: 5 BDTs against specific backgrounds
    int xgtest = -1;
    DMatrixHandle mat_nuNCpi0;
    xgtest = XGDMatrixCreateFromMat(data.data(),1,data.size(),0,&mat_nuNCpi0);
    bst_ulong out_len_nuNCpi0 = 0;
    const float* out_result_nuNCpi0   = NULL;
    xgtest = XGBoosterPredict(booster_nuNCpi0  , mat_nuNCpi0, 0, 317/*319*/, &out_len_nuNCpi0, &out_result_nuNCpi0  );// NB: setting by hand ntree_limit=booster.best_iteration
    _bdt_nuNCpi0   = *out_result_nuNCpi0  ;
    if (fVerbose) std::cout << "bdt_nuNCpi0=" << *out_result_nuNCpi0 << " " << _bdt_nuNCpi0 << std::endl;
    //
    DMatrixHandle mat_numuCCpi0;
    xgtest = XGDMatrixCreateFromMat(data.data(),1,data.size(),0,&mat_numuCCpi0);
    bst_ulong out_len_numuCCpi0 = 0;
    const float* out_result_numuCCpi0 = NULL;
    xgtest = XGBoosterPredict(booster_numuCCpi0, mat_numuCCpi0, 0,  343/*705*/, &out_len_numuCCpi0, &out_result_numuCCpi0);
    _bdt_numuCCpi0 = *out_result_numuCCpi0;
    if (fVerbose) std::cout << "bdt_numuCCpi0=" << *out_result_numuCCpi0 << " " << _bdt_numuCCpi0 << std::endl;
    //
    DMatrixHandle mat_numuCC;
    xgtest = XGDMatrixCreateFromMat(data.data(),1,data.size(),0,&mat_numuCC);
    bst_ulong out_len_numuCC = 0;
    const float* out_result_numuCC    = NULL;
    xgtest = XGBoosterPredict(booster_numuCC   , mat_numuCC, 0, 491/*514*/, &out_len_numuCC, &out_result_numuCC   );
    _bdt_numuCC    = *out_result_numuCC   ;
    if (fVerbose) std::cout << "bdt_numuCC=" << *out_result_numuCC << " " << _bdt_numuCC << std::endl;
    //
    DMatrixHandle mat_ext;
    xgtest = XGDMatrixCreateFromMat(data.data(),1,data.size(),0,&mat_ext);
    bst_ulong out_len_ext = 0;
    const float* out_result_ext       = NULL;
    xgtest = XGBoosterPredict(booster_ext      , mat_ext, 0, 265/*390*/, &out_len_ext, &out_result_ext      );
    _bdt_ext       = *out_result_ext      ;
    if (fVerbose) std::cout << "bdt_ext=" << *out_result_ext << " " << _bdt_ext << std::endl;
    //
    DMatrixHandle mat_cosmic;
    xgtest = XGDMatrixCreateFromMat(data.data(),1,data.size(),0,&mat_cosmic);
    bst_ulong out_len_cosmic = 0;
    const float* out_result_cosmic    = NULL;
    xgtest = XGBoosterPredict(booster_cosmic   , mat_cosmic, 0, 304/*782*/, &out_len_cosmic, &out_result_cosmic   );
    _bdt_cosmic    = *out_result_cosmic   ;
    if (fVerbose) std::cout << "bdt_cosmic=" << *out_result_cosmic << " " << _bdt_cosmic << std::endl;
    //
    // const char** dump_result_cosmic       = NULL;
    // xgtest = XGBoosterDumpModel(booster_cosmic,"",0, &out_len_cosmic, &dump_result_cosmic);
    // for (bst_ulong k=0;k<out_len_cosmic;k++){
    //   printf("%s \n",dump_result_cosmic[k]);
    // }
    //

    //get global BDT prediction
    std::vector<float> data_global;
    data_global.push_back(_bdt_ext      );
    data_global.push_back(_bdt_nuNCpi0  );
    data_global.push_back(_bdt_numuCC   );
    data_global.push_back(_bdt_numuCCpi0);
    data_global.push_back(_bdt_cosmic   );
    DMatrixHandle mat_global;
    xgtest = XGDMatrixCreateFromMat(data_global.data(),1,data_global.size(),0,&mat_global);
    bst_ulong out_len_global = 0;
    const float* out_result_global    = NULL;
    xgtest = XGBoosterPredict(booster_global, mat_global, 0, 115/*60*/, &out_len_global, &out_result_global   );
    _bdt_global    = *out_result_global   ;
    if (fVerbose) std::cout << "bdt_global=" << *out_result_global << " " << _bdt_global << std::endl;

    assert(xgtest==0);

    return;
  }

  void BDT::analyzeEvent(art::Event const &e, bool fData)
  {
    // std::cout << "analyze event" << std::endl;
  }

  void BDT::setBranches(TTree* _tree)
  {
    _tree->Branch("bdt_nuNCpi0"  ,&_bdt_nuNCpi0  ,"bdt_nuNCpi0/F"  );
    _tree->Branch("bdt_numuCCpi0",&_bdt_numuCCpi0,"bdt_numuCCpi0/F");
    _tree->Branch("bdt_numuCC"   ,&_bdt_numuCC   ,"bdt_numuCC/F"   );
    _tree->Branch("bdt_ext"      ,&_bdt_ext      ,"bdt_ext/F"      );
    _tree->Branch("bdt_cosmic"   ,&_bdt_cosmic   ,"bdt_cosmic/F"   );
    _tree->Branch("bdt_global"   ,&_bdt_global   ,"bdt_global/F"   );
    _tree->Branch("pass_antibdt_filter", &_pass_antibdt_filter,"bdt_global/I");
    _mytree = _tree;//not ideal, be careful...
  }
  
  void BDT::resetTTree(TTree* _tree)
  {
    _bdt_nuNCpi0   = 9999.;
    _bdt_numuCCpi0 = 9999.;
    _bdt_numuCC    = 9999.;
    _bdt_ext       = 9999.;
    _bdt_cosmic    = 9999.;
    _bdt_global    = 9999.;
    _pass_antibdt_filter = -9999;
  }

  
  DEFINE_ART_CLASS_TOOL(BDT)
} // namespace analysis

#endif
