////////////////////////////////////////////////////////////////////////
// Class:       CommonOpticalFilter
// Plugin Type: analyzer (art v3_01_02)
// File:        CommonOpticalFilter_module.cc
//
// Generated at Wed Apr 17 16:21:07 2019 by David Caratelli using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"

// save info associated to common optical filter                                                                                                                                       
#include "ubobj/Optical/UbooneOpticalFilter.h"


#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"

class CommonOpticalFilter;


class CommonOpticalFilter : public art::EDAnalyzer {
public:
  explicit CommonOpticalFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CommonOpticalFilter(CommonOpticalFilter const&) = delete;
  CommonOpticalFilter(CommonOpticalFilter&&) = delete;
  CommonOpticalFilter& operator=(CommonOpticalFilter const&) = delete;
  CommonOpticalFilter& operator=(CommonOpticalFilter&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // TTree
  TTree* _tree;
  float _nu_e, _elec_e;
  int   _nelec, _nreco;
  float _vtx_x, _vtx_y, _vtx_z, _vtx_t;
  float  _opfilter_pe_beam, _opfilter_pe_beam_tot, _opfilter_pe_veto, _opfilter_pe_veto_tot;
  int _n_mcs;
  float _mcshr_edep;

  // Declare member data here.

};


CommonOpticalFilter::CommonOpticalFilter(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

  // TTree
  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("opfilter", "opfilter");
  // truth
  _tree->Branch("nu_e"  ,&_nu_e  ,"nu_e/F"  );
  _tree->Branch("vtx_x" ,&_vtx_x ,"vtx_x/F" );
  _tree->Branch("vtx_y" ,&_vtx_y ,"vtx_y/F" );
  _tree->Branch("vtx_z" ,&_vtx_z ,"vtx_z/F" );
  _tree->Branch("vtx_t" ,&_vtx_t ,"vtx_t/F" );
  _tree->Branch("nelec",&_nelec,"nelec/I");
  _tree->Branch("elec_e",&_elec_e,"elec_e/F");
  _tree->Branch("_opfilter_pe_beam",&_opfilter_pe_beam,"opfilter_pe_beam/F");
  _tree->Branch("_opfilter_pe_beam_tot",&_opfilter_pe_beam_tot,"opfilter_pe_beam_tot/F");
  _tree->Branch("_opfilter_pe_veto",&_opfilter_pe_veto,"opfilter_pe_veto/F");
  _tree->Branch("_opfilter_pe_veto_tot",&_opfilter_pe_veto_tot,"opfilter_pe_veto_tot/F");
  _tree->Branch("_n_mcs",&_n_mcs,"n_mcs/I");
  _tree->Branch("_mcshr_edep",&_mcshr_edep,"mcshr_edep/F");

}

void CommonOpticalFilter::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  // load commont-optical-filter output                                                                                                                                                
  art::Handle<uboone::UbooneOpticalFilter> CommonOpticalFilter_h;
  art::InputTag fCommonOpFiltTag("opfiltercommon");
  e.getByLabel(fCommonOpFiltTag, CommonOpticalFilter_h);

  _opfilter_pe_beam     = CommonOpticalFilter_h->PE_Beam();
  _opfilter_pe_beam_tot = CommonOpticalFilter_h->PE_Beam_Total();
  _opfilter_pe_veto     = CommonOpticalFilter_h->PE_Veto();
  _opfilter_pe_veto_tot = CommonOpticalFilter_h->PE_Veto_Total();


  // load MCTruth
  auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth> >("generator");

  auto mct      = mct_h->at(0);
  auto neutrino = mct.GetNeutrino();
  auto nu       = neutrino.Nu();

  _nu_e  = nu.Trajectory().E(0);
  _vtx_x = nu.EndX();
  _vtx_y = nu.EndY();
  _vtx_z = nu.EndZ();
  _vtx_t = nu.T();

  auto const& mcs_h = e.getValidHandle<std::vector<sim::MCShower> >("mcreco");

  _n_mcs = mcs_h->size();
  _mcshr_edep = 0.;

  for (size_t s=0; s < mcs_h->size(); s++) {
    if (mcs_h->at(s).DetProfile().E() > _mcshr_edep) 
      _mcshr_edep = mcs_h->at(s).DetProfile().E();
  }// for all mcshowers

  _tree->Fill();
}

void CommonOpticalFilter::beginJob()
{
  // Implementation of optional member function here.
}

void CommonOpticalFilter::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(CommonOpticalFilter)
