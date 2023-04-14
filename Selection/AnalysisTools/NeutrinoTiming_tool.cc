#ifndef ANALYSIS_NEUTRINOTIMING_CXX
#define ANALYSIS_NEUTRINOTIMING_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       NeutrinoTiming
    // File:        NeutrinoTiming.cc
    //
    //              A basic analysis example
    //
    // Configuration parameters:
    //
    // TBD
    //
    // Created by David Caratelli (dcaratelli@ucsb.edu) on 10/06/2022
    // Code developed by Dante Totani and exported to LArSoft by Erin Yandel
    //
    ////////////////////////////////////////////////////////////////////////

  class NeutrinoTiming : public AnalysisToolBase {

  public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    NeutrinoTiming(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~NeutrinoTiming(){ };
    
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

    void AddDaughters(const art::Ptr<recob::PFParticle>& pfp,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v);
    
  private:

    art::InputTag fPFPproducer;
    art::InputTag fT0producer;
    art::InputTag fFLASHproducer;
    art::InputTag fPMTWFproducer;
    art::InputTag fLogicWFproducer;
    art::InputTag fSpacePointproducer;

    //beam timing plots
    TH1F *H_time;
    TH1F *H_maxH;
    TH1F *H_t0_Beam;
    TH2F *H_TimeVsPh;

    // beam timing variables
    //int _tickSum;
    //int _tickP[32];
    //double _maxSum; 
    //double _timeSum; 
    double _maxP[32]; 
    double _timeP[32];
    double _RWM_T;
    //double _BeamT0;

    float f_shiftoffset;
    bool f_isrun3;
    float f_ccnd1_a;
    float f_ccnd1_b;
    float f_ccnd2_a;
    float f_ccnd2_b;
    float f_ccnd3_a;
    float f_ccnd3_b;
    float f_ccnd3_c;
    float f_ccnd3_d;
    int _run;

    float _interaction_time_modulo; // variable stored in TTree
    float _interaction_time_abs;    // variable stored in TTree
    
    double _nuvtx_x, _nuvtx_y, _nuvtx_z;

    Float_t _evtDeltaTimeNS, _evtTimeNS;
    double  calib[32];

    std::map<unsigned int, unsigned int> _pfpmap;

    void getBeamWF(std::vector<double> *wf_w_03);
    void getPMTwf(std::vector<double> wfsum,
		  std::vector< std::vector<double> > _wf_v);
    //int& tickSum, int tickP[32],
    //double& maxSum, double& timeSum, double maxP[32], double timeP[32]);
    
    void nsbeamtiming(std::vector<std::vector<art::Ptr<recob::SpacePoint> > > spacepoint_v_v);




  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  NeutrinoTiming::NeutrinoTiming(const fhicl::ParameterSet& p)
  {

    fPFPproducer = p.get< art::InputTag >("PFPproducer");
    fT0producer  = p.get< art::InputTag >("T0producer" );
    fPMTWFproducer  = p.get< art::InputTag >("PMTWFproducer" );
    fLogicWFproducer = p.get< art::InputTag >("LogicWFproducer");
    fFLASHproducer  = p.get< art::InputTag >("FLASHproducer" );
    fSpacePointproducer = p.get< art::InputTag >("SpacePointproducer");

    f_shiftoffset = p.get<float>("ShiftOffset", 0);
    f_isrun3 = p.get<bool>("isRun3", false);
    f_ccnd1_a = p.get<float>("ccnd1_a", 0.529594);
    f_ccnd1_b = p.get<float>("ccnd1_b", 7.13804);
    f_ccnd2_a = p.get<float>("ccnd2_a", 0.068752);
    f_ccnd2_b = p.get<float>("ccnd2_b", 2.32023);
    f_ccnd3_a = p.get<float>("ccnd3_a", 0.4697);
    f_ccnd3_b = p.get<float>("ccnd3_b", 0.004233);
    f_ccnd3_c = p.get<float>("ccnd3_c", 0.000001006);
    f_ccnd3_d = p.get<float>("ccnd3_d", -0.195);

    //ns beam timing plots for validation
    art::ServiceHandle<art::TFileService> tfs;
    H_time= tfs->make<TH1F>("H_time","Time PMT",500, 0,6000);
    H_maxH= tfs->make<TH1F>("H_maxH","Max amplitude",800,2000,2100);
    H_t0_Beam= tfs->make<TH1F>("H_t0_Beam","T_0 beam",800,3000,8450);
    H_TimeVsPh= tfs->make<TH2F>("H_TimeVsPh","H_TimeVsPh",  100, -50,50,  100, 0,500);

  }

  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void NeutrinoTiming::configure(fhicl::ParameterSet const & p)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void NeutrinoTiming::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
  {


    _run = e.run();

    const lariov::PmtGainProvider& gain_provider = art::ServiceHandle<lariov::PmtGainService>()->GetProvider();
    for (int i=0; i<32; i++){
      calib[i] = gain_provider.ExtraInfo(i).GetFloatData("amplitude_gain");
      //std::cout << "[NeutrinoTimingDebug] calib[i] for i "  << i << " is " << calib[i] << std::endl;
    }

    // load PMT waveforms
    art::Handle< std::vector< raw::OpDetWaveform > > wf_handle;
    e.getByLabel( fPMTWFproducer, wf_handle );


    if(!wf_handle.isValid()) {
      std::cerr<<"\033[93m[WARNING]\033[00m no raw::OpDetWaveform PMT Waveforms. Skip neutrino timing"<<std::endl;
      return;
    }


    // load logic waveforms
    art::Handle< std::vector< raw::OpDetWaveform > > wf_handle_beam;
    e.getByLabel( fLogicWFproducer, wf_handle_beam );
    
    if(!wf_handle_beam.isValid()) {
      std::cerr<<"\033[93m[WARNING]\033[00m no raw::OpDetWaveform Logic Signal. Skip neutrino timing"<<std::endl;
      return;
    }

    /*
     load waveforms into vectors
    */
    //std::cout << "[NeutrinoTiming] load waveforms into vectors " << std::endl;

    //clear PMT waveform
    auto wfsum = std::vector<double>(1500,0);
    auto _wf_v = std::vector< std::vector<double> >(32,std::vector<double>(1500,0));
    
    for (size_t i=0; i < wf_handle->size(); i++) {
      auto const wf = wf_handle->at(i);
      auto ch = wf.ChannelNumber();
      if (ch >= 32) continue;
      
      for (size_t n=0; n < wf.size(); n++){
        if (n < 1500){
          _wf_v[ch][n] = (wf_handle->at(i))[n];
          wfsum[n] += wf[n];
        }
      }// for all channels
    }// for PMT all waveforms

    //std::cout << "[NeutrinoTiming] calling getPMTwf() " << std::endl;
    getPMTwf(wfsum, _wf_v);
    
    /*
      load Logic waveform
     */
    //std::cout << "[NeutrinoTiming] load logic waveform " << std::endl;

    std::vector<double> *wf_w_03 = new std::vector<double>;
    
    for (size_t i=0; i < wf_handle_beam->size(); i++) {
      auto const wf = wf_handle_beam->at(i);
      auto ch = wf.ChannelNumber();
      //std::cout<<"Channel: "<<ch<<std::endl;
      if (ch != 39) continue;
      
      for (size_t n=0; n < wf.size(); n++){
        if (n < 1500){
          wf_w_03->emplace_back((wf_handle_beam->at(i))[n]);
        }
      }// for all channels
    }// for all waveforms

    getBeamWF(wf_w_03);

    
    // load tracks previously created for which T0 reconstruction is requested                                                                                                                                

    //std::cout << "[NeutrinoTiming] load T0 tagged information " << std::endl;

    art::Handle<std::vector<anab::T0> > t0_h;
    e.getByLabel( fT0producer , t0_h );

    // make sure tracks look good                                                                                                                                                                             
    if(!t0_h.isValid()) {
      std::cerr<<"\033[93m[WARNING]\033[00m no anab::T0 for flash-matching. Skip flash-matching"<<std::endl;
      return;
    }
    
    art::ValidHandle<std::vector<recob::PFParticle>> pfp_h = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    
    art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fPFPproducer);
    auto const& spacepoint_h = e.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointproducer);  
    art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, e, fSpacePointproducer);
    
    
    // figure out which PFP is the neutrino
    //size_t nupfp = 0;
    for (auto pfp : slice_pfp_v) {
      
      if (pfp->IsPrimary() == false) continue;
      
      if ( (pfp->PdgCode() == 12) || (pfp->PdgCode() == 14) ) {
	
	//nupfp = pfp->Self();
	
	auto vtx = pfp.get<recob::Vertex>();
	if (vtx.size() != 1) {
	  std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
	  return;
	}
	std::cout<<"Neutrino! \n";
	// save vertex to array
	TVector3 nuvtx;
	Double_t xyz[3] = {};
	vtx.at(0)->XYZ(xyz);
	nuvtx = TVector3(xyz[0],xyz[1],xyz[2]);
	_nuvtx_x = nuvtx.X();
	_nuvtx_y = nuvtx.Y();
	_nuvtx_z = nuvtx.Z();
	
      }// if neutrino
    }// for all PFPs in slice
    
    
    // now build vectors of PFParticles, space-points, and hits for this slice
    // std::vector<recob::PFParticle> pfp_v;
    std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
    std::vector<std::vector<art::Ptr<recob::SpacePoint> > > spacepoint_v_v;
    std::vector<std::vector<art::Ptr<recob::Hit> > > hit_v_v; 
    
    // fill map: pfparticle Self() -> index/key
    _pfpmap.clear();
    for (unsigned int p=0; p < pfp_h->size(); p++) _pfpmap[pfp_h->at(p).Self()] = p;
    
    // loop through all PFParticles
    for (size_t p=0; p < pfp_h->size(); p++) {
      
      auto const& pfp = pfp_h->at(p);
      const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
      
      // only primary PFPs have a flash-match score
      if (pfp.IsPrimary() == false) continue;
      
      // here only interested in the neutrino slice
      auto PDG = fabs(pfp.PdgCode());
      if ( (PDG != 12) && (PDG != 14) ) continue;
      
      AddDaughters(pfp_ptr, pfp_h, pfp_ptr_v);
      
      // go through these pfparticles and fill info needed for matching
      for (size_t i=0; i < pfp_ptr_v.size(); i++) {    
	auto key = pfp_ptr_v.at(i).key();
	// recob::PFParticle ipfp = *pfp_ptr_v.at(i);
	// pfp_v.push_back(ipfp);  
	//auto const& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
	const std::vector< art::Ptr<recob::SpacePoint> >& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
	std::vector< art::Ptr<recob::Hit> > hit_ptr_v;
	for (size_t sp=0; sp < spacepoint_ptr_v.size(); sp++) {
	  auto const& spkey = spacepoint_ptr_v.at(sp).key();
	  const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = spacepoint_hit_assn_v.at( spkey );
	  for (size_t h=0; h < this_hit_ptr_v.size(); h++) {
	    hit_ptr_v.push_back( this_hit_ptr_v.at( h ) );
	  }// for all hits associated to this spacepoint
	}// fpr all spacepoints  
	spacepoint_v_v.push_back( spacepoint_ptr_v );
	hit_v_v.push_back( hit_ptr_v );
      }// for all pfp pointers
      //flashmatch::SliceCandidate slice(pfp_ptr_v, spacepoint_v_v, hit_v_v, m_chargeToNPhotonsTrack, m_chargeToNPhotonsShower);
      
    }// for all PFPs
    
    //std::cout << "[NeutrinoTiming] run ns beam timing " << std::endl;
    nsbeamtiming(spacepoint_v_v);

    _interaction_time_modulo = _evtDeltaTimeNS;
    _interaction_time_abs = _evtTimeNS;
  
  return;
  }
  
  void NeutrinoTiming::analyzeEvent(art::Event const &e, bool fData)
  {
 
    return;
  }

  void NeutrinoTiming::setBranches(TTree* _tree)
  {
    _tree->Branch("interaction_time_abs",&_interaction_time_abs,"interaction_time_abs/F");
    _tree->Branch("interaction_time_modulo",&_interaction_time_modulo,"interaction_time_modulo/F");
  }
  
  void NeutrinoTiming::resetTTree(TTree* _tree)
  {
    _interaction_time_abs    = std::numeric_limits<float>::lowest();
    _interaction_time_modulo = std::numeric_limits<float>::lowest();
  }

  void NeutrinoTiming::AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v) {
  
    auto daughters = pfp_ptr->Daughters();
  
    pfp_v.push_back(pfp_ptr);
  
    std::cout << "\t PFP w/ PdgCode " << pfp_ptr->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;
  
    for(auto const& daughterid : daughters) {

      if (_pfpmap.find(daughterid) == _pfpmap.end()) {
	std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
	continue;
      }
    
      const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, _pfpmap.at(daughterid) );
    
      AddDaughters(pfp_ptr, pfp_h, pfp_v);
    
    }// for all daughters


  
    return;
  }


  void NeutrinoTiming::getPMTwf(std::vector<double> wfsum,
				std::vector< std::vector<double> > _wf_v)
  //art::Handle< std::vector< raw::OpDetWaveform > > wf_handle, 
  //int& tickSum, int tickP[32],
  //double& maxSum, double& timeSum, double maxP[32], double timeP[32])
  {
    
    //analyze waveforms
    //std::cout << "[NeutrinoTiming:getPMTwf()] check1 " << std::endl;

    
    //=======================================================================================================
    //=======================================================================================================
    double Help_wf_v[32][500];
    double x_wf_v[500], Raw_wf_v[500], Base_wf_v[500], Norm_wf_v[500];
    double maxZ,max0,base;   int basebinmax,tick;
    double tca,tcb,tcc, TT[32], max[32];
    
    int TF,TB,tickF,tickB,FB, Nss,is;
    double  maxZhelp1,maxZhelp2,maxZhelp3, tickFit1,tickFit2;
    
    //Raw waveform saturation
    int saturation=4094;
    //Parameters for discrete saturated WF reconstruction (1st step)
    double Frac[100]={1.0, 0.951931, 0.93356, 0.838637, 0.719408, 0.701042, 0.565673, 0.44655, 0.39447, 0.352336, 0.28716, 0.245364, 0.216771, 0.194888, 0.178976, 0.16844, 0.161732, 0.157486, 0.154762, 0.153136, 0.152468, 0.152299, 0.152147, 0.151456, 0.150069, 0.148089, 0.145565, 0.142516, 0.139369, 0.136467, 0.133871, 0.13141, 0.128926, 0.126237, 0.12313, 0.119611, 0.115946, 0.112453, 0.111706, 0.109225, 0.106281, 0.103499, 0.100926, 0.0985215, 0.0961512, 0.0938219, 0.0917209, 0.0898817, 0.0882937, 0.0868338, 0.0854753, 0.0841353, 0.0827237, 0.081198, 0.0794796, 0.077565, 0.0756475, 0.0738863, 0.0722821, 0.0708473, 0.0695119, 0.0682504, 0.0672023, 0.0663549, 0.0656337, 0.0649918, 0.0645003, 0.0641535, 0.0638046, 0.0633435, 0.0627506, 0.0621379, 0.0615464, 0.0609178, 0.0601846, 0.0592098, 0.0580465, 0.0568861, 0.0559024, 0.0550731, 0.0541904, 0.0532532, 0.0524181, 0.0517606, 0.0512326, 0.0507392, 0.0502093, 0.0495968, 0.0488915, 0.0480173, 0.0470195, 0.0459744, 0.0448855, 0.0437359, 0.0425199, 0.0412832, 0.0400036, 0.038688, 0.0373173, 0.0358925};
    //double Delay[100]={0.0, 0.0686367, 0.0948574, 0.230854, 0.405159, 0.432604, 0.642806, 0.845613, 0.942555, 1.02616, 1.16775, 1.26923, 1.34523, 1.40802, 1.45675, 1.49068, 1.51306, 1.52757, 1.53702, 1.54272, 1.54507, 1.54567, 1.54621, 1.54865, 1.55359, 1.56069, 1.56984, 1.58104, 1.59279, 1.60379, 1.61377, 1.62337, 1.63319, 1.64398, 1.65665, 1.67129, 1.68688, 1.70208, 1.70537, 1.71644, 1.72981, 1.74271, 1.75486, 1.76644, 1.77806, 1.78969, 1.80036, 1.80987, 1.81819, 1.82594, 1.83324, 1.84054, 1.84831, 1.85684, 1.86658, 1.87763, 1.88891, 1.89947, 1.90926, 1.91815, 1.92656, 1.93462, 1.9414, 1.94695, 1.95171, 1.95598, 1.95928, 1.96162, 1.96398, 1.96711, 1.97117, 1.97539, 1.97951, 1.98391, 1.98909, 1.99605, 2.00448, 2.01302, 2.02038, 2.02665, 2.03342, 2.04069, 2.04726, 2.0525, 2.05674, 2.06073, 2.06505, 2.0701, 2.07597, 2.08334, 2.09188, 2.10099, 2.11065, 2.12106, 2.13232, 2.14403, 2.15646, 2.16957, 2.18362, 2.19868 };
    //Fit function parameters for saturated WF reconstruction (2nd step)
    double pLL[8]={3.93256,4.31002,2.44182,5.12491,0.830928,0.231375,50.9081,-2.69014};
    
    //=====================================================================================================
    //=======================================================================================================
    
    std::vector<float> max_pmt_v(32,0);
    std::vector<int>   max_tick_v(32,0);
    
    for(int q=0; q<32; q++){
      for(int i=0; i<500; i++){
	Help_wf_v[q][i]=_wf_v[q][i];
	if (i > 100 && i < 110) 
	  //std::cout << "[NeutrinoTimingDebug] PMT " << q << " @ tick " << i << " has ADC " << Help_wf_v[q][i] << std::endl;
	if (Help_wf_v[q][i] > max_pmt_v[q]) { max_pmt_v[q] = Help_wf_v[q][i]; max_tick_v[q] = i; }
      }
    }

    for (size_t nn=0; nn < 32; nn++) {
      //std::cout << "[NeutrinoTimingDebug] PMT " << nn << " has maximum tick @ " <<  max_tick_v[nn] << " with value " << max_pmt_v[nn] << "." << std::endl;
    }

    _wf_v.clear(); _wf_v.shrink_to_fit();
    
    //std::cout << "[NeutrinoTiming:getPMTwf()] check2 " << std::endl;
    
    //-----------------------------------------------------------------------------------------------------
    for(int q=0; q<32; q++){
      TT[q]=-9999.;
      max[q]=-9999.;
      maxZ=0.;
      max0=0.;
      base=0.;
      tick=0;
      tickB=0;
      tickF=0;
      TF=0;
      TB=0;
      //Getting raw waveform (Raw_wf_v[i]) only for i<500 since the beam window is between 3 and 5 us -> [i>3*64 && i<5*64]
      for(int i=0; i<500; i++){
	x_wf_v[i]=i*1.0;
	Raw_wf_v[i]=Help_wf_v[q][i];
      }
      //Getting raw wf max amplitude and max amp tick
      for(int i=3*64; i<5*64; i++){if(maxZ<Raw_wf_v[i]){maxZ=Raw_wf_v[i]; tick=i;}}
      //Baseline removal
      TH1F *basehelp= new TH1F("basehelp","basehelp",400, 1900,2200);
      basebinmax=0; for(int i=0; i<3*64; i++){basehelp->Fill(Raw_wf_v[i]);}
      basebinmax=basehelp->GetMaximumBin(); base=basehelp->GetXaxis()->GetBinCenter(basebinmax);
      basehelp->Delete();
      //Getting wf max amp after baseline removal (this is proportional to number of Photons in the rising endge)
      //getting wf baseline subtracted and wf baseline subtracted and normalized for the max amp.
      for(int i=0; i<500; i++){max0=maxZ-base;
	Base_wf_v[i]=Raw_wf_v[i]-base; Norm_wf_v[i]=Base_wf_v[i]/max0;}
      //fitting the normalized baseline subtracted wf
      TGraph *gr = new TGraph(500,x_wf_v,Norm_wf_v);
      TF1 *fit = new TF1("fit","[2]*exp(-TMath::Power(([0]-x)/[1],4))",tick-10, tick);
      fit->SetParameters(tick,2,1);  gr->Fit("fit","Q","",tick-10, tick);
      tca=fit->GetParameter(0);  tcb=fit->GetParameter(1);  tcc=fit->GetParameter(2);
      //timing is the risign edge half height
      TT[q]=(tca-abs(tcb*TMath::Power(-log(0.5/tcc),0.25)))/0.064; max[q]=max0;
      //----------------------------------------------
      //check for saturated wf
      //std::cout << "[NeutrinoTiming:getPMTwf()] check3 " << std::endl;
      if(maxZ<=saturation){
	TT[q]=TT[q];
	max[q]=max[q];
	//std::cout << "[NeutrinoTimingDebug] not saturated!" << std::endl;
      }
      else if(maxZ>saturation) {
	
	//std::cout << "[NeutrinoTiming:getPMTwf()] check3.0 " << std::endl;
	//counting the number of ticks above the saturation
	for(int i=3*64; i<7*64; i++){
	  if(TF==0){if(Raw_wf_v[i+1]>4094 && Raw_wf_v[i]<=4094){tickF=i; TF=1;}}
	  if(TB==0){if(Raw_wf_v[i]>4094 && Raw_wf_v[i+1]<=4094){tickB=i; TB=1;}}}
	FB=tickB-tickF;  if(FB>99){FB=99;}
	//amplitude discrete correction
	//std::cout << "[NeutrinoTiming:getPMTwf()] check3.0.0 " << std::endl;
	maxZhelp1=maxZ/Frac[FB]; tick=tickF; Nss=0; is=0;
      for(int i=3*64; i<7*64; i++){if(Raw_wf_v[i]<4095){Nss=Nss+1;}}
      
      //std::cout << "[NeutrinoTiming:getPMTwf()] check3.1 " << std::endl;
      
      double txSS[256],tySS[256],txSS2[256],tySS2[256];
      
      for(int i=3*64; i<7*64; i++){if(Raw_wf_v[i]<4095){txSS[is]=i*1.0; tySS[is]=Raw_wf_v[i]/maxZhelp1; is=is+1;}}
      TGraph *g1 = new TGraph(Nss,txSS,tySS);
      //for (int uu=0; uu < (7*64); uu++)
      //std::cout << "[NeutrinoTimingDebug] value @ tick " << uu << " is " << tySS[uu] << std::endl;
      TF1 *fitS1 = new TF1("fitS1","[9]*(exp(-TMath::Power(([0]-(x-[8]))/[1],4))*0.5*(TMath::Erf(-(x-[8])-[7])+1.0)+([5]+[4]*exp(-TMath::Power(([2]-(x-[8]))/[3],2)))*exp((-(x-[8]))/[6])*0.5*(TMath::Erf([7]+(x-[8]))+1.0))",tick-30, tick+250);
      fitS1->SetParameters(pLL[0],pLL[1],pLL[2],pLL[3],pLL[4],pLL[5],pLL[6],pLL[7],tick,1.);
      for(int i=0; i<8; i++){fitS1->FixParameter(i,pLL[i]);} g1->Fit("fitS1","Q","",tick-30, tick+250);
      tickFit1=fitS1->GetParameter(8); maxZhelp2=fitS1->GetParameter(9);  maxZhelp3=maxZhelp1/maxZhelp2;
      //amplitude fit correction
      //std::cout << "[NeutrinoTiming:getPMTwf()] check3.2 " << std::endl;
      for(int i=0; i<Nss; i++){txSS2[i]=txSS[i]; tySS2[i]=tySS[i]/maxZhelp2;}
      TGraph *g2 = new TGraph(Nss,txSS2,tySS2);
      //std::cout << "[NeutrinoTiming:getPMTwf()] check3.3 " << std::endl;
      TF1 *fitS2 = new TF1("fitS2","exp(-TMath::Power(([0]-(x-[8]))/[1],4))*0.5*(TMath::Erf(-(x-[8])-[7])+1.0)+([5]+[4]*exp(-TMath::Power(([2]-(x-[8]))/[3],2)))*exp((-(x-[8]))/[6])*0.5*(TMath::Erf([7]+(x-[8]))+1.0)",tick-30, tick+250);
      //std::cout << "[NeutrinoTiming:getPMTwf()] check3.4 " << std::endl;
      fitS2->SetParameters(pLL[0],pLL[1],pLL[2],pLL[3],pLL[4],pLL[5],pLL[6],pLL[7],tickFit1);
      for(int i=0; i<8; i++){fitS2->FixParameter(i,pLL[i]);}
      //std::cout << "[NeutrinoTiming:getPMTwf()] check3.5 " << std::endl;
      g2->Fit("fitS2","Q","",tick-30, tick+250);  tickFit2=fitS2->GetParameter(8);
      //timing is the risign edge half height
      TT[q]=tickFit2/0.064; max[q]=maxZhelp3;
    }// if not saturated
    //-------------------------------------------------------------------------------------------------------
    //std::cout << "[NeutrinoTiming:getPMTwf()] check3a " << std::endl;
    H_time->Fill(TT[q]);
    //std::cout << "[NeutrinoTiming:getPMTwf()] check3b " << std::endl;
    }
    //-------------------------------------------------------------------------------------------------------
    
    //std::cout << "[NeutrinoTiming:getPMTwf()] check4 " << std::endl;
    for(int q=0; q<32; q++){
      _maxP[q]=max[q];
      _timeP[q]=TT[q];
      //std::cout << "[NeutrinoTimingDebug] : PMT "  << q << " has maximum @ tick " << _timeP[q] << " with value " << _maxP[q] << std::endl;
    } //only two variables needed
    
  }
  
  
  void NeutrinoTiming::getBeamWF(
				 std::vector<double> *wf_w_03
				 //art::Handle< std::vector< raw::OpDetWaveform > > wf_handle_beam
				 )
  {

    
    //=======================================================================================================
    //=======================================================================================================

    double beamBase,BBmax,wx[500],wy[500],pca,pcb,pcc, TT = -99999.0;
    int Btick,tickMax;
    //=======================================================================================================
    //-------baseline calculation-------------------------------------------
    beamBase=0.; BBmax=0.; Btick=0; TT=-9999.;
    for(int i=0; i<4*64; i++){beamBase=beamBase+wf_w_03->at(i);}
    beamBase=beamBase/(4.0*64.0);
    //baseline subtraction
    for(int i=0; i<500; i++){wx[i]=i*1.0; wy[i]=wf_w_03->at(i)-beamBase;
    //max amplitude
    if(BBmax<wy[i]){BBmax=wy[i]; Btick=i;}}
    H_maxH->Fill(BBmax);
    if(BBmax>2000 && BBmax<2100){
    //wf normalization
    for(int i=0; i<500; i++){wy[i]=wy[i]/BBmax;}
    //wf max check
    tickMax=0;
    for(int i=Btick-20; i<Btick+10; i++){if(wy[i-1]<1 && wy[i]==1){tickMax=i;}}
    //wf fit
    TGraph *gr0 = new TGraph(500,wx,wy);
    TF1 *fit = new TF1("fit","[2]*exp(-TMath::Power((x-[0])/[1],4))",tickMax-6, tickMax);
    fit->SetParameters(tickMax,2,1);     gr0->Fit("fit","Q","",tickMax-6, tickMax);
    pca=fit->GetParameter(0); pcb=fit->GetParameter(1); pcc=fit->GetParameter(2);
    //timing is the risign edge half height
    TT=(pca-abs(pcb*TMath::Power(-log(0.5/pcc),0.25)))/0.064;
    //std::cout << "[NeutrinoTimingDebug] RWM time : " << TT << std::endl;
    H_t0_Beam->Fill(TT);
    _RWM_T = TT;

    }
}

  void NeutrinoTiming::nsbeamtiming(
				     std::vector<std::vector<art::Ptr<recob::SpacePoint> > > spacepoint_v_v
				    )
  {
    /*    
	  int tickSum = _tickSum;
	  int tick[32] = _tickP;
	  double maxSum= _maxSum;
	  double timeSum= _timeSum;
	  double max[32] = _maxP;
	  double time[32] = _tickP;
	  double BeamT0 = _BeamT0;
    */
    //getPMTwf(e,tickSum,tick, maxSum,timeSum,max,time);
    //std::cout<<"tickSum: "<<tickSum<<" maxSum: "<<maxSum<<" timeSum: "<<timeSum<<std::endl;
    //std::cout<<"tick[0]: "<<tick[0]<<" max[0]: "<<max[0]<<" time[0]: "<<time[0]<<std::endl;
    //BeamT0 = getBeamWF(e);
    //std::cout<<"Beam T0: "<<BeamT0<<std::endl;
    

    Float_t x =	_nuvtx_x;
    Float_t y =	_nuvtx_y;
    Float_t z =	_nuvtx_z;

    //std::cout << "[NeutrinoTimingDebug] vertex @ (" << x << ", " << y << ", " << z << ")" << std::endl;
    
    std::vector<float> *sps_x = new std::vector<float>;
    std::vector<float> *sps_y = new std::vector<float>;
    std::vector<float> *sps_z = new std::vector<float>;
    
    for (size_t s=0; s < spacepoint_v_v.size(); s++) {
      auto sps_v = spacepoint_v_v[s];
      for (size_t ss=0; ss < sps_v.size(); ss++) {
	auto sps = sps_v[ss];
	sps_x->push_back( sps->XYZ()[0] );
	sps_y->push_back( sps->XYZ()[1] );
	sps_z->push_back( sps->XYZ()[2] );
      }
    }
    
    //std::cout << "[NeutrinoTimingDebug] there are " << sps_z->size() << " spacepoints" << std::endl;
    
    double max[32],time[32];
    double BeamT0 = -99999.;
    for (int ii=0; ii < 32; ii++){
      time[ii] = _timeP[ii];
      max[ii] = _maxP[ii];
    }
    
    //getPMTwf(e,max,time);
    BeamT0 = _RWM_T;
    
    double PMT0[3]={-11.4545, -28.625, 990.356};  double PMT1[3]={-11.4175, 27.607, 989.712};
    double PMT2[3]={-11.7755, -56.514, 951.865};  double PMT3[3]={-11.6415, 55.313, 951.861};
    double PMT4[3]={-12.0585, -56.309, 911.939};  double PMT5[3]={-11.8345, 55.822, 911.065};
    double PMT6[3]={-12.1765, -0.722, 865.599};   double PMT7[3]={-12.3045, -0.502, 796.208};
    double PMT8[3]={-12.6045, -56.284, 751.905};  double PMT9[3]={-12.5405, 55.625, 751.884};
    double PMT10[3]={-12.6125, -56.408, 711.274}; double PMT11[3]={-12.6615, 55.8, 711.073};
    double PMT12[3]={-12.6245, -0.051, 664.203};  double PMT13[3]={-12.6515, -0.549, 585.284};
    double PMT14[3]={-12.8735, 55.822, 540.929};  double PMT15[3]={-12.6205, -56.205, 540.616};
    double PMT16[3]={-12.5945, -56.323, 500.221}; double PMT17[3]={-12.9835, 55.771, 500.134};
    double PMT18[3]={-12.6185, -0.875, 453.096};  double PMT19[3]={-13.0855, -0.706, 373.839};
    double PMT20[3]={-12.6485, -57.022, 328.341}; double PMT21[3]={-13.1865, 54.693, 328.212};
    double PMT22[3]={-13.4175, 54.646, 287.976};  double PMT23[3]={-13.0075, -56.261, 287.639};
    double PMT24[3]={-13.1505, -0.829, 242.014};  double PMT25[3]={-13.4415, -0.303, 173.743};
    double PMT26[3]={-13.3965, 55.249, 128.354};  double PMT27[3]={-13.2784, -56.203, 128.18};
    double PMT28[3]={-13.2375, -56.615, 87.8695}; double PMT29[3]={-13.5415, 55.249, 87.7605};
    double PMT30[3]={-13.4345, 27.431, 51.1015};  double PMT31[3]={-13.1525, -28.576, 50.4745};
    double PMT[32][3];    for(int j=0; j<3; j++){ PMT[30][j]=PMT30[j]; PMT[31][j]=PMT31[j];
      PMT[0][j]=PMT0[j];   PMT[10][j]=PMT10[j]; PMT[20][j]=PMT20[j]; PMT[1][j]=PMT1[j];   PMT[11][j]=PMT11[j];
      PMT[21][j]=PMT21[j]; PMT[2][j]=PMT2[j];   PMT[12][j]=PMT12[j]; PMT[22][j]=PMT22[j]; PMT[3][j]=PMT3[j];
      PMT[13][j]=PMT13[j]; PMT[23][j]=PMT23[j]; PMT[4][j]=PMT4[j];   PMT[14][j]=PMT14[j]; PMT[24][j]=PMT24[j];
      PMT[5][j]=PMT5[j];   PMT[15][j]=PMT15[j]; PMT[25][j]=PMT25[j]; PMT[6][j]=PMT6[j];   PMT[16][j]=PMT16[j];
      PMT[18][j]=PMT18[j]; PMT[28][j]=PMT28[j]; PMT[9][j]=PMT9[j];   PMT[19][j]=PMT19[j]; PMT[29][j]=PMT29[j];
      PMT[26][j]=PMT26[j]; PMT[7][j]=PMT7[j];   PMT[17][j]=PMT17[j]; PMT[27][j]=PMT27[j]; PMT[8][j]=PMT8[j];}
    double offset[32]={1.03002, -5.18104, -2.11164, -5.99395, -1.25798, 0.633079, 2.87666, 2.21969, 0.885092, 2.35423,
		       -1.63039, -1.83775, -0.859883, 3.4741, 1.84833, 1.58233, -2.71783, 0, 3.18776, 0.982666, 0.728438, 0.280592, -5.27068,
		       -3.27857, -1.41196, 1.59643, 1.41425, -1.62682, -2.55772, 1.49136, -0.522791, 0.974533};
    
    
    //================================================================================================================
    double gap=18.936;
    double MaxLim=2.5;
    std::vector<int> N_pmt;
    double ccnd1, ccnd2,ccnd3;
    double Ph_Tot, RWM_T, nuToF, DPh,DLh, tPhelp,tp;
    //double Med_TT0,Med_TT1,Med_TT2;
    double Med_TT3=-9999.;
    double TT_merged = -9999.;
    //===================================================================================================================
    //===================================================================================================================
    Ph_Tot=0.;
    N_pmt.clear();
    for(int q=0; q<32; q++) {
      max[q]=max[q]/calib[q];
      //std::cout << "[NeutrinoTimingDebug] max[q] for q = " << q << " is " << max[q] << std::endl;
      if((max[q]>MaxLim && q!=17 && q!=28) && (time[q]>3000.0 && time[q]<5000.0)) {
	N_pmt.push_back(q);
	Ph_Tot=Ph_Tot+max[q];
      }
    } 
   //--------------------------------------------------------------------------------------------------------------------
    //std::cout << "[NeutrinoTimingDebug] PMT size : " << N_pmt.size() << std::endl;
    if(N_pmt.size()>2){
      RWM_T=BeamT0;
      nuToF=z*0.033356;
      std::vector<double> timeProp = std::vector<double>(N_pmt.size(),0);
      for(uint i=0; i<N_pmt.size(); i++){
        tp=5000000000.0;
        for(uint j=0; j<sps_x->size(); j++){
          DPh=abs(sqrt(TMath::Power(x-sps_x->at(j),2)+TMath::Power(y-sps_y->at(j),2)+TMath::Power(z-sps_z->at(j),2)));
          DLh=abs(sqrt(TMath::Power(PMT[N_pmt.at(i)][0]-sps_x->at(j),2)+TMath::Power(PMT[N_pmt.at(i)][1]-sps_y->at(j),2)+TMath::Power(PMT[N_pmt.at(i)][2]-sps_z->at(j),2)));
          tPhelp=(DPh*0.033356)+(DLh*0.0746);
          if(tPhelp<tp){tp=tPhelp;}
        }
        timeProp[i]=tp;
      }
      
      double TT3_array[32];
      if(f_isrun3 && _run>17200 && _run<17400){if(RWM_T>5450){ f_shiftoffset=118.3;}}
      float RWM_offset = 5700.0 - f_shiftoffset;
      
      for(uint i=0; i<N_pmt.size(); i++){
	ccnd1= timeProp[i]*(f_ccnd1_a)-(f_ccnd1_b);
	ccnd2= max[N_pmt.at(i)]*(f_ccnd2_a)-(f_ccnd2_b);
	if(Ph_Tot>150){ccnd3=f_ccnd3_a-f_ccnd3_b*Ph_Tot+f_ccnd3_c*Ph_Tot*Ph_Tot;}
	else if (Ph_Tot<=150){ccnd3=f_ccnd3_d;}
	
	//all the corrections
	TT3_array[i]=(time[N_pmt.at(i)])-RWM_T+RWM_offset-nuToF-timeProp[i]-offset[N_pmt.at(i)]+ccnd1+ccnd2+ccnd3;
      }
      //for (int yy=0; yy < 32; yy++) { std::cout << "[NeutrinoTimingDebug] PMT " << yy << " timing " << TT3_array[yy] << std::endl;}
      Med_TT3=TMath::Median((Long64_t)N_pmt.size(),TT3_array);
      //Fill a 2d histogram with  TT3_array[i] vs max[N_pmt.at(i)] this is usefull to check for any errors
      for(uint i=0; i<N_pmt.size(); i++){
	H_TimeVsPh->Fill( TT3_array[i]-Med_TT3, max[N_pmt.at(i)]);
      }
    }// if PMT size > 2

    _evtTimeNS = Med_TT3;
    //std::cout<<"[NeutrinoTimingDebug] evtTimeNS: "<< Med_TT3<<std::endl;
    //Merge Peaks
    double Shift=3166.9;
    double TThelp=Med_TT3-Shift+gap*0.5;

    _evtTimeNS = TThelp;
    
    //std::cout << "[NeutrinoTimingDebug] TThelp : "  << TThelp << std::endl;
    //merge peaks
    if(TThelp>=0 && TThelp<gap*81.0){
      TT_merged=(TThelp-(int((TThelp)/gap))*gap)-gap*0.5;
    }
    else {TT_merged=-9999;}

    _evtDeltaTimeNS = TT_merged;
    
    //std::cout<<"evtDeltaTimeNS: "<< TT_merged<<std::endl;
    
  }
  
  
  
  DEFINE_ART_CLASS_TOOL(NeutrinoTiming)
} // namespace analysis

#endif
