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
    // Created by Giuseppe Cerati (cerati@fnal.gov) on 03/15/2019
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
    TH1F *H_MaxHSum;
    TH1F *H_TickSum;
    TH1F *H_MaxH;
    TH1F *H_Tick;
    TH2F *H_TickSumVsPMT;
    TH1F *H_time;
    TH1F *H_timeS;
    TH1F *H_maxH;
    TH1F *H_tickHalf;
    TH1F *H_MaxDTick;
    TH1F *H_t0_Beam;

    // beam timing variables
    int _tickSum;
    int _tickP[32];
    double _maxSum; 
    double _timeSum; 
    double _maxP[32]; 
    double _timeP[32];
    double _BeamT0;
    
    double _nuvtx_x, _nuvtx_y, _nuvtx_z;

    Float_t _evtDeltaTimeNS, _evtTimeNS;

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
    fLogicWFproducer = p.get< art::InputTag >("pmtreadout:UnspecifiedLogic");
    fFLASHproducer  = p.get< art::InputTag >("FLASHproducer" );
    fSpacePointproducer = p.get< art::InputTag >("SpacePointproducer");

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

    getPMTwf(wfsum, _wf_v);
    
    /*
      load Logic waveform
     */
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

    nsbeamtiming(spacepoint_v_v);
  
  return;
  }
  
  void NeutrinoTiming::analyzeEvent(art::Event const &e, bool fData)
  {
 
    return;
  }

  void NeutrinoTiming::setBranches(TTree* _tree)
  {
  }
  
  void NeutrinoTiming::resetTTree(TTree* _tree)
  {
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
    
    //=======================================================================================================
    //=======================================================================================================
    
    double tx[500],tySum[500],ty[32][500], tyN[32][500],tySN[500];
    double BBs,maxS,base[32],max[32]; int basebinmax,tickS,tick[32],pass[32], Fns;
    double TTS,TT[32],txF[10],tyF[10],tca,tcb,tcc;
    
    //=====================================================================================================
    //=======================================================================================================
    
    BBs=0.; maxS=0.; tickS=0; for(int q=0; q<32; q++){base[q]=0.; max[q]=0.; tick[q]=0; pass[q]=0;}
    for(int j=0; j<500; ++j){tx[j]=j;
      ty[0][j]= _wf_v[0][j]; ty[10][j]= _wf_v[10][j]; ty[20][j]= _wf_v[20][j];
      ty[1][j]= _wf_v[1][j]; ty[11][j]= _wf_v[11][j]; ty[21][j]= _wf_v[21][j];
      ty[2][j]= _wf_v[2][j]; ty[12][j]= _wf_v[12][j]; ty[22][j]= _wf_v[22][j];
      ty[3][j]= _wf_v[3][j]; ty[13][j]= _wf_v[13][j]; ty[23][j]= _wf_v[23][j];
      ty[4][j]= _wf_v[4][j]; ty[14][j]= _wf_v[14][j]; ty[24][j]= _wf_v[24][j];
      ty[5][j]= _wf_v[5][j]; ty[15][j]= _wf_v[15][j]; ty[25][j]= _wf_v[25][j];
      ty[6][j]= _wf_v[6][j]; ty[16][j]= _wf_v[16][j]; ty[26][j]= _wf_v[26][j];
      ty[7][j]= _wf_v[7][j]; ty[17][j]= _wf_v[17][j]; ty[27][j]= _wf_v[27][j];
      ty[8][j]= _wf_v[8][j]; ty[18][j]= _wf_v[18][j]; ty[28][j]= _wf_v[28][j];
      ty[9][j]= _wf_v[9][j]; ty[19][j]= _wf_v[19][j]; ty[29][j]= _wf_v[29][j];
      tySum[j]= wfsum.at(j); ty[30][j]= _wf_v[30][j]; ty[31][j]= _wf_v[31][j];}
    
    _wf_v.clear(); _wf_v.shrink_to_fit(); wfsum.clear(); wfsum.shrink_to_fit();
    
    
    //-----------------------------------------------------------------------------------------------------
    TH1F *basehelpS= new TH1F("basehelpS","basehelpS",300, 65400,65700);
    basebinmax=0; for(int i=0; i<180; i++){basehelpS->Fill(tySum[i]);}
    basebinmax=basehelpS->GetMaximumBin(); BBs=basehelpS->GetXaxis()->GetBinCenter(basebinmax);
    basehelpS->Delete();    for(int j=0; j<500; ++j){tySum[j]=tySum[j]-BBs;}
    for(int q=0; q<32; q++){TH1F *basehelp= new TH1F("basehelp","basehelp",400, 1900,2200);
      basebinmax=0; for(int i=0; i<180; i++){basehelp->Fill(ty[q][i]);}
      basebinmax=basehelp->GetMaximumBin(); base[q]=basehelp->GetXaxis()->GetBinCenter(basebinmax);
      basehelp->Delete(); for(int j=0; j<500; ++j){ty[q][j]=ty[q][j]-base[q];}}
    //----------------------------------------------------------------------------------------------------
    for(int i=180; i<350; i++){if(maxS<tySum[i]){maxS=tySum[i]; tickS=i;}}
    for(int j=0; j<500; ++j){tySN[j]=tySum[j]/maxS;}
    H_TickSum->Fill(tickS); H_MaxHSum->Fill(maxS);
    for(int q=0; q<32; q++){for(int i=180; i<480; i++){if(max[q]<ty[q][i]){max[q]=ty[q][i]; tick[q]=i;}}
      for(int j=0; j<500; ++j){tyN[q][j]=ty[q][j]/max[q];}
      if(max[q]>40 && (tick[q]>tickS-10 && tick[q]<tickS+10)){pass[q]=1;} else {pass[q]=-1;}
      H_Tick->Fill(tick[q]); H_MaxH->Fill(max[q]); if(pass[q]==1){H_TickSumVsPMT->Fill(tickS,tick[q]);}
    }
    //-----------------------------------------------------------------------------------------------------
    Fns=tickS-9;
    for(int i=0; i<10; i++){txF[i]=tx[Fns+i]; tyF[i]=tySN[Fns+i];} TGraph *grSt = new TGraph(10,txF,tyF);
    TF1 *fit = new TF1("fit","[2]*exp(-TMath::Power(([0]-x)/[1],4))",txF[0], txF[9]);
    fit->SetParameters(tx[tickS],2,1);  grSt->Fit("fit","Q","",txF[0],txF[9]);
    tca=fit->GetParameter(0);  tcb=fit->GetParameter(1);  tcc=fit->GetParameter(2);
    TTS=(tca-abs(tcb*TMath::Power(-log(0.5/tcc),0.25)))/0.064;
    H_timeS->Fill(TTS);
    //-------------------------------------------------------------------------------------------------------
    for(int q=0; q<32; q++){TT[q]=-1.;   if(pass[q]==1){Fns=tick[q]-9;
	for(int i=0; i<10; i++){txF[i]=tx[Fns+i]; tyF[i]=tyN[q][Fns+i];}  TGraph *grHt = new TGraph(10,txF,tyF);
	TF1 *fit = new TF1("fit","[2]*exp(-TMath::Power(([0]-x)/[1],4))",txF[0], txF[9]);
	fit->SetParameters(tx[tick[q]],2,1);  grHt->Fit("fit","Q","",txF[0],txF[9]);
	tca=fit->GetParameter(0);  tcb=fit->GetParameter(1);  tcc=fit->GetParameter(2);
	TT[q]=(tca-abs(tcb*TMath::Power(-log(0.5/tcc),0.25)))/0.064;
	H_time->Fill(TT[q]);
      }}
    //-------------------------------------------------------------------------------------------------------
    
    _tickSum=tickS; _maxSum=maxS; _timeSum=TTS;
    for(int q=0; q<32; q++){_tickP[q]=tick[q]; _maxP[q]=max[q]; _timeP[q]=TT[q];}
    
    /*
    std::cout<<"tickSum: "<<tickSum<<" maxSum: "<<maxSum<<" timeSum: "<<timeSum<<std::endl;
    std::cout<<"tickP[0]: "<<tickP[0]<<" maxP[0]: "<<maxP[0]<<" timeP[0]: "<<timeP[0]<<std::endl;
    */
    
  }
  
  
  void NeutrinoTiming::getBeamWF(
				 std::vector<double> *wf_w_03
				   //art::Handle< std::vector< raw::OpDetWaveform > > wf_handle_beam
				   )
  {

    
    //=======================================================================================================
    //=======================================================================================================

    //std::cout<<"get Beam T0"<<std::endl;

    //TMultiGraph *mg0=new TMultiGraph();
    double BBw,BBwQ,BBErr,BBmax,wx[1500],wy[1500],wxh[150],wyh[150],wxhE[150],wyhE[150];
    int Nmc;  double XNmax, pca,pcb,pcc, TT = -99999.0;
    int Nmax, Nrt; double XNmc;
    //=======================================================================================================
    BBw=0.; BBwQ=0.; BBmax=0.; BBErr=0.;
    for(size_t j=0; j<wf_w_03->size(); ++j){wx[j]=j; wy[j]=wf_w_03->at(j);
    if(wx[j]<300){BBw=BBw+wy[j]; BBwQ=BBwQ+TMath::Power(wy[j],2);}}
    BBErr=sqrt((BBwQ/300.0-TMath::Power(BBw/300.0 ,2))/300.0);
    for(size_t j=0; j<wf_w_03->size(); ++j){wy[j]=wy[j]-BBw/300.0; if(BBmax<wy[j]){BBmax=wy[j];}}
    H_maxH->Fill(BBmax);
     if(BBmax>2000 && BBmax<2100){
        for(size_t j=0; j<wf_w_03->size(); ++j){wy[j]=wy[j]/BBmax;}
        //-------------------------------------------------------------------------------------------------------
        for(int i=0; i<150; i++){wxh[i]=wx[i+300]; wyh[i]=wy[i+300]; wxhE[i]=0.; wyhE[i]=BBErr*0.0;}
        TGraphErrors *gr0 = new TGraphErrors(150,wxh,wyh,wxhE,wyhE);
        gr0->SetMarkerStyle(24); gr0->SetLineColor(4); gr0->SetMarkerColor(4);  //if(l<grlimit2){mg0->Add(gr0);}
        //-------------------------------------------------------------------------------------------------------
        Nmc=0; XNmax=0.; Nrt=0; XNmc=0.; Nmax=0;
        for(int i=0; i<80; i++){if(wyh[i-1]<0.5 && wyh[i]>=0.5){XNmc=wxh[i]; Nmc=i;}}
        H_tickHalf->Fill(XNmc);
        for(int i=Nmc; i<Nmc+5; i++){if(wyh[i-1]<1 && wyh[i]==1){Nmax=i; XNmax=wxh[i];}}
        Nrt=Nmax-Nmc;
        H_MaxDTick->Fill(Nrt);
        //-------------------------------------------------------------------------------------------------------
        TF1 *fit = new TF1("fit","[2]*exp(-TMath::Power((x-[0])/[1],4))",XNmax-6, XNmax+0);
        fit->SetParameters(XNmax,2,1);     gr0->Fit("fit","Q","",XNmax-6, XNmax+0);
        pca=fit->GetParameter(0);  pcb=fit->GetParameter(1);  pcc=fit->GetParameter(2);
        TT=(pca-abs(pcb*TMath::Power(-log(0.5/pcc),0.25)))/0.064;
        H_t0_Beam->Fill(TT);

    }

     _BeamT0 = TT;
     //std::cout<<"Beam T0: "<<TT<<std::endl;

     //return TT;
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
    double calib[32];  calib[30]=23.1762; calib[31]=22.7841;
    calib[0]=22.8508;  calib[1]=23.4035;  calib[2]=23.8295;  calib[3]=21.9985;  calib[4]=22.6061;  calib[5]=27.9068;
    calib[6]=23.4256;  calib[7]=26.3513;  calib[8]=24.7808;  calib[9]=24.5895;  calib[10]=25.1391; calib[11]=25.4291;
    calib[12]=24.1004; calib[13]=27.3829; calib[14]=25.4302; calib[15]=25.8601; calib[16]=24.7963; calib[17]=22;
    calib[18]=25.7944; calib[19]=23.2512; calib[20]=25.3285; calib[21]=25.5014; calib[22]=23.5923; calib[23]=23.139;
    calib[24]=25.4383; calib[25]=21.9842; calib[26]=22.5661; calib[27]=24.2682; calib[28]=20.7202; calib[29]=22.8939;
    //double offset[32]={1.03002, -5.18104, -2.11164, -5.99395, -1.25798, 0.633079, 2.87666, 2.21969, 0.885092, 2.35423,
    //  -1.63039, -1.83775, -0.859883, 3.4741, 1.84833, 1.58233, -2.71783, 0, 3.18776, 0.982666, 0.728438, 0.280592, -5.27068,
    //  -3.27857, -1.41196, 1.59643, 1.41425, -1.62682, -2.55772, 1.49136, -0.522791, 0.974533};
    
    
    //================================================================================================================
    double gap=18.936;
    double MaxLim=2.5;
    std::vector<int> N_pmt; //double ccnd1, ccnd2;
    double PhTot, RWM_T, nuToF, DPh,DLh, tPhelp,tp;
    //double Med_TT0,Med_TT1,Med_TT2;
    double Med_TT3=-9999.;
    double TT_merged = -9999.;
    //===================================================================================================================
    //===================================================================================================================
    PhTot=0.;
    N_pmt.clear();
    for(int q=0; q<32; q++){_maxP[q]=_maxP[q]/calib[q];
      if((_maxP[q]>MaxLim && q!=17 && q!=28) && (_timeP[q]>3000.0 && _timeP[q]<5000.0)){N_pmt.push_back(q); PhTot=PhTot+_maxP[q];}}
    //--------------------------------------------------------------------------------------------------------------------
    if(N_pmt.size()>2){
      RWM_T=_BeamT0;///0.064;
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
      for(uint i=0; i<N_pmt.size(); i++){
	//ccnd1= timeProp[i]     *(0.529594)- (7.13804);
	//ccnd2= max[N_pmt.at(i)]*(0.068752)- (2.32023);
	
	//just the time of flight and time prop corrections
	TT3_array[i]=(_timeP[N_pmt.at(i)]/*/0.064*/)-RWM_T+5312.0-nuToF-timeProp[i];
	//all the corrections
	//TT3_array[i]=(time[N_pmt.at(i)]/*/0.064*/)-RWM_T+5312.0-nuToF-timeProp[i]  -offset[N_pmt.at(i)] +ccnd1+ccnd2;
      }
      Med_TT3=TMath::Median((Long64_t)N_pmt.size(),TT3_array);
    }
    _evtTimeNS = Med_TT3;
    std::cout<<"evtTimeNS: "<< Med_TT3<<std::endl;
    //Merge Peaks
    double Shift=3166.1216;
    double TThelp=Med_TT3-Shift+gap*0.5;
    if(TThelp>=0 && TThelp<gap*81.0){TT_merged=(TThelp-(int((TThelp)/gap))*gap)-gap*0.5;}
    else {TT_merged=-9999;}
    
    _evtDeltaTimeNS = TT_merged;
    
    std::cout<<"evtDeltaTimeNS: "<< TT_merged<<std::endl;
    
  }


  
  DEFINE_ART_CLASS_TOOL(NeutrinoTiming)
} // namespace analysis

#endif
