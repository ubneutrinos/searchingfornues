#ifndef ANALYSIS_PMTNOISE_CXX
#define ANALYSIS_PMTNOISE_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Cluster.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"

namespace analysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       PMTNoise
    // File:        PMTNoise.cc
    //
    //             PMT Noise
    //
    // Configuration parameters:
    //
    // TBD
    //
    // Created by S. Berkman (sberkman@fnal.gov) on 03/19/2019
    //
    ////////////////////////////////////////////////////////////////////////

  class PMTNoise : public AnalysisToolBase {
    
  public:
    
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    PMTNoise(const fhicl::ParameterSet& pset);
    
    /**
     *  @brief  Destructor
     */
    ~PMTNoise(){ };
    
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
    
    //  private:

    //    art::InputTag HitInputTag("gaushit");
    //    art::InputTag OpHitCosmicInputTag("ophitCosmic");
    /*
    //optical information                                                           
    int totophits;
    std::vector<int> opchannel;
    std::vector<double> peaktime;
    std::vector<double> width;
    std::vector<double> area;
    std::vector<double> amplitude;
    std::vector<double> pe;
    std::vector<double> ophit_wire;
    std::vector<double> ophit_wire_cm;
    std::vector<double> ophit_time;
    std::vector<int> ophit_plane;

    //hits                                                                          
    std::vector<float> hit_peaktime;
    std::vector<double> hit_wire;
    std::vector<double> hit_wire_cm;
    std::vector<int> hit_plane;

    //hits in neutrino slice                                                        
    std::vector<float> slicehit_peaktime;
    std::vector<double> slicehit_wire;
    std::vector<int> slicehit_plane;
    */
    //noise frac in neutrino slice on plane 1
  private:    
    float frac_slnoise_pl1;

    int nflag_pl1;
    int nnoise_pl1;
    int nslhits_pl1;
    int nslnoise_pl1;
    int nhits_pl1;    
  };
  
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///    PMTNoise(const fhicl::ParameterSet& pset);

  PMTNoise::PMTNoise(const fhicl::ParameterSet& p)
  {
    /*
    fCRTVetoproducer = p.get< art::InputTag > ("CRTVetoproducer",""); // default is no CRT veto
    fCLSproducer = p.get< art::InputTag > ("CLSproducer");
    fMCTproducer = p.get< art::InputTag > ("MCTproducer");
    fBacktrackTag = p.get< art::InputTag > ("BacktrackTag");
    // kinematic thresholds for defining signal
    fProtonThreshold = p.get< float > ("ProtonThreshold");
    */
    /*
    fFlagNoiseTimeLow = p.get< float > ("FlagNoiseTimeLow");
    fFlagNoiseTimeHigh = p.get< float > ("FlagNoiseTimeHigh");
    fFlagNoiseWireLow = p.get< float > ("FlagNoiseWireLow");
    fFlagNoiseWireHigh = p.get< float > ("FlagNoiseWireHigh");

    fNoiseTimeLow = p.get< float > ("NoiseTimeLow");
    fNoiseTimeHigh = p.get< float > ("NoiseTimeHigh");
    fNoiseWireHigh =  p.get< float> ("NoiseWireHigh");
    fNoiseWireLow =  p.get< float > ("fNoiseWireLow");

    fHighPEThresh = p.get< float > ("fHighPEThresh");
    */
    

  
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void PMTNoise::configure(fhicl::ParameterSet const & p)
  {
  }
  
  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void PMTNoise::analyzeEvent(art::Event const& e, bool fData)
  {
    /*
    _evt = e.event();
    _sub = e.subRun();
    _run = e.run();

    return;
    */
  }



  //  void PMTNoise::analyzeSlice(art::Event const& e){
  void PMTNoise::analyzeSlice(art::Event const& e, std::vector<ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected){

      //      ProxySliceColl_t const& slice_proxy = proxy::getCollection<std::vector<recob::Slice> >(e, fCLSproducer, proxy::withAssociated<recob::Hit>(HitInputTag));


    art::InputTag HitInputTag("gaushit");
    art::ValidHandle<std::vector<recob::Hit> > inputHits = e.getValidHandle<std::vector<recob::Hit> >(HitInputTag);
    art::InputTag OpHitCosmicInputTag("ophitCosmic");
    art::ValidHandle<std::vector<recob::OpHit> > inputOpHits = e.getValidHandle<std::vector<recob::OpHit> >(OpHitCosmicInputTag);
    art::ServiceHandle<geo::Geometry> geom;

    art::InputTag PfInputTag("pandora");
    art::ValidHandle<std::vector<recob::PFParticle> > inputPfParticle = e.getValidHandle<std::vector<recob::PFParticle> >(PfInputTag);
    art::ValidHandle<std::vector<recob::Cluster> > inputCluster = e.getValidHandle<std::vector<recob::Cluster> >(PfInputTag);
    auto assocCluster = std::unique_ptr<art::FindManyP<recob::Cluster> >(new art::FindManyP<recob::Cluster>(inputPfParticle, e, PfInputTag));
    auto assocHit = std::unique_ptr<art::FindManyP<recob::Hit> >(new art::FindManyP<recob::Hit>(inputCluster, e, PfInputTag));

    art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(inputPfParticle, e, PfInputTag);

    auto assocSlice = std::unique_ptr<art::FindManyP<recob::Slice> >(new art::FindManyP<recob::Slice>(inputPfParticle, e, PfInputTag));
    art::ValidHandle<std::vector<recob::Slice> > inputSlice = e.getValidHandle<std::vector<recob::Slice> >(PfInputTag);
    auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit> >(new art::FindManyP<recob::Hit>(inputSlice, e, PfInputTag));

    //hit                                                                                                                                              
    std::vector<float> hit_peaktime;
    std::vector<double> hit_wire;
    std::vector<int> hit_plane;

    //hits in neutrino slice                                                                                                                           
    std::vector<float> slicehit_peaktime;
    std::vector<double> slicehit_wire;
    std::vector<int> slicehit_plane;

    //optical                                                                              
    std::vector<int> opchannel;
    std::vector<double> peaktime;
    std::vector<double> pe;
    std::vector< std::vector<double> > ophit_wire;
    std::vector<int> ophit_plane;
    std::vector<double> ophit_time;

    bool ispl1_flag=false;
    bool ispl1_noise=false;
    bool ispl1_slicenoise=false;

    int n_hits_pl1=0;
    int n_flag_pl1=0;
    int n_noise_pl1=0;
    int n_slhits_pl1=0;
    int n_slnoise_pl1=0;


    //loop over optical hits                                                               
    for ( unsigned int ioh=0; ioh<inputOpHits->size(); ioh++) {
      const auto& ophit = inputOpHits->at(ioh);
      std::vector<double> opwiretmp;
      //save if "high" p.e. PMT pulse and within TPC truncated time window
      if(ophit.PE()>450 && ophit.PeakTime()>-400 && ophit.PeakTime()<3200){
	opchannel.push_back(ophit.OpChannel());
	peaktime.push_back(ophit.PeakTime());
	pe.push_back(ophit.PE());
	
	double PMTxyz[3];
	unsigned int opch;
	
	opch=ophit.OpChannel();
	
	geom->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);
	
	for( int ipl=0; ipl<3; ipl++){
	  auto plid=geo::PlaneID(0,0,ipl);
	  auto wire = geom->WireCoordinate(PMTxyz[1],PMTxyz[2],plid);
	  //auto time = PMTxyz[0]; //should be the same +/- ~cm for all PMTs                       
	  //	  std::cout<<"PMT HIT: "<<ophit.PE()<<std::endl;
	  opwiretmp.push_back(wire);
	  //	  ophit_wire.push_back(wire);
	  //ophit_time.push_back(time);
	  //ophit_plane.push_back(ipl);
	}
	ophit_wire.push_back(opwiretmp);
      }
    }//finish ophit loop


    //only look at TPC if there are high p.e. pmt pulses
    if(pe.size()==0){
      frac_slnoise_pl1=0;
    }
   
    if(pe.size()!=0){
      //get TPC hits in the event
      for ( unsigned int ih=0; ih<inputHits->size(); ih++) {
	
	const auto& tpchit = inputHits->at(ih);
	auto hit_w=tpchit.WireID().Wire;
	int hit_pl=tpchit.WireID().Plane;
	
	//hit_peaktime.push_back(tpchit.PeakTime());
	//hit_wire.push_back(tpchit.WireID().Wire);
	//hit_plane.push_back(tpchit.WireID().Plane);
	
	//convert hit time to coordinates of TPC time, relative to trigger
	float hit_time_pmtcoord=(tpchit.PeakTime()/2)-400;
	
	//count hits in boxes on plane 1
	if(hit_pl==1){
	  n_hits_pl1=n_hits_pl1+1;  
	  //loop over high p.e. PMT hits
	  for(unsigned int i=0; i<pe.size(); i++){
	    //cut on position in flag region
	    if(ophit_wire[i][hit_pl]-hit_w>-15 && ophit_wire[i][hit_pl]-hit_w<15 && ispl1_flag==false){
	      //cut on time
	      if(peaktime[i]-hit_time_pmtcoord>-1 && peaktime[i]-hit_time_pmtcoord<0.5){
		n_flag_pl1=n_flag_pl1+1;//count number of hits in flag region
		ispl1_flag=true;
	      }
	    }
	    //cut on noise box
	    if(ophit_wire[i][hit_pl]-hit_w>-50 && ophit_wire[i][hit_pl]-hit_w<50 && ispl1_noise==false){
	      //cut on time                                                                                            
	      if(peaktime[i]-hit_time_pmtcoord>-2 && peaktime[i]-hit_time_pmtcoord<2){
		n_noise_pl1=n_noise_pl1+1;//count number of hits in flag region                                          
		ispl1_noise=true;
	      }
	    }
	  }// PMT loop
	}//plane 1 only
	ispl1_flag=false;
	ispl1_noise=false; 
      }// hit loop 
     

      
      for (unsigned int inpf=0; inpf<inputPfParticle->size(); ++inpf) {
	art::Ptr<recob::PFParticle> npfp(inputPfParticle,inpf);
	bool isTheNeutrino = false;
	
	const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(inpf));
	if (!pfParticleMetadataList.empty()) {
	  for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j) {
	    const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
	    
	    const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
	    for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	      if (it->first=="IsNeutrino" && it->second==1) isTheNeutrino = true;
	    }
	  }
	}
	if (npfp->IsPrimary()==false) continue;
	//  recob::Slices are unique for the neutrino slice and for clear cosmics. Ambiguous cosmics may have >1 primary PFParticle per slice 
	//(i.e. the slice is not unique in terms of primary PFParticles)                           
	auto slices = assocSlice->at(npfp.key());
	if (slices.size()!=1) {
	  std::cout << "WRONG!!! n slices = " << slices.size() << std::endl;
	}
	if(isTheNeutrino==true){
	  auto slicehit = assocSliceHit->at(slices[0].key());
          for ( unsigned int isl=0; isl<slicehit.size(); isl++) {
	    const auto& slhit = slicehit.at(isl);
	    //now see if there is noise in the neutrino slice
	    // slicehit_peaktime.push_back(slhit->PeakTime());
	    //slicehit_wire.push_back(slhit->WireID().Wire);
	    //slicehit_plane.push_back(slhit->WireID().Plane);
	    double slhit_time_pmtcoord=(slhit->PeakTime()/2)-400;
	    int slhit_pl=slhit->WireID().Plane;
	    double slhit_wire=slhit->WireID().Wire;
	    double frac_flag_to_noise=(double)n_flag_pl1/(double)n_noise_pl1;
	    if(slhit_pl==1){
	      //count hits in slice on plane 1
	      n_slhits_pl1=n_slhits_pl1+1;
	      //loop over high p.e. PMT hits                                                                                      
	      for(unsigned int i=0; i<pe.size(); i++){
		//cut on noise box                                                                                                             
		//at least 2% of the noise hits in the event need to be in the flag box to count slice noise
		if(ophit_wire[i][slhit_pl]-slhit_wire>-50 && ophit_wire[i][slhit_pl]-slhit_wire<50 && ispl1_slicenoise==false && frac_flag_to_noise>0.02){
		  //cut on time                                                                                                                
		  if(peaktime[i]-slhit_time_pmtcoord>-2 && peaktime[i]-slhit_time_pmtcoord<2){
		    n_slnoise_pl1=n_slnoise_pl1+1;//count number of hits in noise region                     
		    ispl1_slicenoise=true;
		  }// if hit is in PMT time noise box
		}// if hit is in the noise box
	      }//for each large PE PMT Hit
	    }// if pl1
	    ispl1_slicenoise=false;
	  }//for each hit in the slice
	}//istheneutrinoslice
      }// slice-pfp assoc loop

      std::cout<<"slicenoise: "<<n_slnoise_pl1<<" nslhits: "<<n_slhits_pl1<<std::endl;
      
      frac_slnoise_pl1=(float)n_slnoise_pl1/(float)n_slhits_pl1;
      
      std::cout<<"SLICE NOISE FRACTION: "<<frac_slnoise_pl1<<std::endl;
      
    }// only look through TPC noise if there is one high PE pulse
    
    nflag_pl1=n_flag_pl1;
    nnoise_pl1=n_noise_pl1;
    nslhits_pl1=n_slhits_pl1;
    nslnoise_pl1=n_slnoise_pl1;
    nhits_pl1=n_hits_pl1;

    n_hits_pl1=0;
    n_flag_pl1=0;
    n_noise_pl1=0;
    n_slhits_pl1=0;
    n_slnoise_pl1=0;

  }//end of function

  void PMTNoise::setBranches(TTree* _tree) 
  {
    /*
    _tree->Branch("opchannel",&opchannel);
    _tree->Branch("peaktime",&peaktime);
    _tree->Branch("width",&width);
    _tree->Branch("area",&area);
    _tree->Branch("amplitude",&amplitude);
    _tree->Branch("pe",&pe);
    _tree->Branch("ophit_wire",&ophit_wire);
    _tree->Branch("ophit_wire_cm",&ophit_wire_cm);
    _tree->Branch("ophit_plane",&ophit_plane);
    _tree->Branch("ophit_time",&ophit_time);

    _tree->Branch("hit_peaktime",&hit_peaktime);
    _tree->Branch("hit_wire",&hit_wire);
    _tree->Branch("hit_wire_cm",&hit_wire_cm);
    _tree->Branch("hit_plane",&hit_plane);

    _tree->Branch("slicehit_peaktime",&slicehit_peaktime);
    _tree->Branch("slicehit_wire",&slicehit_wire);
    _tree->Branch("slicehit_plane",&slicehit_plane);
    */

    _tree->Branch("nflag_pl1",&nflag_pl1);
    _tree->Branch("nnoise_pl1",&nnoise_pl1);
    _tree->Branch("nslhits_pl1",&nslhits_pl1);
    _tree->Branch("nslnoise_pl1",&nslnoise_pl1);
    _tree->Branch("nhits_pl1",&nhits_pl1);
    _tree->Branch("frac_slnoise_pl1",&frac_slnoise_pl1);

  }

  void PMTNoise::resetTTree(TTree* _tree)
  {

    /*
    opchannel.clear();
    peaktime.clear();
    width.clear();
    area.clear();
    amplitude.clear();
    pe.clear();
    ophit_wire.clear();
    ophit_wire_cm.clear();
    ophit_plane.clear();
    ophit_time.clear();

    hit_peaktime.clear();
    hit_wire.clear();
    hit_wire_cm.clear();
    hit_plane.clear();

    slicehit_peaktime.clear();
    slicehit_wire.clear();
    slicehit_plane.clear();
    */
    frac_slnoise_pl1=0;
    nflag_pl1=0;
    nnoise_pl1=0;
    nslhits_pl1=0;
    nslnoise_pl1=0;
    nhits_pl1=0;
}

  DEFINE_ART_CLASS_TOOL(PMTNoise)
} // namespace analysis

#endif
