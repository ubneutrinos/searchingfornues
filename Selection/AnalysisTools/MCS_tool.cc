#ifndef ANALYSIS_DEFAULTANALYSIS_CXX
#define ANALYSIS_DEFAULTANALYSIS_CXX

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include "TROOT.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TError.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TStyle.h"
#include "TString.h"
#include "TVector3.h"
#include "TCanvas.h"

#include "AnalysisToolBase.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "lardataobj/RecoBase/SpacePoint.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/Containment.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/ProximityClustering.h"

#include "canvas/Persistency/Common/TriggerResults.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

//#include "/grid/fermiapp/products/larsoft/eigen/v3_3_4a/include/eigen3/Eigen/Dense" //Needed on uboonegpvm
#include <Eigen/Dense>

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"

struct PCAResults {
  TVector3 centroid;
  std::pair<TVector3,TVector3> endPoints;
  float length;
  TVector3 eVals;
  TVector3 eRatios;
  std::vector<TVector3> eVecs;
  int n_pts, pdg;
};

struct SpcPoint{
  Double_t x, y, z, dvtx, dshrvtx;
  Int_t n_pts, pdg;
};

struct PCAASResult{
  float AS1C, AS2C, AS3C;
};

typedef std::vector<SpcPoint> SpcPointCloud;

// Driver function to sort the vector elements 
// by second element of pairs 
/*
bool sortbysec(const std::vector<SpcPoint> &a, 
               const std::vector<SpcPoint> &b) 
{ 
    return (a.dvtx < b.dvtx); 
} 
*/
using namespace std;

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       MCS
// File:        MCS.cc
//
//              A basic analysis example
//
// Configuration parameters:
//
// TBD
//
// Created by Ivan Caro Terrazas (icaro@colostate.edu) on 12/06/2019
//
////////////////////////////////////////////////////////////////////////

class MCS : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  MCS(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~MCS(){};

  // provide for initialization
  void configure(fhicl::ParameterSet const &pset);

  /**
     * @brief Analysis function
     */
  void analyzeEvent(art::Event const &e, bool fData) override;

  /**
     * @brief Analyze slice
     */
  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

  /**
     * @brief Save truth info for event associated to neutrino
     */
  void SaveTruth(art::Event const &e);

  /**
     * @brief set branches for TTree
     */
  void setBranches(TTree *_tree) override;

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree) override;

  //void LoadPointCloud(PCAPointCloud &points, const WCClst &ord_trk);

  PCAResults DoPCA(const SpcPointCloud &points);
  float CalcMed(vector<float> scores);
  double pythagoras(double x1, double x2, double y1, double y2, double z1, double z2);
  TVector3 CrossProd(TVector3 u, TVector3 v);
  double linepointdist(TVector3 lStart, TVector3 lEnd, TVector3 point);
  double dotProd(TVector3 v1, TVector3 v2);
  double DoDeltaRMS(SpcPointCloud cloud, TVector3 StartP, TVector3 EndP);
  float DoMCSCalc(SpcPointCloud cloud, TVector3 StartP, TVector3 EndP);
  float DoMomASCalc(SpcPointCloud cloud, TVector3 vtx);
  PCAASResult DoPCAASCalc(SpcPointCloud cloud, float nuvX, float nuvY, float nuvZ);
  float DoCylFracCalc(SpcPointCloud cloud, TVector3 StartP, TVector3 EndP, float rad);
  float DoDeltaMed(SpcPointCloud cloud, TVector3 LStart, TVector3 LEnd);
  void Project3Dto2D(const TVector3& pt3d, const int& pl,
         const float& wire2cm, const float& time2cm,
         float& wirecm, float& timecm);

private:
  art::InputTag fHproducer;
  art::InputTag fCLSproducer;
  art::InputTag fSLCproducer; // slice associated to PFP
  art::InputTag fMCRproducer;
  art::InputTag fPFPproducer;
  art::InputTag fSpacePointproducer;

  float fPCADistThreshold;
  float fTrkShrscore;  /**< Threshold on the Pandora track score (default 0.5) */

  //std::string fPFPproducer, fSpacePointproducer;

  //float _mcsrms;
  int slpdg;
  int slnhits;  // number of hits in slice
  std::vector<int> pfpdg;          // PDG code of pfp in slice

  unsigned int _n_pfps_fpca;
  float _reco_nu_vtx_x, _reco_nu_vtx_y, _reco_nu_vtx_z;
  float _reco_nu_vtx_sce_x, _reco_nu_vtx_sce_y, _reco_nu_vtx_sce_z;

  float _reco_nu_vtx_sce_time_plane0, _reco_nu_vtx_sce_wire_plane0;
  float _reco_nu_vtx_sce_time_plane1, _reco_nu_vtx_sce_wire_plane1;
  float _reco_nu_vtx_sce_time_plane2, _reco_nu_vtx_sce_wire_plane2;


  std::vector<float> _X_SpcPts_v;
  std::vector<float> _Y_SpcPts_v;
  std::vector<float> _Z_SpcPts_v;

  std::vector<float> _trkshrscore_v;
  unsigned int _shr_hits_max;        /**< Number of hits of the leading shower */
  size_t _shr_pfp_id; /**< Index of the leading shower in the PFParticle vector */
  unsigned int _hits_outfv;      /**< Number of hits of PFParticles outside the fiducial volume */
  int _n_showers_contained;
  double _shr_score_MCStool;
  int _n_shrSpcPts; /**< Number of spacepoints of the leading shower */

  float fFidvolZstart; /**< Fiducial volume distance from the start of the TPC on the z axis (default 10 cm) */
  float fFidvolZend;   /**< Fiducial volume distance from the end of the TPC on the z axis (default 50 cm) */
  float fFidvolYstart; /**< Fiducial volume distance from the bottom of the TPC on the y axis (default 15 cm) */
  float fFidvolYend;   /**< Fiducial volume distance from the top of the TPC on the y axis (default 15 cm) */
  float fFidvolXstart; /**< Fiducial volume distance from the start of the TPC on the x axis (default 10 cm) */
  float fFidvolXend;   /**< Fiducial volume distance from the end of the TPC on the x axis (default 10 cm) */

  TVector3 _shrvtx;
  
  float _shrPCALen;
  float _shrPCAAS;

  float _shrPCA_1Cr, _shrPCA_2Cr, _shrPCA_3Cr;
  float _shrPCA_1Ce, _shrPCA_2Ce, _shrPCA_3Ce;
  float _shrMCSMom;
  float _shrPCA1CAS; float _shrPCA2CAS; float _shrPCA3CAS;
  float _shrPCA_1Cr2h; float _shrPCA_2Cr2h; float _shrPCA_3Cr2h;
  float _shrMCSMom2h;

  float _shrPCA_1Cr1h; float _shrPCA_2Cr1h; float _shrPCA_3Cr1h;
  float _shrMCSMom1h;

  float _DeltaMed, _DeltaMed1h, _DeltaMed2h;
  float _DeltaRMS, _DeltaRMS1h, _DeltaRMS2h;
  float _CylFrac, _CylFrac1h, _CylFrac2h;

  float _CylFrac_1cm, _CylFrac1h_1cm, _CylFrac2h_1cm;
  float _CylFrac_2cm, _CylFrac1h_2cm, _CylFrac2h_2cm;
  float _CylFrac_3cm, _CylFrac1h_3cm, _CylFrac2h_3cm;
  float _CylFrac_4cm, _CylFrac1h_4cm, _CylFrac2h_4cm;
  float _CylFrac_5cm, _CylFrac1h_5cm, _CylFrac2h_5cm;

  std::vector<float> _PCAWin_1Cr_5cm, _PCAWin_2Cr_5cm, _PCAWin_3Cr_5cm; 
  std::vector<float> _PCAWin_dist_5cm;
  std::vector<int> _PCAWin_npts_5cm;
  float _shrStart_5cm; float _shrStartMCS_5cm; float _shrMCSAS_5cm; 
  float _shrPCA1CAS_5cm; float _shrPCA2CAS_5cm; float _shrPCA3CAS_5cm;

  
  float _shrPCA1CMed_5cm; /**< Median 1Cr value of all the 5cm windows of the showers*/
  //float _shrPCA2CMed; float _shrPCA3CMed;

  std::vector<float> _PCAWin_1Cr_2_5cm, _PCAWin_2Cr_2_5cm, _PCAWin_3Cr_2_5cm; 
  std::vector<float> _PCAWin_dist_2_5cm;
  std::vector<int> _PCAWin_npts_2_5cm;
  float _shrStart_2_5cm; float _shrStartMCS_2_5cm; float _shrMCSAS_2_5cm; 
  float _shrPCA1CAS_2_5cm; float _shrPCA2CAS_2_5cm; float _shrPCA3CAS_2_5cm;
  float _shrPCA1CMed_2_5cm; /**< Median 1Cr value of all the 5cm windows of the showers*/
  //float _shrPCA2CMed; float _shrPCA3CMed;

};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
MCS::MCS(const fhicl::ParameterSet &p)
{
  fCLSproducer = p.get<art::InputTag>("CLSproducer");
  fHproducer = p.get<art::InputTag>("Hproducer");
  fSLCproducer = p.get<art::InputTag>("SLCproducer");
  fMCRproducer = p.get<art::InputTag>("MCRproducer");
  fSpacePointproducer  = p.get< art::InputTag >("SpacePointproducer");
  fPFPproducer  = p.get<art::InputTag>("PFPproducer");

  fPCADistThreshold = p.get<float>("PCADistThreshold", 10); //5 cm from vertex
  fTrkShrscore = p.get<float>("TrkShrscore", 0.5);

  fFidvolZstart = p.get<float>("FidvolZstart");
  fFidvolZend = p.get<float>("FidvolZend");
  fFidvolYstart = p.get<float>("FidvolYstart");
  fFidvolYend = p.get<float>("FidvolYend");
  fFidvolXstart = p.get<float>("FidvolXstart");
  fFidvolXend = p.get<float>("FidvolXend");

}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void MCS::configure(fhicl::ParameterSet const &p)
{

}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void MCS::analyzeEvent(art::Event const &e, bool fData)
{
  std::cout << "[MCS::analyzeEvent] Run: " << e.run() << ", SubRun: " << e.subRun() << ", Event: " << e.event() << std::endl;
  art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
  //evnhits = inputHits->size();
}

void MCS::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected){
  //std::cout << "IVAN START Event: " << e.event() << ", from run = " << e.run() << std::endl;
/*
  //TEST
  SpcPoint PStart; PStart.x = 0; PStart.y = 0; PStart.z = 0; PStart.dvtx = 0; PStart.dshrvtx = 0; PStart.pdg = 0; PStart.n_pts = 0;
  SpcPoint PEnd; PEnd.x = 1; PEnd.y = 1; PEnd.z = 0; PEnd.dvtx = 0; PEnd.dshrvtx = 0; PEnd.pdg = 0; PEnd.n_pts = 0;
  SpcPoint PTest; PTest.x = 2; PTest.y = 4; PTest.z = 0; PTest.dvtx = 0; PTest.dshrvtx = 0; PTest.pdg = 0; PTest.n_pts = 0;

  TVector3 PStart_v(PStart.x,PStart.y,PStart.z); TVector3 PEnd_v(PEnd.x,PEnd.y,PEnd.z); TVector3 PTest_v(PTest.x,PTest.y,PTest.z);
  TVector3 PTest2_v(0,1,0);
  float TestDist = linepointdist(PStart_v, PEnd_v, PTest_v);
  float TestDist2 = linepointdist(PStart_v, PEnd_v, PTest2_v);

  cout << "IVAN The test distance1 is " << TestDist <<" and test distance2 is " << TestDist2 << endl;
*/
//////////END Test Area//////////
//  ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
//                                                                                        proxy::withAssociated<recob::Hit>(fCLSproducer));
  ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
                                                                                        proxy::withAssociated<recob::Hit>(fCLSproducer));
  art::ServiceHandle<geo::Geometry> geom;
  //auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  // somehow proxies don't work for the slice-hit association, so go back to old assns
  art::ValidHandle<std::vector<recob::Slice>> inputSlice = e.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
  auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSlice, e, fSLCproducer));

  // Build larpandora info:
  lar_pandora::LArPandoraHelper larpandora;
  lar_pandora::PFParticleVector pfparticles;
  lar_pandora::PFParticleMap particleMap;
  larpandora.CollectPFParticles(e, "pandora", pfparticles);
  larpandora.BuildPFParticleMap(pfparticles, particleMap);

  // load backtrack information
  std::vector<searchingfornues::BtPart> btparts_v;
  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
  //float _wire2cm = geom->WirePitch(0,0,0);
  //float _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );


    // START checking if vertex is in the fiducial volume
    double nu_vtx[3] = {};
    TVector3 nuvtx;
    for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
    {
        auto PDG = fabs(slice_pfp_v[i_pfp]->PdgCode());

        if (PDG == 12 || PDG == 14)
        {

            auto vtx = slice_pfp_v[i_pfp].get<recob::Vertex>();
            if (vtx.size() != 1)
            {
                std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
            }
            else
            {
                vtx.at(0)->XYZ(nu_vtx);
                if (!searchingfornues::isFiducial(nu_vtx,fFidvolXstart,fFidvolYstart,fFidvolZstart,fFidvolXend,fFidvolYend,fFidvolZend))
                {
                    return;
                }
                nuvtx.SetXYZ(nu_vtx[0], nu_vtx[1], nu_vtx[2]);
            }

            break;
        }
    }// for all PFParticles
    // DONE checking if vertex is in the fiducial volume

  //Find leading shower id _shr_pfp_id
  for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++){
    auto const &pfp_pxy = slice_pfp_v.at(i_pfp);
    auto PDG = fabs(pfp_pxy->PdgCode());
    if (PDG == 12 || PDG == 14) continue;

    auto trkshrscore = searchingfornues::GetTrackShowerScore(pfp_pxy);
    if (trkshrscore < fTrkShrscore){
      unsigned int shr_hits = 0;
      for (const auto &shr : pfp_pxy.get<recob::Shower>()){
        double shr_vertex[3] = {shr->ShowerStart().X(), shr->ShowerStart().Y(), shr->ShowerStart().Z()};
        auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();
        if (!searchingfornues::isFiducial(shr_vertex,fFidvolXstart,fFidvolYstart,fFidvolZstart,fFidvolXend,fFidvolYend,fFidvolZend)){
          //int ipfp = i_pfp;
          for (auto ass_clus : clus_pxy_v){
            const auto &clus = clus_proxy[ass_clus.key()];
            auto clus_hit_v = clus.get<recob::Hit>();
            _hits_outfv += clus_hit_v.size();
          }
          continue;
        }
        _n_showers_contained++;
        for (auto ass_clus : clus_pxy_v){
          const auto &clus = clus_proxy[ass_clus.key()];
          auto clus_hit_v = clus.get<recob::Hit>();
          shr_hits += clus_hit_v.size();
        }
        if (shr_hits > _shr_hits_max){
          _shrvtx = TVector3(shr_vertex[0],shr_vertex[1],shr_vertex[2]);
          _shr_pfp_id = i_pfp;
          _shr_hits_max = shr_hits;
          _shr_score_MCStool = trkshrscore;
        }
      }
    }
  }

  unsigned ipfp = 0;
  
  for (auto pfp : slice_pfp_v){

    if (pfp->IsPrimary()){ // Get Neutrino Vertex //
      
      slpdg = pfp->PdgCode();
      auto slice_pxy_v = pfp.get<recob::Slice>();
      if (slice_pxy_v.size() != 1)
      {
        std::cout << "WRONG!!! n slices = " << slice_pxy_v.size() << " " << __FILE__ << " " << __LINE__ << std::endl;
        return;
      }
      auto slicehits = assocSliceHit->at(slice_pxy_v[0].key());
      slnhits = slicehits.size();

      // grab vertex
      double xyz[3] = {};

      auto vtx = pfp.get<recob::Vertex>();
      if (vtx.size() == 1)
      {
        // save vertex to array
        vtx.at(0)->XYZ(xyz);
        auto nuvtx = TVector3(xyz[0], xyz[1], xyz[2]);

        _reco_nu_vtx_x = nuvtx.X();
        _reco_nu_vtx_y = nuvtx.Y();
        _reco_nu_vtx_z = nuvtx.Z();

        float _reco_nu_vtx_sce[3];
        searchingfornues::ApplySCECorrectionXYZ(_reco_nu_vtx_x, _reco_nu_vtx_y, _reco_nu_vtx_z, _reco_nu_vtx_sce);
        _reco_nu_vtx_sce_x = _reco_nu_vtx_sce[0];
        _reco_nu_vtx_sce_y = _reco_nu_vtx_sce[1];
        _reco_nu_vtx_sce_z = _reco_nu_vtx_sce[2];
      }else{
        std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
      }
      
      continue;
    } // if neutrino PFParticle -> Found the Neutrino Vertex

    TVector3 vtx3d(_reco_nu_vtx_sce_x, _reco_nu_vtx_sce_y, _reco_nu_vtx_sce_z);
    ++ipfp;
    if (_shr_pfp_id == 0) continue;
    if (ipfp != _shr_pfp_id) continue;
    auto trkshrscore = searchingfornues::GetTrackShowerScore(pfp);
    _trkshrscore_v.push_back(trkshrscore);
    auto spcpts_v = pfp.get<recob::SpacePoint>();
    
    SpcPointCloud PFPcloud; //cloud of PFP particle
    //Store PFP particle space points
    for (size_t s=0; s < spcpts_v.size(); s++) {
      SpcPoint point;
      auto sps = spcpts_v[s];
      auto spsx = sps->XYZ()[0]; auto spsy = sps->XYZ()[1]; auto spsz = sps->XYZ()[2]; 
      _X_SpcPts_v.push_back(spsx); _Y_SpcPts_v.push_back(spsy); _Z_SpcPts_v.push_back(spsz);
      double dist = pythagoras(spsx, _reco_nu_vtx_sce_x, spsy, _reco_nu_vtx_sce_y, spsz, _reco_nu_vtx_sce_z);
      double shrdist = pythagoras(spsx, _shrvtx(0), spsy, _shrvtx(1),spsz, _shrvtx(2));

      point.x = spsx; point.y = spsy; point.z = spsz; point.n_pts = spcpts_v.size(); point.pdg = pfp->PdgCode(); 
      point.dvtx = dist; point.dshrvtx = shrdist;
      PFPcloud.push_back(point);
    }// for all spacepoints

    //Sort w.r.t. vtx
    sort( PFPcloud.begin( ), PFPcloud.end( ), [ ]( const SpcPoint& lhs, const SpcPoint& rhs ){return lhs.dshrvtx < rhs.dshrvtx;});
    
    _n_shrSpcPts = PFPcloud.size();

    PCAResults PFP_PCAres = DoPCA(PFPcloud);
    TVector3 PFPcentroid = PFP_PCAres.centroid;
    _shrMCSMom = DoMCSCalc(PFPcloud, _shrvtx, PFPcentroid);
    
    float dShrCent = pythagoras(PFPcentroid(0), _shrvtx(0), PFPcentroid(1), _shrvtx(1), PFPcentroid(2), _shrvtx(2));
    SpcPointCloud PFPcloud2h, PFPcloud1h;

    //cout << "IVAN: evt.sub.evt: " << e.run() << "." << e.subRun() << "." << e.event() << ", dShrCent = " << dShrCent << endl;
    for (unsigned i = 0; i < PFPcloud.size(); ++i){
      if(PFPcloud[i].dshrvtx > dShrCent){
        //cout << "IVAN 2h evt.sub.evt: " << e.run() << "." << e.subRun() << "." << e.event() << ",PFPcloud[" << i << "].dshrvtx = " << PFPcloud[i].dshrvtx << endl;
        PFPcloud2h.push_back(PFPcloud[i]);
      }else{
        PFPcloud1h.push_back(PFPcloud[i]);
        //cout << "IVAN 1h evt.sub.evt: " << e.run() << "." << e.subRun() << "." << e.event() << ",PFPcloud[" << i << "].dshrvtx = " << PFPcloud[i].dshrvtx << endl;
      }
    }

    PCAResults PFP_PCAres2h = DoPCA(PFPcloud2h);
    _shrPCA_1Cr2h = PFP_PCAres2h.eRatios(0); _shrPCA_2Cr2h = PFP_PCAres2h.eRatios(1); _shrPCA_3Cr2h = PFP_PCAres2h.eRatios(2);
    _shrMCSMom2h = DoMCSCalc(PFPcloud2h, _shrvtx, PFPcentroid);

    PCAResults PFP_PCAres1h = DoPCA(PFPcloud1h);
    _shrPCA_1Cr1h = PFP_PCAres1h.eRatios(0); _shrPCA_2Cr1h = PFP_PCAres1h.eRatios(1); _shrPCA_3Cr1h = PFP_PCAres1h.eRatios(2);
    _shrMCSMom1h = DoMCSCalc(PFPcloud1h, _shrvtx, PFPcentroid);

    PCAASResult shrPCAAS = DoPCAASCalc(PFPcloud, _shrvtx(0), _shrvtx(1), _shrvtx(2));
    _shrPCA1CAS = shrPCAAS.AS1C; _shrPCA2CAS = shrPCAAS.AS2C; _shrPCA3CAS = shrPCAAS.AS3C;

    TVector3 shr_eRatios(PFP_PCAres.eRatios(0), PFP_PCAres.eRatios(1), PFP_PCAres.eRatios(2));
    _shrPCA_1Cr = shr_eRatios(0); _shrPCA_2Cr = shr_eRatios(1); _shrPCA_3Cr = shr_eRatios(2);
    _shrPCA_1Ce = PFP_PCAres.eVals(0); _shrPCA_2Ce = PFP_PCAres.eVals(1); _shrPCA_3Ce = PFP_PCAres.eVals(2);

    float PFPLen = PFP_PCAres.length; _shrPCALen = PFPLen;
    TVector3 cent1h = PFP_PCAres1h.centroid;
    //TVector3 cent2h = PFP_PCAres2h.centroid;

    //Determine First point of shower
    //TVector3 shrPointA = PFP_PCAres.endPoints.first;
    //TVector3 shrPointB = PFP_PCAres.endPoints.second;
    //double dShrVtxA = pythagoras(shrPointA(0), _shrvtx(0), shrPointA(1), _shrvtx(1), shrPointA(2), _shrvtx(2));
    //double dShrVtxB = pythagoras(shrPointB(0), _shrvtx(0), shrPointB(1), _shrvtx(1), shrPointB(2), _shrvtx(2));
    TVector3 shrFirst(PFPcloud[0].x,PFPcloud[0].y,PFPcloud[0].z);
  
    //Do Delta Median Calculations
    _DeltaMed = DoDeltaMed(PFPcloud, shrFirst, PFPcentroid);
    _DeltaMed1h = DoDeltaMed(PFPcloud1h, shrFirst, PFPcentroid);
    _DeltaMed2h = DoDeltaMed(PFPcloud2h, shrFirst, PFPcentroid);

    _DeltaRMS = DoDeltaRMS(PFPcloud, shrFirst, PFPcentroid);
    _DeltaRMS1h = DoDeltaRMS(PFPcloud1h, shrFirst, PFPcentroid);
    _DeltaRMS2h = DoDeltaRMS(PFPcloud2h, shrFirst, PFPcentroid);

    _CylFrac_1cm = DoCylFracCalc(PFPcloud, shrFirst, PFPcentroid, 1.);
    _CylFrac1h_1cm = DoCylFracCalc(PFPcloud1h, shrFirst, PFPcentroid, 1.);
    _CylFrac2h_1cm = DoCylFracCalc(PFPcloud2h, shrFirst, PFPcentroid, 1.);

    _CylFrac_2cm = DoCylFracCalc(PFPcloud, shrFirst, PFPcentroid, 2.);
    _CylFrac1h_2cm = DoCylFracCalc(PFPcloud1h, shrFirst, PFPcentroid, 2.);
    _CylFrac2h_2cm = DoCylFracCalc(PFPcloud2h, shrFirst, PFPcentroid, 2.);

    _CylFrac_3cm = DoCylFracCalc(PFPcloud, shrFirst, PFPcentroid, 3.);
    _CylFrac1h_3cm = DoCylFracCalc(PFPcloud1h, shrFirst, PFPcentroid, 3.);
    _CylFrac2h_3cm = DoCylFracCalc(PFPcloud2h, shrFirst, PFPcentroid, 3.);

    _CylFrac_4cm = DoCylFracCalc(PFPcloud, shrFirst, PFPcentroid, 4.);
    _CylFrac1h_4cm = DoCylFracCalc(PFPcloud1h, shrFirst, PFPcentroid, 4.);
    _CylFrac2h_4cm = DoCylFracCalc(PFPcloud2h, shrFirst, PFPcentroid, 4.);

    _CylFrac_5cm = DoCylFracCalc(PFPcloud, shrFirst, PFPcentroid, 5.);
    _CylFrac1h_5cm = DoCylFracCalc(PFPcloud1h, shrFirst, PFPcentroid, 5.);
    _CylFrac2h_5cm = DoCylFracCalc(PFPcloud2h, shrFirst, PFPcentroid, 5.);



    
    //Showering Start finder
    double winSize = 5.;// 5 cm
    std::vector<std::vector<unsigned>> windows_v;
    std::vector<unsigned> window;
    double newWinSize = winSize;
    for (unsigned i = 0; i < PFPcloud.size(); ++i){
      if(PFPcloud.at(i).dshrvtx < newWinSize){
        window.push_back(i);
      }else{
        windows_v.push_back(window);
        window.clear();
        newWinSize += winSize;
        window.push_back(i);
      }
    }
    windows_v.push_back(window);
    //double old1Cr = 1;
    SpcPointCloud tempCloud;
    
    //Calculate the window PCA 1C
    for (unsigned i = 0; i < windows_v.size(); i++){
      float lastdist;
      bool smallWindow = false;
      if(windows_v[i].size() <= 5) smallWindow = true; 
      for (unsigned j = 0; j < windows_v[i].size(); ++j){
        unsigned i_win = windows_v[i][j];
        tempCloud.push_back(PFPcloud[i_win]);
        lastdist = PFPcloud[i_win].dshrvtx;
      }
      if(smallWindow && tempCloud.size() < 5) continue;
      PCAResults tempPCAResults = DoPCA(tempCloud);
      int npts_5cm = tempCloud.size();
      tempCloud.clear();
      float new1Cr = tempPCAResults.eRatios(0);
      float new2Cr = tempPCAResults.eRatios(1);
      float new3Cr = tempPCAResults.eRatios(2);
      if(isnan(new1Cr)) new1Cr = -9999;
      if(isnan(new2Cr)) new2Cr = -9999;
      if(isnan(new3Cr)) new3Cr = -9999;
      _PCAWin_1Cr_5cm.push_back(new1Cr);
      _PCAWin_2Cr_5cm.push_back(new2Cr);
      _PCAWin_3Cr_5cm.push_back(new3Cr);
      _PCAWin_dist_5cm.push_back(lastdist);
      _PCAWin_npts_5cm.push_back(npts_5cm);
    }
    _shrPCA1CMed_5cm = CalcMed(_PCAWin_1Cr_5cm); //

    //Find Showering start
    //cout << "IVAN: evt.sub.evt: " << e.run() << "." << e.subRun() << "." << e.event() << endl;
    for (unsigned i = 0; i < _PCAWin_1Cr_5cm.size(); ++i){
      //cout << "IVAN: evt.sub.evt: " << e.run() << "." << e.subRun() << "." << e.event() << "_PCAWin_1Cr_5cm[" << i << "] = " << _PCAWin_1Cr_5cm[i] << ", _shrStart_5cm = " << _PCAWin_dist_5cm[i] << endl;
      if (_PCAWin_1Cr_5cm[i] < 0.9 && _PCAWin_1Cr_5cm[i] >= 0 && i == 0){
        _shrStart_5cm = 0;
       break;
      }else if(_PCAWin_1Cr_5cm[i] < 0.9 && _PCAWin_1Cr_5cm[i] >= 0){
        _shrStart_5cm = _PCAWin_dist_5cm[i];
        break;
      }
    }
    //cout << "IVAN: evt.sub.evt: "<< e.run() << "." << e.subRun() << "." << e.event() << ", _shrStart_5cm chosen = " << _shrStart_5cm << endl;
    //Calculate MCS up to showering start
    SpcPointCloud shrStartcloud_5cm;
    for (unsigned i = 0; i < PFPcloud.size(); ++i){
      if (PFPcloud[i].dshrvtx < _shrStart_5cm){
        shrStartcloud_5cm.push_back(PFPcloud[i]);
      }
    }
    _shrStartMCS_5cm = DoMCSCalc(shrStartcloud_5cm, _shrvtx, PFPcentroid);

    //Calculate the MCS AS up to the showering start
    _shrMCSAS_5cm = DoMomASCalc(shrStartcloud_5cm, _shrvtx);
    PCAASResult PCAAS_5cm = DoPCAASCalc(shrStartcloud_5cm, _shrvtx(0), _shrvtx(1), _shrvtx(2));
    _shrPCA1CAS_5cm = PCAAS_5cm.AS1C; _shrPCA2CAS_5cm = PCAAS_5cm.AS2C; _shrPCA3CAS_5cm = PCAAS_5cm.AS3C;

    if(isnan(_shrStartMCS_5cm) || isinf(_shrStartMCS_5cm)) _shrStartMCS_5cm = -9999;
    if(isnan(_shrMCSAS_5cm) || isinf(_shrMCSAS_5cm)) _shrMCSAS_5cm = -9999;
    if(isnan(_shrPCA1CAS_5cm) || isinf(_shrPCA1CAS_5cm)) _shrPCA1CAS_5cm = -9999;
    if(isnan(_shrPCA2CAS_5cm) || isinf(_shrPCA2CAS_5cm)) _shrPCA2CAS_5cm = -9999;
    if(isnan(_shrPCA3CAS_5cm) || isinf(_shrPCA3CAS_5cm)) _shrPCA3CAS_5cm = -9999;

    //2.5 cm study///
    double winSize_2_5 = 2.5;// 5 cm
    std::vector<std::vector<unsigned>> windows_2_5_v;
    std::vector<unsigned> window_2_5;
    double newWinSize_2_5 = winSize_2_5;
    for (unsigned i = 0; i < PFPcloud.size(); ++i){
      if(PFPcloud.at(i).dshrvtx < newWinSize_2_5){
        window_2_5.push_back(i);
      }else{
        windows_2_5_v.push_back(window_2_5);
        window_2_5.clear();
        newWinSize_2_5 += winSize_2_5;
        window_2_5.push_back(i);
      }
    }
    windows_2_5_v.push_back(window_2_5);
    //double old1Cr = 1;
    SpcPointCloud tempCloud_2_5;
    
    //Calculate the window PCA 1C
    for (unsigned i = 0; i < windows_2_5_v.size(); i++){
      float lastdist_2_5;
      bool smallWindow_2_5 = false;
      if(windows_2_5_v[i].size() <= 5) smallWindow_2_5 = true; 
      for (unsigned j = 0; j < windows_2_5_v[i].size(); ++j){
        unsigned i_win = windows_2_5_v[i][j];
        tempCloud_2_5.push_back(PFPcloud[i_win]);
        lastdist_2_5 = PFPcloud[i_win].dshrvtx;
      }
      if(smallWindow_2_5 && tempCloud_2_5.size() < 5) continue;
      PCAResults tempPCAResults_2_5 = DoPCA(tempCloud_2_5);
      int npts_2_5cm = tempCloud_2_5.size();
      tempCloud_2_5.clear();
      float new1Cr_2_5 = tempPCAResults_2_5.eRatios(0);
      float new2Cr_2_5 = tempPCAResults_2_5.eRatios(1);
      float new3Cr_2_5 = tempPCAResults_2_5.eRatios(2);
      if(isnan(new1Cr_2_5)) new1Cr_2_5 = -9999;
      if(isnan(new2Cr_2_5)) new2Cr_2_5 = -9999;
      if(isnan(new3Cr_2_5)) new3Cr_2_5 = -9999;
      _PCAWin_1Cr_2_5cm.push_back(new1Cr_2_5);
      _PCAWin_2Cr_2_5cm.push_back(new2Cr_2_5);
      _PCAWin_3Cr_2_5cm.push_back(new3Cr_2_5);
      _PCAWin_dist_2_5cm.push_back(lastdist_2_5);
      _PCAWin_npts_2_5cm.push_back(npts_2_5cm);
    }
    _shrPCA1CMed_2_5cm = CalcMed(_PCAWin_1Cr_2_5cm); //


    //Find Showering start
    for (unsigned i = 0; i < _PCAWin_1Cr_2_5cm.size(); ++i){
      if (_PCAWin_1Cr_2_5cm[i] < 0.9 && _PCAWin_1Cr_2_5cm[i] >= 0 && i == 0){
        _shrStart_2_5cm = 0;
       break;
      }else if(_PCAWin_1Cr_2_5cm[i] < 0.9 && _PCAWin_1Cr_2_5cm[i] >= 0){
        _shrStart_2_5cm = _PCAWin_dist_2_5cm[i];
        break;
      }
    }
    //Calculate MCS up to showering start
    SpcPointCloud shrStartcloud_2_5cm;
    for (unsigned i = 0; i < PFPcloud.size(); ++i){
      if (PFPcloud[i].dshrvtx < _shrStart_2_5cm){
        shrStartcloud_2_5cm.push_back(PFPcloud[i]);
      }
    }
    _shrStartMCS_2_5cm = DoMCSCalc(shrStartcloud_2_5cm, _shrvtx, PFPcentroid);

    //Calculate the MCS AS up to the showering start
    _shrMCSAS_2_5cm = DoMomASCalc(shrStartcloud_2_5cm, _shrvtx);
    PCAASResult PCAAS_2_5cm = DoPCAASCalc(shrStartcloud_2_5cm, _shrvtx(0), _shrvtx(1), _shrvtx(2));
    _shrPCA1CAS_2_5cm = PCAAS_2_5cm.AS1C; _shrPCA2CAS_2_5cm = PCAAS_2_5cm.AS2C; _shrPCA3CAS_2_5cm = PCAAS_2_5cm.AS3C;

    if(isnan(_shrStartMCS_2_5cm) || isinf(_shrStartMCS_2_5cm)) _shrStartMCS_2_5cm = -9999;
    if(isnan(_shrMCSAS_2_5cm) || isinf(_shrMCSAS_2_5cm)) _shrMCSAS_2_5cm = -9999;
    if(isnan(_shrPCA1CAS_2_5cm) || isinf(_shrPCA1CAS_2_5cm)) _shrPCA1CAS_2_5cm = -9999;
    if(isnan(_shrPCA2CAS_2_5cm) || isinf(_shrPCA2CAS_2_5cm)) _shrPCA2CAS_2_5cm = -9999;
    if(isnan(_shrPCA3CAS_2_5cm) || isinf(_shrPCA3CAS_2_5cm)) _shrPCA3CAS_2_5cm = -9999;

  }//End of looping through pfparticles in neutrino slice//
  
}

void MCS::setBranches(TTree *_tree){

  _tree->Branch("X_SpcPts_v", "std::vector< float >", &_X_SpcPts_v);
  _tree->Branch("Y_SpcPts_v", "std::vector< float >", &_Y_SpcPts_v);
  _tree->Branch("Z_SpcPts_v", "std::vector< float >", &_Z_SpcPts_v);

  _tree->Branch("shr_id_MCStool", &_shr_pfp_id, "shr_pfp_id/i");
  _tree->Branch("shr_hits_max_MCStool", &_shr_hits_max, "shr_hits_max_MCStool/i");
  _tree->Branch("n_showers_contained_MCStool", &_n_showers_contained, "n_showers_contained_MCStool/i");

  _tree->Branch("trkshrscore_v", "std::vector< float >", &_trkshrscore_v);

  _tree->Branch("shrPCA_1Cr", &_shrPCA_1Cr, "shrPCA_1Cr/f");
  _tree->Branch("shrPCA_2Cr", &_shrPCA_2Cr, "shrPCA_2Cr/f");
  _tree->Branch("shrPCA_3Cr", &_shrPCA_3Cr, "shrPCA_3Cr/f");

  _tree->Branch("shrPCA_1Ce", &_shrPCA_1Ce, "shrPCA_1Ce/f");
  _tree->Branch("shrPCA_2Ce", &_shrPCA_2Ce, "shrPCA_2Ce/f");
  _tree->Branch("shrPCA_3Ce", &_shrPCA_3Ce, "shrPCA_3Ce/f");

  _tree->Branch("shrPCA1CAS", &_shrPCA1CAS, "shrPCA1CAS/f");
  _tree->Branch("shrPCA2CAS", &_shrPCA2CAS, "shrPCA2CAS/f");
  _tree->Branch("shrPCA3CAS", &_shrPCA3CAS, "shrPCA3CAS/f");

  _tree->Branch("shrPCA_1Cr2h", &_shrPCA_1Cr2h, "shrPCA_1Cr2h/f");//new
  _tree->Branch("shrPCA_2Cr2h", &_shrPCA_2Cr2h, "shrPCA_2Cr2h/f");//new
  _tree->Branch("shrPCA_3Cr2h", &_shrPCA_3Cr2h, "shrPCA_3Cr2h/f");//new

  _tree->Branch("shrPCA_1Cr1h", &_shrPCA_1Cr1h, "shrPCA_1Cr1h/f");//new
  _tree->Branch("shrPCA_2Cr1h", &_shrPCA_2Cr1h, "shrPCA_2Cr1h/f");//new
  _tree->Branch("shrPCA_3Cr1h", &_shrPCA_3Cr1h, "shrPCA_3Cr1h/f");//new

  _tree->Branch("shrMCSMom", &_shrMCSMom, "shrMCSMom/f");
  _tree->Branch("shrMCSMom1h", &_shrMCSMom1h, "shrMCSMom1h/f");
  _tree->Branch("shrMCSMom2h", &_shrMCSMom2h, "shrMCSMom2h/f");

  _tree->Branch("shrPCALen", &_shrPCALen, "shrPCALen/f");
  _tree->Branch("n_shrSpcPts", &_n_shrSpcPts, "n_shrSpcPts/i");

  _tree->Branch("PCAWin_1Cr_5cm", "std::vector<float>", &_PCAWin_1Cr_5cm);
  _tree->Branch("PCAWin_2Cr_5cm", "std::vector<float>", &_PCAWin_2Cr_5cm);
  _tree->Branch("PCAWin_3Cr_5cm", "std::vector<float>", &_PCAWin_3Cr_5cm);
  _tree->Branch("PCAWin_dist_5cm", "std::vector<float>", &_PCAWin_dist_5cm);
  _tree->Branch("PCAWin_npts_5cm", "std::vector<int>", &_PCAWin_npts_5cm);
  
  _tree->Branch("shrStart_5cm", &_shrStart_5cm, "shrStart_5cm/f");
  _tree->Branch("shrStartMCS_5cm", &_shrStartMCS_5cm, "shrStartMCS_5cm/f");
  _tree->Branch("shrMCSAS_5cm",&_shrMCSAS_5cm,"shrMCSAS_5cm/f");
  _tree->Branch("shrPCA1CAS_5cm", &_shrPCA1CAS_5cm, "shrPCA1CAS_5cm/f");
  _tree->Branch("shrPCA2CAS_5cm", &_shrPCA2CAS_5cm, "shrPCA2CAS_5cm/f");
  _tree->Branch("shrPCA3CAS_5cm", &_shrPCA3CAS_5cm, "shrPCA3CAS_5cm/f");
  _tree->Branch("shrPCA1CMed_5cm", &_shrPCA1CMed_5cm, "_shrPCA1CMed_5cm/f");


  _tree->Branch("PCAWin_1Cr_2_5cm", "std::vector<float>", &_PCAWin_1Cr_2_5cm);
  _tree->Branch("PCAWin_2Cr_2_5cm", "std::vector<float>", &_PCAWin_2Cr_2_5cm);
  _tree->Branch("PCAWin_3Cr_2_5cm", "std::vector<float>", &_PCAWin_3Cr_2_5cm);
  _tree->Branch("PCAWin_dist_2_5cm", "std::vector<float>", &_PCAWin_dist_2_5cm);
  _tree->Branch("PCAWin_npts_2_5cm", "std::vector<int>", &_PCAWin_npts_2_5cm);
  
  _tree->Branch("shrStart_2_5cm", &_shrStart_2_5cm, "shrStart_2_5cm/f");
  _tree->Branch("shrStartMCS_2_5cm", &_shrStartMCS_2_5cm, "shrStartMCS_2_5cm/f");
  _tree->Branch("shrMCSAS_2_5cm",&_shrMCSAS_2_5cm,"shrMCSAS_2_5cm/f");
  _tree->Branch("shrPCA1CAS_2_5cm", &_shrPCA1CAS_2_5cm, "shrPCA1CAS_2_5cm/f");
  _tree->Branch("shrPCA2CAS_2_5cm", &_shrPCA2CAS_2_5cm, "shrPCA2CAS_2_5cm/f");
  _tree->Branch("shrPCA3CAS_2_5cm", &_shrPCA3CAS_2_5cm, "shrPCA3CAS_2_5cm/f");
  _tree->Branch("shrPCA1CMed_2_5cm", &_shrPCA1CMed_2_5cm, "_shrPCA1CMed_2_5cm/f");

  _tree->Branch("DeltaMed", &_DeltaMed, "DeltaMed/f");
  _tree->Branch("DeltaMed1h", &_DeltaMed1h, "DeltaMed1h/f");
  _tree->Branch("DeltaMed2h", &_DeltaMed2h, "DeltaMed2h/f");

  _tree->Branch("DeltaRMS", &_DeltaRMS, "DeltaRMS/f");
  _tree->Branch("DeltaRMS1h", &_DeltaRMS1h, "DeltaRMS1h/f");
  _tree->Branch("DeltaRMS2h", &_DeltaRMS2h, "DeltaRMS2h/f");

  _tree->Branch("CylFrac_1cm", &_CylFrac_1cm, "CylFrac_1cm/f");
  _tree->Branch("CylFrac1h_1cm", &_CylFrac1h_1cm, "CylFrac1h_1cm/f");
  _tree->Branch("CylFrac2h_1cm", &_CylFrac2h_1cm, "CylFrac2h_1cm/f");

  _tree->Branch("CylFrac_2cm", &_CylFrac_2cm, "CylFrac_2cm/f");
  _tree->Branch("CylFrac1h_2cm", &_CylFrac1h_2cm, "CylFrac1h_2cm/f");
  _tree->Branch("CylFrac2h_2cm", &_CylFrac2h_2cm, "CylFrac2h_2cm/f");

  _tree->Branch("CylFrac_3cm", &_CylFrac_3cm, "CylFrac_3cm/f");
  _tree->Branch("CylFrac1h_3cm", &_CylFrac1h_3cm, "CylFrac1h_3cm/f");
  _tree->Branch("CylFrac2h_3cm", &_CylFrac2h_3cm, "CylFrac2h_3cm/f");

  _tree->Branch("CylFrac_4cm", &_CylFrac_4cm, "CylFrac_4cm/f");
  _tree->Branch("CylFrac1h_4cm", &_CylFrac1h_4cm, "CylFrac1h_4cm/f");
  _tree->Branch("CylFrac2h_4cm", &_CylFrac2h_4cm, "CylFrac2h_4cm/f");

  _tree->Branch("CylFrac_5cm", &_CylFrac_5cm, "CylFrac_5cm/f");
  _tree->Branch("CylFrac1h_5cm", &_CylFrac1h_5cm, "CylFrac1h_5cm/f");
  _tree->Branch("CylFrac2h_5cm", &_CylFrac2h_5cm, "CylFrac2h_5cm/f");

}

void MCS::resetTTree(TTree *_tree)
{
  _reco_nu_vtx_x = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_y = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_z = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_sce_x = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_sce_y = std::numeric_limits<float>::lowest();
  _reco_nu_vtx_sce_z = std::numeric_limits<float>::lowest();

  _X_SpcPts_v.clear(); _Y_SpcPts_v.clear(); _Z_SpcPts_v.clear();

  _shr_pfp_id = 0; _shr_hits_max = 0; _n_showers_contained = 0;  _hits_outfv = 0;

  _trkshrscore_v.clear(); 

  _shrPCA_1Cr = -1; _shrPCA_2Cr = -1; _shrPCA_3Cr = -1;
  _shrPCA_1Ce = -1; _shrPCA_2Ce = -1; _shrPCA_3Ce = -1;
  _shrMCSMom = -1;
  _shrPCA1CAS = -9999; _shrPCA2CAS = -9999; _shrPCA3CAS = -9999;
  _shrMCSMom2h = -1; _shrMCSMom1h = -1;

  _shrPCA_1Cr2h = -1; _shrPCA_2Cr2h = -1; _shrPCA_3Cr2h = -1;
  _shrPCA_1Cr1h = -1; _shrPCA_2Cr1h = -1; _shrPCA_3Cr1h = -1;

  _shrPCALen = -1; _n_shrSpcPts = -1;

  _DeltaMed = -1; _DeltaMed1h = -1; _DeltaMed2h = -1;
  _DeltaMed = -1; _DeltaRMS1h = -1; _DeltaRMS2h = -1;

  _CylFrac_1cm = -1; _CylFrac1h_1cm = -1; _CylFrac2h_1cm = -1;
  _CylFrac_2cm = -1; _CylFrac1h_2cm = -1; _CylFrac2h_2cm = -1;
  _CylFrac_3cm = -1; _CylFrac1h_3cm = -1; _CylFrac2h_3cm = -1;
  _CylFrac_4cm = -1; _CylFrac1h_4cm = -1; _CylFrac2h_4cm = -1;
  _CylFrac_5cm = -1; _CylFrac1h_5cm = -1; _CylFrac2h_5cm = -1;


  _PCAWin_1Cr_5cm.clear(); _PCAWin_2Cr_5cm.clear(); _PCAWin_3Cr_5cm.clear();
  _PCAWin_dist_5cm.clear(); _PCAWin_npts_5cm.clear();
  _shrStart_5cm = -1; _shrStartMCS_5cm = -1; _shrMCSAS_5cm = -9999; 
  _shrPCA1CAS_5cm = -9999; _shrPCA2CAS_5cm = -9999; _shrPCA3CAS_5cm = -9999;
  _shrPCA1CMed_5cm = -1;

  _PCAWin_1Cr_2_5cm.clear(); _PCAWin_2Cr_2_5cm.clear(); _PCAWin_3Cr_2_5cm.clear();
  _PCAWin_dist_2_5cm.clear(); _PCAWin_npts_2_5cm.clear();
  _shrStart_2_5cm = -1; _shrStartMCS_2_5cm = -1; _shrMCSAS_2_5cm = -9999; 
  _shrPCA1CAS_2_5cm = -9999; _shrPCA2CAS_2_5cm = -9999; _shrPCA3CAS_2_5cm = -9999;
  _shrPCA1CMed_2_5cm = -1;
}

PCAResults MCS::DoPCA(const SpcPointCloud &points) {
    TVector3 outputCentroid;
    std::pair<TVector3,TVector3> outputEndPoints;
    float outputLength;
    TVector3 outputEigenValues;
    TVector3 outputEigenRatios;
    std::vector<TVector3> outputEigenVecs;
    float meanPosition[3] = {0., 0., 0.};
    unsigned int nThreeDHits = 0;
    int outputN_pts = points.size();
    int outputPDG = points[0].pdg;
    for (unsigned int i = 0; i < points.size(); i++) {
      meanPosition[0] += points[i].x;
      meanPosition[1] += points[i].y;
      meanPosition[2] += points[i].z;
      ++nThreeDHits;
    }
    if (nThreeDHits == 0) {
      PCAResults results;
      return results;
    }
    const float nThreeDHitsAsFloat(static_cast<float>(nThreeDHits));
    meanPosition[0] /= nThreeDHitsAsFloat;
    meanPosition[1] /= nThreeDHitsAsFloat;
    meanPosition[2] /= nThreeDHitsAsFloat;
    outputCentroid = TVector3(meanPosition[0], meanPosition[1], meanPosition[2]);
    float xi2 = 0.0;
    float xiyi = 0.0;
    float xizi = 0.0;
    float yi2 = 0.0;
    float yizi = 0.0;
    float zi2 = 0.0;
    float weightSum = 0.0;
    for (unsigned int i = 0; i < points.size(); i++) {
      const float weight(1.);
      const float x((points[i].x - meanPosition[0]) * weight);
      const float y((points[i].y - meanPosition[1]) * weight);
      const float z((points[i].z - meanPosition[2]) * weight);
      xi2 += x * x;
      xiyi += x * y;
      xizi += x * z;
      yi2 += y * y;
      yizi += y * z;
      zi2 += z * z;
      weightSum += weight * weight;
    }

    Eigen::Matrix3f sig;

    sig <<  xi2, xiyi, xizi,
        xiyi, yi2, yizi,
        xizi, yizi, zi2;

    sig *= 1.0 / weightSum;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(sig);

    typedef std::pair<float,size_t> EigenValColPair;
    typedef std::vector<EigenValColPair> EigenValColVector;

    EigenValColVector eigenValColVector;
    const auto &resultEigenMat(eigenMat.eigenvalues());
    eigenValColVector.emplace_back(resultEigenMat(0), 0);
    eigenValColVector.emplace_back(resultEigenMat(1), 1);
    eigenValColVector.emplace_back(resultEigenMat(2), 2);

    std::sort(eigenValColVector.begin(), eigenValColVector.end(), [](const EigenValColPair &left, const EigenValColPair &right){return left.first > right.first;} );

    outputEigenValues = TVector3(eigenValColVector.at(0).first, eigenValColVector.at(1).first, eigenValColVector.at(2).first);
    double eigTot = eigenValColVector.at(0).first + eigenValColVector.at(1).first + eigenValColVector.at(2).first;
    outputEigenRatios = TVector3(eigenValColVector.at(0).first/eigTot, eigenValColVector.at(1).first/eigTot, eigenValColVector.at(2).first/eigTot);

    const Eigen::Matrix3f &eigenVecs(eigenMat.eigenvectors());

    for (const EigenValColPair &pair : eigenValColVector) {
      outputEigenVecs.emplace_back(eigenVecs(0, pair.second), eigenVecs(1, pair.second), eigenVecs(2, pair.second));
    }

    PCAResults results;

    Eigen::ParametrizedLine<float,3> priAxis(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)),Eigen::Vector3f(outputEigenVecs[0](0),outputEigenVecs[0](1),outputEigenVecs[0](2)));

    Eigen::Vector3f endPoint1(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));
    Eigen::Vector3f endPoint2(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));

    Eigen::Vector3f testPoint;
    Eigen::Vector3f projTestPoint;
    float maxDist1 = -1.0;
    float maxDist2 = -1.0;
    float dist;
    float dotP;
    for (unsigned int i = 0; i < points.size(); i++) {
      testPoint = Eigen::Vector3f(points[i].x,points[i].y,points[i].z);
      projTestPoint = priAxis.projection(testPoint);
      dist = sqrt(pow(projTestPoint(0)-outputCentroid(0),2.0)+pow(projTestPoint(1)-outputCentroid(1),2.0)+pow(projTestPoint(2)-outputCentroid(2),2.0));
      dotP = (projTestPoint(0)-outputCentroid(0))*outputEigenVecs[0](0) + (projTestPoint(1)-outputCentroid(1))*outputEigenVecs[0](1) + (projTestPoint(2)-outputCentroid(2))*outputEigenVecs[0](2);


      if ((dotP < 0.0) && (dist > maxDist1)) {
        endPoint1 = projTestPoint;
        maxDist1 = dist;
      }
      else if ((dotP > 0.0) && (dist > maxDist2)) {
        endPoint2 = projTestPoint;
        maxDist2 = dist;
      }
    }
    outputEndPoints.first = TVector3(endPoint1(0),endPoint1(1),endPoint1(2));
    outputEndPoints.second = TVector3(endPoint2(0),endPoint2(1),endPoint2(2));
    outputLength = sqrt(pow(endPoint2(0)-endPoint1(0),2.0)+pow(endPoint2(1)-endPoint1(1),2.0)+pow(endPoint2(2)-endPoint1(2),2.0));
    results.centroid = outputCentroid;
    results.endPoints = outputEndPoints;
    results.length = outputLength;
    results.eVals = outputEigenValues;
    results.eRatios = outputEigenRatios;
    results.eVecs = outputEigenVecs;
    results.n_pts = outputN_pts;
    results.pdg = outputPDG;
    return results;
  }

  double MCS::DoDeltaRMS(SpcPointCloud cloud, TVector3 StartP, TVector3 EndP){
    double drms = 0.0;
    if(cloud.size() < 2) return drms;
    //TVector3 StartP(cloud.at(0).x, cloud.at(0).y, cloud.at(0).z);
    //TVector3 EndP(cloud.back().x, cloud.back().y, cloud.back().z);
    
    for (unsigned i = 0; i < cloud.size(); ++i){
      TVector3 P(cloud.at(i).x, cloud.at(i).y, cloud.at(i).z);
      double delta = linepointdist(StartP, EndP, P);
      drms += pow(delta, 2.);
    }
    drms = sqrt(drms/((double) cloud.size()));
    if(isnan(drms) || isinf(drms)) drms = -1;
    return drms;
  }

  double MCS::pythagoras(double x1, double x2, double y1, double y2, double z1, double z2){
    return std::sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
  }

  double MCS::dotProd(TVector3 v1, TVector3 v2){
    double dotP;
    dotP = v1(0)*v2(0) + v1(1)*v2(1) + v1(2)*v2(2);
    return dotP;
  }
  TVector3 MCS::CrossProd(TVector3 u, TVector3 v){
    double x,y,z;
    x = u(1)*v(2) - u(2)*v(1);
    y = -(u(0)*v(2) - u(2)*v(0));
    z = u(0)*v(1) - u(1)*v(0);
    TVector3 result(x,y,z);
    return result;
  }

  double MCS::linepointdist(TVector3 lStart, TVector3 lEnd, TVector3 point){
    // Return minimum distance between line segment lStart-lEnd and point
    
    TVector3 line_v(lEnd(0) - lStart(0), lEnd(1) - lStart(1), lEnd(2) - lStart(2));
    TVector3 StartPoint_v(lStart(0) - point(0), lStart(1) - point(1), lStart(2) - point(2));

    double lmod = pythagoras(line_v(0), 0, line_v(1), 0, line_v(2), 0);
    /*
    TVector3 line_uv(line_v(0)/lmod, line_v(1)/lmod, line_v(2)/lmod);

    double dP_lp = dotProd(point,line_v); 
    TVector3 projPoint((dP_lp/lmod)*line_uv(0),(dP_lp/lmod)*line_uv(1), (dP_lp/lmod)*line_uv(2));

    double dist = pythagoras(point(0), projPoint(0), point(1),projPoint(1), point(2), projPoint(2));
    */
    TVector3 cross = CrossProd(line_v, StartPoint_v);
    double numerator = pythagoras(cross(0),0,cross(1),0,cross(2),0);
    double dist = numerator / lmod;
    return dist;
  }

  float MCS::DoMCSCalc(SpcPointCloud cloud, TVector3 StartP, TVector3 EndP){
    PCAResults PCAcloud = DoPCA(cloud);
    double X = PCAcloud.length;
    double drms = DoDeltaRMS(cloud, StartP, EndP);
    float X_0 = 14.;
    //double X = pythagoras(StartP(0),EndP(0),StartP(1),EndP(1),StartP(2),EndP(2))
    float mcs = (13.6/(4.*sqrt(3.)))*(X*sqrt(X/X_0))/drms;
    if(isnan(mcs) || isinf(mcs)) mcs = -1;
    return mcs;
  }

  float MCS::DoCylFracCalc(SpcPointCloud cloud, TVector3 StartP, TVector3 EndP, float rad){
    float cylfrac = 0.0;
    for(unsigned i = 0; i < cloud.size(); ++i){
      TVector3 point(cloud[i].x, cloud[i].y, cloud[i].z);
      double delta = linepointdist(StartP, EndP, point);
      if (delta <= rad) cylfrac += 1;
    }
    cylfrac = cylfrac/((float) cloud.size());
    if(isnan(cylfrac) || isinf(cylfrac)) cylfrac = -1;
    return cylfrac;
  }

  float MCS::DoMomASCalc(SpcPointCloud cloud, TVector3 vtx){
    SpcPointCloud cloud1h, cloud2h;
    //float X_0 = 14.;
    PCAResults PCACloud = DoPCA(cloud);
    double vcDist = pythagoras(PCACloud.centroid(0), vtx(0), PCACloud.centroid(1), vtx(1), PCACloud.centroid(2), vtx(2));
    for(unsigned i = 0; i < cloud.size(); ++i){
      if(cloud[i].dshrvtx < vcDist){
        cloud1h.push_back(cloud[i]);
      }else{
        cloud2h.push_back(cloud[i]);
      }
    }
    PCAResults PCACloud1h = DoPCA(cloud1h);
    PCAResults PCACloud2h = DoPCA(cloud2h);
    if(cloud2h.size() == 0) return -1;
    sort( cloud2h.begin(), cloud2h.end(), [ ]( const SpcPoint& lhs, const SpcPoint& rhs ){return lhs.dshrvtx < rhs.dshrvtx;});
    TVector3 FirstP2h(cloud2h[0].x, cloud2h[0].y, cloud2h[0].z);
    float mcs1h = DoMCSCalc(cloud1h, vtx, PCACloud1h.centroid); float mcs2h = DoMCSCalc(cloud2h, FirstP2h, PCACloud2h.centroid);
    float mcsAS = (mcs1h - mcs2h)/(mcs1h + mcs2h);
    if(isnan(mcsAS) || isinf(mcsAS)) mcsAS = -1;
    return mcsAS;
  }

PCAASResult MCS::DoPCAASCalc(SpcPointCloud cloud, float nuvX, float nuvY, float nuvZ){
  SpcPointCloud cloud1h, cloud2h;
    PCAResults PCACloud = DoPCA(cloud);
    double vcDist = pythagoras(PCACloud.centroid(0), nuvX, PCACloud.centroid(1), nuvY, PCACloud.centroid(2), nuvZ);
    for(unsigned i = 0; i < cloud.size(); ++i){
      if(cloud[i].dshrvtx < vcDist){
        cloud1h.push_back(cloud[i]);
      }else{
        cloud2h.push_back(cloud[i]);
      }
    }
    PCAResults pcaCloud1h = DoPCA(cloud1h); PCAResults pcaCloud2h = DoPCA(cloud2h);
    float pca1C_1h = pcaCloud1h.eRatios(0); float pca1C_2h = pcaCloud2h.eRatios(0); 
    float pca1C_AS = (pca1C_1h - pca1C_2h)/(pca1C_1h + pca1C_2h);
    
    float pca2C_1h = pcaCloud1h.eRatios(1); float pca2C_2h = pcaCloud2h.eRatios(1); 
    float pca2C_AS = (pca2C_1h - pca2C_2h)/(pca2C_1h + pca2C_2h);

    float pca3C_1h = pcaCloud1h.eRatios(2); float pca3C_2h = pcaCloud2h.eRatios(2); 
    float pca3C_AS = (pca3C_1h - pca3C_2h)/(pca3C_1h + pca3C_2h);
    if(isnan(pca1C_AS) || isinf(pca1C_AS)) pca1C_AS = -9999;
    if(isnan(pca2C_AS) || isinf(pca2C_AS)) pca2C_AS = -9999;
    if(isnan(pca3C_AS) || isinf(pca3C_AS)) pca3C_AS = -9999;
    PCAASResult results;
    results.AS1C = pca1C_AS; results.AS2C = pca2C_AS; results.AS3C = pca3C_AS;
    return results;
}

float MCS::CalcMed(vector<float> scores){
  size_t size = scores.size();
  if (size == 0){
    return 0;  // Undefined, really.
  }
  else{
    sort(scores.begin(), scores.end());
    if (size % 2 == 0){
      return (scores[size / 2 - 1] + scores[size / 2]) / 2.;
    }else{
      return scores[size / 2];
    }
  }
}


float MCS::DoDeltaMed(SpcPointCloud cloud, TVector3 LStart, TVector3 LEnd){
  std::vector<float> deltas_v;
  for (unsigned i = 0; i < cloud.size(); ++i){
    TVector3 point(cloud[i].x,cloud[i].y,cloud[i].z);
    float delta = linepointdist(LStart, LEnd, point);
    deltas_v.push_back(delta);
  }
  float median = CalcMed(deltas_v);
  if (isnan(median)||isinf(median)) median = -1.;
  return median;
}  

void MCS::Project3Dto2D(const TVector3& pt3d, const int& pl,
  const float& wire2cm, const float& time2cm,
  float& wirecm, float& timecm) {

  auto const* geom = ::lar::providerFrom<geo::Geometry>();
    
  wirecm = geom->WireCoordinate(pt3d[1],pt3d[2],geo::PlaneID(0,0,pl)) * wire2cm;
  timecm = pt3d[0];

  return;
}

DEFINE_ART_CLASS_TOOL(MCS)
} // namespace analysis

#endif
