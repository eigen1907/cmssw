#ifndef MuonIdentification_MuonIdProducer_h
#define MuonIdentification_MuonIdProducer_h

// -*- C++ -*-
//
// Package:    MuonIdentification
// Class:      MuonIdProducer
//
/*

 Description: reco::Muon producer that can fill various information:
              - track-segment matching
              - energy deposition
              - muon isolation
              - muon hypothesis compatibility (calorimeter)
              Acceptable inputs:
              - reco::TrackCollection
              - reco::MuonCollection
              - reco::MuonTrackLinksCollection
*/
//
// Original Author:  Dmytro Kovalskyi
//
//

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
// #include "Utilities/Timing/interface/TimerStack.h"

#include "RecoMuon/MuonIdentification/interface/MuonTimingFiller.h"
#include "RecoMuon/MuonIdentification/interface/MuonCaloCompatibility.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractor.h"
#include "RecoMuon/MuonIdentification/interface/MuonShowerDigiFiller.h"

// RPC-Muon stuffs
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/MuonReco/interface/MuonRPCHitMatch.h"
#include "DataFormats/MuonReco/interface/MuonGEMHitMatch.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"

#include "RecoMuon/MuonIdentification/interface/MuonIdTruthInfo.h"
#include "RecoMuon/MuonIdentification/interface/MuonArbitrationMethods.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment2D.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0RecHitCollection.h"

#include "TTree.h"
#include <chrono> 

class MuonMesh;
class MuonKinkFinder;

class MuonIdProducer : public edm::stream::EDProducer<> {
public:
  typedef reco::Muon::MuonTrackType TrackType;

  explicit MuonIdProducer(const edm::ParameterSet&);

  ~MuonIdProducer() override;

  void produce(edm::Event&, const edm::EventSetup&) override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override;

  static double sectorPhi(const DetId& id);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void fillMuonId(edm::Event&,
                  const edm::EventSetup&,
                  reco::Muon&,
                  TrackDetectorAssociator::Direction direction = TrackDetectorAssociator::InsideOut);
  void fillArbitrationInfo(reco::MuonCollection*, unsigned int muonType = reco::Muon::TrackerMuon);
  void arbitrateMuons(reco::MuonCollection*, reco::CaloMuonCollection*);
  void fillMuonIsolation(edm::Event&,
                         const edm::EventSetup&,
                         reco::Muon& aMuon,
                         reco::IsoDeposit& trackDep,
                         reco::IsoDeposit& ecalDep,
                         reco::IsoDeposit& hcalDep,
                         reco::IsoDeposit& hoDep,
                         reco::IsoDeposit& jetDep);
  void fillGlbQuality(edm::Event&, const edm::EventSetup&, reco::Muon& aMuon);
  void fillTrackerKink(reco::Muon& aMuon);
  void init(edm::Event&, const edm::EventSetup&);

  // make a muon based on a track ref
  reco::Muon makeMuon(edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::TrackRef& track, TrackType type);
  // make a global muon based on the links object
  reco::Muon makeMuon(const reco::MuonTrackLinks& links);

  // make a muon based on track (p4)
  reco::Muon makeMuon(const reco::Track& track);

  reco::CaloMuon makeCaloMuon(const reco::Muon&);

  // check if a silicon track satisfies the trackerMuon requirements
  bool isGoodTrack(const reco::Track& track);

  bool isGoodTrackerMuon(const reco::Muon& muon);
  bool isGoodCaloMuon(const reco::CaloMuon& muon);
  bool isGoodRPCMuon(const reco::Muon& muon);
  bool isGoodGEMMuon(const reco::Muon& muon);
  bool isGoodME0Muon(const reco::Muon& muon);

  // check number of common DetIds for a given trackerMuon and a stand alone
  // muon track
  int overlap(const reco::Muon& muon, const reco::Track& track);

  unsigned int chamberId(const DetId&);

  double phiOfMuonInteractionRegion(const reco::Muon& muon) const;

  bool checkLinks(const reco::MuonTrackLinks*) const;
  inline bool approxEqual(const double a, const double b, const double tol = 1E-3) const {
    return std::abs(a - b) < tol;
  }

  /// get the segment matches of the appropriate type
  std::vector<reco::MuonSegmentMatch>* getSegmentMatches(reco::MuonChamberMatch& chamber, unsigned int muonType) const {
    if (muonType == reco::Muon::TrackerMuon)
      return &chamber.segmentMatches;
    else if (muonType == reco::Muon::ME0Muon)
      return &chamber.me0Matches;
    else if (muonType == reco::Muon::GEMMuon)
      return &chamber.gemMatches;
    else
      throw cms::Exception("getSegmentMatches called with unsupported muonType");
  }

  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters parameters_;

  struct ICTypes {
    enum ICTypeKey { INNER_TRACKS, OUTER_TRACKS, LINKS, MUONS, TEV_FIRSTHIT, TEV_PICKY, TEV_DYT, NONE };

    static ICTypeKey toKey(const std::string& s) {
      if (s == "inner tracks")
        return INNER_TRACKS;
      else if (s == "outer tracks")
        return OUTER_TRACKS;
      else if (s == "links")
        return LINKS;
      else if (s == "muons")
        return MUONS;
      else if (s == "tev firstHit")
        return TEV_FIRSTHIT;
      else if (s == "tev picky")
        return TEV_PICKY;
      else if (s == "tev dyt")
        return TEV_DYT;

      throw cms::Exception("FatalError") << "Unknown input collection type: " << s;
    }

    static std::string toStr(const ICTypeKey k) {
      switch (k) {
        case INNER_TRACKS:
          return "inner tracks";
        case OUTER_TRACKS:
          return "outer tracks";
        case LINKS:
          return "links";
        case MUONS:
          return "muons";
        case TEV_FIRSTHIT:
          return "tev firstHit";
        case TEV_PICKY:
          return "tev picky";
        case TEV_DYT:
          return "tev dyt";
        default:
          throw cms::Exception("FatalError") << "Unknown input collection type";
      }
      return "";
    }
  };
  std::vector<edm::InputTag> inputCollectionLabels_;
  std::vector<ICTypes::ICTypeKey> inputCollectionTypes_;

  std::unique_ptr<MuonTimingFiller> theTimingFiller_;

  std::unique_ptr<MuonShowerDigiFiller> theShowerDigiFiller_;

  // selections
  double minPt_;
  double minP_;
  double minPCaloMuon_;
  int minNumberOfMatches_;
  double maxAbsEta_;
  bool addExtraSoftMuons_;

  // matching
  double maxAbsDx_;
  double maxAbsPullX2_;
  double maxAbsDy_;
  double maxAbsPullY2_;

  // what information to fill
  bool fillCaloCompatibility_;
  bool fillEnergy_;
  bool storeCrossedHcalRecHits_;
  bool fillMatching_;
  bool fillShowerDigis_;
  bool fillIsolation_;
  bool writeIsoDeposits_;
  double ptThresholdToFillCandidateP4WithGlobalFit_;
  double sigmaThresholdToFillCandidateP4WithGlobalFit_;

  bool arbitrateTrackerMuons_;

  bool debugWithTruthMatching_;

  edm::Handle<reco::TrackCollection> innerTrackCollectionHandle_;
  edm::Handle<reco::TrackCollection> outerTrackCollectionHandle_;
  edm::Handle<reco::TrackCollection> outerTrackSecondaryCollectionHandle_;
  edm::Handle<reco::MuonCollection> muonCollectionHandle_;
  edm::Handle<reco::MuonTrackLinksCollection> linkCollectionHandle_;
  edm::Handle<reco::TrackToTrackMap> tpfmsCollectionHandle_;
  edm::Handle<reco::TrackToTrackMap> pickyCollectionHandle_;
  edm::Handle<reco::TrackToTrackMap> dytCollectionHandle_;
  edm::Handle<reco::VertexCollection> pvHandle_;

  edm::EDGetTokenT<reco::TrackCollection> innerTrackCollectionToken_;
  edm::EDGetTokenT<reco::TrackCollection> outerTrackCollectionToken_;
  edm::EDGetTokenT<reco::TrackCollection> outerTrackSecondaryCollectionToken_;
  edm::EDGetTokenT<reco::MuonCollection> muonCollectionToken_;
  edm::EDGetTokenT<reco::MuonTrackLinksCollection> linkCollectionToken_;
  edm::EDGetTokenT<reco::TrackToTrackMap> tpfmsCollectionToken_;
  edm::EDGetTokenT<reco::TrackToTrackMap> pickyCollectionToken_;
  edm::EDGetTokenT<reco::TrackToTrackMap> dytCollectionToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;

  edm::EDGetTokenT<edm::ValueMap<reco::MuonQuality> > glbQualToken_;
  
  edm::EDGetTokenT<RPCRecHitCollection> rpcHitToken_;
  edm::EDGetTokenT<GEMRecHitCollection> gemHitToken_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentToken_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegmentToken_;
  edm::EDGetTokenT<GEMSegmentCollection> gemSegmentToken_;

  edm::Handle<RPCRecHitCollection> rpcHitHandle_;
  edm::Handle<GEMRecHitCollection> gemHitHandle_;
  edm::Handle<DTRecSegment4DCollection> dtSegmentHandle_;
  edm::Handle<CSCSegmentCollection> cscSegmentHandle_;
  edm::Handle<GEMSegmentCollection> gemSegmentHandle_;

  edm::Handle<edm::ValueMap<reco::MuonQuality> > glbQualHandle_;

  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> geomTokenRun_;
  const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatorToken_;
  edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> globalGeomToken_;

  MuonCaloCompatibility muonCaloCompatibility_;
  std::unique_ptr<reco::isodeposit::IsoDepositExtractor> muIsoExtractorCalo_;
  std::unique_ptr<reco::isodeposit::IsoDepositExtractor> muIsoExtractorTrack_;
  std::unique_ptr<reco::isodeposit::IsoDepositExtractor> muIsoExtractorJet_;
  std::string trackDepositName_;
  std::string ecalDepositName_;
  std::string hcalDepositName_;
  std::string hoDepositName_;
  std::string jetDepositName_;

  bool fillGlobalTrackQuality_;
  bool fillGlobalTrackRefits_;
  edm::InputTag globalTrackQualityInputTag_;

  edm::InputTag pvInputTag_;
  bool selectHighPurity_;
  bool fillTrackerKink_;
  std::unique_ptr<MuonKinkFinder> trackerKinkFinder_;

  double caloCut_;

  bool arbClean_;
  std::unique_ptr<MuonMesh> meshAlgo_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemgeomToken_;
  const GEMGeometry* gemgeom;
  double GEM_edgecut_;


  struct DeepMuonRecoStruct {
    std::vector<unsigned int>    trackKey;
    std::vector<double> trackVx;
    std::vector<double> trackVy;
    std::vector<double> trackVz;
    std::vector<double> trackPx;
    std::vector<double> trackPy;
    std::vector<double> trackPz;
    std::vector<double> trackQOverP;
    std::vector<double> trackLambda;
    std::vector<double> trackPhi;
    std::vector<double> trackDxy;
    std::vector<double> trackDsz;
    std::vector<double> trackQOverPErr;
    std::vector<double> trackLambdaErr;
    std::vector<double> trackPhiErr;
    std::vector<double> trackDxyErr;
    std::vector<double> trackDszErr;
    std::vector<int>    trackCharge;
    std::vector<double> trackChi2;
    std::vector<double> trackNdof;
    std::vector<unsigned int> trackIsSplit;
    std::vector<unsigned int> trackIsNormalTrkMuon;
    std::vector<unsigned int> trackIsInOutTrkMuon;
    std::vector<unsigned int> trackIsOutInTrkMuon;
    std::vector<unsigned int> trackIsTrkMuon;
    std::vector<unsigned int> trackIsTrkMuonWithArb;

    std::vector<unsigned int> rpcHitRawId;
    std::vector<float>        rpcHitPosX;
    std::vector<float>        rpcHitPosY;
    std::vector<float>        rpcHitPosZ;
    std::vector<float>        rpcHitPosXErr;
    std::vector<float>        rpcHitPosYErr;
    std::vector<int>          rpcHitCls;
    std::vector<int>          rpcHitBX;

    std::vector<unsigned int> gemHitRawId;
    std::vector<float>        gemHitPosX;
    std::vector<float>        gemHitPosY;
    std::vector<float>        gemHitPosZ;
    std::vector<float>        gemHitPosXErr;
    std::vector<float>        gemHitPosYErr;
    std::vector<int>          gemHitCls;
    std::vector<int>          gemHitBX;

    std::vector<unsigned int> dtSegRawId;
    std::vector<float>        dtSegPosX;
    std::vector<float>        dtSegPosY;
    std::vector<float>        dtSegPosZ;
    std::vector<float>        dtSegPosXErr;
    std::vector<float>        dtSegPosYErr;
    std::vector<float>        dtSegDirX;
    std::vector<float>        dtSegDirY;
    std::vector<float>        dtSegDirZ;
    std::vector<float>        dtSegDirXErr;
    std::vector<float>        dtSegDirYErr;
    std::vector<float>        dtSegChi2;
    std::vector<int>          dtSegNdof;
    
    std::vector<unsigned int> cscSegRawId;
    std::vector<float>        cscSegPosX;
    std::vector<float>        cscSegPosY;
    std::vector<float>        cscSegPosZ;
    std::vector<float>        cscSegPosXErr;
    std::vector<float>        cscSegPosYErr;
    std::vector<float>        cscSegDirX;
    std::vector<float>        cscSegDirY;
    std::vector<float>        cscSegDirZ;
    std::vector<float>        cscSegDirXErr;
    std::vector<float>        cscSegDirYErr;
    std::vector<float>        cscSegChi2;
    std::vector<int>          cscSegNdof;
    std::vector<float>        cscSegTime;

    void setBranches(TTree* tree) {
      tree->Branch("track_key",                              &trackKey);
      tree->Branch("track_vx",                               &trackVx);
      tree->Branch("track_vy",                               &trackVy);
      tree->Branch("track_vz",                               &trackVz);
      tree->Branch("track_px",                               &trackPx);
      tree->Branch("track_py",                               &trackPy);
      tree->Branch("track_pz",                               &trackPz);
      tree->Branch("track_qoverp",                           &trackQOverP);
      tree->Branch("track_lamda",                            &trackLambda);
      tree->Branch("track_phi",                              &trackPhi);
      tree->Branch("track_dxy",                              &trackDxy);
      tree->Branch("track_dsz",                              &trackDsz);
      tree->Branch("track_qoverp_err",                     &trackQOverPErr);
      tree->Branch("track_lambda_err",                     &trackLambdaErr);
      tree->Branch("track_phi_err",                        &trackPhiErr);
      tree->Branch("track_dxy_err",                        &trackDxyErr);
      tree->Branch("track_dsz_err",                        &trackDszErr);
      tree->Branch("track_charge",                           &trackCharge);
      tree->Branch("track_chi2",                             &trackChi2);
      tree->Branch("track_ndof",                             &trackNdof);
      tree->Branch("track_is_split",                         &trackIsSplit);
      tree->Branch("track_is_normal_trk_muon",           &trackIsNormalTrkMuon);
      tree->Branch("track_is_inout_trk_muon",            &trackIsInOutTrkMuon);
      tree->Branch("track_is_outin_trk_muon",            &trackIsOutInTrkMuon);
      tree->Branch("track_is_trk_muon",                  &trackIsTrkMuon);
      tree->Branch("track_is_trk_muon_with_arb",          &trackIsTrkMuonWithArb);

      tree->Branch("rpc_hit_rawid",     &rpcHitRawId);
      tree->Branch("rpc_hit_pos_x",     &rpcHitPosX);
      tree->Branch("rpc_hit_pos_y",     &rpcHitPosY);
      tree->Branch("rpc_hit_pos_z",     &rpcHitPosZ);
      tree->Branch("rpc_hit_pos_x_err", &rpcHitPosXErr);
      tree->Branch("rpc_hit_pos_y_err", &rpcHitPosYErr);
      tree->Branch("rpc_hit_cls",  &rpcHitCls);
      tree->Branch("rpc_hit_bx",        &rpcHitBX);

      tree->Branch("gem_hit_rawid",     &gemHitRawId);
      tree->Branch("gem_hit_pos_x",     &gemHitPosX);
      tree->Branch("gem_hit_pos_y",     &gemHitPosY);
      tree->Branch("gem_hit_pos_z",     &gemHitPosZ);
      tree->Branch("gem_hit_pos_x_err", &gemHitPosXErr);
      tree->Branch("gem_hit_pos_y_err", &gemHitPosYErr);
      tree->Branch("gem_hit_cls",  &gemHitCls);
      tree->Branch("gem_hit_bx",        &gemHitBX);

      tree->Branch("dt_seg_rawid",     &dtSegRawId);
      tree->Branch("dt_seg_pos_x",     &dtSegPosX);
      tree->Branch("dt_seg_pos_y",     &dtSegPosY);
      tree->Branch("dt_seg_pos_z",     &dtSegPosZ);
      tree->Branch("dt_seg_pos_x_err", &dtSegPosXErr);
      tree->Branch("dt_seg_pos_y_err", &dtSegPosYErr);
      tree->Branch("dt_seg_dir_x",     &dtSegDirX);
      tree->Branch("dt_seg_dir_y",     &dtSegDirY);
      tree->Branch("dt_seg_dir_z",     &dtSegDirZ);
      tree->Branch("dt_seg_dir_x_err", &dtSegDirXErr);
      tree->Branch("dt_seg_dir_y_err", &dtSegDirYErr);
      tree->Branch("dt_seg_chi2",      &dtSegChi2);
      tree->Branch("dt_seg_ndof",       &dtSegNdof);

      tree->Branch("csc_seg_rawid",     &cscSegRawId);
      tree->Branch("csc_seg_pos_x",     &cscSegPosX);
      tree->Branch("csc_seg_pos_y",     &cscSegPosY);
      tree->Branch("csc_seg_pos_z",     &cscSegPosZ);
      tree->Branch("csc_seg_pos_x_err", &cscSegPosXErr);
      tree->Branch("csc_seg_pos_y_err", &cscSegPosYErr);
      tree->Branch("csc_seg_dir_x",     &cscSegDirX);
      tree->Branch("csc_seg_dir_y",     &cscSegDirY);
      tree->Branch("csc_seg_dir_z",     &cscSegDirZ);
      tree->Branch("csc_seg_dir_x_err", &cscSegDirXErr);
      tree->Branch("csc_seg_dir_y_err", &cscSegDirYErr);
      tree->Branch("csc_seg_chi2",      &cscSegChi2);
      tree->Branch("csc_seg_ndof",       &cscSegNdof);
      tree->Branch("csc_seg_time",      &cscSegTime);
    }

    void clear() {
      trackKey.clear();
      trackVx.clear();
      trackVy.clear();
      trackVz.clear();
      trackPx.clear();
      trackPy.clear();
      trackPz.clear();
      trackQOverP.clear();
      trackLambda.clear();
      trackPhi.clear();
      trackDxy.clear();
      trackDsz.clear();
      trackQOverPErr.clear();
      trackLambdaErr.clear();
      trackPhiErr.clear();
      trackDxyErr.clear();
      trackDszErr.clear();
      trackCharge.clear();
      trackChi2.clear();
      trackNdof.clear();
      trackIsSplit.clear();
      trackIsNormalTrkMuon.clear();
      trackIsInOutTrkMuon.clear();
      trackIsOutInTrkMuon.clear();
      trackIsTrkMuon.clear();
      trackIsTrkMuonWithArb.clear();

      rpcHitRawId.clear();
      rpcHitPosX.clear();
      rpcHitPosY.clear();
      rpcHitPosZ.clear();
      rpcHitPosXErr.clear();
      rpcHitPosYErr.clear();
      rpcHitCls.clear();
      rpcHitBX.clear();

      gemHitRawId.clear();
      gemHitPosX.clear();
      gemHitPosY.clear();
      gemHitPosZ.clear();
      gemHitPosXErr.clear();
      gemHitPosYErr.clear();
      gemHitCls.clear();
      gemHitBX.clear();

      dtSegRawId.clear();
      dtSegPosX.clear();
      dtSegPosY.clear();
      dtSegPosZ.clear();
      dtSegPosXErr.clear();
      dtSegPosYErr.clear();
      dtSegDirX.clear();
      dtSegDirY.clear();
      dtSegDirZ.clear();
      dtSegDirXErr.clear();
      dtSegDirYErr.clear();
      dtSegChi2.clear();
      dtSegNdof.clear();

      cscSegRawId.clear();
      cscSegPosX.clear();
      cscSegPosY.clear();
      cscSegPosZ.clear();
      cscSegPosXErr.clear();
      cscSegPosYErr.clear();
      cscSegDirX.clear();
      cscSegDirY.clear();
      cscSegDirZ.clear();
      cscSegDirXErr.clear();
      cscSegDirYErr.clear();
      cscSegChi2.clear();
      cscSegNdof.clear();
      cscSegTime.clear();
    }
  };
  
  edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeomToken_;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord>  dtGeomToken_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeomToken_;
  DeepMuonRecoStruct DeepMuonRecoStruct_;
  TTree* TTree_;
};
#endif
