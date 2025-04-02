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

  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentToken_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegmentToken_;
  edm::EDGetTokenT<GEMSegmentCollection> gemSegmentToken_;
  edm::EDGetTokenT<RPCRecHitCollection> rpcHitToken_;
  edm::EDGetTokenT<GEMRecHitCollection> gemHitToken_;
  edm::EDGetTokenT<edm::ValueMap<reco::MuonQuality> > glbQualToken_;

  edm::Handle<DTRecSegment4DCollection> dtSegmentHandle_;
  edm::Handle<CSCSegmentCollection> cscSegmentHandle_;
  edm::Handle<GEMSegmentCollection> gemSegmentHandle_;
  edm::Handle<RPCRecHitCollection> rpcHitHandle_;
  edm::Handle<GEMRecHitCollection> gemHitHandle_;
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

  struct ProfileInfo {
    int nCallMuonIdProducer;
    int nInnerTrack;
    int nGoodInnerTrack;
    int nCallFillMuonIdTrackerMuon;
    int nCallFillMuonIdStandAloneMuon;
    int nCallAssociate;
    double timeMuonIdProducer;
    double timeAssociate;
  };

  TTree* profileTree_;
  ProfileInfo profileInfo_;

  struct TrackerMuonInfo {
    std::vector<double> trackPt;
    std::vector<double> trackP;
    std::vector<double> trackEta;
    std::vector<double> trackPhi;

    std::vector<bool>   isGoodTrackerMuon;
    std::vector<bool>   isGoodRPCMuon;
    std::vector<bool>   isGoodGEMMuon;
    std::vector<bool>   isGoodME0Muon;
    std::vector<bool>   isOutputMuon;
  };

  TTree* trackerMuonTree_;
  TrackerMuonInfo trackerMuonInfo_;
  
  struct MuonHitSegInfo {
    std::vector<unsigned int>   rpcRecHitRawId;
    std::vector<float>          rpcRecHitPosX;
    std::vector<float>          rpcRecHitPosY;
    std::vector<float>          rpcRecHitPosZ;

    std::vector<unsigned int>   gemRecHitRawId;
    std::vector<float>          gemRecHitPosX;
    std::vector<float>          gemRecHitPosY;
    std::vector<float>          gemRecHitPosZ;

    std::vector<unsigned int>   dtSegmentRawId;
    std::vector<float>          dtSegmentPosX;
    std::vector<float>          dtSegmentPosY;
    std::vector<float>          dtSegmentPosZ;
    std::vector<float>          dtSegmentDirX;
    std::vector<float>          dtSegmentDirY;
    std::vector<float>          dtSegmentDirZ;

    std::vector<unsigned int>   cscSegmentRawId;
    std::vector<float>          cscSegmentPosX;
    std::vector<float>          cscSegmentPosY;
    std::vector<float>          cscSegmentPosZ;
    std::vector<float>          cscSegmentDirX;
    std::vector<float>          cscSegmentDirY;
    std::vector<float>          cscSegmentDirZ;

    std::vector<unsigned int>   gemSegmentRawId;
    std::vector<float>          gemSegmentPosX;
    std::vector<float>          gemSegmentPosY;
    std::vector<float>          gemSegmentPosZ;
    std::vector<float>          gemSegmentDirX;
    std::vector<float>          gemSegmentDirY;
    std::vector<float>          gemSegmentDirZ;
  };

  TTree* muonHitSegTree_;
  MuonHitSegInfo muonHitSegInfo_;
  
  edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeomToken_;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord>  dtGeomToken_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord>  gemGeomToken_;
};
#endif