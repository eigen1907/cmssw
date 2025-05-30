#include "Validation/RecoMuon/plugins/MuonTrackValidator.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "SimTracker/TrackAssociation/interface/TrackingParticleIP.h"

#include "TMath.h"
#include <optional>

using namespace std;
using namespace edm;

void MuonTrackValidator::bookHistograms(DQMEDAnalyzer::DQMStore::IBooker& ibooker,
                                        edm::Run const&,
                                        edm::EventSetup const& setup) {
  // No booking if there are no tracks in a collection
  if (label.empty()) {
    return;
  }

  const auto minColl = -0.5;
  const auto maxColl = label.size() - 0.5;
  const auto nintColl = label.size();
  edm::LogVerbatim("MuonTrackValidator") << "Label size: " << label.size() << "\n"
                                         << "Associator size: " << associators.size() << "\n"
                                         << "HistoParameters size: " << histoParameters.size();

  auto binLabels = [&](dqm::reco::MonitorElement* me) -> dqm::reco::MonitorElement* {
    for (size_t i = 0; i < label.size(); ++i) {
      std::string labelName =
          label[i].instance().empty() ? label[i].label() : label[i].label() + "_" + label[i].instance();
      me->setBinLabel(i + 1, labelName);
    }
    me->disableAlphanumeric();
    return me;
  };

  size_t outerLoopSize = UseAssociators ? associators.size() : 1;
  for (unsigned int ww = 0; ww < outerLoopSize; ++ww) {
    if (doSummaryPlots_) {
      ibooker.setCurrentFolder(dirName_);
      h_assoc_coll.push_back(binLabels(ibooker.book1D("num_asso_SimToReco_coll",
                                                      "N of associated (SimToReco) muons vs muon collection",
                                                      nintColl,
                                                      minColl,
                                                      maxColl)));
      h_simul_coll.push_back(binLabels(
          ibooker.book1D("num_simul_coll", "N of simulated muons vs muon collection", nintColl, minColl, maxColl)));
      h_reco_coll.push_back(
          binLabels(ibooker.book1D("num_reco_coll", "N of reco muons vs muon collection", nintColl, minColl, maxColl)));
      h_assoc2_coll.push_back(binLabels(ibooker.book1D("num_asso_RecoToSim_coll",
                                                       "N of associated (recoToSim) muons vs muon collection",
                                                       nintColl,
                                                       minColl,
                                                       maxColl)));
    }
    for (unsigned int www = 0; www < label.size(); ++www) {
      edm::LogVerbatim("MuonTrackValidator")
          << " Label " << www << ", Associator: " << associators[www] << '\n'
          << "Track Label: " << label[www].label() << ", instance: " << label[www].instance()
          << ", process name: " << label[www].process() << "\n";

      ibooker.cd();
      InputTag algo = label[www];
      string dirName = dirName_;

      auto setBinLogX = [&](TH1* th1) {
        if (histoParameters[www].useLogPt) {
          BinLogX(th1);
        }
      };

      if (!algo.process().empty())
        dirName += algo.process() + "_";
      if (!algo.label().empty())
        dirName += algo.label();
      if (!algo.instance().empty())
        dirName += ("_" + algo.instance());
      if (dirName.find("Tracks") < dirName.length()) {
        dirName.replace(dirName.find("Tracks"), 6, "Trks");
      }
      if (dirName.find("UpdatedAtVtx") < dirName.length()) {
        dirName.replace(dirName.find("UpdatedAtVtx"), 12, "UpdAtVtx");
      }
      if (associators[www].find("tpToTkmuTrackAssociation") != std::string::npos) {
        dirName += "_TkAsso";
      }
      std::replace(dirName.begin(), dirName.end(), ':', '_');
      ibooker.setCurrentFolder(dirName);

      h_tracks.push_back(ibooker.book1D("Ntracks",
                                        "Number of reconstructed tracks",
                                        histoParameters[www].nintNTracks,
                                        histoParameters[www].minNTracks,
                                        histoParameters[www].maxNTracks));
      h_fakes.push_back(ibooker.book1D("Nfakes",
                                       "Number of fake reco tracks",
                                       histoParameters[www].nintFTracks,
                                       histoParameters[www].minFTracks,
                                       histoParameters[www].maxFTracks));
      h_charge.push_back(ibooker.book1D("Ncharge", "track charge", 3, -1.5, 1.5));

      h_recoeta.push_back(ibooker.book1D("num_reco_eta",
                                         "N of reco track vs eta",
                                         histoParameters[www].nintEta,
                                         histoParameters[www].minEta,
                                         histoParameters[www].maxEta));
      h_assoceta.push_back(ibooker.book1D("num_assoSimToReco_eta",
                                          "N of associated tracks (simToReco) vs eta",
                                          histoParameters[www].nintEta,
                                          histoParameters[www].minEta,
                                          histoParameters[www].maxEta));
      h_assoc2eta.push_back(ibooker.book1D("num_assoRecoToSim_eta",
                                           "N of associated (recoToSim) tracks vs eta",
                                           histoParameters[www].nintEta,
                                           histoParameters[www].minEta,
                                           histoParameters[www].maxEta));
      h_simuleta.push_back(ibooker.book1D("num_simul_eta",
                                          "N of simulated tracks vs eta",
                                          histoParameters[www].nintEta,
                                          histoParameters[www].minEta,
                                          histoParameters[www].maxEta));
      h_misideta.push_back(ibooker.book1D("num_chargemisid_eta",
                                          "N of associated (simToReco) tracks with charge misID vs eta",
                                          histoParameters[www].nintEta,
                                          histoParameters[www].minEta,
                                          histoParameters[www].maxEta));

      h_recopT.push_back(ibooker.book1D("num_reco_pT",
                                        "N of reco track vs pT",
                                        histoParameters[www].nintPt,
                                        histoParameters[www].minPt,
                                        histoParameters[www].maxPt,
                                        setBinLogX));
      h_assocpT.push_back(ibooker.book1D("num_assoSimToReco_pT",
                                         "N of associated tracks (simToReco) vs pT",
                                         histoParameters[www].nintPt,
                                         histoParameters[www].minPt,
                                         histoParameters[www].maxPt,
                                         setBinLogX));
      h_assoc2pT.push_back(ibooker.book1D("num_assoRecoToSim_pT",
                                          "N of associated (recoToSim) tracks vs pT",
                                          histoParameters[www].nintPt,
                                          histoParameters[www].minPt,
                                          histoParameters[www].maxPt,
                                          setBinLogX));
      h_simulpT.push_back(ibooker.book1D("num_simul_pT",
                                         "N of simulated tracks vs pT",
                                         histoParameters[www].nintPt,
                                         histoParameters[www].minPt,
                                         histoParameters[www].maxPt,
                                         setBinLogX));
      h_misidpT.push_back(ibooker.book1D("num_chargemisid_pT",
                                         "N of associated (simToReco) tracks with charge misID vs pT",
                                         histoParameters[www].nintPt,
                                         histoParameters[www].minPt,
                                         histoParameters[www].maxPt,
                                         setBinLogX));

      h_assocpTB.push_back(ibooker.book1D("num_assoSimToReco_pT_barrel",
                                          "N of associated tracks (simToReco) vs pT - BARREL",
                                          histoParameters[www].nintPt,
                                          histoParameters[www].minPt,
                                          histoParameters[www].maxPt,
                                          setBinLogX));
      h_simulpTB.push_back(ibooker.book1D("num_simul_pT_barrel",
                                          "N of simulated tracks vs pT - BARREL",
                                          histoParameters[www].nintPt,
                                          histoParameters[www].minPt,
                                          histoParameters[www].maxPt,
                                          setBinLogX));

      h_assocpTO.push_back(ibooker.book1D("num_assoSimToReco_pT_overlap",
                                          "N of associated tracks (simToReco) vs pT - OVERLAP",
                                          histoParameters[www].nintPt,
                                          histoParameters[www].minPt,
                                          histoParameters[www].maxPt,
                                          setBinLogX));
      h_simulpTO.push_back(ibooker.book1D("num_simul_pT_overlap",
                                          "N of simulated tracks vs pT - OVERLAP",
                                          histoParameters[www].nintPt,
                                          histoParameters[www].minPt,
                                          histoParameters[www].maxPt,
                                          setBinLogX));

      h_assocpTE.push_back(ibooker.book1D("num_assoSimToReco_pT_endcap",
                                          "N of associated tracks (simToReco) vs pT - ENCAP",
                                          histoParameters[www].nintPt,
                                          histoParameters[www].minPt,
                                          histoParameters[www].maxPt,
                                          setBinLogX));
      h_simulpTE.push_back(ibooker.book1D("num_simul_pT_endcap",
                                          "N of simulated tracks vs pT - ENDCAP",
                                          histoParameters[www].nintPt,
                                          histoParameters[www].minPt,
                                          histoParameters[www].maxPt,
                                          setBinLogX));

      h_recophi.push_back(ibooker.book1D("num_reco_phi",
                                         "N of reco track vs phi",
                                         histoParameters[www].nintPhi,
                                         histoParameters[www].minPhi,
                                         histoParameters[www].maxPhi));
      h_assocphi.push_back(ibooker.book1D("num_assoSimToReco_phi",
                                          "N of associated tracks (simToReco) vs phi",
                                          histoParameters[www].nintPhi,
                                          histoParameters[www].minPhi,
                                          histoParameters[www].maxPhi));
      h_assoc2phi.push_back(ibooker.book1D("num_assoRecoToSim_phi",
                                           "N of associated (recoToSim) tracks vs phi",
                                           histoParameters[www].nintPhi,
                                           histoParameters[www].minPhi,
                                           histoParameters[www].maxPhi));
      h_simulphi.push_back(ibooker.book1D("num_simul_phi",
                                          "N of simulated tracks vs phi",
                                          histoParameters[www].nintPhi,
                                          histoParameters[www].minPhi,
                                          histoParameters[www].maxPhi));
      h_misidphi.push_back(ibooker.book1D("num_chargemisid_phi",
                                          "N of associated (simToReco) tracks with charge misID vs phi",
                                          histoParameters[www].nintPhi,
                                          histoParameters[www].minPhi,
                                          histoParameters[www].maxPhi));

      h_assocphiB.push_back(ibooker.book1D("num_assoSimToReco_phi_barrel",
                                           "N of associated tracks (simToReco) vs phi - BARREL",
                                           histoParameters[www].nintPhi,
                                           histoParameters[www].minPhi,
                                           histoParameters[www].maxPhi));
      h_simulphiB.push_back(ibooker.book1D("num_simul_phi_barrel",
                                           "N of simulated track vs phi - BARREL",
                                           histoParameters[www].nintPhi,
                                           histoParameters[www].minPhi,
                                           histoParameters[www].maxPhi));

      h_assocphiO.push_back(ibooker.book1D("num_assoSimToReco_phi_overlap",
                                           "N of associated tracks (simToReco) vs phi - OVERLAP",
                                           histoParameters[www].nintPhi,
                                           histoParameters[www].minPhi,
                                           histoParameters[www].maxPhi));
      h_simulphiO.push_back(ibooker.book1D("num_simul_phi_overlap",
                                           "N of simulated track vs phi - OVERLAP",
                                           histoParameters[www].nintPhi,
                                           histoParameters[www].minPhi,
                                           histoParameters[www].maxPhi));

      h_assocphiE.push_back(ibooker.book1D("num_assoSimToReco_phi_endcap",
                                           "N of associated tracks (simToReco) vs phi - ENDCAP",
                                           histoParameters[www].nintPhi,
                                           histoParameters[www].minPhi,
                                           histoParameters[www].maxPhi));
      h_simulphiE.push_back(ibooker.book1D("num_simul_phi_endcap",
                                           "N of simulated track vs phi - ENDCAP",
                                           histoParameters[www].nintPhi,
                                           histoParameters[www].minPhi,
                                           histoParameters[www].maxPhi));

      h_recohit.push_back(ibooker.book1D("num_reco_hit",
                                         "N of reco tracks vs N SimHits",
                                         histoParameters[www].nintNHit,
                                         histoParameters[www].minNHit,
                                         histoParameters[www].maxNHit));
      h_assochit.push_back(ibooker.book1D("num_assoSimToReco_hit",
                                          "N of associated tracks (simToReco) vs N SimHits",
                                          histoParameters[www].nintNHit,
                                          histoParameters[www].minNHit,
                                          histoParameters[www].maxNHit));
      h_assoc2hit.push_back(ibooker.book1D("num_assoRecoToSim_hit",
                                           "N of associated (recoToSim) tracks vs N Rechits",
                                           histoParameters[www].nintNHit,
                                           histoParameters[www].minNHit,
                                           histoParameters[www].maxNHit));
      h_simulhit.push_back(ibooker.book1D("num_simul_hit",
                                          "N of simulated tracks vs N SimHits",
                                          histoParameters[www].nintNHit,
                                          histoParameters[www].minNHit,
                                          histoParameters[www].maxNHit));
      h_misidhit.push_back(ibooker.book1D("num_chargemisid_hit",
                                          "N of associated (recoToSim) tracks with charge misID vs N RecHits",
                                          histoParameters[www].nintNHit,
                                          histoParameters[www].minNHit,
                                          histoParameters[www].maxNHit));

      h_recodR.push_back(ibooker.book1D("num_reco_dR",
                                        "N of reco track vs dR",
                                        histoParameters[www].nintdR,
                                        histoParameters[www].mindR,
                                        histoParameters[www].maxdR));
      h_assocdR.push_back(ibooker.book1D("num_assoSimToReco_dR",
                                         "N of associated tracks (simToReco) vs dR",
                                         histoParameters[www].nintdR,
                                         histoParameters[www].mindR,
                                         histoParameters[www].maxdR));
      h_assoc2dR.push_back(ibooker.book1D("num_assoRecoToSim_dR",
                                          "N of associated (recoToSim) tracks vs dR",
                                          histoParameters[www].nintdR,
                                          histoParameters[www].mindR,
                                          histoParameters[www].maxdR));
      h_simuldR.push_back(ibooker.book1D("num_simul_dR",
                                         "N of simulated tracks vs dR",
                                         histoParameters[www].nintdR,
                                         histoParameters[www].mindR,
                                         histoParameters[www].maxdR));
      h_misiddR.push_back(ibooker.book1D("num_chargemisid_dR",
                                         "N of associated (simToReco) tracks with charge misID vs dR",
                                         histoParameters[www].nintdR,
                                         histoParameters[www].mindR,
                                         histoParameters[www].maxdR));

      h_recodxy.push_back(ibooker.book1D("num_reco_dxy",
                                         "N of reco track vs dxy",
                                         histoParameters[www].nintDxy,
                                         histoParameters[www].minDxy,
                                         histoParameters[www].maxDxy));
      h_assocdxy.push_back(ibooker.book1D("num_assoSimToReco_dxy",
                                          "N of associated tracks (simToReco) vs dxy",
                                          histoParameters[www].nintDxy,
                                          histoParameters[www].minDxy,
                                          histoParameters[www].maxDxy));
      h_assoc2dxy.push_back(ibooker.book1D("num_assoRecoToSim_dxy",
                                           "N of associated (recoToSim) tracks vs dxy",
                                           histoParameters[www].nintDxy,
                                           histoParameters[www].minDxy,
                                           histoParameters[www].maxDxy));
      h_simuldxy.push_back(ibooker.book1D("num_simul_dxy",
                                          "N of simulated tracks vs dxy",
                                          histoParameters[www].nintDxy,
                                          histoParameters[www].minDxy,
                                          histoParameters[www].maxDxy));
      h_misiddxy.push_back(ibooker.book1D("num_chargemisid_dxy",
                                          "N of associated (simToReco) tracks with charge misID vs dxy",
                                          histoParameters[www].nintDxy,
                                          histoParameters[www].minDxy,
                                          histoParameters[www].maxDxy));
      h_recodz.push_back(ibooker.book1D("num_reco_dz",
                                        "N of reco track vs dz",
                                        histoParameters[www].nintDz,
                                        histoParameters[www].minDz,
                                        histoParameters[www].maxDz));
      h_assocdz.push_back(ibooker.book1D("num_assoSimToReco_dz",
                                         "N of associated tracks (simToReco) vs dz",
                                         histoParameters[www].nintDz,
                                         histoParameters[www].minDz,
                                         histoParameters[www].maxDz));
      h_assoc2dz.push_back(ibooker.book1D("num_assoRecoToSim_dz",
                                          "N of associated (recoToSim) tracks vs dz",
                                          histoParameters[www].nintDz,
                                          histoParameters[www].minDz,
                                          histoParameters[www].maxDz));
      h_simuldz.push_back(ibooker.book1D("num_simul_dz",
                                         "N of simulated tracks vs dz",
                                         histoParameters[www].nintDz,
                                         histoParameters[www].minDz,
                                         histoParameters[www].maxDz));
      h_misiddz.push_back(ibooker.book1D("num_chargemisid_dz",
                                         "N of associated (simToReco) tracks with charge misID vs dz",
                                         histoParameters[www].nintDz,
                                         histoParameters[www].minDz,
                                         histoParameters[www].maxDz));

      h_assocRpos.push_back(ibooker.book1D("num_assoSimToReco_Rpos",
                                           "N of associated tracks (simToReco) vs Radius",
                                           histoParameters[www].nintRpos,
                                           histoParameters[www].minRpos,
                                           histoParameters[www].maxRpos));
      h_simulRpos.push_back(ibooker.book1D("num_simul_Rpos",
                                           "N of simulated tracks vs Radius",
                                           histoParameters[www].nintRpos,
                                           histoParameters[www].minRpos,
                                           histoParameters[www].maxRpos));

      h_assocZpos.push_back(ibooker.book1D("num_assoSimToReco_Zpos",
                                           "N of associated tracks (simToReco) vs Z",
                                           histoParameters[www].nintZpos,
                                           histoParameters[www].minZpos,
                                           histoParameters[www].maxZpos));
      h_simulZpos.push_back(ibooker.book1D("num_simul_Zpos",
                                           "N of simulated tracks vs Z",
                                           histoParameters[www].nintZpos,
                                           histoParameters[www].minZpos,
                                           histoParameters[www].maxZpos));

      h_recopu.push_back(ibooker.book1D("num_reco_pu",
                                        "N of reco track vs pu",
                                        histoParameters[www].nintPU,
                                        histoParameters[www].minPU,
                                        histoParameters[www].maxPU));
      h_assocpu.push_back(ibooker.book1D("num_assoSimToReco_pu",
                                         "N of associated tracks (simToReco) vs pu",
                                         histoParameters[www].nintPU,
                                         histoParameters[www].minPU,
                                         histoParameters[www].maxPU));
      h_assoc2pu.push_back(ibooker.book1D("num_assoRecoToSim_pu",
                                          "N of associated (recoToSim) tracks vs pu",
                                          histoParameters[www].nintPU,
                                          histoParameters[www].minPU,
                                          histoParameters[www].maxPU));
      h_simulpu.push_back(ibooker.book1D("num_simul_pu",
                                         "N of simulated tracks vs pu",
                                         histoParameters[www].nintPU,
                                         histoParameters[www].minPU,
                                         histoParameters[www].maxPU));
      h_misidpu.push_back(ibooker.book1D("num_chargemisid_pu",
                                         "N of associated (simToReco) charge misIDed tracks vs pu",
                                         histoParameters[www].nintPU,
                                         histoParameters[www].minPU,
                                         histoParameters[www].maxPU));

      h_nchi2.push_back(ibooker.book1D("chi2", "Track normalized #chi^{2}", 80, 0., 20.));
      h_nchi2_prob.push_back(ibooker.book1D("chi2prob", "Probability of track normalized #chi^{2}", 100, 0., 1.));

      chi2_vs_nhits.push_back(ibooker.book2D("chi2_vs_nhits",
                                             "#chi^{2} vs nhits",
                                             histoParameters[www].nintNHit,
                                             histoParameters[www].minNHit,
                                             histoParameters[www].maxNHit,
                                             20,
                                             0.,
                                             10.));
      chi2_vs_eta.push_back(ibooker.book2D("chi2_vs_eta",
                                           "chi2_vs_eta",
                                           histoParameters[www].nintEta,
                                           histoParameters[www].minEta,
                                           histoParameters[www].maxEta,
                                           40,
                                           0.,
                                           20.));
      chi2_vs_phi.push_back(ibooker.book2D("chi2_vs_phi",
                                           "#chi^{2} vs #phi",
                                           histoParameters[www].nintPhi,
                                           histoParameters[www].minPhi,
                                           histoParameters[www].maxPhi,
                                           40,
                                           0.,
                                           20.));

      h_nhits.push_back(ibooker.book1D("nhits",
                                       "Number of hits per track",
                                       histoParameters[www].nintNHit,
                                       histoParameters[www].minNHit,
                                       histoParameters[www].maxNHit));
      nhits_vs_eta.push_back(ibooker.book2D("nhits_vs_eta",
                                            "Number of Hits vs eta",
                                            histoParameters[www].nintEta,
                                            histoParameters[www].minEta,
                                            histoParameters[www].maxEta,
                                            histoParameters[www].nintNHit,
                                            histoParameters[www].minNHit,
                                            histoParameters[www].maxNHit));
      nhits_vs_phi.push_back(ibooker.book2D("nhits_vs_phi",
                                            "#hits vs #phi",
                                            histoParameters[www].nintPhi,
                                            histoParameters[www].minPhi,
                                            histoParameters[www].maxPhi,
                                            histoParameters[www].nintNHit,
                                            histoParameters[www].minNHit,
                                            histoParameters[www].maxNHit));

      if (histoParameters[www].do_MUOhitsPlots) {
        nDThits_vs_eta.push_back(ibooker.book2D("nDThits_vs_eta",
                                                "Number of DT hits vs eta",
                                                histoParameters[www].nintEta,
                                                histoParameters[www].minEta,
                                                histoParameters[www].maxEta,
                                                histoParameters[www].nintDTHit,
                                                histoParameters[www].minDTHit,
                                                histoParameters[www].maxDTHit));
        nCSChits_vs_eta.push_back(ibooker.book2D("nCSChits_vs_eta",
                                                 "Number of CSC hits vs eta",
                                                 histoParameters[www].nintEta,
                                                 histoParameters[www].minEta,
                                                 histoParameters[www].maxEta,
                                                 histoParameters[www].nintCSCHit,
                                                 histoParameters[www].minCSCHit,
                                                 histoParameters[www].maxCSCHit));
        nRPChits_vs_eta.push_back(ibooker.book2D("nRPChits_vs_eta",
                                                 "Number of RPC hits vs eta",
                                                 histoParameters[www].nintEta,
                                                 histoParameters[www].minEta,
                                                 histoParameters[www].maxEta,
                                                 histoParameters[www].nintRPCHit,
                                                 histoParameters[www].minRPCHit,
                                                 histoParameters[www].maxRPCHit));
        if (useGEMs_)
          nGEMhits_vs_eta.push_back(ibooker.book2D("nGEMhits_vs_eta",
                                                   "Number of GEM hits vs eta",
                                                   histoParameters[www].nintEta,
                                                   histoParameters[www].minEta,
                                                   histoParameters[www].maxEta,
                                                   histoParameters[www].nintNHit,
                                                   histoParameters[www].minNHit,
                                                   histoParameters[www].maxNHit));
        if (useME0_)
          nME0hits_vs_eta.push_back(ibooker.book2D("nME0hits_vs_eta",
                                                   "Number of ME0 hits vs eta",
                                                   histoParameters[www].nintEta,
                                                   histoParameters[www].minEta,
                                                   histoParameters[www].maxEta,
                                                   histoParameters[www].nintNHit,
                                                   histoParameters[www].minNHit,
                                                   histoParameters[www].maxNHit));
      }

      if (histoParameters[www].do_TRKhitsPlots) {
        nTRK_LayersWithMeas_vs_eta.push_back(ibooker.book2D("nTRK_LayersWithMeas_vs_eta",
                                                            "# TRK Layers with measurement vs eta",
                                                            histoParameters[www].nintEta,
                                                            histoParameters[www].minEta,
                                                            histoParameters[www].maxEta,
                                                            histoParameters[www].nintLayers,
                                                            histoParameters[www].minLayers,
                                                            histoParameters[www].maxLayers));
        nPixel_LayersWithMeas_vs_eta.push_back(ibooker.book2D("nPixel_LayersWithMeas_vs_eta",
                                                              "Number of Pixel Layers with measurement vs eta",
                                                              histoParameters[www].nintEta,
                                                              histoParameters[www].minEta,
                                                              histoParameters[www].maxEta,
                                                              histoParameters[www].nintPixels,
                                                              histoParameters[www].minPixels,
                                                              histoParameters[www].maxPixels));
        h_nmisslayers_inner.push_back(ibooker.book1D("nTRK_misslayers_inner",
                                                     "Number of missing inner TRK layers",
                                                     histoParameters[www].nintLayers,
                                                     histoParameters[www].minLayers,
                                                     histoParameters[www].maxLayers));
        h_nmisslayers_outer.push_back(ibooker.book1D("nTRK_misslayers_outer",
                                                     "Number of missing outer TRK layers",
                                                     histoParameters[www].nintLayers,
                                                     histoParameters[www].minLayers,
                                                     histoParameters[www].maxLayers));
        h_nlosthits.push_back(ibooker.book1D("nlosthits", "Number of lost hits per track", 6, -0.5, 5.5));
        nlosthits_vs_eta.push_back(ibooker.book2D("nlosthits_vs_eta",
                                                  "Number of lost hits per track vs eta",
                                                  histoParameters[www].nintEta,
                                                  histoParameters[www].minEta,
                                                  histoParameters[www].maxEta,
                                                  6,
                                                  -0.5,
                                                  5.5));
      }

      ptres_vs_eta.push_back(ibooker.book2D("ptres_vs_eta",
                                            "p_{T} Relative Residual vs #eta",
                                            histoParameters[www].nintEta,
                                            histoParameters[www].minEta,
                                            histoParameters[www].maxEta,
                                            histoParameters[www].ptRes_nbin,
                                            histoParameters[www].ptRes_rangeMin,
                                            histoParameters[www].ptRes_rangeMax));
      ptres_vs_phi.push_back(ibooker.book2D("ptres_vs_phi",
                                            "p_{T} Relative Residual vs #phi",
                                            histoParameters[www].nintPhi,
                                            histoParameters[www].minPhi,
                                            histoParameters[www].maxPhi,
                                            histoParameters[www].ptRes_nbin,
                                            histoParameters[www].ptRes_rangeMin,
                                            histoParameters[www].ptRes_rangeMax));
      ptres_vs_pt.push_back(ibooker.book2D("ptres_vs_pt",
                                           "p_{T} Relative Residual vs p_{T}",
                                           histoParameters[www].nintPt,
                                           histoParameters[www].minPt,
                                           histoParameters[www].maxPt,
                                           histoParameters[www].ptRes_nbin,
                                           histoParameters[www].ptRes_rangeMin,
                                           histoParameters[www].ptRes_rangeMax,
                                           setBinLogX));
      h_ptpull.push_back(ibooker.book1D("ptpull", "p_{T} Pull", 100, -10., 10.));
      ptpull_vs_eta.push_back(ibooker.book2D("ptpull_vs_eta",
                                             "p_{T} Pull vs #eta",
                                             histoParameters[www].nintEta,
                                             histoParameters[www].minEta,
                                             histoParameters[www].maxEta,
                                             100,
                                             -10.,
                                             10.));
      ptpull_vs_phi.push_back(ibooker.book2D("ptpull_vs_phi",
                                             "p_{T} Pull vs #phi",
                                             histoParameters[www].nintPhi,
                                             histoParameters[www].minPhi,
                                             histoParameters[www].maxPhi,
                                             100,
                                             -10.,
                                             10.));
      h_qoverppull.push_back(ibooker.book1D("qoverppull", "q/p Pull", 100, -10., 10.));

      h_etaRes.push_back(ibooker.book1D("etaRes",
                                        "#eta residual",
                                        histoParameters[www].etaRes_nbin,
                                        histoParameters[www].etaRes_rangeMin,
                                        histoParameters[www].etaRes_rangeMax));
      etares_vs_eta.push_back(ibooker.book2D("etares_vs_eta",
                                             "#eta Residual vs #eta",
                                             histoParameters[www].nintEta,
                                             histoParameters[www].minEta,
                                             histoParameters[www].maxEta,
                                             histoParameters[www].etaRes_nbin,
                                             histoParameters[www].etaRes_rangeMin,
                                             histoParameters[www].etaRes_rangeMax));

      thetaCotres_vs_eta.push_back(ibooker.book2D("thetaCotres_vs_eta",
                                                  "cot(#theta) Residual vs #eta",
                                                  histoParameters[www].nintEta,
                                                  histoParameters[www].minEta,
                                                  histoParameters[www].maxEta,
                                                  histoParameters[www].cotThetaRes_nbin,
                                                  histoParameters[www].cotThetaRes_rangeMin,
                                                  histoParameters[www].cotThetaRes_rangeMax));
      thetaCotres_vs_pt.push_back(ibooker.book2D("thetaCotres_vs_pt",
                                                 "cot(#theta) Residual vs p_{T}",
                                                 histoParameters[www].nintPt,
                                                 histoParameters[www].minPt,
                                                 histoParameters[www].maxPt,
                                                 histoParameters[www].cotThetaRes_nbin,
                                                 histoParameters[www].cotThetaRes_rangeMin,
                                                 histoParameters[www].cotThetaRes_rangeMax,
                                                 setBinLogX));
      h_thetapull.push_back(ibooker.book1D("thetapull", "#theta Pull", 100, -10., 10.));
      thetapull_vs_eta.push_back(ibooker.book2D("thetapull_vs_eta",
                                                "#theta Pull vs #eta",
                                                histoParameters[www].nintEta,
                                                histoParameters[www].minEta,
                                                histoParameters[www].maxEta,
                                                100,
                                                -10,
                                                10));
      thetapull_vs_phi.push_back(ibooker.book2D("thetapull_vs_phi",
                                                "#theta Pull vs #phi",
                                                histoParameters[www].nintPhi,
                                                histoParameters[www].minPhi,
                                                histoParameters[www].maxPhi,
                                                100,
                                                -10,
                                                10));

      phires_vs_eta.push_back(ibooker.book2D("phires_vs_eta",
                                             "#phi Residual vs #eta",
                                             histoParameters[www].nintEta,
                                             histoParameters[www].minEta,
                                             histoParameters[www].maxEta,
                                             histoParameters[www].phiRes_nbin,
                                             histoParameters[www].phiRes_rangeMin,
                                             histoParameters[www].phiRes_rangeMax));
      phires_vs_pt.push_back(ibooker.book2D("phires_vs_pt",
                                            "#phi Residual vs p_{T}",
                                            histoParameters[www].nintPt,
                                            histoParameters[www].minPt,
                                            histoParameters[www].maxPt,
                                            histoParameters[www].phiRes_nbin,
                                            histoParameters[www].phiRes_rangeMin,
                                            histoParameters[www].phiRes_rangeMax,
                                            setBinLogX));
      phires_vs_phi.push_back(ibooker.book2D("phires_vs_phi",
                                             "#phi Residual vs #phi",
                                             histoParameters[www].nintPhi,
                                             histoParameters[www].minPhi,
                                             histoParameters[www].maxPhi,
                                             histoParameters[www].phiRes_nbin,
                                             histoParameters[www].phiRes_rangeMin,
                                             histoParameters[www].phiRes_rangeMax));
      h_phipull.push_back(ibooker.book1D("phipull", "#phi Pull", 100, -10., 10.));
      phipull_vs_eta.push_back(ibooker.book2D("phipull_vs_eta",
                                              "#phi Pull vs #eta",
                                              histoParameters[www].nintEta,
                                              histoParameters[www].minEta,
                                              histoParameters[www].maxEta,
                                              100,
                                              -10,
                                              10));
      phipull_vs_phi.push_back(ibooker.book2D("phipull_vs_phi",
                                              "#phi Pull vs #phi",
                                              histoParameters[www].nintPhi,
                                              histoParameters[www].minPhi,
                                              histoParameters[www].maxPhi,
                                              100,
                                              -10,
                                              10));

      dxyres_vs_eta.push_back(ibooker.book2D("dxyres_vs_eta",
                                             "dxy Residual vs #eta",
                                             histoParameters[www].nintEta,
                                             histoParameters[www].minEta,
                                             histoParameters[www].maxEta,
                                             histoParameters[www].dxyRes_nbin,
                                             histoParameters[www].dxyRes_rangeMin,
                                             histoParameters[www].dxyRes_rangeMax));
      dxyres_vs_pt.push_back(ibooker.book2D("dxyres_vs_pt",
                                            "dxy Residual vs p_{T}",
                                            histoParameters[www].nintPt,
                                            histoParameters[www].minPt,
                                            histoParameters[www].maxPt,
                                            histoParameters[www].dxyRes_nbin,
                                            histoParameters[www].dxyRes_rangeMin,
                                            histoParameters[www].dxyRes_rangeMax,
                                            setBinLogX));
      h_dxypull.push_back(ibooker.book1D("dxypull", "dxy Pull", 100, -10., 10.));
      dxypull_vs_eta.push_back(ibooker.book2D("dxypull_vs_eta",
                                              "dxy Pull vs #eta",
                                              histoParameters[www].nintEta,
                                              histoParameters[www].minEta,
                                              histoParameters[www].maxEta,
                                              100,
                                              -10,
                                              10));

      dzres_vs_eta.push_back(ibooker.book2D("dzres_vs_eta",
                                            "dz Residual vs #eta",
                                            histoParameters[www].nintEta,
                                            histoParameters[www].minEta,
                                            histoParameters[www].maxEta,
                                            histoParameters[www].dzRes_nbin,
                                            histoParameters[www].dzRes_rangeMin,
                                            histoParameters[www].dzRes_rangeMax));
      dzres_vs_pt.push_back(ibooker.book2D("dzres_vs_pt",
                                           "dz Residual vs p_{T}",
                                           histoParameters[www].nintPt,
                                           histoParameters[www].minPt,
                                           histoParameters[www].maxPt,
                                           histoParameters[www].dzRes_nbin,
                                           histoParameters[www].dzRes_rangeMin,
                                           histoParameters[www].dzRes_rangeMax,
                                           setBinLogX));
      h_dzpull.push_back(ibooker.book1D("dzpull", "dz Pull", 100, -10., 10.));
      dzpull_vs_eta.push_back(ibooker.book2D("dzpull_vs_eta",
                                             "dz Pull vs #eta",
                                             histoParameters[www].nintEta,
                                             histoParameters[www].minEta,
                                             histoParameters[www].maxEta,
                                             100,
                                             -10,
                                             10));

      nRecHits_vs_nSimHits.push_back(ibooker.book2D("nRecHits_vs_nSimHits",
                                                    "nRecHits vs nSimHits",
                                                    histoParameters[www].nintNHit,
                                                    histoParameters[www].minNHit,
                                                    histoParameters[www].maxNHit,
                                                    histoParameters[www].nintNHit,
                                                    histoParameters[www].minNHit,
                                                    histoParameters[www].maxNHit));

      if (MABH) {
        h_PurityVsQuality.push_back(
            ibooker.book2D("PurityVsQuality", "Purity vs Quality (MABH)", 20, 0.01, 1.01, 20, 0.01, 1.01));
      }
      if (UseAssociators) {
        if (associators[ww] == "trackAssociatorByChi2") {
          h_assochi2.push_back(ibooker.book1D("assocChi2", "track association #chi^{2}", 1000, 0., 100.));
          h_assochi2_prob.push_back(
              ibooker.book1D("assocChi2_prob", "probability of association #chi^{2}", 100, 0., 1.));
        } else if (associators[ww] == "trackAssociatorByHits") {
          h_assocFraction.push_back(ibooker.book1D("assocFraction", "fraction of shared hits", 22, 0., 1.1));
          h_assocSharedHit.push_back(ibooker.book1D("assocSharedHit", "number of shared hits", 41, -0.5, 40.5));
        }
      }

    }  //for (unsigned int www=0;www<label.size();www++)
  }  //for (unsigned int ww=0;ww<associators.size();ww++)
}

void MuonTrackValidator::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  using namespace reco;

  edm::LogInfo("MuonTrackValidator") << "\n===================================================="
                                     << "\n"
                                     << "Analyzing new event"
                                     << "\n"
                                     << "====================================================\n"
                                     << "\n";

  edm::Handle<std::vector<PileupSummaryInfo> > puinfoH;
  int PU_NumInteractions(-1);

  if (parametersDefiner == "LhcParametersDefinerForTP") {
    // PileupSummaryInfo is contained only in collision events
    event.getByToken(pileupinfo_Token, puinfoH);
    for (std::vector<PileupSummaryInfo>::const_iterator puInfoIt = puinfoH->begin(); puInfoIt != puinfoH->end();
         ++puInfoIt) {
      if (puInfoIt->getBunchCrossing() == 0) {
        PU_NumInteractions = puInfoIt->getPU_NumInteractions();
        break;
      }
    }

  } else if (parametersDefiner == "CosmicParametersDefinerForTP") {
    edm::Handle<SimHitTPAssociationProducer::SimHitTPAssociationList> simHitsTPAssoc;
    //warning: make sure the TP collection used in the map is the same used here
    event.getByToken(_simHitTpMapTag, simHitsTPAssoc);
    cosmicParametersDefinerTP_->initEvent(simHitsTPAssoc);
    cosmictpSelector.initEvent(simHitsTPAssoc);
  }

  TrackingParticleRefVector TPrefV;
  const TrackingParticleRefVector* ptr_TPrefV = nullptr;
  edm::Handle<TrackingParticleCollection> TPCollection_H;
  edm::Handle<TrackingParticleRefVector> TPCollectionRefVector_H;

  if (label_tp_refvector) {
    event.getByToken(tp_refvector_Token, TPCollectionRefVector_H);
    ptr_TPrefV = TPCollectionRefVector_H.product();
  } else {
    event.getByToken(tp_Token, TPCollection_H);
    size_t nTP = TPCollection_H->size();
    for (size_t i = 0; i < nTP; ++i) {
      TPrefV.push_back(TrackingParticleRef(TPCollection_H, i));
    }
    ptr_TPrefV = &TPrefV;
  }
  TrackingParticleRefVector const& tPC = *ptr_TPrefV;

  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  bool bs_Available = event.getByToken(bsSrc_Token, recoBeamSpotHandle);
  reco::BeamSpot bs;
  if (bs_Available)
    bs = *recoBeamSpotHandle;
  edm::LogVerbatim("MuonTrackValidator") << bs;

  std::vector<const reco::TrackToTrackingParticleAssociator*> associator;
  if (UseAssociators) {
    edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
    for (auto const& associatorName : associators) {
      event.getByLabel(associatorName, theAssociator);
      associator.push_back(theAssociator.product());
    }
  }

  int w = 0;
  unsigned int muonPlotIndex = 0;
  unsigned int trackerPlotIndex = 0;
  const size_t outerLoopSize = UseAssociators ? associators.size() : 1;
  for (unsigned int ww = 0; ww < outerLoopSize; ++ww) {
    for (unsigned int www = 0; www < label.size(); ++www, ++w) {  // increment w here to account for continue
      bool muonIndexIncremented = false;
      bool trackerIndexIncremented = false;
      //
      //get collections from the event
      //
      edm::Handle<edm::View<Track> > trackCollection;
      unsigned int trackCollectionSize = 0;

      reco::RecoToSimCollection recSimColl;
      reco::SimToRecoCollection simRecColl;

      // account for missing track collections (HLT case)
      if (!event.getByToken(track_Collection_Token[www], trackCollection) && ignoremissingtkcollection_) {
        recSimColl.post_insert();
        simRecColl.post_insert();
      }

      //associate tracks to TrackingParticles
      else {
        trackCollectionSize = trackCollection->size();

        if (UseAssociators) {
          edm::LogVerbatim("MuonTrackValidator")
              << "Analyzing " << label[www].process() << ":" << label[www].label() << ":" << label[www].instance()
              << " with " << associators[ww].c_str() << "\n";

          LogTrace("MuonTrackValidator") << "Calling associateRecoToSim method"
                                         << "\n";
          recSimColl = associator[ww]->associateRecoToSim(trackCollection, TPCollection_H);
          LogTrace("MuonTrackValidator") << "Calling associateSimToReco method"
                                         << "\n";
          simRecColl = associator[ww]->associateSimToReco(trackCollection, TPCollection_H);
        } else {
          edm::LogVerbatim("MuonTrackValidator")
              << "Analyzing " << label[www].process() << ":" << label[www].label() << ":" << label[www].instance()
              << " with " << associatormap[www].process() << ":" << associatormap[www].label() << ":"
              << associatormap[www].instance() << "\n";

          Handle<reco::SimToRecoCollection> simtorecoCollectionH;
          event.getByToken(simToRecoCollection_Token[www], simtorecoCollectionH);
          simRecColl = *simtorecoCollectionH.product();

          Handle<reco::RecoToSimCollection> recotosimCollectionH;
          event.getByToken(recoToSimCollection_Token[www], recotosimCollectionH);
          recSimColl = *recotosimCollectionH.product();
        }
      }

      //
      //fill simulation histograms
      //
      edm::LogVerbatim("MuonTrackValidator") << "\n# of TrackingParticles: " << tPC.size() << "\n";
      int ats = 0;
      int st = 0;

      std::optional<double> dR = [&]() -> std::optional<double> {
        if (tpSelector.isSignalOnly() && tPC.size() == 2) {
          const auto& tp1 = *tPC[0];
          const auto& tp2 = *tPC[1];
          return reco::deltaR(tp1.momentum(), tp2.momentum());
        }
        return std::nullopt;
      }();

      for (size_t i = 0; i < tPC.size(); i++) {
        bool TP_is_matched = false;
        bool isChargeOK = true;
        double quality = 0.;

        const TrackingParticleRef& tpr = tPC[i];
        const TrackingParticle& tp = *tpr;

        TrackingParticle::Vector momentumTP;
        TrackingParticle::Point vertexTP;
        double dxySim = 0;
        double dzSim = 0;

        //If the TrackingParticle is collision-like, get the momentum and vertex at production state
        //and the impact parameters w.r.t. PCA
        if (parametersDefiner == "LhcParametersDefinerForTP") {
          LogTrace("MuonTrackValidator") << "TrackingParticle " << i;
          if (!tpSelector(tp))
            continue;
          momentumTP = tp.momentum();
          vertexTP = tp.vertex();
          TrackingParticle::Vector momentum = lhcParametersDefinerTP_->momentum(event, setup, tpr);
          TrackingParticle::Point vertex = lhcParametersDefinerTP_->vertex(event, setup, tpr);
          dxySim = TrackingParticleIP::dxy(vertex, momentum, bs.position());
          dzSim = TrackingParticleIP::dz(vertex, momentum, bs.position());
        }
        //for cosmics get the momentum and vertex at PCA
        else if (parametersDefiner == "CosmicParametersDefinerForTP") {
          edm::LogVerbatim("MuonTrackValidator") << "TrackingParticle " << i;
          if (!cosmictpSelector(tpr, &bs, event, setup))
            continue;
          momentumTP = cosmicParametersDefinerTP_->momentum(event, setup, tpr);
          vertexTP = cosmicParametersDefinerTP_->vertex(event, setup, tpr);
          dxySim = TrackingParticleIP::dxy(vertexTP, momentumTP, bs.position());
          dzSim = TrackingParticleIP::dz(vertexTP, momentumTP, bs.position());
        }

        // Number of counted SimHits depend on the selection of tracker and muon detectors (via cfg parameters)
        int nSimHits = 0;
        if (histoParameters[www].usetracker && histoParameters[www].usemuon) {
          nSimHits = tpr.get()->numberOfHits();
        } else if (!histoParameters[www].usetracker && histoParameters[www].usemuon) {
          nSimHits = tpr.get()->numberOfHits() - tpr.get()->numberOfTrackerHits();
        } else if (histoParameters[www].usetracker && !histoParameters[www].usemuon) {
          nSimHits = tpr.get()->numberOfTrackerHits();
        }

        edm::LogVerbatim("MuonTrackValidator") << "--------------------Selected TrackingParticle #" << tpr.key()
                                               << "  (N counted simhits = " << nSimHits << ")";
        edm::LogVerbatim("MuonTrackValidator")
            << "momentumTP: pt = " << sqrt(momentumTP.perp2()) << ", pz = " << momentumTP.z()
            << ", \t vertexTP: radius = " << sqrt(vertexTP.perp2()) << ",  z = " << vertexTP.z();
        st++;

        double TPeta = momentumTP.eta();
        double xTPeta = getEta(TPeta);
        double absEta = fabs(TPeta);
        double TPpt = sqrt(momentumTP.perp2());
        double xTPpt = getPt(TPpt);
        double TPphi = momentumTP.phi();
        double TPrpos = sqrt(vertexTP.perp2());
        double TPzpos = vertexTP.z();

        int assoc_recoTrack_NValidHits = 0;
        if (simRecColl.find(tpr) != simRecColl.end()) {
          auto const& rt = simRecColl[tpr];
          if (!rt.empty()) {
            RefToBase<Track> assoc_recoTrack = rt.begin()->first;
            TP_is_matched = true;
            ats++;
            if (assoc_recoTrack->charge() != tpr->charge())
              isChargeOK = false;
            quality = rt.begin()->second;
            assoc_recoTrack_NValidHits = assoc_recoTrack->numberOfValidHits();
            edm::LogVerbatim("MuonTrackValidator") << "-----------------------------associated to Track #"
                                                   << assoc_recoTrack.key() << " with quality:" << quality << "\n";
          }
        } else {
          edm::LogVerbatim("MuonTrackValidator")
              << "TrackingParticle #" << tpr.key() << " with pt,eta,phi: " << sqrt(momentumTP.perp2()) << " , "
              << momentumTP.eta() << " , " << momentumTP.phi() << " , "
              << " NOT associated to any reco::Track"
              << "\n";
        }

        // histos for efficiency vs eta
        fillPlotNoFlow(h_simuleta[w], xTPeta);
        if (TP_is_matched) {
          fillPlotNoFlow(h_assoceta[w], xTPeta);
          if (!isChargeOK)
            fillPlotNoFlow(h_misideta[w], xTPeta);
        }

        // histos for efficiency vs phi
        fillPlotNoFlow(h_simulphi[w], TPphi);
        if (absEta < 0.9) {
          fillPlotNoFlow(h_simulphiB[w], TPphi);
        } else if (absEta < 1.2) {
          fillPlotNoFlow(h_simulphiO[w], TPphi);
        } else {
          fillPlotNoFlow(h_simulphiE[w], TPphi);
        }
        if (TP_is_matched) {
          fillPlotNoFlow(h_assocphi[w], TPphi);
          if (absEta < 0.9) {
            fillPlotNoFlow(h_assocphiB[w], TPphi);
          } else if (absEta < 1.2) {
            fillPlotNoFlow(h_assocphiO[w], TPphi);
          } else {
            fillPlotNoFlow(h_assocphiE[w], TPphi);
          }
          if (!isChargeOK)
            fillPlotNoFlow(h_misidphi[w], TPphi);
        }

        // histos for efficiency vs dR
        if (dR) {
          fillPlotNoFlow(h_simuldR[w], *dR);
          if (TP_is_matched) {
            fillPlotNoFlow(h_assocdR[w], *dR);
            if (!isChargeOK)
              fillPlotNoFlow(h_misiddR[w], *dR);
          }
        }

        // histos for efficiency vs pT
        fillPlotNoFlow(h_simulpT[w], xTPpt);
        if (absEta < 0.9) {
          fillPlotNoFlow(h_simulpTB[w], xTPpt);
        } else if (absEta < 1.2) {
          fillPlotNoFlow(h_simulpTO[w], xTPpt);
        } else {
          fillPlotNoFlow(h_simulpTE[w], xTPpt);
        }
        if (TP_is_matched) {
          fillPlotNoFlow(h_assocpT[w], xTPpt);
          if (absEta < 0.9) {
            fillPlotNoFlow(h_assocpTB[w], xTPpt);
          } else if (absEta < 1.2) {
            fillPlotNoFlow(h_assocpTO[w], xTPpt);
          } else {
            fillPlotNoFlow(h_assocpTE[w], xTPpt);
          }
          if (!isChargeOK)
            fillPlotNoFlow(h_misidpT[w], xTPpt);
        }

        // histos for efficiency vs dxy
        fillPlotNoFlow(h_simuldxy[w], dxySim);
        if (TP_is_matched) {
          fillPlotNoFlow(h_assocdxy[w], dxySim);
          if (!isChargeOK)
            fillPlotNoFlow(h_misiddxy[w], dxySim);
        }

        // histos for efficiency vs dz
        fillPlotNoFlow(h_simuldz[w], dzSim);
        if (TP_is_matched) {
          fillPlotNoFlow(h_assocdz[w], dzSim);
          if (!isChargeOK)
            fillPlotNoFlow(h_misiddz[w], dzSim);
        }

        // histos for efficiency vs Radius
        fillPlotNoFlow(h_simulRpos[w], TPrpos);
        if (TP_is_matched)
          fillPlotNoFlow(h_assocRpos[w], TPrpos);

        // histos for efficiency vs z position
        fillPlotNoFlow(h_simulZpos[w], TPzpos);
        if (TP_is_matched)
          fillPlotNoFlow(h_assocZpos[w], TPzpos);

        // histos for efficiency vs Number of Hits
        fillPlotNoFlow(h_simulhit[w], nSimHits);
        if (TP_is_matched) {
          fillPlotNoFlow(h_assochit[w], nSimHits);
          nRecHits_vs_nSimHits[w]->Fill(nSimHits, assoc_recoTrack_NValidHits);

          // charge misid is more useful w.r.t. nRecHits (filled after)
          //if (!isChargeOK) fillPlotNoFlow(h_misidhit[w], nSimHits);
        }

        // histos for efficiency vs PileUp
        fillPlotNoFlow(h_simulpu[w], PU_NumInteractions);
        if (TP_is_matched) {
          fillPlotNoFlow(h_assocpu[w], PU_NumInteractions);
          if (!isChargeOK)
            fillPlotNoFlow(h_misidpu[w], PU_NumInteractions);
        }

        if (doSummaryPlots_) {
          fillPlotNoFlow(h_simul_coll[ww], www);
          if (TP_is_matched) {
            fillPlotNoFlow(h_assoc_coll[ww], www);
          }
        }

      }  // End for (size_t i = 0; i < tPCeff.size(); i++) {

      //
      //fill reconstructed track histograms
      //
      edm::LogVerbatim("MuonTrackValidator")
          << "\n# of reco::Tracks with " << label[www].process() << ":" << label[www].label() << ":"
          << label[www].instance() << ": " << trackCollectionSize << "\n";

      int at = 0;
      int rT = 0;
      for (edm::View<Track>::size_type i = 0; i < trackCollectionSize; ++i) {
        bool Track_is_matched = false;
        bool isChargeOK = true;
        RefToBase<Track> track(trackCollection, i);
        int nRecHits = track->numberOfValidHits();
        rT++;

        std::vector<std::pair<TrackingParticleRef, double> > tp;
        TrackingParticleRef tpr;

        // new logic (bidirectional)
        if (BiDirectional_RecoToSim_association) {
          edm::LogVerbatim("MuonTrackValidator") << "----------------------------------------Track #" << track.key()
                                                 << " (N valid rechits = " << nRecHits << ")";

          if (recSimColl.find(track) != recSimColl.end()) {
            tp = recSimColl[track];
            if (!tp.empty()) {
              tpr = tp.begin()->first;
              // RtS and StR must associate the same pair !
              if (simRecColl.find(tpr) != simRecColl.end()) {
                auto const& assoc_track_checkback = simRecColl[tpr].begin()->first;

                if (assoc_track_checkback.key() == track.key()) {
                  Track_is_matched = true;
                  at++;
                  if (track->charge() != tpr->charge())
                    isChargeOK = false;
                  double Purity = tp.begin()->second;
                  double Quality = simRecColl[tpr].begin()->second;
                  edm::LogVerbatim("MuonTrackValidator")
                      << "with pt=" << track->pt() << " associated with purity:" << Purity << " to TrackingParticle #"
                      << tpr.key() << "\n";
                  if (MABH)
                    h_PurityVsQuality[w]->Fill(Quality, Purity);
                }
              }
            }
          }

          if (!Track_is_matched)
            edm::LogVerbatim("MuonTrackValidator")
                << "with pt=" << track->pt() << " NOT associated to any TrackingParticle"
                << "\n";
        }
        // old logic, valid for cosmics 2-legs reco (bugged for collision scenario)
        else {
          if (recSimColl.find(track) != recSimColl.end()) {
            tp = recSimColl[track];
            if (!tp.empty()) {
              tpr = tp.begin()->first;
              Track_is_matched = true;
              at++;
              if (track->charge() != tpr->charge())
                isChargeOK = false;
              edm::LogVerbatim("MuonTrackValidator") << "reco::Track #" << track.key() << " with pt=" << track->pt()
                                                     << " associated with quality:" << tp.begin()->second << "\n";
            }
          } else {
            edm::LogVerbatim("MuonTrackValidator") << "reco::Track #" << track.key() << " with pt=" << track->pt()
                                                   << " NOT associated to any TrackingParticle"
                                                   << "\n";
          }
        }

        double etaRec = track->eta();
        double xetaRec = getEta(etaRec);

        double ptRec = track->pt();
        double xptRec = getPt(ptRec);

        double qoverpRec = track->qoverp();
        double phiRec = track->phi();
        double thetaRec = track->theta();
        double dxyRec = track->dxy(bs.position());
        double dzRec = track->dz(bs.position());

        double qoverpError = track->qoverpError();
        double ptError = track->ptError();
        double thetaError = track->thetaError();
        double phiError = track->phiError();
        double dxyError = track->dxyError();
        double dzError = track->dzError();

        // histos for coll
        if (doSummaryPlots_) {
          fillPlotNoFlow(h_reco_coll[ww], www);
          if (Track_is_matched) {
            fillPlotNoFlow(h_assoc2_coll[ww], www);
          }
        }

        // histos for fake rate vs eta
        fillPlotNoFlow(h_recoeta[w], xetaRec);
        if (Track_is_matched) {
          fillPlotNoFlow(h_assoc2eta[w], xetaRec);
        }

        // histos for fake rate vs phi
        fillPlotNoFlow(h_recophi[w], phiRec);
        if (Track_is_matched) {
          fillPlotNoFlow(h_assoc2phi[w], phiRec);
        }

        // histos for fake rate vs pT
        fillPlotNoFlow(h_recopT[w], xptRec);
        if (Track_is_matched) {
          fillPlotNoFlow(h_assoc2pT[w], xptRec);
        }

        if (dR) {
          fillPlotNoFlow(h_recodR[w], *dR);
          if (Track_is_matched) {
            fillPlotNoFlow(h_assoc2dR[w], *dR);
          }
        }

        // histos for fake rate vs dxy
        fillPlotNoFlow(h_recodxy[w], dxyRec);
        if (Track_is_matched) {
          fillPlotNoFlow(h_assoc2dxy[w], dxyRec);
        }

        // histos for fake rate vs dz
        fillPlotNoFlow(h_recodz[w], dzRec);
        if (Track_is_matched) {
          fillPlotNoFlow(h_assoc2dz[w], dzRec);
        }

        // histos for fake rate vs Number of RecHits in track
        fillPlotNoFlow(h_recohit[w], nRecHits);
        if (Track_is_matched) {
          fillPlotNoFlow(h_assoc2hit[w], nRecHits);
          // charge misid w.r.t. nRecHits
          if (!isChargeOK)
            fillPlotNoFlow(h_misidhit[w], nRecHits);
        }

        // histos for fake rate vs Number of PU interactions
        fillPlotNoFlow(h_recopu[w], PU_NumInteractions);
        if (Track_is_matched) {
          fillPlotNoFlow(h_assoc2pu[w], PU_NumInteractions);
        }

        // Fill other histos
        TrackingParticle* tpp = const_cast<TrackingParticle*>(tpr.get());
        // TrackingParticle parameters at point of closest approach to the beamline
        TrackingParticle::Vector momentumTP;
        TrackingParticle::Point vertexTP;

        if (parametersDefiner == "LhcParametersDefinerForTP") {
          // following reco plots are made only from tracks associated to selected signal TPs
          if (!(Track_is_matched && tpSelector(*tpp)))
            continue;
          else {
            momentumTP = lhcParametersDefinerTP_->momentum(event, setup, tpr);
            vertexTP = lhcParametersDefinerTP_->vertex(event, setup, tpr);
          }
        } else if (parametersDefiner == "CosmicParametersDefinerForTP") {
          // following reco plots are made only from tracks associated to selected signal TPs
          if (!(Track_is_matched && cosmictpSelector(tpr, &bs, event, setup)))
            continue;
          else {
            momentumTP = cosmicParametersDefinerTP_->momentum(event, setup, tpr);
            vertexTP = cosmicParametersDefinerTP_->vertex(event, setup, tpr);
          }
        }

        if (UseAssociators) {
          if (associators[ww] == "trackAssociatorByChi2") {
            //association chi2
            double assocChi2 = -tp.begin()->second;  //in association map is stored -chi2
            h_assochi2[www]->Fill(assocChi2);
            h_assochi2_prob[www]->Fill(TMath::Prob((assocChi2) * 5, 5));
          } else if (associators[ww] == "trackAssociatorByHits") {
            double fraction = tp.begin()->second;
            h_assocFraction[www]->Fill(fraction);
            h_assocSharedHit[www]->Fill(fraction * nRecHits);
          }
        }

        h_charge[w]->Fill(track->charge());

        // Hits
        h_nhits[w]->Fill(nRecHits);
        nhits_vs_eta[w]->Fill(xetaRec, nRecHits);
        nhits_vs_phi[w]->Fill(phiRec, nRecHits);

        if (histoParameters[www].do_MUOhitsPlots) {
          nDThits_vs_eta[muonPlotIndex]->Fill(xetaRec, track->hitPattern().numberOfValidMuonDTHits());
          nCSChits_vs_eta[muonPlotIndex]->Fill(xetaRec, track->hitPattern().numberOfValidMuonCSCHits());
          nRPChits_vs_eta[muonPlotIndex]->Fill(xetaRec, track->hitPattern().numberOfValidMuonRPCHits());
          if (useGEMs_) {
            nGEMhits_vs_eta[muonPlotIndex]->Fill(xetaRec, track->hitPattern().numberOfValidMuonGEMHits());
          }
          if (useME0_) {
            nME0hits_vs_eta[muonPlotIndex]->Fill(xetaRec, track->hitPattern().numberOfValidMuonME0Hits());
          }
        }

        if (histoParameters[www].do_TRKhitsPlots) {
          nTRK_LayersWithMeas_vs_eta[trackerPlotIndex]->Fill(xetaRec,
                                                             track->hitPattern().trackerLayersWithMeasurement());
          nPixel_LayersWithMeas_vs_eta[trackerPlotIndex]->Fill(xetaRec,
                                                               track->hitPattern().pixelLayersWithMeasurement());
          h_nlosthits[trackerPlotIndex]->Fill(track->numberOfLostHits());
          h_nmisslayers_inner[trackerPlotIndex]->Fill(
              track->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
          h_nmisslayers_outer[trackerPlotIndex]->Fill(
              track->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS));
          nlosthits_vs_eta[trackerPlotIndex]->Fill(xetaRec, track->numberOfLostHits());
        }

        // normalized chi2
        h_nchi2[w]->Fill(track->normalizedChi2());
        h_nchi2_prob[w]->Fill(TMath::Prob(track->chi2(), (int)track->ndof()));
        chi2_vs_nhits[w]->Fill(nRecHits, track->normalizedChi2());
        chi2_vs_eta[w]->Fill(xetaRec, track->normalizedChi2());
        chi2_vs_phi[w]->Fill(phiRec, track->normalizedChi2());

        double ptSim = sqrt(momentumTP.perp2());
        double xptSim = getPt(ptSim);
        double qoverpSim = tpr->charge() / sqrt(momentumTP.x() * momentumTP.x() + momentumTP.y() * momentumTP.y() +
                                                momentumTP.z() * momentumTP.z());
        double etaSim = momentumTP.eta();
        double thetaSim = momentumTP.theta();
        double phiSim = momentumTP.phi();
        double dxySim = TrackingParticleIP::dxy(vertexTP, momentumTP, bs.position());
        double dzSim = TrackingParticleIP::dz(vertexTP, momentumTP, bs.position());

        double etares = etaRec - etaSim;
        double ptRelRes = (ptRec - ptSim) / ptSim;  // relative residual -> resolution
        double ptPull = (ptRec - ptSim) / ptError;
        double qoverpPull = (qoverpRec - qoverpSim) / qoverpError;
        double thetaPull = (thetaRec - thetaSim) / thetaError;
        double phiDiff = reco::deltaPhi(phiRec, phiSim);
        double phiPull = phiDiff / phiError;
        double dxyPull = (dxyRec - dxySim) / dxyError;
        double dzPull = (dzRec - dzSim) / dzError;

        h_etaRes[w]->Fill(etares);
        etares_vs_eta[w]->Fill(xetaRec, etares);

        ptres_vs_eta[w]->Fill(xetaRec, ptRelRes);
        ptres_vs_pt[w]->Fill(xptSim, ptRelRes);
        ptres_vs_phi[w]->Fill(phiRec, ptRelRes);
        h_ptpull[w]->Fill(ptPull);
        ptpull_vs_eta[w]->Fill(xetaRec, ptPull);
        ptpull_vs_phi[w]->Fill(phiRec, ptPull);
        h_qoverppull[w]->Fill(qoverpPull);

        thetaCotres_vs_eta[w]->Fill(xetaRec, cos(thetaRec) / sin(thetaRec) - cos(thetaSim) / sin(thetaSim));
        thetaCotres_vs_pt[w]->Fill(xptSim, cos(thetaRec) / sin(thetaRec) - cos(thetaSim) / sin(thetaSim));
        h_thetapull[w]->Fill(thetaPull);
        thetapull_vs_eta[w]->Fill(xetaRec, thetaPull);
        thetapull_vs_phi[w]->Fill(phiRec, thetaPull);

        phires_vs_eta[w]->Fill(xetaRec, phiDiff);
        phires_vs_pt[w]->Fill(xptSim, phiDiff);
        phires_vs_phi[w]->Fill(phiRec, phiDiff);
        h_phipull[w]->Fill(phiPull);
        phipull_vs_eta[w]->Fill(xetaRec, phiPull);
        phipull_vs_phi[w]->Fill(phiRec, phiPull);

        dxyres_vs_eta[w]->Fill(xetaRec, dxyRec - dxySim);
        dxyres_vs_pt[w]->Fill(xptSim, dxyRec - dxySim);
        h_dxypull[w]->Fill(dxyPull);
        dxypull_vs_eta[w]->Fill(xetaRec, dxyPull);

        dzres_vs_eta[w]->Fill(xetaRec, dzRec - dzSim);
        dzres_vs_pt[w]->Fill(xptSim, dzRec - dzSim);
        h_dzpull[w]->Fill(dzPull);
        dzpull_vs_eta[w]->Fill(xetaRec, dzPull);

        double contrib_Qoverp = qoverpPull * qoverpPull / 5;
        double contrib_dxy = dxyPull * dxyPull / 5;
        double contrib_dz = dzPull * dzPull / 5;
        double contrib_theta = thetaPull * thetaPull / 5;
        double contrib_phi = phiPull * phiPull / 5;
        double assoChi2 = contrib_Qoverp + contrib_dxy + contrib_dz + contrib_theta + contrib_phi;

        LogTrace("MuonTrackValidator") << "normalized Chi2 (track 5-dofs matching) = " << assoChi2 << "\n"
                                       << "\t contrib_Qoverp = " << contrib_Qoverp << "\n"
                                       << "\t contrib_theta = " << contrib_theta << "\n"
                                       << "\t contrib_phi = " << contrib_phi << "\n"
                                       << "\t contrib_dxy = " << contrib_dxy << "\n"
                                       << "\t contrib_dz = " << contrib_dz << "\n";

        LogTrace("MuonTrackValidator") << "ptRec = " << ptRec << "\n"
                                       << "etaRec = " << etaRec << "\n"
                                       << "qoverpRec = " << qoverpRec << "\n"
                                       << "thetaRec = " << thetaRec << "\n"
                                       << "phiRec = " << phiRec << "\n"
                                       << "dxyRec = " << dxyRec << "\n"
                                       << "dzRec = " << dzRec << "\n"
                                       << ""
                                       << "\n"
                                       << "qoverpError = " << qoverpError << "\n"
                                       << "thetaError = " << thetaError << "\n"
                                       << "phiError = " << phiError << "\n"
                                       << "dxyError = " << dxyError << "\n"
                                       << "dzError = " << dzError << "\n"
                                       << ""
                                       << "\n"
                                       << "ptSim = " << ptSim << "\n"
                                       << "etaSim = " << etaSim << "\n"
                                       << "qoverpSim = " << qoverpSim << "\n"
                                       << "thetaSim = " << thetaSim << "\n"
                                       << "phiSim = " << phiSim << "\n"
                                       << "dxySim = " << dxySim << "\n"
                                       << "dzSim = " << dzSim << "\n";
      }  // End of for(edm::View<Track>::size_type i=0; i<trackCollectionSize; ++i) {

      if (histoParameters[www].do_MUOhitsPlots && !muonIndexIncremented) {
        ++muonPlotIndex;
        muonIndexIncremented = true;
      }
      if (histoParameters[www].do_TRKhitsPlots && !trackerIndexIncremented) {
        ++trackerPlotIndex;
        trackerIndexIncremented = true;
      }

      h_tracks[w]->Fill(at);
      h_fakes[w]->Fill(rT - at);
      edm::LogVerbatim("MuonTrackValidator") << "Total Simulated: " << st << "\n"
                                             << "Total Associated (simToReco): " << ats << "\n"
                                             << "Total Reconstructed: " << rT << "\n"
                                             << "Total Associated (recoToSim): " << at << "\n"
                                             << "Total Fakes: " << rT - at << "\n";
    }  // End of for (unsigned int www=0;www<label.size();www++){
  }  //END of for (unsigned int ww=0;ww<associators.size();ww++){
}
