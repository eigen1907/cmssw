import FWCore.ParameterSet.Config as cms

process = cms.Process("DeepMuonReco")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.Geometry.GeometryExtendedRun4D121Reco_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:recosim.root')
    fileNames = cms.untracked.vstring('file:/home/joshin/workspace-gate/Store/store-hdfs/DeepMuonReco/CMSSW_16_0_0_pre2/Run4D121/tenmu-pileup/step4_recosim/output_0.root')
)

process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

process.tpClusterProducer.pixelSimLinkSrc = cms.InputTag("mixData", "Pixel")
process.tpClusterProducer.phase2OTSimLinkSrc = cms.InputTag("mixData", "Phase2OTDigiSimLink")
process.tpClusterProducer.stripSimLinkSrc = cms.InputTag("")
process.tpClusterProducer.pixelClusterSrc = cms.InputTag("siPixelClusters")
process.tpClusterProducer.phase2OTClusterSrc = cms.InputTag("siPhase2Clusters")

process.load("DeepMuonReco.Ntuplizer.deepMuonRecoNtuplizer_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output.root")
)

process.p = cms.Path(
    process.tpClusterProducer *
    process.quickTrackAssociatorByHits *
    process.deepMuonRecoNtuplizer
)
