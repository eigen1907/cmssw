import FWCore.ParameterSet.Config as cms

deepMuonRecoNtuplizer = cms.EDAnalyzer('DeepMuonRecoNtuplizer',
    muons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"), 
    trackingParticles = cms.InputTag("mix", "MergedTrackTruth"),
    associator = cms.InputTag("quickTrackAssociatorByHits"),

    rpcRecHits = cms.InputTag("rpcRecHits"),
    gemRecHits = cms.InputTag("gemRecHits"),
    dtSegments = cms.InputTag("dt4DSegments"), 
    cscSegments = cms.InputTag("cscSegments")
)