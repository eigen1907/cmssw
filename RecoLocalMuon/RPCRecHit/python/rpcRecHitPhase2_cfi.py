import FWCore.ParameterSet.Config as cms

rpcRecHitPhase2 = cms.EDProducer(
    "RPCRecHitPhase2Producer",
    rpcDigiPhase2Label = cms.InputTag("simMuonRPCDigis"),
    irpcDigiLabel      = cms.InputTag("simMuoniRPCDigis"),
    useIRPC            = cms.bool(False),
)