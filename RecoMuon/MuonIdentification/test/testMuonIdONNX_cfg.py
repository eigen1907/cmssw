"""Compare CPU/CUDA scores on one completed offline RECO event."""

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.setDefault("maxEvents", 1)
options.setDefault("outputFile", "muon_id_onnx_scores.root")
options.register("modelPath", "", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Absolute ONNX model path")
options.parseArguments()
if not options.inputFiles or not options.modelPath:
    raise ValueError("Pass inputFiles=file:/path/to/offline-reco.root and modelPath=/path/to/model.onnx")

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process("ONNXSCORES", Phase2C17I13M9)
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.Geometry.GeometryExtended2026D110Reco_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, "140X_mcRun4_realistic_v4", "")
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(1),
)
from RecoMuon.MuonIdentification.muonIdONNX_cfi import muonIdONNX

process.cpuScores = muonIdONNX.clone(modelPath=options.modelPath, backend="cpu")
process.cudaScores = muonIdONNX.clone(modelPath=options.modelPath, backend="cuda")
process.scores = cms.Path(process.cpuScores + process.cudaScores)
process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName=cms.untracked.string(options.outputFile),
    outputCommands=cms.untracked.vstring("drop *", "keep *_cpuScores_*_ONNXSCORES",
                                        "keep *_cudaScores_*_ONNXSCORES"),
)
process.writeOutput = cms.EndPath(process.output)
