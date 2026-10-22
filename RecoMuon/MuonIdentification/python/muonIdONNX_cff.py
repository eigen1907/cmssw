"""Opt in to ONNX scoring before the standard muons1stStep producer."""

import math

import FWCore.ParameterSet.Config as cms


def customizeMuonIdONNX(process, modelPath, backend="cpu", applyFilter=False, threshold=0.38735):
    if backend not in ("cpu", "cuda"):
        raise ValueError("MuonIdONNX backend must be 'cpu' or 'cuda'")
    if not modelPath:
        raise ValueError("MuonIdONNX requires a modelPath")
    if not math.isfinite(threshold) or not 0.0 <= threshold <= 1.0:
        raise ValueError("MuonIdONNX threshold must be a probability between 0 and 1")
    if not hasattr(process, "muons1stStep") or not hasattr(process, "muonIdProducerTask"):
        raise ValueError("Load the standard reconstruction before customizeMuonIdONNX")

    if not hasattr(process, "muonIdONNX"):
        process.load("RecoMuon.MuonIdentification.muonIdONNX_cfi")
    process.muonIdONNX.modelPath = cms.string(str(modelPath))
    process.muonIdONNX.backend = cms.string(backend)
    process.muonIdProducerTask.add(process.muonIdONNX)
    process.muons1stStep.deepMuonScores = cms.InputTag("muonIdONNX")
    process.muons1stStep.applyDeepMuonFilter = cms.bool(applyFilter)
    process.muons1stStep.deepMuonThreshold = cms.double(threshold)
    return process
