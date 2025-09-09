import FWCore.ParameterSet.Config as cms

#
# This cfi should be included to build the muon DetLayers.
#
MuonDetLayerGeometryESProducer = cms.ESProducer(
    "MuonDetLayerGeometryESProducer",
    rpcUseLegacyIsFront = cms.bool(False)
)

from Configuration.Eras.Modifier_run2_common_cff import run2_common
from Configuration.Eras.Modifier_run3_common_cff import run3_common

run2_common.toModify(MuonDetLayerGeometryESProducer, rpcUseLegacyIsFront = True)
run3_common.toModify(MuonDetLayerGeometryESProducer, rpcUseLegacyIsFront = True)

from Configuration.Eras.Modifier_pp_on_PbPb_run3_2025_cff import pp_on_PbPb_run3_2025
pp_on_PbPb_run3_2025.toModify(MuonDetLayerGeometryESProducer, rpcUseLegacyIsFront = False)
