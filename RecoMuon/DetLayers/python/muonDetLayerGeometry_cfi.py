import FWCore.ParameterSet.Config as cms

#
# This cfi should be included to build the muon DetLayers.
#
MuonDetLayerGeometryESProducer = cms.ESProducer(
    "MuonDetLayerGeometryESProducer",
    rpcUseLegacyIsFront = cms.bool(False)
)

from Configuration.Eras.Modifier_run1_cff import run1
from Configuration.Eras.Modifier_run2_common_cff import run2_common
from Configuration.Eras.Modifier_run3_common_cff import run3_common

run1.toModify(MuonDetLayerGeometryESProducer, rpcUseLegacyIsFront = True)
run2_common.toModify(MuonDetLayerGeometryESProducer, rpcUseLegacyIsFront = True)
run3_common.toModify(MuonDetLayerGeometryESProducer, rpcUseLegacyIsFront = True)

from Configuration.Eras.python.Era_Run3_pp_on_PbPb_2025_cff import Run3_pp_on_PbPb_2025
Run3_2025_pp_on_PbPb.toModify(MuonDetLayerGeometryESProducer, rpcUseLegacyIsFront = False)
