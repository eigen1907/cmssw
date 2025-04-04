import FWCore.ParameterSet.Config as cms

# Reconstruction geometry services
#  Tracking Geometry
from Geometry.CommonTopologies.globalTrackingGeometry_cfi import *

#Tracker
from RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi import *
# TrackerAdditionalParametersPerDet contains only default values, needed for consistency with Phase 2
from Geometry.TrackerGeometryBuilder.TrackerAdditionalParametersPerDet_cfi import *
from Geometry.TrackerGeometryBuilder.trackerParameters_cfi import *
from Geometry.TrackerNumberingBuilder.trackerTopology_cfi import *

#Muon
from Geometry.MuonNumbering.muonNumberingInitialization_cfi import *
from RecoMuon.DetLayers.muonDetLayerGeometry_cfi import *

#  Alignment
from Geometry.TrackerGeometryBuilder.idealForDigiTrackerGeometry_cff import *
from Geometry.CSCGeometryBuilder.idealForDigiCscGeometry_cff import *
from Geometry.DTGeometryBuilder.idealForDigiDtGeometry_cff import *

#  Calorimeters
from Geometry.CaloEventSetup.CaloTopology_cfi import *
from Geometry.CaloEventSetup.CaloGeometry_cff import *
from Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi import *
from Geometry.EcalMapping.EcalMapping_cfi import *
from Geometry.EcalMapping.EcalMappingRecord_cfi import *
from Geometry.EcalCommonData.ecalSimulationParameters_cff import *
from Geometry.HcalCommonData.hcalDDConstants_cff import *
from Geometry.HcalEventSetup.hcalTopologyIdeal_cfi import *
from Geometry.ForwardGeometry.zdcTopologyEP_cfi import *
from Geometry.ForwardGeometry.ForwardGeometry_cfi import *
