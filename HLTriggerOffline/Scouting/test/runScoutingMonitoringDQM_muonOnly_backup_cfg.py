import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("DQMServices.Components.DQMStoreStats_cfi")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.Services_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_dataRun3_Prompt_v4')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.debugModules = cms.untracked.vstring('DQMGenericClient')
# Enable LogInfo
process.MessageLogger.cerr = cms.untracked.PSet(
    # threshold = cms.untracked.string('ERROR'),
    WARNING = cms.untracked.PSet(
        limit = cms.untracked.int32(10000)
    ),
 )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/442/00000/234b60ac-a73e-470c-8274-cb1ef449c53f.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/420/00000/8de25231-3ffe-49d2-a96a-4c713cff5991.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/425/00000/767f561d-9957-403c-b722-8ccb9031affc.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/416/00000/50c5ef90-3737-4ba2-98ba-aa33ad0caa82.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/420/00000/ca22fd3d-b378-409d-89a9-7866cf6228d8.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/420/00000/780a6ec0-5061-4ffd-b86a-1a73aef0588a.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/420/00000/5aaef7fd-ab63-497c-ac23-ca1f98d99ab8.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/420/00000/86cc8673-48be-4a38-a7aa-311d43dbb3e7.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/433/00000/e6522372-1987-42c7-8495-c77605c8081d.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/416/00000/450d2326-8353-4751-8783-b3bcd5b9721c.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/416/00000/b6979a30-03de-4182-8c40-f3b11cd98c3b.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/442/00000/8c0d3bbc-4ee6-4b88-a918-83406f032528.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/442/00000/9bbb5a81-dae7-42db-beba-f732b60194be.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/425/00000/ecb0947b-a475-40b3-9455-d71c77658b0b.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/441/00000/a5717237-180d-4428-9be5-6942e268fd07.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/456/00000/66b984d5-bcea-4015-b71b-17cab38c4e80.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/454/00000/44fbb3a0-4396-423e-a5b4-a6fee3b8587e.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/456/00000/0b274bfd-f312-4aaf-a9ee-48a9516f6457.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/456/00000/99cc3453-fc79-42ad-bdf3-6fdcdc141ae1.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/456/00000/9c08fb35-2853-4737-9541-3e17433077f6.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/456/00000/e2ac4031-5865-499d-9cb1-fb32c05ab0ec.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/456/00000/696d5ed0-d967-4a58-b890-38f200546a7e.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/456/00000/5c16eb47-17d6-4bd3-8087-df06d4138866.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/470/00000/bedfe11a-8cf6-4b89-800a-3173d00131eb.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/470/00000/0ab5a672-58a0-4c55-b912-02282ac5f011.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/470/00000/366e9dc4-9d6c-475b-91f0-e1e62ce2c3ae.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/454/00000/061075e8-a52a-460b-b3eb-c98eeb0485b0.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/470/00000/abe1acc3-a5e1-41e4-a690-b27837bdd97c.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/470/00000/6bbc191e-6585-43e7-a084-7c5956a50bc6.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/470/00000/11590d22-6c4a-4a05-89d2-0f0dc9859f39.root',
        'root://cms-xrd-global.cern.ch///store/data/Run2024C/ScoutingPFMonitor/MINIAOD/PromptReco-v1/000/379/530/00000/521074ea-46d4-42ab-80db-560b0ae366f9.root'
    )
)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 

#dataFormat = DataFormat.MiniAOD
#switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Winter22_122X_V1_cff']

#add them to the VID producer
#for idmod in my_id_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#process.load("DQMOffline.Scouting.ScoutingMonitoring_cfi")
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )
process.load("HLTriggerOffline.Scouting.ScoutingMuonTagProbeAnalyzer_cfi")
process.load("HLTriggerOffline.Scouting.ScoutingMuonTriggerAnalyzer_cfi")
process.load("HLTriggerOffline.Scouting.ScoutingMuonMonitoring_Client_cff")
#process.load("DQMOffline.Scouting.PatMuonTagProbeAnalyzer_cfi")
process.DQMStore = cms.Service("DQMStore")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load("DQMServices.FileIO.DQMFileSaverOnline_cfi")
process.dqmSaver.tag = 'SCOUTMONIT'

process.options = cms.untracked.PSet(numberOfThreads = cms.untracked.uint32(1))

#process.dqmSaver.convention = 'Offline'
#process.content = cms.EDAnalyzer("EventContentAnalyzer")
#process.allPath = cms.Path(process.scoutingMonitoringTagProbeMuon * process.scoutingMonitoringTriggerMuon * process.muonEfficiency)
process.allPath = cms.Path(process.scoutingMonitoringTriggerMuon*process.muonTriggerEfficiency)
process.p = cms.EndPath(process.dqmSaver)
