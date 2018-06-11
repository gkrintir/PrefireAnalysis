import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing()
options.register(
    "is2016",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "2016 vs 2017+ label name for trigger objects (selectedPatTrigger vs. slimmedPatTrigger)"
)
options.parseArguments()

process = cms.Process("TEST")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring([
        #'/store/data/Run2017F/SingleElectron/MINIAOD/31Mar2018-v1/100000/64244276-8937-E811-AD7C-0CC47A4C8EBA.root',
            #'file:/eos/cms/store/data/Run2017G/HighEGJet/MINIAOD/17Nov2017-v2/100000/48BCD29C-C022-E811-A80C-0242AC130002.root'
          'file:/eos/cms/store/data/Run2017G/HighEGJet/AOD/17Nov2017-v2/100000/FADB26FE-4322-E811-8E1A-0CC47A0AD48A.root',

    ])
)

process.prefireVetoFilter = cms.EDFilter("TriggerRulePrefireVetoFilter",
    l1AcceptRecordLabel = cms.InputTag("scalersRawToDigi"),
)

# Latest 2016 is cutBasedElectronID-Summer16-80X-V1-medium but close enough
eleID = 'cutBasedElectronID-Spring15-25ns-V1-standalone-medium' if options.is2016 else 'cutBasedElectronID-Fall17-94X-V1-medium'

process.ntuple = cms.EDAnalyzer("PrefiringZAna",
    #electronSrc = cms.InputTag("slimmedElectrons"),
    #vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    #tagElectronCut = cms.string("pt > 30 && abs(eta) < 2.5 && electronID('%s') && userInt('HLT_Ele32_WPTight_Gsf')" % eleID),
    l1egSrc = cms.InputTag("caloStage2Digis:EGamma"),
    #triggerObjects = cms.InputTag("selectedPatTrigger" if options.is2016 else "slimmedPatTrigger"),
    #triggerPrescales = cms.InputTag("patTrigger"),
    processName = cms.string("HLT"),
    loadTriggersFromHLT = cms.untracked.bool(False),
    triggerNames = cms.vstring("HLT_HIEle20_WPLoose_Gsf_v*"),
    #triggerNames = cms.vstring("HLT_HIDoublePhoton15_Eta3p1ForPPRef_Mass50to1000_v*"),
    triggerResults = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    gsfElectronLabel = cms.InputTag("gedGsfElectrons")
)

process.skimPath = cms.Path(process.prefireVetoFilter+process.ntuple)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
    closeFileFast = cms.untracked.bool(True)
)

