import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("TEST")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring([
        #'/store/data/Run2017F/JetHT/MINIAOD/31Mar2018-v1/30000/54B9283C-C637-E811-A002-B496910A0554.root',
        
            'file:/eos/cms/store/data/Run2017G/HighEGJet/AOD/17Nov2017-v2/100000/FADB26FE-4322-E811-8E1A-0CC47A0AD48A.root',


    ])
)

process.prefireVetoFilter = cms.EDFilter("TriggerRulePrefireVetoFilter",
    l1AcceptRecordLabel = cms.InputTag("scalersRawToDigi"),
)

process.ntuple = cms.EDAnalyzer("PrefiringJetAna",
    jetSrc = cms.InputTag("ak4PFJets"),
    trackTag = cms.InputTag("generalTracks"),
    tagJetCut = cms.string("pt > 30 && userInt('looseJetId')"),
    l1egSrc = cms.InputTag("caloStage2Digis:EGamma"),
    l1GtSrc = cms.InputTag("gtStage2Digis"),
)

process.skimPath = cms.Path(process.prefireVetoFilter+process.ntuple)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
    closeFileFast = cms.untracked.bool(True)
)

