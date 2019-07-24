import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource", fileNames = 
	cms.untracked.vstring('file:/afs/cern.ch/work/r/rymuelle/public/EDAnalyzerMuons/test_miniAOD.root'))

process.demo = cms.EDAnalyzer('Low_Occupancy_Analyzer', muons = cms.InputTag('ctfWithMaterialTracks'))

process.TFileService = cms.Service("TFileService", fileName = cms.string('TightMuonProfile.root'))

process.p = cms.Path(process.demo)
