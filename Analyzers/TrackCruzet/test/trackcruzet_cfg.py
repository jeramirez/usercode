import FWCore.ParameterSet.Config as cms

#from Analyzers.TrackCruzet.Commissioning08-PromptReco-v2 import source

process = cms.Process("TrkC")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# geometry
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")

process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

process.load("Analyzers.TrackCruzet.trackcruzet_cfi")
#process.load("Analyzers.TrackCruzet.Commissioning08-PromptReco-v2")
#process.source = source
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#       '/store/data/Commissioning08/Cosmics/RECO/v1/000/069/522/0213BB34-89AB-DD11-A612-000423D99EEE.root',
#       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/085/104DC578-26A0-DD11-8422-001617C3B69C.root',
#       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/114/0E1D6F21-3DA0-DD11-9EC5-001617C3B70E.root',
#
       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/539/024F12F4-58A2-DD11-9C2E-000423D952C0.root',
       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/541/060E174C-89A2-DD11-96AA-001D09F23F2A.root',
       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/544/00061DC9-9BA2-DD11-BD60-001D09F24EC0.root',
       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/548/0A70494C-AAA2-DD11-BE41-001D09F241B9.root',
       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/557/0C1A23A3-B9A2-DD11-9115-0016177CA778.root',
       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/573/027EC372-B8A2-DD11-A2A2-000423D98EC4.root',
       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/645/084B809F-E7A2-DD11-B560-001D09F2910A.root',
#       
       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/838/001E79A3-6DA5-DD11-977D-001617E30CC8.root',
#
       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/647/0032035E-08A3-DD11-8F70-001617C3B6CC.root',
       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/810/00092FF2-A6A4-DD11-BC1F-000423D944F8.root',

#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/818/0037F1C0-DDA4-DD11-B560-000423D99B3E.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/838/001E79A3-6DA5-DD11-977D-001617E30CC8.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/068/124/0A8114BB-00A7-DD11-B974-001617C3B6DC.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/068/141/005C47C4-61A7-DD11-828D-000423D174FE.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/068/264/042413BC-9AA7-DD11-9BDD-000423D986A8.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/068/273/06FED6C9-C4A7-DD11-BB66-000423D98E6C.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/068/276/002A93BA-CEA7-DD11-8E98-000423D6C8E6.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/068/279/002D32CB-F8A7-DD11-AD55-001617C3B76A.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/068/288/00ACA63C-48A8-DD11-B89C-000423D6B2D8.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/068/665/00C70214-03A9-DD11-9778-0019B9F709A4.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/068/926/022ED031-E6A8-DD11-ACAF-000423D98E30.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/068/949/061BF1B8-01A9-DD11-9F8B-000423D951D4.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/068/958/001366F1-16A9-DD11-A3DC-001617E30F58.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/069/046/00335B56-B0A9-DD11-88E2-000423D98930.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/069/269/0229CBE3-B4AA-DD11-8845-000423D991F0.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/069/276/02AFB52E-6DAA-DD11-8047-001617C3B69C.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/069/382/00ACFB84-48AB-DD11-A5E8-001617E30CA4.root',
#        '/store/data/Commissioning08/Cosmics/RECO/v1/000/069/396/00D0E3CF-5AAB-DD11-AF40-000423D98DB4.root',
        '/store/data/Commissioning08/Cosmics/RECO/v1/000/069/522/0213BB34-89AB-DD11-A612-000423D99EEE.root'

))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(757992))

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
#
process.GoodRunFilter = cms.EDFilter("ReadGoodRunFilter",
    listGoodRun = cms.string('GoodRunList.txt')
)
process.FiltroMuons = cms.Sequence(process.GoodRunFilter*process.trackcruzet)

#
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('CRAFT_Commissioning08-PromptReco-v2_m0_all_theta.root')
)

process.p = cms.Path(process.FiltroMuons)
#process.p = cms.Path(process.trackcruzet)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000000000
process.GlobalTag.connect = "frontier://PromptProd/CMS_COND_21X_GLOBALTAG"
process.GlobalTag.globaltag = "CRAFT_V2P::All"
process.prefer("GlobalTag")
process.trackcruzet.ctfTrackCollection = 'ctfWithMaterialTracksP5'
process.trackcruzet.rsTrackCollection = 'rsWithMaterialTracksP5'
process.trackcruzet.cosmicTFTrackCollection = 'cosmictrackfinderP5'
process.trackcruzet.ascciFileName = 'CRAFT_Commissioning08-PromptReco-v2.txt'


