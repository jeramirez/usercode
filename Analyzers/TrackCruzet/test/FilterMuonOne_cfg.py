import FWCore.ParameterSet.Config as cms

#from Analyzers.TrackCruzet.CRAFT_V1P_TrackerPointing_v3 import source

process = cms.Process("FilterTEST")
# geometry
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")

# Conditions (Global Tag is used here):
# Magnetic fiuld: force mag field to be 0.0 tesla
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.GlobalRuns.ForceZeroTeslaField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.ReconstructionCosmics_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagator_cfi")
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

#process.load("Analyzers.TrackCruzet.CRAFT_V1P_TrackerPointing_v3")
#process.source = source
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/data/Commissioning08/Cosmics/RECO/v1/000/069/522/0213BB34-89AB-DD11-A612-000423D99EEE.root',
#       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/085/104DC578-26A0-DD11-8422-001617C3B69C.root',
#       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/114/0E1D6F21-3DA0-DD11-9EC5-001617C3B70E.root',
#
       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/810/00092FF2-A6A4-DD11-BC1F-000423D944F8.root',
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
       '/store/data/Commissioning08/Cosmics/RECO/v1/000/067/647/0032035E-08A3-DD11-8F70-001617C3B6CC.root'

    )
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


#process.prefer("GlobalTag")

process.GoodRunFilter = cms.EDFilter("ReadGoodRunFilter",
    listGoodRun = cms.string('GoodRunList.txt')
)

process.PositionFilter = cms.EDFilter("CruzetMuonFilter",
    MuonONEFilter = cms.string('Muons_More_ONE_Filter.txt'),
    position_x = cms.double(300.0)
    position_y = cms.double(300.0)
    position_z = cms.double(400.0),
    nmuon  = cms.uint32(0),
    MuCollectLabel = cms.string('cosmicMuons'),
)

process.SiTrkFilter = cms.EDFilter("CruzetMuonPointingFilter",
    maxZ = cms.double(230.0),
    radius = cms.double(110.0),
    PropagatorName = cms.string('SteppingHelixPropagatorAny'),
    MuCollectLabel = cms.string('cosmicMuons'),
    PointigFilter = cms.string('SiTrkPointigFilter.txt')
)

process.PixelFilter = cms.EDFilter("CruzetMuonPointingFilter",
    maxZ = cms.double(46.5),
    radius = cms.double(15.0),
    PropagatorName = cms.string('SteppingHelixPropagatorAny'),
    MuCollectLabel = cms.string('cosmicMuons'),
    PointigFilter = cms.string('PixelPointigFilter.txt')
)

process.IPFilter = cms.EDFilter("CruzetMuonPointingFilter",
    maxZ = cms.double(20.0),
    radius = cms.double(4.0),
    PropagatorName = cms.string('SteppingHelixPropagatorAny'),
    MuCollectLabel = cms.string('cosmicMuons'),
    PointigFilter = cms.string('IPPointigFilter.txt')
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.out = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('FiltroMuonsPath')
    ),
    fileName = cms.untracked.string('CRAFT_data_muons_more_One.root')
)

process.GlobalTag.connect = "frontier://PromptProd/CMS_COND_21X_GLOBALTAG"
process.GlobalTag.globaltag = "CRAFT_V2P::All"
process.prefer("GlobalTag")
process.FiltroMuons = cms.Sequence(process.GoodRunFilter*process.PositionFilter*process.SiTrkFilter*process.PixelFilter*process.IPFilter)
process.FiltroMuonsPath = cms.Path(process.FiltroMuons)
process.o = cms.EndPath(process.out)


