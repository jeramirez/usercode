import FWCore.ParameterSet.Config as cms

#from Configuration.StandardSequences.MagneticField_cff import *
from Configuration.StandardSequences.MagneticField_38T_cff import *
#from Configuration.GlobalRuns.ForceZeroTeslaField_cff import *
from Configuration.StandardSequences.Geometry_cff import *
from Configuration.StandardSequences.ReconstructionCosmics_cff import *
from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagator_cfi import *
from TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff import *


trackcruzet = cms.EDFilter("TrackCruzet",
    cosmicTFTrackCollection = cms.InputTag("cosmicTrackFinder"),
    rsTrackCollection = cms.InputTag("rsWithMaterialTracks"),
    ctfTrackCollection = cms.InputTag("ctfWithMaterialTracks"),
    MuCollectLabel = cms.string('cosmicMuons'),
    PropagatorName = cms.string('SteppingHelixPropagatorAny'),
#    PropagatorName = cms.string('SinglePointAlongTrajectory'),
#    PropagatorName = cms.string('RungeKuttaTrackerPropagato'),
#    PropagatorName = cms.string('RKTrackerPropagator'),
#    PropagatorName = cms.string('PropagatorWithMaterial'),
#    PropagatorName = cms.string('SteppingHelixPropagatorAlong'),
#    PropagatorName = cms.string('SteppingHelixPropagatorOpposite'),
#    PropagatorName = cms.string('SteppingHelixPropagator'),
    ascciFileName = cms.string('EnentCross.txt'),
                           
    radius_IP = cms.double(4.0),
    maxZ_IP = cms.double(20.0),
    radius_Pixel = cms.double(15.0),
    maxZ_Pixel = cms.double(46.5),
    radius_SiTrk = cms.double(110.0),
    maxZ_SiTrk = cms.double(230.0),
    nmuon  = cms.uint32(0),

    position_x = cms.double(300.0),
    position_y = cms.double(300.0),
    position_z = cms.double(400.0),

    days  = cms.uint32(21),
    minRun = cms.uint32(67539),
    maxRun = cms.uint32(67810),

    date_year  = cms.uint32(2008),
    date_month = cms.uint32(10),
    date_day   = cms.uint32(23),
    date_hour  = cms.uint32(0),
    date_min   = cms.uint32(0),
    date_sec   = cms.uint32(0),
                           
)



