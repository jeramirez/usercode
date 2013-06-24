import FWCore.ParameterSet.Config as cms

hlt_d0_phi_analyzer = cms.EDAnalyzer("HLTBeamSpotAnalyzer",
    BSAnalyzerParameters = cms.PSet(
        RunAllFitters = cms.bool(False), ## run only default fitter

        WriteToDB = cms.bool(False), ## do not write results to DB

        MaximumNtracks = cms.int32(500), ## disable for the moment 
       
        minSHit = cms.int32(0), ## Number of Hit in Strip (8)
	
        TrackCollection = cms.untracked.string('hltPixelTracks'),
#        TrackCollection = cms.untracked.string('hltBLifetimeRegionalCtfWithMaterialTracks'),
 #       TrackCollection = cms.untracked.string('generalTracks'),
        InputBeamWidth = cms.untracked.double(-1.0), ## if -1 use the value calculated by the analyzer

        MinimumPt = cms.double(2.0) ## Gev/c

    ),
    OutputFileName = cms.untracked.string('hlt_analyze_d0_phi.root')
)


