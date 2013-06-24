import FWCore.ParameterSet.Config as cms

looseV0Candidates = cms.EDProducer("V0Producer",
    # InputTag that tells which TrackCollection to use for vertexing
    trackRecoAlgorithm = cms.InputTag('generalTracks'),

    # These bools decide whether or not to reconstruct
    #  specific V0 particles
    selectKshorts = cms.bool(True),
    selectLambdas = cms.bool(True),

    useSmoothing = cms.bool(True),
    # Select tracks using TrackBase::TrackQuality.
    # Select ALL tracks by leaving this vstring empty, which
    #   is equivalent to using 'loose'
    #trackQualities = cms.vstring('highPurity', 'goodIterative'),
    trackQualities = cms.vstring('loose'),

#    storeSmoothedTracksInRecoVertex = cms.bool(False),
    # The next parameters are cut values
    # Track quality cuts
    #   Normalized track Chi2:
    tkChi2Cut = cms.double(15.0),
    #   Number of valid hits on track:
    tkNhitsCut = cms.int32(6),

    # Vertex cuts
    vtxChi2Cut = cms.double(7.0),
    collinearityCut = cms.double(0.02),
    #  Setting this one to zero; significance cut is sufficient
    rVtxCut = cms.double(0.0),
#    vtxSignificanceCut = cms.double(22.0),
#    vtxSignificanceCut = cms.double(15.0),
    vtxSignificance2DCut = cms.double(5.0),
#   Not cut (UNUSED)
    vtxSignificance3DCut = cms.double(0.0),
    kShortMassCut = cms.double(0.07),
    lambdaMassCut = cms.double(0.05),
    impactParameterSigCut = cms.double(0.5),
    mPiPiCut = cms.double(1.),
    tkDCACut = cms.double(1.),

    # We check if either track has a hit inside (radially) the vertex position
    #  minus this number times the sigma of the vertex fit
    #  NOTE: Set this to -1 to disable this cut, which MUST be done
    #  if you want to run V0Producer on the AOD track collection!
    innerHitPosCut = cms.double(4.),

    vertexFitter = cms.InputTag('KalmanVertexFitter')

)

