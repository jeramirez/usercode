import FWCore.ParameterSet.Config as cms

looseCascade = cms.EDProducer('CascadeProducer',
    v0Algo           = cms.string('looseV0Candidates'),
    v0DecayName      = cms.untracked.string('Lambda'),
    trackingAlgo     = cms.InputTag('generalTracks'),
    VtxFromFit         = cms.bool(False),
    MassFromFit        = cms.bool(False),
    XiDCACut           = cms.double(2.0),
    mLambdaPiCut       = cms.double(1.4),
    mLambdaKCut        = cms.double(1.8),
    LambdaMassWidthCut = cms.double(0.008),
    refitPrimary        = cms.bool(True),
    vtxSignificance2DCut= cms.double(0.0),
    vtxSignificance3DCut= cms.double(3.0),
    IP3DCut             = cms.double(3.0)
)

looseCasFilter = cms.EDFilter("VeeCountFilter",
    algo = cms.string('looseCascade'),
    name = cms.string('Cascade'),
    minNumber = cms.uint32(1)
)

looseOmeFilter = cms.EDFilter("VeeCountFilter",
    algo = cms.string('looseCascade'),
    name = cms.string('Omega'),
    minNumber = cms.uint32(1)
)
