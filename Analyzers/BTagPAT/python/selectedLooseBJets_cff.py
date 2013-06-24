import FWCore.ParameterSet.Config as cms

from Analyzers.BTagPAT.Jetpt30Selector_cfi import *
from Analyzers.BTagPAT.selectedLooseBJetsTC2_cfi import *
from Analyzers.BTagPAT.selectedLooseBJetsTC3_cfi import *
from Analyzers.BTagPAT.selectedLooseBJetsTP_cfi import *
from Analyzers.BTagPAT.selectedLooseBJetsSSV_cfi import *
from Analyzers.BTagPAT.selectedLooseBJetsCSV_cfi import *
from Analyzers.BTagPAT.selectedLooseBJetsMSV_cfi import *
from Analyzers.BTagPAT.selectedLooseBJetsSET_cfi import *
from Analyzers.BTagPAT.selectedLooseBJetsSMT_cfi import *

#individual taggers with a loose cut
LooseBJetsTC2 = cms.Sequence(selectedLayer1JetsPt30*selectedLooseBJetsTC2)
LooseBJetsTC3 = cms.Sequence(selectedLayer1JetsPt30*selectedLooseBJetsTC3)
LooseBJetsTP  = cms.Sequence(selectedLayer1JetsPt30*selectedLooseBJetsTP)
LooseBJetsSSV = cms.Sequence(selectedLayer1JetsPt30*selectedLooseBJetsSSV)
LooseBJetsCSV = cms.Sequence(selectedLayer1JetsPt30*selectedLooseBJetsCSV)
LooseBJetsMSV = cms.Sequence(selectedLayer1JetsPt30*selectedLooseBJetsMSV)
LooseBJetsSET = cms.Sequence(selectedLayer1JetsPt30*selectedLooseBJetsSET)
LooseBJetsSMT = cms.Sequence(selectedLayer1JetsPt30*selectedLooseBJetsSMT)

#default tagger for loose 'trackCountingHighEffBJetTags'
selectedLooseBJets = cms.Sequence(selectedLayer1JetsPt30*selectedLooseBJetsTC2)


