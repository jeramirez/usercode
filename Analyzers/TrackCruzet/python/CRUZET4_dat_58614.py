import FWCore.ParameterSet.Config as cms
#dataset /Cosmics/CRUZET4_v1_CRZT210_V1_TrackerPointing_v1/RECO :: run 58614

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( ( 
       '/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/728B06EF-BF72-DD11-8472-001731AF67BD.root',
       '/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/8421AD2E-4073-DD11-96F2-001731AF66A5.root',
       '/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0000/E4B83F57-3873-DD11-87D7-00304876A147.root',
       '/store/data/CRUZET4_v1/Cosmics/RECO/CRZT210_V1_TrackerPointing_v1/0001/50748CD9-1E73-DD11-8C95-001A92971B04.root') );


source.secFiles.extend( (
       '/store/data/Commissioning08/Cosmics/RECO/CRUZET4_v1/000/058/614/DC60A34D-6371-DD11-916A-001617E30F56.root',
       '/store/data/Commissioning08/Cosmics/RECO/CRUZET4_v1/000/058/614/8055E0CF-6071-DD11-BEF9-00161757BF42.root',
       '/store/data/Commissioning08/Cosmics/RECO/CRUZET4_v1/000/058/614/6E232C18-6671-DD11-9E5A-000423D94534.root',
       '/store/data/Commissioning08/Cosmics/RECO/CRUZET4_v1/000/058/614/4650EF3E-6871-DD11-9122-001D09F290BF.root',
       '/store/data/Commissioning08/Cosmics/RECO/CRUZET4_v1/000/058/614/A4DEE891-6171-DD11-992B-001617C3B76E.root',
       '/store/data/Commissioning08/Cosmics/RECO/CRUZET4_v1/000/058/614/2E2D355A-6571-DD11-AD58-000423D98A44.root',
       '/store/data/Commissioning08/Cosmics/RECO/CRUZET4_v1/000/058/614/A8FC8853-6F71-DD11-818A-000423D992A4.root',
       '/store/data/Commissioning08/Cosmics/RECO/CRUZET4_v1/000/058/614/D6BE1ECA-6D71-DD11-A2D9-000423D98950.root',
       '/store/data/Commissioning08/Cosmics/RECO/CRUZET4_v1/000/058/614/B2D80290-6971-DD11-9125-0019B9F730D2.root',
       '/store/data/Commissioning08/Cosmics/RECO/CRUZET4_v1/000/058/614/92F43381-7371-DD11-84C3-000423D98E54.root',
       '/store/data/Commissioning08/Cosmics/RECO/CRUZET4_v1/000/058/614/640744EC-6871-DD11-89BB-001D09F2525D.root') );
