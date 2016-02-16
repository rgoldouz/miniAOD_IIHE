# This is the standard pset for the IIHETree ntuple maker.  You need to change this pset
# to tell CMSSW which samples you are using.  You can do this in a semi-automated way
# by using sed to replace the field "###JOBTYPE###".  You can also change the
# pt_threshold of leptons, but be careful not to go above 20 if you want to perform the
# scale factor study.

#job_type = 'MC_50ns'
#job_type = 'MC_25ns'
#job_type = 'data_Run2015BC'
#job_type = 'data_Run2015D'
#job_type = 'data_Run2015_Rereco'
# Run2 Fall15 MiniAOD v2 campaign (76X version 2) 76X releases 7_6_3= or later; Global tags: 76X_mcRun2_asymptotic_v12 (25ns) 
job_type = '76X_mc_Run2'
###JOBTYPE###

# Only keep leptons with at least 25 GeV of pT or ET
pt_threshold = 15

# Not actually used yet, but there are plans for this in the future.
miniAOD = False

##########################################################################################
#                                      Global tags                                       #
##########################################################################################
if job_type == 'MC_50ns':
    globalTag = 'MCRUN2_74_V9A::All' # 50ns asymptotic
elif job_type == 'MC_25ns':
    globalTag = 'MCRUN2_74_V9::All'  # 25ns asymptotic
elif job_type == 'data_Run2015BC':
    globalTag = '74X_dataRun2_Prompt_v1'
elif job_type == 'data_Run2015D':
    globalTag = '74X_dataRun2_Prompt_v4'
elif job_type == '76X_mc_Run2':
    globalTag = '76X_mcRun2_asymptotic_v12'
elif job_type == 'data_Run2015_Rereco':
    globalTag = '76X_dataRun2_v15'
##########################################################################################
#                                  Start the sequences                                   #
##########################################################################################
import FWCore.ParameterSet.Config as cms

process = cms.Process("IIHEAnalysis")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.EventContent.EventContent_cff")

if job_type=='data_Run2015BC' or job_type=='data_Run2015D':
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
else:
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

##########################################################################################
#                                         Files                                          #
##########################################################################################
readFiles = cms.untracked.vstring()
secFiles  = cms.untracked.vstring()
process.source = cms.Source("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend([
'file:/afs/cern.ch/work/r/rgoldouz/IIHETree/miniAOD/CMSSW_7_6_3/src/UserCode/IIHETree/test/ZToEE_NNPDF30_13TeV-powheg_M_2300_3500_76X_mcRun2_asymptotic_v12-v1.root'
])

filename_out = 'outfile.root'
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string(filename_out) )
process.TFileService = cms.Service("TFileService", fileName = cms.string(          filename_out) )

##########################################################################################
#                                     Main options                                       #
##########################################################################################
process.GlobalTag.globaltag = globalTag
print "Global Tag is ", process.GlobalTag.globaltag

process.options = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

# Taken from Sherif's code
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

##########################################################################################
#                                   IIHETree options                                     #
##########################################################################################
process.load("UserCode.IIHETree.IIHETree_cfi")

# Only save some triggers.
process.IIHEAnalysis.TriggerResults = cms.InputTag('TriggerResults', '', 'HLT')
process.IIHEAnalysis.triggerEvent = cms.InputTag('selectedPatTrigger')
triggers = 'singleElectron;doubleElectron'
process.IIHEAnalysis.triggers = cms.untracked.string(triggers)

#process.IIHEAnalysis.triggers = cms.untracked.string('doubleElectron')

process.IIHEAnalysis.globalTag = cms.string(globalTag)

# Collections.
process.IIHEAnalysis.photonCollection    = cms.InputTag('slimmedPhotons'        )
process.IIHEAnalysis.electronCollection  = cms.InputTag('slimmedElectrons')
process.IIHEAnalysis.muonCollection      = cms.InputTag('slimmedMuons'          )
process.IIHEAnalysis.METCollection      = cms.InputTag('slimmedMETs'          )
process.IIHEAnalysis.superClusterCollection = cms.InputTag('reducedEgamma', 'reducedSuperClusters','PAT')
process.IIHEAnalysis.reducedBarrelRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
process.IIHEAnalysis.reducedEndcapRecHitCollection = cms.InputTag('reducedEcalRecHitsEE')
process.IIHEAnalysis.eventRho = cms.InputTag('fixedGridRhoAll')
process.IIHEAnalysis.ebReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEBRecHits")
process.IIHEAnalysis.eeReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEERecHits")
process.IIHEAnalysis.generatorLabel = cms.InputTag("generator")
process.IIHEAnalysis.PileUpSummaryInfo = cms.untracked.InputTag('addPileupInfo')
process.IIHEAnalysis.genParticleSrc = cms.InputTag("prunedGenParticles")

if miniAOD:
    process.IIHEAnalysis.primaryVertex = cms.InputTag('offlineSlimmedPrimaryVertices')
    process.IIHEAnalysis.photonCollection    = cms.InputTag('slimmedPhotons'  )
    process.IIHEAnalysis.electronCollection  = cms.InputTag('slimmedElectrons')
    process.IIHEAnalysis.muonCollection      = cms.InputTag('slimmedMuons'    )

# Trigger matching stuff.  0.5 should be sufficient.
process.IIHEAnalysis.muon_triggerDeltaRThreshold = cms.untracked.double(0.5)
process.IIHEAnalysis.HEEP_triggerDeltaRThreshold = cms.untracked.double(0.5)

# In the absence of high ET electrons, only save events with really high Z candidates.
process.IIHEAnalysis.ZBosonZMassAcceptLower    = cms.untracked.double(850)
# Don't bother with J/psi or Upsilon, they will only weigh us down!
process.IIHEAnalysis.ZBosonJPsiAcceptMassLower = cms.untracked.double(1e6)
process.IIHEAnalysis.ZBosonJPsiAcceptMassUpper = cms.untracked.double(1e6)
process.IIHEAnalysis.ZBosonUpsAcceptMassLower  = cms.untracked.double(1e6)
process.IIHEAnalysis.ZBosonUpsAcceptMassUpper  = cms.untracked.double(1e6)

# But make sure we save Z bosons from 50 GeV and up.
process.IIHEAnalysis.ZBosonZMassLowerCuttoff   = cms.untracked.double( 50)
process.IIHEAnalysis.ZBosonDeltaRCut           = cms.untracked.double(1e-3)

# Only save Z->ee, Z->em.
process.IIHEAnalysis.ZBosonEtThreshold = cms.untracked.double(pt_threshold)
process.IIHEAnalysis.ZBosonSaveZee  = cms.untracked.bool(True )
process.IIHEAnalysis.ZBosonSaveZmm  = cms.untracked.bool(True )
process.IIHEAnalysis.ZBosonSaveZem  = cms.untracked.bool(True )
process.IIHEAnalysis.ZBosonSaveZeeg = cms.untracked.bool(False)
process.IIHEAnalysis.ZBosonSaveZmmg = cms.untracked.bool(False)

process.IIHEAnalysis.electrons_ETThreshold = cms.untracked.double(pt_threshold)
process.IIHEAnalysis.muon_pTThreshold      = cms.untracked.double(pt_threshold)

process.IIHEAnalysis.LeptonsAccept_pTThreshold = cms.untracked.double(pt_threshold)
# Require at least two leptons...
process.IIHEAnalysis.LeptonsAccept_nLeptons    = cms.untracked.double(2)
# ...at least one of which is an electron.
process.IIHEAnalysis.LeptonsAccept_nElectrons  = cms.untracked.double(1)


process.IIHEAnalysis.includeLeptonsAcceptModule  = cms.untracked.bool(True)
process.IIHEAnalysis.includeTriggerModule        = cms.untracked.bool(False)
process.IIHEAnalysis.includeEventModule          = cms.untracked.bool(True)
process.IIHEAnalysis.includeMCTruthModule        = cms.untracked.bool(True)
process.IIHEAnalysis.includeVertexModule         = cms.untracked.bool(True) 
process.IIHEAnalysis.includeSuperClusterModule   = cms.untracked.bool(True)
process.IIHEAnalysis.includePhotonModule         = cms.untracked.bool(True)
process.IIHEAnalysis.includeElectronModule       = cms.untracked.bool(True)
process.IIHEAnalysis.includeMuonModule           = cms.untracked.bool(True)
process.IIHEAnalysis.includeMETModule            = cms.untracked.bool(True)
process.IIHEAnalysis.includeHEEPModule           = cms.untracked.bool(True)
process.IIHEAnalysis.includeZBosonModule         = cms.untracked.bool(True)
process.IIHEAnalysis.includeAutoAcceptEventModule= cms.untracked.bool(True)
process.IIHEAnalysis.includeTracksModule         = cms.untracked.bool(True)


#process.IIHEAnalysis.includeTriggerModule         = cms.untracked.bool(True )
#process.IIHEAnalysis.includeMCTruthModule         = cms.untracked.bool(('MC' in job_type))
#process.IIHEAnalysis.includeAutoAcceptEventModule = cms.untracked.bool(False)
#process.IIHEAnalysis.PileUpSummaryInfo = cms.untracked.InputTag('addPileupInfo')
#process.IIHEAnalysis.genParticleSrc = cms.InputTag("prunedGenParticles")


process.IIHEAnalysis.debug = cms.bool(False)

##########################################################################################
#                            Woohoo!  We're ready to start!                              #
##########################################################################################
#process.p1 = cms.Path(process.kt6PFJetsForIsolation+process.IIHEAnalysis)
process.p1 = cms.Path(process.IIHEAnalysis)
