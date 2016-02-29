# In order to run the code for MC
#[cmsRun IIHE.py DataProcessing='mc' dataset='RunIIFall15MiniAODv2' sample='TT_TuneCUETP8M1_13TeV-powheg-pythia8' address='MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000' file='F8D2CEAA-C5D1-E511-9895-001E675A6C2A.root'  ]
#root://eoscms//cms/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000/F8D2CEAA-C5D1-E511-9895-001E675A6C2A.root

# In order to run the code for DATA
#[cmsRun IIHE.py DataProcessing='mc' dataset='RunIIFall15MiniAODv2' sample='TT_TuneCUETP8M1_13TeV-powheg-pythia8' address='MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000' file='F8D2CEAA-C5D1-E511-9895-001E675A6C2A.root'  ]
#root://eoscms//cms/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F6E918C9-87A6-E511-B3D3-0CC47A4D76B2.root

import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts
import copy
import os

options = opts.VarParsing ('analysis')
options.register('sample',
                 'ZToEE_NNPDF30_13TeV-powheg_M_2300_3500', 
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Sample to analyze')
options.register('address',
                 '',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'address of sample in eos like: MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000')
options.register('file',
                 'abcd.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'file to analyze')
options.register('DataProcessing',
                 "mc",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Data processing types. Options are:mc,data')
options.register('dataset',
                 "SingleElectron",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'datasets to analyze: SingleElectron, DoubleEG')
options.parseArguments()


##########################################################################################
#                                      Global tags                                       #
##########################################################################################
if options.DataProcessing == "mc":
  globalTag = '76X_mcRun2_asymptotic_v12'
if options.DataProcessing == "data":
  globalTag = '76X_dataRun2_v15'

##########################################################################################
#                                  Start the sequences                                   #
##########################################################################################

process = cms.Process("IIHEAnalysis")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.GlobalTag.globaltag = globalTag
print "Global Tag is ", process.GlobalTag.globaltag

process.options = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound') )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

##########################################################################################
#                                         Files                                          #
##########################################################################################
if options.DataProcessing == "mc":
  path = 'root://eoscms//eos/cms/store/'+ options.DataProcessing + '/' + options.dataset  + '/' + options.sample + '/' + options.address + '/' + options.file

if options.DataProcessing == "data":
  path = 'root://eoscms//cms/eos/store/'+ options.DataProcessing + '/' + options.sample + '/' + options.address + '/' + options.file

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring())

process.source.fileNames.append( path )
#readFiles.extend([
#'file:/pnfs/iihe/cms/store/user/rgoldouz/ZToEE_NNPDF30_13TeV-powheg_M_2300_3500/C84500E1-6BB8-E511-A9D6-002590E505FE.root'
#])

filename_out = "output_%s" % (options.sample + '_' + options.file)
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('file:/tmp/' + filename_out) )
process.TFileService = cms.Service("TFileService", fileName = cms.string(        'file:/tmp/' +  filename_out) )

##########################################################################################
#                                   IIHETree options                                     #
##########################################################################################
process.load("UserCode.IIHETree.IIHETree_cfi")

pt_threshold = 15

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
process.IIHEAnalysis.eventRho = cms.InputTag('fixedGridRhoFastjetAll')
process.IIHEAnalysis.ebReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEBRecHits")
process.IIHEAnalysis.eeReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEERecHits")
process.IIHEAnalysis.esReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedESRecHits")
process.IIHEAnalysis.generatorLabel = cms.InputTag("generator")
process.IIHEAnalysis.PileUpSummaryInfo = cms.untracked.InputTag('slimmedAddPileupInfo')
process.IIHEAnalysis.genParticleSrc = cms.InputTag("prunedGenParticles")

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
process.IIHEAnalysis.includeTriggerModule        = cms.untracked.bool(True)
process.IIHEAnalysis.includeEventModule          = cms.untracked.bool(True)
process.IIHEAnalysis.includeVertexModule         = cms.untracked.bool(True) 
process.IIHEAnalysis.includeSuperClusterModule   = cms.untracked.bool(True)
process.IIHEAnalysis.includePhotonModule         = cms.untracked.bool(True)
process.IIHEAnalysis.includeElectronModule       = cms.untracked.bool(True)
process.IIHEAnalysis.includeMuonModule           = cms.untracked.bool(True)
process.IIHEAnalysis.includeMETModule            = cms.untracked.bool(True)
process.IIHEAnalysis.includeHEEPModule           = cms.untracked.bool(True)
process.IIHEAnalysis.includeZBosonModule         = cms.untracked.bool(True)
process.IIHEAnalysis.includeTracksModule         = cms.untracked.bool(True)


process.IIHEAnalysis.includeMCTruthModule         = cms.untracked.bool(('mc' in options.DataProcessing))
#change it to true if you want to save all events
process.IIHEAnalysis.includeAutoAcceptEventModule= cms.untracked.bool(False)

process.IIHEAnalysis.debug = cms.bool(False)

##########################################################################################
#                            Woohoo!  We're ready to start!                              #
##########################################################################################
#process.p1 = cms.Path(process.kt6PFJetsForIsolation+process.IIHEAnalysis)
process.p1 = cms.Path(process.IIHEAnalysis)
