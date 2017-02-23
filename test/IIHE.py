# In order to run the code for MC on lxplus
#[cmsRun IIHE.py DataProcessing="mc" dataset="RunIIFall15MiniAODv2" sample="TT_TuneCUETP8M1_13TeV-powheg-pythia8" address="MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000" file="F8D2CEAA-C5D1-E511-9895-001E675A6C2A.root"  ]
#root://eoscms//cms/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000/F8D2CEAA-C5D1-E511-9895-001E675A6C2A.root

# In order to run the code for DATA on lxplus
#[cmsRun IIHE.py DataProcessing="data" dataset="Run2015D" sample="SingleElectron" address="MINIAOD/16Dec2015-v1/20000" file="001E76A5-D3A6-E511-BC32-008CFA05E874.root"  ]
#root://eoscms//cms/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F6E918C9-87A6-E511-B3D3-0CC47A4D76B2.root

import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts
import copy
import os

options = opts.VarParsing ("analysis")
options.register("sample",
                 "", 
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "Sample to analyze")
options.register("address",
                 "",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "address of sample in eos like: MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000")
options.register("file",
                 "",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "file to analyze")
options.register("DataProcessing",
                 "data",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "Data processing types. Options are:mc,rerecodata,promptdata")
options.register("dataset",
                 "",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "datasets to analyze: SingleElectron, DoubleEG")
options.parseArguments()


##########################################################################################
#                                      Global tags                                       #
##########################################################################################
if options.DataProcessing == "mc2016":
  globalTag = "80X_mcRun2_asymptotic_2016_TrancheIV_v6"
if options.DataProcessing == "rerecodata":
  globalTag = "80X_dataRun2_2016SeptRepro_v7"
if options.DataProcessing == "promptdata":
  globalTag = "80X_dataRun2_Prompt_v16"

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

process.options = cms.untracked.PSet( SkipEvent = cms.untracked.vstring("ProductNotFound") )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.MessageLogger = cms.Service("MessageLogger")

##########################################################################################
#                                         Files                                          #
##########################################################################################
if options.DataProcessing == "mc":
  path = "root://eoscms//eos/cms/store/"+ options.DataProcessing + "/" + options.dataset  + "/" + options.sample + "/" + options.address + "/" + options.file

if options.DataProcessing == "data":
  path = "root://eoscms//eos/cms/store/"+ options.DataProcessing + "/" + options.dataset  + "/" + options.sample + "/" + options.address + "/" + options.file

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring())

process.source.fileNames.append( "file:MC_MINIAOD2.root" )
#process.source.fileNames.append( "file:03Feb2017data.root" )
###
filename_out = "outfile.root"
if options.DataProcessing == "mc":
#  filename_out = "file:/tmp/output_%s" % (options.sample + "_" + options.file)
  filename_out = "outfile.root"
if options.DataProcessing == "data":
  filename_out = "outfile.root"

#filename_out = "outfile.root"
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string(filename_out) )
process.TFileService = cms.Service("TFileService", fileName = cms.string(filename_out) )

##########################################################################################
#                                   IIHETree options                                     #
##########################################################################################
#from TrkIsoCorr.CorrectedElectronTrkisoProducers.CorrectedElectronTrkisoProducers_cfi import *
#process.CorrectedEle = CorrectedElectronTrkiso.clone()


process.load("UserCode.IIHETree.IIHETree_cfi")

triggers = "singleElectron;doubleElectron;singleMuon;singlePhoton;singleElectronSingleMuon"
process.IIHEAnalysis.triggers = cms.untracked.string(triggers)
process.IIHEAnalysis.globalTag = cms.string(globalTag)


#Track isolation correction
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoEgamma.ElectronIdentification.heepIdVarValueMapProducer_cfi")
#EGamma VID
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff','RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff']
#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

# Collections for DATA and MC.
process.IIHEAnalysis.triggerResultsCollectionHLT                 = cms.InputTag("TriggerResults"        ,""                           ,"HLT")
process.IIHEAnalysis.triggerResultsCollectionPAT                 = cms.InputTag("TriggerResults"        ,""                           ,"PAT")
process.IIHEAnalysis.triggerObjectStandAloneCollection           = cms.InputTag("selectedPatTrigger"                                        )
process.IIHEAnalysis.patTriggerCollection                        = cms.InputTag("patTrigger"                                                )
process.IIHEAnalysis.triggerEvent                                = cms.InputTag("hltTriggerSummaryAOD"  ,""                           ,"PAT")
process.IIHEAnalysis.photonCollection                            = cms.InputTag("slimmedPhotons"                                            )
process.IIHEAnalysis.electronCollection                          = cms.InputTag("slimmedElectrons"                                          )
process.IIHEAnalysis.muonCollection                              = cms.InputTag("slimmedMuons"                                              )
process.IIHEAnalysis.METCollection                               = cms.InputTag("slimmedMETs"                                               )
process.IIHEAnalysis.JetCollection                               = cms.InputTag("slimmedJets"                                               )
process.IIHEAnalysis.tauCollection                               = cms.InputTag("slimmedTaus"                                               )
process.IIHEAnalysis.superClusterCollection                      = cms.InputTag("reducedEgamma"         , "reducedSuperClusters"            )
process.IIHEAnalysis.eventRho                                    = cms.InputTag("fixedGridRhoFastjetAll"                                    )
process.IIHEAnalysis.ebReducedRecHitCollection                   = cms.InputTag("reducedEgamma"         , "reducedEBRecHits"                )
process.IIHEAnalysis.eeReducedRecHitCollection                   = cms.InputTag("reducedEgamma"         , "reducedEERecHits"                )
process.IIHEAnalysis.esReducedRecHitCollection                   = cms.InputTag("reducedEgamma"         , "reducedESRecHits"                )
process.IIHEAnalysis.PileUpSummaryInfo                           = cms.InputTag("slimmedAddPileupInfo"                                      )

# VID output
process.IIHEAnalysis.eleTrkPtIsoLabel                            = cms.InputTag("heepIDVarValueMaps"    ,"eleTrkPtIso"       ,"IIHEAnalysis")
process.IIHEAnalysis.VIDVeto                                     = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto")
process.IIHEAnalysis.VIDLoose                                    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose")
process.IIHEAnalysis.VIDMedium                                   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium")
process.IIHEAnalysis.VIDTight                                    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight")
process.IIHEAnalysis.VIDmvaEleIDwp90                             = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90")
process.IIHEAnalysis.VIDmvaEleIDwp80                             = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80")
process.IIHEAnalysis.VIDHEEP7                                    = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70")

# Collections for MC only.
process.IIHEAnalysis.generatorLabel                              = cms.InputTag("generator"                                                 )
process.IIHEAnalysis.genParticleSrc                              = cms.InputTag("prunedGenParticles"                                        )
# Collections for DATA only.
process.IIHEAnalysis.electronsBeforeGSFixCollection              = cms.InputTag("slimmedElectronsBeforeGSFix"                               )
process.IIHEAnalysis.particleFlowEGammaGSFixedCollection         = cms.InputTag("particleFlowEGammaGSFixed", "dupECALClusters"              )
process.IIHEAnalysis.ecalMultiAndGSGlobalRecHitEBCollection      = cms.InputTag("ecalMultiAndGSGlobalRecHitEB","dupESClusters"        ,"PAT")
process.IIHEAnalysis.METsMuEGCleanCollection                     = cms.InputTag("slimmedMETsMuEGClean"                                      )
process.IIHEAnalysis.discardedMuonCollection                     = cms.InputTag("packedPFCandidatesDiscarded"                               )


# Trigger matching stuff.  0.5 should be sufficient.
process.IIHEAnalysis.muon_triggerDeltaRThreshold = cms.untracked.double(0.5)
process.IIHEAnalysis.HEEP_triggerDeltaRThreshold = cms.untracked.double(0.5)

# In the absence of high ET electrons, only save events with really high Z candidates.
process.IIHEAnalysis.ZBosonZMassAcceptLower    = cms.untracked.double(850)
# Don"t bother with J/psi or Upsilon, they will only weigh us down!
process.IIHEAnalysis.ZBosonJPsiAcceptMassLower = cms.untracked.double(1e6)
process.IIHEAnalysis.ZBosonJPsiAcceptMassUpper = cms.untracked.double(1e6)
process.IIHEAnalysis.ZBosonUpsAcceptMassLower  = cms.untracked.double(1e6)
process.IIHEAnalysis.ZBosonUpsAcceptMassUpper  = cms.untracked.double(1e6)

# But make sure we save Z bosons from 50 GeV and up.
process.IIHEAnalysis.ZBosonZMassLowerCuttoff   = cms.untracked.double( 50)
process.IIHEAnalysis.ZBosonDeltaRCut           = cms.untracked.double(1e-3)

# Only save Z->ee, Z->em.
process.IIHEAnalysis.ZBosonEtThreshold = cms.untracked.double(15)
process.IIHEAnalysis.ZBosonSaveZee  = cms.untracked.bool(True )
process.IIHEAnalysis.ZBosonSaveZmm  = cms.untracked.bool(True )
process.IIHEAnalysis.ZBosonSaveZem  = cms.untracked.bool(True )
process.IIHEAnalysis.ZBosonSaveZeeg = cms.untracked.bool(False)
process.IIHEAnalysis.ZBosonSaveZmmg = cms.untracked.bool(False)


# Set pt or mass thresholds for the truth module here
# Setting thresholds reduces the size of the output files significantly
process.IIHEAnalysis.MCTruth_ptThreshold = cms.untracked.double(10.0)
process.IIHEAnalysis.MCTruth_mThreshold  = cms.untracked.double(20.0)
process.IIHEAnalysis.MCTruth_DeltaROverlapThreshold = cms.untracked.double(0.001)
# IMPORTANT         ****SKIM OBJECT****
process.IIHEAnalysis.electronPtThreshold  = cms.untracked.double(15)
process.IIHEAnalysis.muonPtThreshold      = cms.untracked.double(15)
process.IIHEAnalysis.photonPtThreshold    = cms.untracked.double(15)
process.IIHEAnalysis.jetPtThreshold       = cms.untracked.double(25)
process.IIHEAnalysis.tauPtTThreshold      = cms.untracked.double(15)

# IMPORTANT         ****SKIM EVENT****
process.IIHEAnalysis.leptonsAcceptPtThreshold    = cms.untracked.double(20)
process.IIHEAnalysis.leptonsAccept_nEle           = cms.untracked.int32(2)
process.IIHEAnalysis.leptonsAccept_nEleMu         = cms.untracked.int32(2)
process.IIHEAnalysis.leptonsAccept_nEleTau        = cms.untracked.int32(2)
process.IIHEAnalysis.leptonsAccept_nMu            = cms.untracked.int32(999)
process.IIHEAnalysis.leptonsAccept_nMuTau         = cms.untracked.int32(2)
process.IIHEAnalysis.leptonsAccept_nTau           = cms.untracked.int32(999)
#***********************************************************************

process.IIHEAnalysis.includeLeptonsAcceptModule  = cms.untracked.bool(True)
process.IIHEAnalysis.includeTriggerModule        = cms.untracked.bool(True)
process.IIHEAnalysis.includeEventModule          = cms.untracked.bool(True)
process.IIHEAnalysis.includeVertexModule         = cms.untracked.bool(True) 
process.IIHEAnalysis.includePhotonModule         = cms.untracked.bool(False)
process.IIHEAnalysis.includeElectronModule       = cms.untracked.bool(True)
process.IIHEAnalysis.includeMuonModule           = cms.untracked.bool(True)
process.IIHEAnalysis.includeMETModule            = cms.untracked.bool(True)
process.IIHEAnalysis.includeJetModule            = cms.untracked.bool(True)
process.IIHEAnalysis.includeTauModule            = cms.untracked.bool(True)
process.IIHEAnalysis.includeZBosonModule         = cms.untracked.bool(False)
process.IIHEAnalysis.includeSuperClusterModule   = cms.untracked.bool(False)
process.IIHEAnalysis.includeTracksModule         = cms.untracked.bool(False)

if ("data" in options.DataProcessing): 
    process.IIHEAnalysis.includeMCTruthModule         = cms.untracked.bool(False)
if ("mc" in options.DataProcessing):
    process.IIHEAnalysis.includeDataModule            = cms.untracked.bool(False)

#change it to true if you want to save all events
process.IIHEAnalysis.includeAutoAcceptEventModule= cms.untracked.bool(False)

process.IIHEAnalysis.debug = cms.bool(False)

##########################################################################################
#                            Woohoo!  We"re ready to start!                              #
##########################################################################################
#process.p1 = cms.Path(process.kt6PFJetsForIsolation+process.IIHEAnalysis)

#process.out = cms.OutputModule(
#    "PoolOutputModule",
#    fileName = cms.untracked.string("test.root")
#)

process.p1 = cms.Path(
    process.egmGsfElectronIDSequence *
    process.heepIDVarValueMaps    *
    process.IIHEAnalysis 
)

#process.outpath = cms.EndPath(process.out)
