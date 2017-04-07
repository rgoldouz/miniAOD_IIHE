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
globalTag = "80"
if options.DataProcessing == "mc2016":
  globalTag = "80X_mcRun2_asymptotic_2016_TrancheIV_v8"
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
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')

process.GlobalTag.globaltag = globalTag
print "Global Tag is ", process.GlobalTag.globaltag
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

##########################################################################################
#                                         Files                                          #
##########################################################################################
if options.DataProcessing == "mc":
  path = "root://eoscms//eos/cms/store/"+ options.DataProcessing + "/" + options.dataset  + "/" + options.sample + "/" + options.address + "/" + options.file

if options.DataProcessing == "data":
  path = "root://eoscms//eos/cms/store/"+ options.DataProcessing + "/" + options.dataset  + "/" + options.sample + "/" + options.address + "/" + options.file

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
#    eventsToProcess = cms.untracked.VEventRange('1:26459:5269847')

)

#process.source.fileNames.append( "file:36CDAE89-B3BE-E611-B022-0025905B8604.root" )
#process.source.fileNames.append("file:pickevents_1.root" )
process.source.fileNames.append( "file:03Feb2017data.root" )
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

#Track isolation correction value for HEEP v7
process.load("RecoEgamma.ElectronIdentification.heepIdVarValueMapProducer_cfi")
#EGamma VID for various working points
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = ["RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff",
                 "RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff"]
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#MET corrections and uncertainties
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData= "data" in options.DataProcessing
                           )
#electron 80 energy regression
from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)
process.load("EgammaAnalysis.ElectronTools.regressionApplication_cff")

#Electron energy scale and smearing
process.load("Configuration.StandardSequences.Services_cff")
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                    calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                    engineName = cms.untracked.string("TRandom3")),
                                                  )
process.load("EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi")
process.calibratedPatElectrons.isMC = cms.bool("mc" in options.DataProcessing)
process.calibratedPatElectrons.electrons = cms.InputTag("slimmedElectrons80","","IIHEAnalysis")
process.calibratedPatElectrons.isSynchronization = cms.bool(False)
# Compatibility with VID 
process.selectedElectrons80 = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("calibratedPatElectrons","","IIHEAnalysis"),
    cut = cms.string("pt>5 && abs(eta)")
)

# Bad Charged Hadron and Bad Muon Filters from MiniAOD
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

##########################################################################################
#                            MY analysis input!                              #
##########################################################################################
process.load("UserCode.IIHETree.IIHETree_cfi")
process.IIHEAnalysis.globalTag = cms.string(globalTag)
process.IIHEAnalysis.isData  = cms.untracked.bool("data" in options.DataProcessing)
process.IIHEAnalysis.isMC    = cms.untracked.bool("mc" in options.DataProcessing)
#process.IIHEAnalysis.isMC    = cms.untracked.bool(False)
#****Collections added before the analysis
# VID output
process.IIHEAnalysis.eleTrkPtIsoLabel                            = cms.InputTag("heepIDVarValueMaps"    ,"eleTrkPtIso"       ,"IIHEAnalysis" )
process.IIHEAnalysis.VIDVeto                                     = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"  )
process.IIHEAnalysis.VIDLoose                                    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose" )
process.IIHEAnalysis.VIDMedium                                   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium")
process.IIHEAnalysis.VIDTight                                    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight" )
process.IIHEAnalysis.VIDmvaEleIDwp90                             = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90" )
process.IIHEAnalysis.VIDmvaEleIDwp80                             = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80" )
process.IIHEAnalysis.VIDHEEP7                                    = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"                   )
# Collections for DATA only.
process.IIHEAnalysis.particleFlowEGammaGSFixedCollection         = cms.InputTag("particleFlowEGammaGSFixed", "dupECALClusters"              )
process.IIHEAnalysis.ecalMultiAndGSGlobalRecHitEBCollection      = cms.InputTag("ecalMultiAndGSGlobalRecHitEB","dupESClusters"        ,"PAT")
process.IIHEAnalysis.METsMuEGCleanCollection                     = cms.InputTag("slimmedMETsMuEGClean"                                      )
process.IIHEAnalysis.discardedMuonCollection                     = cms.InputTag("packedPFCandidatesDiscarded"                               )

#use 80 regression + scale/smearing for electron
process.IIHEAnalysis.electronCollection80    = cms.InputTag("selectedElectrons80","","IIHEAnalysis")

#jet smeared collection
process.IIHEAnalysis.JetCollection                   = cms.InputTag("basicJetsForMet" ,"","IIHEAnalysis")
process.IIHEAnalysis.JetCollectionSmeared            = cms.InputTag("patSmearedJets"               ,"","IIHEAnalysis")
process.IIHEAnalysis.JetCollectionEnUp               = cms.InputTag("shiftedPatJetEnUp"            ,"","IIHEAnalysis")
process.IIHEAnalysis.JetCollectionEnDown             = cms.InputTag("shiftedPatJetEnDown"          ,"","IIHEAnalysis")
process.IIHEAnalysis.JetCollectionSmearedJetResUp    = cms.InputTag("shiftedPatSmearedJetResUp"    ,"","IIHEAnalysis")
process.IIHEAnalysis.JetCollectionSmearedJetResDown  = cms.InputTag("shiftedPatSmearedJetResDown"  ,"","IIHEAnalysis")

#MET collections
process.IIHEAnalysis.patPFMetCollection                        = cms.InputTag("patPFMet"                  , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1Collection                      = cms.InputTag("patPFMetT1"                , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1JetEnDownCollection             = cms.InputTag("patPFMetT1JetEnDown"       , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1JetEnUpCollection               = cms.InputTag("patPFMetT1JetEnUp"         , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1SmearCollection                 = cms.InputTag("patPFMetT1Smear"           , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1SmearJetEnDownCollection        = cms.InputTag("patPFMetT1SmearJetEnDown"  , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1SmearJetEnUpCollection          = cms.InputTag("patPFMetT1SmearJetEnUp"    , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1SmearJetResDownCollection       = cms.InputTag("patPFMetT1SmearJetResDown" , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1SmearJetResUpCollection         = cms.InputTag("patPFMetT1SmearJetResUp"   , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1TxyCollection                   = cms.InputTag("patPFMetT1Txy"             , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetFinalCollection                   = cms.InputTag("slimmedMETs"               , ""                ,"IIHEAnalysis"  )

process.IIHEAnalysis.includeLeptonsAcceptModule  = cms.untracked.bool(True)
process.IIHEAnalysis.includeTriggerModule        = cms.untracked.bool(True)
process.IIHEAnalysis.includeEventModule          = cms.untracked.bool(True)
process.IIHEAnalysis.includeVertexModule         = cms.untracked.bool(True)
process.IIHEAnalysis.includeElectronModule       = cms.untracked.bool(True)
process.IIHEAnalysis.includeMuonModule           = cms.untracked.bool(True)
process.IIHEAnalysis.includeMETModule            = cms.untracked.bool(True)
process.IIHEAnalysis.includeJetModule            = cms.untracked.bool(True)
process.IIHEAnalysis.includeTauModule            = cms.untracked.bool(True)
process.IIHEAnalysis.includeMCTruthModule        = cms.untracked.bool("mc" in options.DataProcessing)
process.IIHEAnalysis.includeDataModule            = cms.untracked.bool("data" in options.DataProcessing)

##########################################################################################
#                            Woohoo!  We"re ready to start!                              #
##########################################################################################
#process.p1 = cms.Path(process.kt6PFJetsForIsolation+process.IIHEAnalysis)
#process.out = cms.OutputModule(
#    "PoolOutputModule",
#    fileName = cms.untracked.string("test.root")
#    )

process.p1 = cms.Path(
    process.regressionApplication     *
    process.calibratedPatElectrons    *
    process.selectedElectrons80       *
    process.egmGsfElectronIDSequence  * 
    process.heepIDVarValueMaps        *
    process.BadPFMuonFilter           *
    process.BadChargedCandidateFilter *
    process.fullPatMetSequence        *
    process.IIHEAnalysis 
    )

#process.outpath = cms.EndPath(process.out)

