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
globalTag = "80X_mcRun2_asymptotic_2016_TrancheIV_v6"
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
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')

process.GlobalTag.globaltag = globalTag
print "Global Tag is ", process.GlobalTag.globaltag
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
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
    fileNames = cms.untracked.vstring())

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
process.load("UserCode.IIHETree.IIHETree_cfi")
process.IIHEAnalysis.globalTag = cms.string(globalTag)
process.IIHEAnalysis.isData  = cms.untracked.bool("data" in options.DataProcessing)
process.IIHEAnalysis.isMC    = cms.untracked.bool("mc" in options.DataProcessing)

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
process.IIHEAnalysis.includeMCTruthModule        = cms.untracked.bool("mc" in options.DataProcessing)
process.IIHEAnalysis.includeDataModule            = cms.untracked.bool("data" in options.DataProcessing)



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

process.load("RecoEgamma.ElectronIdentification.ElectronIDValueMapProducer_cfi")
process.egmGsfElectronIDs.physicsObjectSrc            = cms.InputTag("selectedElectrons80","","IIHEAnalysis")
process.electronIDValueMapProducer.srcMiniAOD         = cms.InputTag("selectedElectrons80","","IIHEAnalysis")
process.electronRegressionValueMapProducer.srcMiniAOD = cms.InputTag("selectedElectrons80","","IIHEAnalysis")
process.electronMVAValueMapProducer.srcMiniAOD        = cms.InputTag("selectedElectrons80","","IIHEAnalysis")


##########################################################################################
#                            Woohoo!  We"re ready to start!                              #
##########################################################################################
#process.p1 = cms.Path(process.kt6PFJetsForIsolation+process.IIHEAnalysis)
process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string("test.root")
    )

process.p1 = cms.Path(
    process.fullPatMetSequence       *
    process.regressionApplication    *
    process.calibratedPatElectrons   *
    process.selectedElectrons80      *
    process.egmGsfElectronIDSequence * 
    process.heepIDVarValueMaps       *
    process.IIHEAnalysis 
    )

process.outpath = cms.EndPath(process.out)

