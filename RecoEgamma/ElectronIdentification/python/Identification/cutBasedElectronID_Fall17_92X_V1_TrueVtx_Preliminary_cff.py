from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry

import FWCore.ParameterSet.Config as cms

# Common functions and classes for ID definition are imported here:
from RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_tools \
    import ( EleWorkingPoint_V3,
             IsolationCutInputs_V2,
             configureVIDCutBasedEleID_V3 )

#
# The ID cuts below are optimized IDs on Spring16 simulation with 80X-based production
# The cut values are taken from the twiki:
#       https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
#       (where they may not stay, if a newer version of cuts becomes available for these
#        conditions)
# See also the presentation explaining these working points:
#        https://indico.cern.ch/event/662751/contributions/2778044/attachments/1562080/2459801/171121_egamma_workshop.pdf
#
# First, define cut values
#

# Veto working point Barrel and Endcap
idName = "cutBasedElectronID-Fall17-92X-V1-Preliminary-veto"
WP_Veto_EB = EleWorkingPoint_V3(
    idName                         = idName  , # idName
    full5x5_sigmaIEtaIEtaCut       = 0.0128  , # full5x5_sigmaIEtaIEtaCut
    dEtaInSeedCut                  = 0.00523 , # dEtaInSeedCut
    dPhiInCut                      = 0.159   , # dPhiInCut
    hOverECut                      = 0.247   , # hOverECut
    relCombIsolationWithEALowPtCut = 0.168   , # relCombIsolationWithEALowPtCut
    relCombIsolationWithEAHighPtCut= 0.168   , # relCombIsolationWithEAHighPtCut
    absEInverseMinusPInverseCut    = 0.193   , # absEInverseMinusPInverseCut
    # conversion veto cut needs no parameters, so not mentioned
    missingHitsCut                 = 2          # missingHitsCut
    )

WP_Veto_EE = EleWorkingPoint_V3(
    idName                         = idName  , # idName
    full5x5_sigmaIEtaIEtaCut       = 0.0445  , # full5x5_sigmaIEtaIEtaCut
    dEtaInSeedCut                  = 0.00984 , # dEtaInSeedCut
    dPhiInCut                      = 0.157   , # dPhiInCut
    hOverECut                      = 0.0982   , # hOverECut
    relCombIsolationWithEALowPtCut = 0.185   , # relCombIsolationWithEALowPtCut
    relCombIsolationWithEAHighPtCut= 0.185   , # relCombIsolationWithEAHighPtCut
    absEInverseMinusPInverseCut    = 0.0962   , # absEInverseMinusPInverseCut
    # conversion veto cut needs no parameters, so not mentioned
    missingHitsCut                 = 3          # missingHitsCut
    )

# Loose working point Barrel and Endcap
idName = "cutBasedElectronID-Fall17-92X-V1-Preliminary-loose"
WP_Loose_EB = EleWorkingPoint_V3(
    idName                         = idName  , # idName
    full5x5_sigmaIEtaIEtaCut       = 0.0105  , # full5x5_sigmaIEtaIEtaCut
    dEtaInSeedCut                  = 0.00387 , # dEtaInSeedCut
    dPhiInCut                      = 0.0716   , # dPhiInCut
    hOverECut                      = 0.236   , # hOverECut
    relCombIsolationWithEALowPtCut = 0.133  , # relCombIsolationWithEALowPtCut
    relCombIsolationWithEAHighPtCut= 0.133  , # relCombIsolationWithEAHighPtCut
    absEInverseMinusPInverseCut    = 0.129   , # absEInverseMinusPInverseCut
    # conversion veto cut needs no parameters, so not mentioned
    missingHitsCut                 = 1          # missingHitsCut
    )

WP_Loose_EE = EleWorkingPoint_V3(
    idName                         = idName  , # idName
    full5x5_sigmaIEtaIEtaCut       = 0.0356  , # full5x5_sigmaIEtaIEtaCut
    dEtaInSeedCut                  = 0.0072 , # dEtaInSeedCut
    dPhiInCut                      = 0.147   , # dPhiInCut
    hOverECut                      = 0.0801   , # hOverECut
    relCombIsolationWithEALowPtCut = 0.146   , # relCombIsolationWithEALowPtCut
    relCombIsolationWithEAHighPtCut= 0.146   , # relCombIsolationWithEAHighPtCut
    absEInverseMinusPInverseCut    = 0.0875   , # absEInverseMinusPInverseCut
    # conversion veto cut needs no parameters, so not mentioned
    missingHitsCut                 = 1         # missingHitsCut
    )

# Medium working point Barrel and Endcap
idName = "cutBasedElectronID-Fall17-92X-V1-Preliminary-medium"
WP_Medium_EB = EleWorkingPoint_V3(
    idName                         = idName , # idName
    full5x5_sigmaIEtaIEtaCut       = 0.0105, # full5x5_sigmaIEtaIEtaCut
    dEtaInSeedCut                  = 0.00365, # dEtaInSeedCut
    dPhiInCut                      = 0.0588  , # dPhiInCut
    hOverECut                      = 0.0859  , # hOverECut
    relCombIsolationWithEALowPtCut = 0.0718 , # relCombIsolationWithEALowPtCut
    relCombIsolationWithEAHighPtCut= 0.0718 , # relCombIsolationWithEAHighPtCut
    absEInverseMinusPInverseCut    = 0.0327  , # absEInverseMinusPInverseCut
    # conversion veto cut needs no parameters, so not mentioned
    missingHitsCut                 = 1          # missingHitsCut
    )

WP_Medium_EE = EleWorkingPoint_V3(
    idName                         = idName  , # idName
    full5x5_sigmaIEtaIEtaCut       = 0.0309  , # full5x5_sigmaIEtaIEtaCut
    dEtaInSeedCut                  = 0.00625 , # dEtaInSeedCut
    dPhiInCut                      = 0.0355  , # dPhiInCut
    hOverECut                      = 0.0604  , # hOverECut
    relCombIsolationWithEALowPtCut = 0.143  , # relCombIsolationWithEALowPtCut
    relCombIsolationWithEAHighPtCut= 0.143  , # relCombIsolationWithEAHighPtCut
    absEInverseMinusPInverseCut    = 0.0335   , # absEInverseMinusPInverseCut
    # conversion veto cut needs no parameters, so not mentioned
    missingHitsCut                 = 1          # missingHitsCut
    )

# Tight working point Barrel and Endcap
idName = "cutBasedElectronID-Fall17-92X-V1-Preliminary-tight"
WP_Tight_EB = EleWorkingPoint_V3(
    idName                         = idName    , # idName
    full5x5_sigmaIEtaIEtaCut       = 0.0104   , # full5x5_sigmaIEtaIEtaCut
    dEtaInSeedCut                  = 0.00353   , # dEtaInSeedCut
    dPhiInCut                      = 0.0499    , # dPhiInCut
    hOverECut                      = 0.0833    , # hOverECut
    relCombIsolationWithEALowPtCut = 0.0361    , # relCombIsolationWithEALowPtCut
    relCombIsolationWithEAHighPtCut= 0.0361    , # relCombIsolationWithEAHighPtCut
    absEInverseMinusPInverseCut    = 0.0278    , # absEInverseMinusPInverseCut
    # conversion veto cut needs no parameters, so not mentioned
    missingHitsCut                 = 1          # missingHitsCut
    )

WP_Tight_EE = EleWorkingPoint_V3(
    idName                         = idName  , # idName
    full5x5_sigmaIEtaIEtaCut       = 0.0305  , # full5x5_sigmaIEtaIEtaCut
    dEtaInSeedCut                  = 0.00567 , # dEtaInSeedCut
    dPhiInCut                      = 0.0165  , # dPhiInCut
    hOverECut                      = 0.0543  , # hOverECut
    relCombIsolationWithEALowPtCut = 0.094  , # relCombIsolationWithEALowPtCut
    relCombIsolationWithEAHighPtCut= 0.094  , # relCombIsolationWithEAHighPtCut
    absEInverseMinusPInverseCut    = 0.0158 , # absEInverseMinusPInverseCut
    # conversion veto cut needs no parameters, so not mentioned
    missingHitsCut                 = 1          # missingHitsCut
    )

# Second, define what effective areas to use for pile-up correction
isoInputs = IsolationCutInputs_V2(
    # phoIsolationEffAreas
    "RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt"
)


#
# Set up VID configuration for all cuts and working points
#

cutBasedElectronID_Fall17_92X_V1_Preliminary_veto = configureVIDCutBasedEleID_V3(WP_Veto_EB, WP_Veto_EE, isoInputs)
cutBasedElectronID_Fall17_92X_V1_Preliminary_loose = configureVIDCutBasedEleID_V3(WP_Loose_EB, WP_Loose_EE, isoInputs)
cutBasedElectronID_Fall17_92X_V1_Preliminary_medium = configureVIDCutBasedEleID_V3(WP_Medium_EB, WP_Medium_EE, isoInputs)
cutBasedElectronID_Fall17_92X_V1_Preliminary_tight = configureVIDCutBasedEleID_V3(WP_Tight_EB, WP_Tight_EE, isoInputs)


# The MD5 sum numbers below reflect the exact set of cut variables
# and values above. If anything changes, one has to 
# 1) comment out the lines below about the registry, 
# 2) run "calculateMD5 <this file name> <one of the VID config names just above>
# 3) update the MD5 sum strings below and uncomment the lines again.
#

#central_id_registry.register(cutBasedElectronID_Fall17_92X_V1_Preliminary_veto.idName,
#                             '0025c1841da1ab64a08d703ded72409b')
#central_id_registry.register(cutBasedElectronID_Fall17_92X_V1_Preliminary_loose.idName,
#                             'c1c4c739f1ba0791d40168c123183475')
#central_id_registry.register(cutBasedElectronID_Fall17_92X_V1_Preliminary_medium.idName,
#                             '71b43f74a27d2fd3d27416afd22e8692')
#central_id_registry.register(cutBasedElectronID_Fall17_92X_V1_Preliminary_tight.idName,
#                             'ca2a9db2976d80ba2c13f9bfccdc32f2')


### for now until we have a database...
cutBasedElectronID_Fall17_92X_V1_Preliminary_veto.isPOGApproved = cms.untracked.bool(False)
cutBasedElectronID_Fall17_92X_V1_Preliminary_loose.isPOGApproved = cms.untracked.bool(False)
cutBasedElectronID_Fall17_92X_V1_Preliminary_medium.isPOGApproved = cms.untracked.bool(False)
cutBasedElectronID_Fall17_92X_V1_Preliminary_tight.isPOGApproved = cms.untracked.bool(False)
