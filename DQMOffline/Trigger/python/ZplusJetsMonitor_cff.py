import FWCore.ParameterSet.Config as cms

from DQMOffline.Trigger.ZplusJetsMonitor_cfi import hltZJetsmonitoring

#---- define trigger paths and module paths ---#

DiMuonMass8_Prommonitoring = hltZJetsmonitoring.clone(
    FolderName = 'HLT/JME/ZplusJets/dimuMass8',
    # pathname = 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v9' 
)

HLTZplusJetmonitoring = cms.Sequence(
    DiMuonMass8_Prommonitoring
        )
