import FWCore.ParameterSet.Config as cms

from DQMOffline.Trigger.zjetsmonitoring_cfi import zjetsmonitoring

hltZJetsmonitoring = zjetsmonitoring.clone(
    FolderName = 'HLT/JME/ZplusJets',

    ptcut = 20.,
    ## TriggerResultsLabel = add this ,
    ## triggerEventObject = add this ,
    ## pathname = add this ,
    ## pathmodule = add this 
)

