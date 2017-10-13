from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()


version  = '{version}'  # PARAMETER
suffix   = '{suffix}'   # PARAMETER for ext datasets of MC (they don't differ in primary name)
dtag     = '{dtag}'     # PARAMETER
dset     = '{dset}'     # PARAMETER
LumiMask = '{LumiMask}' # PARAMETER
config_file = '{config_file}' # PARAMETER

request_tag = 'Ntupler_' + version + '_' + dtag + suffix # apparently it's the only place to distinguish ext datasets of MC for crab

config.General.requestName = request_tag
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = config_file

config.Data.inputDataset = dset
config.Data.inputDBS = 'global'
if not LumiMask:
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 10
else:
    config.Data.splitting = 'LumiBased'
    config.Data.unitsPerJob = 20
    config.Data.lumiMask = LumiMask

config.Data.outLFNDirBase = '/store/user/%s/%s/' % (getUsernameFromSiteDB(), version)
config.Data.publication = False
config.Data.outputDatasetTag = request_tag

config.Site.storageSite = 'T2_PT_NCG_Lisbon'

