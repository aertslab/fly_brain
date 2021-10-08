do_me <- function(darsets,runName){
    packageVersion("RcisTarget") # >1.5.2 required
    packageVersion("BiocParallel") # >1.5.2 required
    source('resources/runRcisTarget_dm6_withBg.R')
    source('resources/Rcistarget_helpfunctions.R')
    source('resources/convertToTargetRegions.R')
    
    load(darsets)
    regionSets <- darSets; rm(darSets)
    regionSets <- regionSets[lengths(regionSets)>=10]
    cbind(lengths(regionSets))
    length(regionSets)

    
    nCores=10

    # with bg:
    load("resources/topicRegions_all_development_230topics.RData", verbose = T)
    allRegionsInTopics <- unique(unlist(topicRegions))
    length(allRegionsInTopics)

    #motifEnrichment <- runRcisTarget_dm6(regionSets, background=allRegionsInTopics, aucMaxRank=round(0.01 * length(allRegionsInTopics)))
    #saveRDS(motifEnrichment, file=paste0(runName,"_motifEnrichment_withBackground_aucMax01.Rds"))

    motifEnrichment <- runRcisTarget_dm6(regionSets, background=allRegionsInTopics, aucMaxRank=round(0.05 * length(allRegionsInTopics)))
    saveRDS(motifEnrichment, file=paste0(runName,"_motifEnrichment_withBackground_aucMax05.Rds"))

    motifEnrichment <- runRcisTarget_dm6(regionSets, background=allRegionsInTopics, aucMaxRank=round(0.10 * length(allRegionsInTopics)))
    saveRDS(motifEnrichment, file=paste0(runName,"_motifEnrichment_withBackground_aucMax10.Rds"))

    motifEnrichment <- runRcisTarget_dm6(regionSets, background=allRegionsInTopics, aucMaxRank=round(0.005 * length(allRegionsInTopics)))
    saveRDS(motifEnrichment, file=paste0(runName,"_motifEnrichment_withBackground_aucMax005.Rds"))
    
    ## No bg
    nRegions <- feather::feather_metadata("/staging/leuven/res_00001/databases/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc9nr/region_based/dm6-regions-11species.mc9nr.feather")$dim[2]
    motifEnrichment <- runRcisTarget_dm6(regionSets, aucMaxRank=round(0.01 * nRegions), nesThreshold=3)
    saveRDS(motifEnrichment, file=paste0(runName,"_motifEnrichment_noBg_aucMax01.Rds"))
    motifEnrichment <- runRcisTarget_dm6(regionSets, aucMaxRank=round(0.05 * nRegions), nesThreshold=3)
    saveRDS(motifEnrichment, file=paste0(runName,"_motifEnrichment_noBg_aucMax05.Rds"))
    motifEnrichment <- runRcisTarget_dm6(regionSets, aucMaxRank=round(0.001 * nRegions), nesThreshold=3)
    saveRDS(motifEnrichment, file=paste0(runName,"_motifEnrichment_noBg_aucMax001.Rds"))
    motifEnrichment <- runRcisTarget_dm6(regionSets, aucMaxRank=round(0.005 * nRegions), nesThreshold=3)
    saveRDS(motifEnrichment, file=paste0(runName,"_motifEnrichment_noBg_aucMax005.Rds"))

    # ChIP
    nRegions <- feather::feather_metadata("/staging/leuven/stg_00002/lcb/saibar/Projects/epiSCENIC/2019-06_ChipDB/Chip_dbs_v2_20190627/encode_modERN_20190621__ChIP_seq.max.feather")$dim[2]; nRegions
    chipEnrichment <- runRcisTarget_dm6(regionSets,
                                        aucMaxRank=round(0.03 * nRegions),
                                        dbPath = "resources/encode_modERN_20190621__ChIP_seq.max.feather",
                                        motifAnnotation = "resources/encode_modERN_20190621_dm6_annotation.tbl")
    saveRDS(chipEnrichment, file=paste0(runName,"_chipEnrichment_aucMax03.Rds"))
    chipEnrichment <- runRcisTarget_dm6(regionSets,
                                        aucMaxRank=round(0.01 * nRegions),
                                        dbPath = "resources/encode_modERN_20190621__ChIP_seq.max.feather",
                                        motifAnnotation = "resources/encode_modERN_20190621_dm6_annotation.tbl")
    saveRDS(chipEnrichment, file=paste0(runName,"_chipEnrichment_aucMax01.Rds"))

}