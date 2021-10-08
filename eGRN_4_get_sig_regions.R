get_sig_regions <- function(runName,dbType,dbPath,bgType,background=NULL){
    source('resources/runRcisTarget_dm6_withBg.R')
    source('resources/Rcistarget_helpfunctions.R')
    source('resources/convertToTargetRegions.R')
    
    dbBg <- paste(dbType, bgType, sep="_"); dbBg # sufix for output file

    # regionSets:
    load(paste0("dars_",runName,"_darSets.RData")); length(darSets)
    regionSets <- darSets
    regionSets <- regionSets[lengths(regionSets)>=10]
    length(regionSets)
    summary(lengths(regionSets))

    #### The rest should be the same for the 3 analyses: ---
    #### Setup RcisTarget analysis (with/without background)  ----
    library(RcisTarget)

    # Convert regions...
    regionSets_db <- convertToDbRegions(regionSets, dbPath)
    allRegionsInSets <- unique(unlist(regionSets_db))

    if(!is.null(background)) 
    {
      background <- c(background, unique(unlist(regionSets)))
      background_db <- convertToDbRegions(list(background), dbPath)[[1]]
      allRegionsInSets <- unique(c(allRegionsInSets, background_db))
    }

    # Load DB
    mDb <- importRankings(dbFile=dbPath, columns=allRegionsInSets)
    if(!is.null(background)) 
    {
      message(format(Sys.time(), "%H:%M:%S"), " Re-ranking...")
      mDb <- reRank(mDb)
      message("Using background with ", getNumColsInDB(mDb), " genes/regions. Set aucMaxRank accordingly...")
    }

    ### Subset motif enrichment accordingly...
    met <- readRDS(paste0("motifs_",runName,"_enrichmentTable_1.Rds"))
    met <- met[which(met$me_BG==bgType),]; nrow(met)
    met <- met[which(met$me_DB==dbType),]; nrow(met) # only motifs
    if(nrow(met)==0) stop("met subset")

    ## TO DO: Add an extra filter for the expressed TFs can be added, we will not need the rest...

    met$setSize <- lengths(regionSets_db)[met$regionSet]
    summary(met$setSize)
    # met <- met[which(met$setSize>=10),]
    nrow(met)
    # 

    ### Run with the corresponding rocThr   ----
    aucThrs <- c("auc05"=0.05,
                  "auc01"=0.01,
                  "auc005"=0.005,
                  "auc001"=0.001) ## Check which ones are used in the final version, not all might be used...


    signifRegionsList <- list()
    signifRegionsStats <- list()
    regionSets <- sort(as.character(unique(met$regionSet)))
    # regionSets <- grep("KC - G", regionSets, value=T)
    regionSets <- regionSets[which(!regionSets %in% names(signifRegionsList))] # do not repeat...
    for(regionSetName in regionSets) #[401:800]) # [801:1154]) # Run in three pieces -in parallel-, takes very long
    {
      message(regionSetName)
      metSS <- met[met$regionSet %in% regionSetName,]
      signifRegionsList[[regionSetName]] <- list()
      signifRegionsStats[[regionSetName]] <- list()
      rocTrhs <- aucThrs[which(names(aucThrs) %in% as.character(unique(metSS$me_rocThr)))] # only the ones needed for this Set
      for(i in seq_along(rocTrhs)){
        aucThr <- rocTrhs[i]
        maxRank <- aucThr * mDb@nColsInDB; maxRank

        # Motifs significant at this ROCtrh:
        motifNames <- as.character(unlist(metSS[which(metSS$me_rocThr %in% names(aucThr)),"motif"])); length(motifNames)
        if(length(motifNames)>0)
        {
          # pdf(paste0("motifEnrichment_DAR/plotsROC_",dbBg,"/ROC_",make.names(regionSetName),"_aucTr",aucThr,"_", dbBg,".pdf"))
          tmp <- getSignificantGenes(regionSets_db[[regionSetName]],
                                     signifRankingNames=motifNames,
                                     rankings=mDb,
                                     genesFormat="geneList",
                                     plotCurve=FALSE, nCores=1,
                                     maxRank=maxRank,
                                     method="iCisTarget")
          # dev.off()

          signifRegionsStats[[regionSetName]][[names(aucThr)]] <- tmp$enrStats
          signifRegionsList[[regionSetName]][[names(aucThr)]] <- lapply(tmp$enrichedGenes, function(r) gsub("dmel_r6.02__","",r)) 
          rm(tmp)
        }else{ print("skipped")}
      }

      tmp <- signifRegionsList[[regionSetName]]
      save(tmp, file=paste0("tmp/",runName,"__", dbBg,"__",make.names(regionSetName),".RData"))
      rm(tmp)
    }
    save(signifRegionsStats, file=paste0(runName,"_signifRegionsStats_",dbBg,".RData")) 
                                    
    files_list = list.files(path = paste0("tmp/"),pattern=runName)
    files_list = files_list[grepl(x = files_list, pattern = dbBg)]
    signifRegionsList <- list()

    for(file in files_list){
        regionSetName <- gsub(paste0(".*",dbBg,"__"),"",file)
        regionSetName <- gsub(".RData","",regionSetName)
        load(paste0("tmp/",file))
        signifRegionsList[[regionSetName]] <- tmp
        rm(tmp)
    }
    save(signifRegionsList, file=paste0(runName,"_signifRegionsList_",dbBg,".RData"))
    head(names(signifRegionsList))

    aucThrs <- c("auc05"=0.05,
                  "auc01"=0.01,
                  "auc005"=0.005,
                  "auc001"=0.001) ## Check which ones are used in the final version, not all might be used...
    aucThrs <- names(aucThrs)
    signifRegionsList_byAUC <- list()
    for(aucTrh in aucThrs)
    {
      signifRegionsList_byAUC[[aucTrh]] <- list()
      for(rs in names(signifRegionsList))
      {
        signifRegionsList_byAUC[[aucTrh]][[rs]] <- signifRegionsList[[rs]][[aucTrh]]
      }
    }

    signifRegionsList.df <-  lapply(signifRegionsList_byAUC,function(x)  reshape2::melt(x))
    save(signifRegionsList.df, file=paste0(runName,"_signifRegionsList_",dbBg,".df.RData"))                                    
}
