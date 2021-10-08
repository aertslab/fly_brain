library(RcisTarget) 
### For version 1.7, load:
# source('/ddn1/vol1/staging/leuven/stg_00002/lcb/saibar/Projects/aux_scripts/cisTarget/convertToTargetRegions.R')
# source('/ddn1/vol1/staging/leuven/stg_00002/lcb/saibar/Projects/aux_scripts/cisTopic/keepUniquePairs.R')
###  (not needed on >=1.11)

library(BiocParallel)
library(GenomicRanges)
# Stats annotation only to TFs on 13 nov 2020

runRcisTarget_dm6 <- function(regionSets, 
                              background=NULL,
                              dbPath="/staging/leuven/res_00001/databases/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc9nr/region_based/dm6-regions-11species.mc9nr.feather",
                              motifAnnotation="/staging/leuven/res_00001/databases/cistarget/motif2tf/motifs-v9-nr.flybase-m0.001-o0.0.tbl",
                              nCores=20,
                              tfFileName="/staging/leuven/stg_00002/lcb/saibar/Projects/aux_resources/TF_lists/dmel/tfs_annotation_dmel_mcv89chip.RData",
                              nesThreshold=3,
                              aucMaxRank=5000,
                              spltChr="__",
                              rnkIndexCol="features") # or motifs or...
{
  if(require(BiocParallel, quietly=TRUE)) # Without this, it is extremely slow!!
  {
    mine <- SnowParam(workers=20, RNGseed=123, type="FORK", stop.on.error=FALSE, progressbar=TRUE) # all platforms
    register(mine, default=TRUE); bpparam(); rm(mine)
  }
  
  ##### Run info:
  message("db: \t", dbPath)
  if((length(motifAnnotation)==1) & all(grepl(".tbl",motifAnnotation))) {
    message("annotation: \t", motifAnnotation)
    motifAnnotation <- RcisTarget::importAnnotations(motifAnnotation)
  }
  idColumn <- colnames(motifAnnotation)[1]
  
  # nested list?
  if(is.list(regionSets[[1]]))
  {
    names(regionSets) <- NULL 
    regionSets <- unlist(regionSets, recursive=FALSE)
    message("First geneset name: \t", names(regionSets)[1])
  }
  
  ##### Convert regions...
  dbRegionsLoc <- getDbRegionsLoc(dbPath)
  regionSets_db <- bplapply(regionSets, function(x) convertToTargetRegions(queryRegions=x, targetRegions=dbRegionsLoc, verbose=FALSE))
  # regionSets_db <- convertToDbRegions(regionSets, dbPath) # replaced nov 25
  
  allRegionsInTopics <- unique(unlist(regionSets_db))
  message("# regions in db in all topics: ", length(allRegionsInTopics))
  print(summary(lengths(regionSets_db)))
  
  if(!is.null(background)) 
  {
    background <- c(background, unique(unlist(regionSets)))
    background_db <- convertToTargetRegions(queryRegions=background, targetRegions=dbRegionsLoc, verbose=FALSE)
    allRegionsInTopics <- unique(c(allRegionsInTopics, background_db))
  }
  
  ##### Run RcisTarget
  message(format(Sys.time(), "%H:%M:%S"), " Importing rankings")
  motifRankings <- importRankings(dbFile=dbPath, columns=allRegionsInTopics, indexCol=rnkIndexCol)
  if(!is.null(background)) 
  {
    message(format(Sys.time(), "%H:%M:%S"), " Re-ranking...")
    motifRankings <- reRank(motifRankings)
    message("Using background with ", getNumColsInDB(motifRankings), " genes/regions. Set aucMaxRank accordingly...")
  }
  
  message(format(Sys.time(), "%H:%M:%S"), " Calculating AUC")
  auc_topics <- calcAUC(regionSets_db, motifRankings, nCores=nCores, aucMaxRank=aucMaxRank)
  message("AUC dims: ", nrow(auc_topics), "x", ncol(auc_topics))
  
  #####
  # Annotation to TFs
  # Default annotation (but all confidences in same column)
  message(format(Sys.time(), "%H:%M:%S"), " Adding annotation")
  allAnnots <- c("directAnnotation", "inferredBy_Orthology", "inferredBy_MotifSimilarity","inferredBy_MotifSimilarity_n_Orthology")
  allAnnots <- allAnnots[which(allAnnots %in% motifAnnotation$annotationSource)]
  motifEnrichmentTable <- addMotifAnnotation(auc_topics, motifAnnot=motifAnnotation, nesThreshold=nesThreshold,
                                             motifAnnot_highConfCat=allAnnots, 
                                             motifAnnot_lowConfCat=NULL,
                                             idColumn=idColumn)
  colnames(motifEnrichmentTable)[which(colnames(motifEnrichmentTable)=="TF_highConf")] <- "any_annot"
  if(grepl(" @ ", motifEnrichmentTable$geneSet[1])) {
    motifEnrichmentTable <- data.frame(data.frame(do.call(rbind, strsplit(motifEnrichmentTable$geneSet,  " @ ", fixed=T))), motifEnrichmentTable)
    colnames(motifEnrichmentTable)[1:2] <- c("topic","binTr")
    motifEnrichmentTable <- motifEnrichmentTable[,-3]
  }
  head(motifEnrichmentTable)
  
  ### Add TF-only as list:
  if(!is.null(tfFileName))
  {
    message(format(Sys.time(), "%H:%M:%S"), " Discarding non-TFs")
    # Remove annotations to non-TFs from DB: 
    load(tfFileName)
    tfs <- tfs_annotation[which(tfs_annotation$tfAnnotationPriority <= 1),]; nrow(tfs)
    motifAnnotation <- motifAnnotation[which(as.character(motifAnnotation$TF) %in% rownames(tfs)),]
  }
  
  # Add to table:
  TF_annot <- sapply(getMotifAnnotation(motifEnrichmentTable[[idColumn]], motifAnnot=motifAnnotation, annotCats=allAnnots,idColumn=idColumn, returnFormat="list"), paste, collapse="; ") # TODO: Could be saved as list instead of "character"
  motifEnrichmentTable$TF_annot <- ""
  motifEnrichmentTable$TF_annot <- TF_annot[motifEnrichmentTable[[idColumn]]]
  motifEnrichmentTable$TF_annot[is.na(motifEnrichmentTable$TF_annot)] <- ""
  
  message(format(Sys.time(), "%H:%M:%S"), " Finished")
  return(list(motifEnrichmentTable=motifEnrichmentTable, auc=auc_topics))
}
