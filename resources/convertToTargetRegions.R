# Last update 25 nov 2020: To get region locations for human/mouse
# Deprecated : convertToDbRegions

library(RcisTarget)
library(GenomicRanges)
library(BiocParallel)
# library(rtracklayer)

## for human/mouse need to use something like this:
# featherFilePath <- "/staging/leuven/res_00001/databases/cistarget/databases/all/hg19-regions-9species.all_regions.mc9nr.feather"
# dbRegionsLoc <- rtracklayer::import.bed("/staging/leuven/stg_00002/lcb/icistarget/data/regions/homo_sapiens/hg19/refseq-r45/hg19__refseq_r45__ClusteredUniformDHS_all_merge_cleaned2_features_rm-insul_rm-exons2_extend.regionid-location.bed")
## For drosophila it is enough with (the region location is on the ID): 
# featherFilePath <- "/staging/leuven/res_00001/databases/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc9nr/region_based/dm6-regions-11species.mc9nr.feather"
# dbRegionsLoc <- getDbRegionsLoc(featherFilePath)
getDbRegionsLoc <- function(featherFilePath)
{
  dbCols <- getColumnNames(featherFilePath)[-1] # skip ID
  dbCols <- setNames(dbCols, dbCols)
  
  if(grepl("__", dbCols[[1]])) dbCols <- sapply(strsplit(dbCols,"__"), function(x) x[[2]])
  dbRegions <- GRanges(unname(dbCols), name=names(dbCols))
  dbRegions <- sortSeqlevels(dbRegions)
  return(dbRegions)
}

# targetRegions <- import.bed(con=dbRegionsLocationFile)
# replaces getIctRegions_0.1.3.R
## To apply on a list of regionSets: 
# applyFun=bplapply
# regionSets_db <- applyFun(regionSets, function(x) convertToTargetRegions(queryRegions=x, targetRegions=dbRegionsLoc, minOverlap=minOverlap, returnCorrespondence=returnCorrespondence, verbose=FALSE))

convertToTargetRegions <- function(queryRegions, targetRegions, minOverlap=0.4, overlapType="any", returnCorrespondence=FALSE, verbose=TRUE)
{
  ### Check types
  if(!is.numeric(minOverlap) || (minOverlap<0 && minOverlap>=1)) stop("minOverlap should be a number between 0 and 1 (percentage of overlap between the regions).")
  if(!isClass(targetRegions, "GRanges")) targetRegions <- GRanges(targetRegions)
  if(!isClass(queryRegions, "GRanges")) queryRegions <- GRanges(queryRegions)
  ###
  
  seqlvls <- intersect(seqlevels(queryRegions), seqlevels(targetRegions))
  queryRegions <- keepSeqlevels(queryRegions, seqlvls, pruning.mode = "coarse")
  overlapHits <- findOverlaps(queryRegions, targetRegions,
                              minoverlap=1,  
                              type=overlapType, select="all", ignore.strand=TRUE)
  
  if(minOverlap>0)
  {
    # In i-cisTarget, the default is 40% minimum overlap. Both ways: It takes the maximum percentage (of the peak or the ict region)
    # To reproduce those results:
    overlaps <- pintersect(queryRegions[queryHits(overlapHits)], targetRegions[subjectHits(overlapHits)],drop.nohit.ranges=TRUE)
    percentOverlapDb <- width(overlaps) / width(targetRegions[subjectHits(overlapHits)])
    percentOverlapQuery <- width(overlaps) / width(queryRegions[queryHits(overlapHits)])
    maxOverlap <- apply(cbind(percentOverlapDb, percentOverlapQuery), 1, max)
    overlapHits <- overlapHits[maxOverlap > 0.4]
  }
  
  # Get regions names
  dbHits <- targetRegions[subjectHits(overlapHits)]
  if("name" %in% colnames(elementMetadata(dbHits)))
  {
    dbHits <- as.character(dbHits$name)
  }else{
    dbHits <- as.character(dbHits)
  }
  ret <- unique(dbHits)
  if(verbose) message(paste("Number of regions selected: ", length(ret)))
  
  if(returnCorrespondence) 
  {
    queryHits <- queryRegions[queryHits(overlapHits)]
    ret <- cbind(query=unname(as.character(queryHits)), db=unname(dbHits))
  }
  
  return(ret)
}
