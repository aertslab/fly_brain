library(GenomicRanges)
# library(rtracklayer)

# Retrieves the DB regions that intersect with the query regions
# Loads the DB regions & calls convertToTargetRegions for each regionSet
# regionSets: List of GRanges or character vectors
convertToDbRegions <- function(regionSets, featherFilePath, minOverlap=0.4, returnCorrespondence=FALSE, applyFun=bplapply, suppressWarning=FALSE, verbose=TRUE)
{
  ### Check types
  if(!is.numeric(minOverlap) || (minOverlap<0 && minOverlap>=1)) stop("minOverlap should be a number between 0 and 1 (percentage of overlap between the regions).")
  ###

  # Get DB regions:
  dbCols <- getColumnNames(featherFilePath)
  dbCols <- setNames(dbCols, dbCols)
  if(grepl("__",dbCols[[1]])) dbCols <- sapply(strsplit(dbCols,"__"), function(x) x[[2]])
  dbRegions <- GRanges(unname(dbCols), name=names(dbCols))
  dbRegions <- sortSeqlevels(dbRegions)
  if(verbose) message(paste("# regions in the database:", length(dbRegions)))

  # Query regions to list if needed
  if(!is.list(regionSets)){
    regionSets <- list(regionSets)
    if(!suppressWarning) warning("Single region set converted to list.")
  }

  # Convert
  ret <- applyFun(regionSets, function(x) convertToTargetRegions(queryRegions=x, targetRegions=dbRegions, minOverlap=minOverlap, returnCorrespondence=returnCorrespondence, verbose=FALSE))
  # if(!rs_isList & !returnCorrespondence) ret <- unlist(ret, recursive=FALSE) # does not work well for returnCorrespondence

  return(ret)
}


# targetRegions <- import.bed(con=dbRegionsLocationFile)
# replaces getIctRegions_0.1.3.R
convertToTargetRegions <- function(queryRegions, targetRegions, minOverlap=0.4, overlapType="any", returnCorrespondence=FALSE, verbose=TRUE)
{
  ### Check types
  if(!is.numeric(minOverlap) || (minOverlap<0 && minOverlap>=1)) stop("minOverlap should be a number between 0 and 1 (percentage of overlap between the regions).")
  if(!isClass(targetRegions, "GRanges")) targetRegions <- GRanges(targetRegions)
  if(!isClass(queryRegions, "GRanges")) queryRegions <- GRanges(queryRegions)
  ###

  seqlvls <- intersect(seqlevels(queryRegions), seqlevels(targetRegions))
  queryRegions <- keepSeqlevels(queryRegions, seqlvls)
  overlapHits <- findOverlaps(queryRegions, targetRegions,
                              minoverlap=1,
                              type=overlapType, select="all", ignore.strand=TRUE)

  if(minOverlap>0)
  {
    # In i-cisTarget, the default is 40% minimum overlap. Both ways: It takes the maximum percentage (of the peak or the ict region)
    # To reproduce those results:
    overlaps <- pintersect(queryRegions[queryHits(overlapHits)], targetRegions[subjectHits(overlapHits)])
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