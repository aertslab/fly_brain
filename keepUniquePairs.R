library(dplyr)

# For each pair Topic-Motif keep best threshold (NES)
keepUniquePairs <- function(met, col1="topic", col2="motif", colScore="NES", colsToKeep=colnames(met), verbose=TRUE)
{
  # met <- metSS[1:20,]
  met <- as_tibble(met)
  nOriginal <- nrow(met)
  met <-  met %>% 
    dplyr::select(unique(c(col1, col2, colScore, colsToKeep))) %>% 
    distinct()
  nrow(met)
  
  # table(paste(met[,col1], met[,col2], sep=" <<===>> "))
  # unique(met[,c("motif","cellType","NES")])
  
  met <- met %>% 
    group_by_at(vars(col1, col2)) %>%
      top_n(!!as.symbol(colScore), n=1) %>%
      arrange(desc(!!as.symbol(colScore))) %>%
      group_by_at(vars(col1, col2, colScore)) %>%
        dplyr::slice(1)
  nrow(met)

  met <-  met %>% 
    dplyr::select(colsToKeep) %>% 
    distinct()

  if(verbose) message("from ", nOriginal, " to ", nrow(met))
  invisible(met)
}


simplifyByTF <- function(met, simplifyCol="cellType", motifCol="motif", tfsCol="TF_expressed",  colScore="NES", isCharacter=TRUE, selectTopMotif=TRUE, randIfTies=TRUE)
{
  # met <- metSS[1:30000,]
  met <- as_tibble(met)
  nrow(met)
  
  met <- met %>%
    dplyr::select(c(simplifyCol, motifCol, tfsCol, colScore)) %>% 
    distinct()
  nrow(met)
  
  allSets <- unique(as.character(unlist(met[,simplifyCol]))); length(allSets)
  tables2show <- list()
  for(setName in allSets)
  {
    metSubset <- met %>%
      filter(!!as.symbol(simplifyCol) %in% setName)
    
    allMotifs <- unique(as.character(unlist(metSubset[,motifCol]))); length(allMotifs)
    if(nrow(metSubset) != length(allMotifs)) stop("The input table should already be simplified by 'simplifyCol'")
    rm(allMotifs)
    
    # split TFs (within cell type, it might take into account expression)
    tfsList <- setNames(unlist(as.list(unname(metSubset[,tfsCol])), recursive=F), 
                        as.character(unlist(metSubset[,motifCol])))
    if(isCharacter) tfsList <- lapply(tfsList, function(x) strsplit(x, "; ")[[1]])
    tfs.df <- reshape2::melt(tfsList)
    colnames(tfs.df) <- c("tf", motifCol)
    tfCol <- "tf"
    
    metSubset <- metSubset %>% 
      left_join(tfs.df, by=motifCol) %>% 
      dplyr::select(c(simplifyCol, tfCol, motifCol,  colScore)) %>%
      distinct()
    nrow(metSubset)
    
    # Keeps all motifs (just selects the best NES for each TF-cellType-motif)
    metSubset <- metSubset %>% 
      group_by_at(vars(simplifyCol, motifCol, tfCol)) %>%
      top_n(!!as.symbol(colScore), n=1) %>%
      distinct() 
    
    if(selectTopMotif){
      # Selects the top motif (for each CellType & TF)
      metSubset <- metSubset %>% 
        group_by_at(vars(simplifyCol, tfCol)) %>%
        top_n(!!as.symbol(colScore), n=1) %>%
        distinct() 
      if(randIfTies)
      {
        metSubset <- metSubset %>%
          arrange(desc(!!as.symbol(colScore))) %>%
          group_by_at(vars(simplifyCol, tfCol, colScore)) %>%
          dplyr::slice(1) 
      }
    }
    
    nrow(metSubset)
    tables2show[[setName]] <- metSubset
    
    rm(metSubset)
  }
  table2show <- data.table::rbindlist(tables2show)
  invisible(table2show)
}



# ### Version without dplyr
# # For each pair Topic-Motif keep best threshold (NES)
# keepUniquePairs <- function(met, col1="topic", col2="motif", colScore="NES", colsToKeep=colnames(met), verbose=TRUE)
# {
#   met <- as.data.frame(met)
#   met[,col1] <- as.character(met[,col1])
# 
#   met$uniqueCombos <- paste(met[,col1], met[,col2], sep=" <<===>> ")
#   uniqueCombos <- unique(met$uniqueCombos)
# 
#   table2show <- data.table::rbindlist(lapply(uniqueCombos, function(x){
#     x <- met[which(met$uniqueCombos %in% x),]
#     x$n <- nrow(x)
#     x[which.max(abs(unlist(x[,colScore]))),,drop=FALSE]
#   }))
#   table2show[,colScore] <- as.numeric(as.character(unlist(table2show[,colScore, with=F])))
#   table2show <- table2show[, ..colsToKeep]
# 
#   if(verbose) message("from ", nrow(met), " to ", nrow(table2show))
#   invisible(table2show)
# }


# simplifyByTF <- function(met, simplifyCol="cellType", motifCol="motif", tfsCol="TF_expressed", isCharacter=TRUE)
# {
#   allSets <- unique(as.character(unlist(met[,simplifyCol,with=F]))); length(allSets)
#   
#   tables2show <- list()
#   for(setName in allSets)
#   {
#     metSubset <- met[which(met[,simplifyCol,with=F]==setName),]
#     allMotifs <- unique(as.character(unlist(metSubset[,motifCol,with=F]))); length(allMotifs)
#     if(nrow(metSubset) != length(allMotifs)) stop("The input table should already be simplified by 'simplifyCol'")
#     
#     tfsList <- setNames(unlist(as.list(unname(metSubset[,tfsCol,with=F])), recursive=F), 
#                         as.character(unlist(metSubset[,motifCol,with=F])))
#     if(isCharacter) tfsList <- lapply(tfsList, function(x) strsplit(x, "; ")[[1]])
#     tfs.df <- reshape2::melt(tfsList)
#     colnames(tfs.df) <- c("tf", "motif")
#     allTFs <- unique(as.character(tfs.df[,"tf"])); length(allTFs)
#     tables2show[[setName]] <- data.table::rbindlist(lapply(allTFs, function(tf){
#       # message(tf)
#       motifs <- unique(tfs.df[which(tfs.df==tf),"motif"])
#       x <- metSubset[which(unlist(metSubset[,motifCol,with=F]) %in% motifs),]
#       x <- x[which.max(x$NES),,drop=FALSE]
#       
#       cbind(TF=tf, x[,c(simplifyCol, motifCol, "NES"), with=F])
#     }))
#     rm(allMotifs)
#     rm(allTFs)
#   }
#   table2show <- data.table::rbindlist(tables2show)
#   invisible(table2show)
# }