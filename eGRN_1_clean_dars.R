# Load DAR analyses
dars <- list()

dars[["vsClosest"]] <- readRDS(paste0("dars_",runName,"_treeContrasts_wGenes_renamedToTypes.Rds"))
dars[["vsRest"]] <- readRDS(paste0("dars_",runName,"_wGenes.Rds"))

dars <- dars[sort(names(dars))]
cbind(sapply(dars, nrow))

cbind(sapply(dars, function(x) length(unique(x$cluster))))

             #### Rename clusters -----
sort(unique(unlist(lapply(dars, function(x) unique(x$cluster)))))
for(i in seq_along(dars))
{
  dars[[i]]$cluster <- gsub("__vsClosest", "",  dars[[i]]$cluster)
}
sort(unique(unlist(lapply(dars, function(x) unique(x$cluster)))))

                          ## Save analysis settings, and add to regionSet name ----
for(i in seq_along(dars))
{
  colnames(dars[[i]])[colnames(dars[[i]])=="cluster"] <- "cellType" 
  dars[[i]]$DAR_contrast <- ""
  dars[[i]]$DAR_within <- ""
  dars[[i]]$regionSet <- ""
}
sort(unique(unlist(lapply(dars, function(x) unique(x$cellType)))))

for(i in grep("vsRest", names(dars), value=T)) dars[[i]]$DAR_contrast <- "vsRest"
for(i in grep("vsClosest", names(dars), value=T)) dars[[i]]$DAR_contrast <- "vsClosestCluster"
sort(unique(unlist(lapply(dars, function(x) unique(x$DAR_contrast)))))


for(i in seq_along(dars))
  dars[[i]]$regionSet <- paste0(dars[[i]]$cellType,"--", dars[[i]]$clResolution, "_", gsub("Cluster", "", dars[[i]]$DAR_contrast), "_w", dars[[i]]$DAR_within)

sort(unique(unlist(lapply(dars, function(x) unique(x$regionSet)))))

tmp <- table(unlist(sapply(dars, function(x) unique(x$regionSet))))
if(any(tmp>1)) stop("clusterIDs are not unique!") # Should be none
rm(tmp)

                           #### Prepare contrast sets -----
dars_DW <- lapply(dars, function(x) x[x$avg_logFC < 0,])
dars <- lapply(dars, function(x) x[x$avg_logFC>0,])

# Exclusive: 
dars_exclusive <- lapply(dars, function(darSS)
{
  cntRegion <- table(darSS$region)
  darSS <- darSS[which(darSS$region %in% names(cntRegion)[which(cntRegion==1)]),]
  darSS$regionSet <- paste0(darSS$regionSet, "___exclusive")
  darSS
})
# sort(sapply(dars, function(x) length(unique(x$regionSet))))
# Remove those which make no sense...

dars_DW <- dars_DW[sapply(dars_DW, nrow)>0]
dars_DW <- lapply(dars_DW, function(darSS)
{
  darSS$regionSet <- paste0(darSS$regionSet, "___DW")
  darSS
})

darsMerged <- c(dars_DW, dars_exclusive, dars)
darsMerged <- do.call(rbind, darsMerged)
save(darsMerged, file=paste0("dars_",runName,"_Merged.df.RData"))

# Region sets ----
darSets <- split(darsMerged$region, factor(darsMerged$regionSet))
darSets <- darSets[lengths(darSets)>=10]
save(darSets, file=paste0("dars_",runName,"_darSets.RData"))

# Format a bit better the DARs table:
rownames(darsMerged) <- NULL
head(darsMerged)
colnames(darsMerged)[colnames(darsMerged)=="geneNearby"] <- "closestGene"

darsMerged$DAR_type <- "Up"
darsMerged$DAR_type[grep("___exclusive", darsMerged$regionSet)] <- "Up (exclusive)"
darsMerged$DAR_type[grep("___DW", darsMerged$regionSet)] <- "Down"

darsMerged <- darsMerged[,c("cellType", "DAR_contrast", "DAR_within","DAR_type", "region", "closestGene", "p_val", "p_val_adj", "avg_logFC", "pct.1", "pct.2","pctDiff", "regionSet")]

darsMerged$cellType <- factor(darsMerged$cellType)
darsMerged$DAR_contrast <- factor(darsMerged$DAR_contrast)
darsMerged$DAR_within <- factor(darsMerged$DAR_within)
darsMerged$DAR_type <- factor(darsMerged$DAR_type)
darsMerged$regionSet <- factor(darsMerged$regionSet)

darsMerged$p_val <- signif(darsMerged$p_val, 2)
darsMerged$p_val_adj <- signif(darsMerged$p_val_adj, 2)
darsMerged$avg_logFC <- signif(darsMerged$avg_logFC, 2)
darsMerged$pct.1 <- signif(darsMerged$pct.1, 2)
darsMerged$pct.2 <- signif(darsMerged$pct.2, 2)
darsMerged$pctDiff <- signif(darsMerged$pctDiff, 2)

darsMerged <- darsMerged[order(-darsMerged$avg_logFC),]
rownames(darsMerged) <- NULL

library(data.table)
darsMerged <- data.table(darsMerged)
saveRDS(darsMerged, file=paste0("dars_",runName,"_darsMerged-Formatted.Rds")) 
