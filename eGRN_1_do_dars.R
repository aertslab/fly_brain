do_dars<-function(predMatrix_i,cellData,regionData,runName,identName='cellType',maxCellsPerIdent=1000,nCores=10) {
################ loading libraries ###############################
library(Seurat)
packageVersion("Seurat")
################ Setup paralelization ###############################
library(future)
library(future.apply)
# check the current active plan
# plan()
# plan("multiprocess", workers = nCores)
# options(future.globals.maxSize = 180*1000*1024^2) # number of GB
plan("sequential")
# Paralellizes: NormalizeData, ScaleData, FindMarkers, FindClusters, FindIntegrationAnchors

############## Seurat #############################
message("Creating seurat object...")
seuratObject <- CreateSeuratObject(predMatrix_i, meta.data=cellData[colnames(predMatrix_i),], assay="ATAC")
rm(predMatrix_i)

message("Normalizing...")
seuratObject <- NormalizeData(seuratObject) # Default (&carmens): normalization.method="LogNormalize", scale.factor=1e4
seuratObject <- FindVariableFeatures(seuratObject, selection.method='mean.var.plot', mean.cutoff=c(0.0125, 3), dispersion.cutoff=c(0.5, Inf))
variableFeatures <- VariableFeatures(seuratObject)
save(variableFeatures, file=paste0("dars_",runName,"_variableFeatures.RData"))

# Remove read depth variation
# Scales and centers features in the dataset.
# If variables are provided in vars.to.regress,
#    they are individually regressed against each feautre,
#    and the resulting residuals are then scaled and centered.
seuratObject <- ScaleData(seuratObject,
                          features=rownames(seuratObject), # Vector of features names to scale/center. Default: all features
                          vars.to.regress=NULL) # e.g. batch correction

############# Markers: vs closest cluster in Tree ################
seuratObject <- SetIdent(seuratObject, value=identName)

#### Calculate tree
seuratObject <- RunPCA(seuratObject)
# Build a phylogenetic tree, and rename/reorder cluster names according to their position on the tree
# See help for details on tree building strategy
# This gives closely related clusters similar cluster IDs, which is occasionally useful for visualization later on
# Assigned cluster will be placed in the 'tree.ident' field of nbt@data.info, and also stored in nbt@ident
seuratObject <- BuildClusterTree(seuratObject, reorder=F, reorder.numeric=F, dims=1:20)
png(paste0("dars_",runName,"_clusterTree.png"),height=700, width=500, type="cairo")
Seurat::PlotClusterTree(seuratObject)
dev.off()

clTree <- seuratObject@tools$BuildClusterTree
save(clTree, file=paste0("dars_",runName,"_clusterTree.RData"))

contrast <- list()
for(nNode in rev(sort(unique(seuratObject@tools$BuildClusterTree$edge[,1]))))
{
  nodeChildren <- seuratObject@tools$BuildClusterTree$edge[which(seuratObject@tools$BuildClusterTree$edge[,1] %in% nNode),2]  
  
  c1 <- seuratObject@tools$BuildClusterTree$tip.label[nodeChildren[1]]
  if(is.na(c1)) c1 <- paste(contrast[[as.character(nodeChildren[1])]], collapse=" + ")
  
  c2 <- seuratObject@tools$BuildClusterTree$tip.label[nodeChildren[2]]
  if(is.na(c2)) c2 <- paste(contrast[[as.character(nodeChildren[2])]], collapse=" + ")
  
  contrast[[as.character(nNode)]] <- c(c1, c2)
}
save(clTree, file=paste0("dars_",runName,"_clusterTree_leafNames.RData"))


#### Calculate markers in tree
library(dplyr)
clusters <- names(contrast)
FCtr <- 0.20
maxCellsTargetCluster <- 5000
maxCellsRest <- NULL # total cells (e.g. randomly from all, ignore clusters)
maxCellsRestCls <- 300 # total: 200*ncusters; NULL to ignore
maxCellsRestCls * length(unique(seuratObject@active.ident))

message("Calculating markers in tree:")
dir.create("tmp")
darsList <- list()
for(cl in clusters[which(!clusters %in% names(darsList))])
{
  message(cl)
  # would probably be better with lower logFC... but it takes too long...
  # wilcox is also faster than MAST
  darsList[[cl]] <- tryCatch(
    FindMarkers(seuratObject, ident.1='clustertree', ident.2=cl, logfc.threshold=FCtr, only.pos=FALSE, max.cells.per.ident=maxCellsPerIdent)
    , error = function(e) { print(e); return(data.frame())})
  dars <- darsList[[cl]]
  saveRDS(dars, file=paste0("tmp/tmp_",runName,"_treeNode", make.names(cl),"_wilcox.Rds"))
  message("nDARs: ", nrow(dars))
  
  # Add gene & contrast name:
  if(nrow( darsList[[cl]])>0)
  {
    colNams <- colnames(darsList[[cl]])
    darsList[[cl]]$region <- rownames(darsList[[cl]])
    darsList[[cl]]$geneNearby <- regionData[rownames(darsList[[cl]]),"SYMBOL"]
    darsList[[cl]][which(is.na(darsList[[cl]]$geneNearby)),"geneNearby"] <- ""
    darsList[[cl]]$cluster <- paste0("treeNode_", rep(cl, nrow(darsList[[cl]])))
    darsList[[cl]] <- darsList[[cl]][,unique(c("cluster","region", colNams,"geneNearby"))]
    rownames(darsList[[cl]]) <- NULL 
  }
}
saveRDS(darsList, file=paste0("dars_",runName,"_treeContrasts.Rds"))

dars <- do.call(rbind, darsList)
dars$pctDiff <- dars$pct.1-dars$pct.2
# dars <- dars[order(-(dars$pct.1*dars$power)),]# or myAUC; plot(dars$power, dars$myAUC)
dars <- dars[order(-(dars$avg_logFC)),]# or myAUC; plot(dars$power, dars$myAUC)
rownames(dars) <- NULL
saveRDS(dars, file=paste0("dars_",runName,"_treeContrasts_wGenes.Rds"))

# Replace node by cluster name
# dars <- readRDS(paste0("dars_",runName,"_treeContrasts_wGenes.Rds"))
names(contrast) <- paste0("treeNode_", names(contrast))
contrast <- do.call(rbind,contrast)
contrast[unique(c(grep(" + ", contrast[,1], fixed=T, invert=T), grep(" + ", contrast[,2], fixed=T, invert=T))),]

darsS <- split(dars, dars$cluster)
darsS <- darsS[intersect(names(darsS), rownames(contrast))]
sapply(darsS, nrow)
darsList <- list()
for(cl in names(darsS))
{
  # print(cl)
  ## DARs up:
  if((!grepl(" + ", contrast[cl,1], fixed=T)) & (contrast[cl,1]!="-")) # just one cluster
  {
    clName <- paste0(contrast[cl,1], "__vsClosest")
    darsList[[clName]] <- darsS[[cl]][which(darsS[[cl]]$avg_logFC>0),]; darsS[[clName]]
    if(nrow(darsList[[clName]])>0) darsList[[clName]]$cluster <- clName
  }

  ## DARs DW --> Up in the other cluster
  if((!grepl(" + ", contrast[cl,2], fixed=T)) & (contrast[cl,2]!="-")) # just one cluster
  {
    clName <- paste0(contrast[cl,2], "__vsClosest")
    darsList[[clName]] <- darsS[[cl]][which(darsS[[cl]]$avg_logFC < 0),];darsList[[clName]]
    if(nrow(darsList[[clName]])>0){
      darsList[[clName]]$cluster <- clName

      ## invert:
      darsList[[clName]]$avg_logFC <- -darsList[[clName]]$avg_logFC
      darsList[[clName]]$pctDiff <- -darsList[[clName]]$pctDiff

      tmp <- darsList[[clName]]$pct.1
      darsList[[clName]]$pct.1 <- darsList[[clName]]$pct.2
      darsList[[clName]]$pct.2 <- tmp
    }
  }

  darsList[[cl]] <- NULL
}
dars <- do.call(rbind, darsList)
rownames(dars) <- NULL
dars$pctDiff <- dars$pct.1-dars$pct.2

dars <- dars[order(-(dars$avg_logFC)),]
saveRDS(dars, file=paste0("dars_",runName,"_treeContrasts_wGenes_renamedToTypes.Rds"))
cbind(table(dars$cluster))

######## Markers: each vs rest ##############
seuratObject <- SetIdent(seuratObject, value=identName)

# @Danie: These numbers of cells are to subsample the clusters 
# (to have them a bit more ballanced and run faster, they are not necesssarily the best settings, just what worked for me...)
FCtr <- 0.20
maxCellsTargetCluster <- 5000
maxCellsRest <- NULL # total cells (e.g. randomly from all, ignore clusters)
maxCellsRestCls <- 300 # total: 200*ncusters; NULL to ignore
maxCellsRestCls * length(unique(seuratObject@active.ident))

library(dplyr)
message(identName)
clusters <- as.character(unique(seuratObject@active.ident))
clusters <- clusters[which(!clusters %in% c("-","other","ignore", "unclear", "NotAssigned"))]
clusters <- clusters[which(!is.na(clusters))]

message("Calculating markers:")
dir.create("tmp")
darsList <- list()
for(cl in clusters[which(!clusters %in% names(darsList))])
{
  message(cl)
  
  seuratObject <- SetIdent(seuratObject, value=identName) # reset (we will modify some values)
  levels(seuratObject@active.ident) <- c(levels(seuratObject@active.ident), "other", "ignore")
  
  set.seed(123)
  # 1. Subset current group if far too big, and add the rest of cells to "ignore"
  cellsCl <- names(seuratObject@active.ident[which(seuratObject@active.ident==cl)])
  if(maxCellsTargetCluster < length(cellsCl)) cellsCl <- sample(cellsCl, maxCellsTargetCluster)
  
  # 2. Select the cells for the contrast trying to keep similar numbers per cluster
  cellsRest <- names(seuratObject@active.ident[which(!seuratObject@active.ident==cl)])
  if(!is.null(maxCellsRest)) cellsRest <- sample(cellsRest, maxCellsRest)
  if(!is.null(maxCellsRestCls)){
    seuratObject@meta.data$row.names <- rownames(seuratObject@meta.data)
    cellsRest <- seuratObject@meta.data[cellsRest,] %>% group_by(!!as.name(identName)) %>% sample_n(maxCellsRestCls, replace=T) %>% distinct %>% select(row.names)
    cellsRest <- cellsRest$row.names
  }
  
  seuratObject@active.ident[names(seuratObject@active.ident)] <- "ignore" #ignore all others
  seuratObject@active.ident[cellsRest] <- "other"
  seuratObject@active.ident[cellsCl] <- cl
  print(cbind(table(factor(seuratObject@active.ident))))
  
  # would probably be better with lower logFC... but it takes too long...
  # wilcox is also faster than MAST
  darsList[[cl]] <- tryCatch(FindMarkers(object=seuratObject, logfc.threshold=FCtr, only.pos=FALSE, ident.1=cl, ident.2="other", test.use="wilcox")
                             , error = function(e) { print(e); return(data.frame())})
  dars <- darsList[[cl]]
  saveRDS(dars, file=paste0("tmp/tmp_",runName,"_", make.names(cl),"_wilcox.Rds"))
  message("nDARs: ", nrow(dars))
  
  # Add gene & contrast name:
  colNams <- colnames(darsList[[cl]])
  darsList[[cl]]$region <- rownames(darsList[[cl]])
  darsList[[cl]]$geneNearby <- regionData[rownames(darsList[[cl]]),"SYMBOL"]
  darsList[[cl]][which(is.na(darsList[[cl]]$geneNearby)),"geneNearby"] <- ""
  darsList[[cl]]$cluster <- rep(cl, nrow(darsList[[cl]]))
  darsList[[cl]] <- darsList[[cl]][,unique(c("cluster","region", colNams,"geneNearby"))]
  rownames(darsList[[cl]]) <- NULL
}
saveRDS(darsList, file=paste0("dars_",runName,".Rds"))

dars <- do.call(rbind, darsList)
dars$pctDiff <- dars$pct.1-dars$pct.2
# dars <- dars[order(-(dars$pct.1*dars$power)),]# or myAUC; plot(dars$power, dars$myAUC)
dars <- dars[order(-(dars$avg_logFC)),]# or myAUC; plot(dars$power, dars$myAUC)
rownames(dars) <- NULL
saveRDS(dars, file=paste0("dars_",runName,"_wGenes.Rds"))
    

    

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

}


