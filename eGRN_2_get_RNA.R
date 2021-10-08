sce <- readRDS("~/jjans/data/Desplan/GSE142787_P15.rds")
rnaMeta <- read.table("~/jjans/analysis/development/development_atlas/Figures/Desplan/Neset_predictions_ownmeancentered_P15.tsv",sep="\t",header=T,row.names=1)
rnaMeta$annotation <- gsub("T4-5c/d","T5",rnaMeta$annotation)
rnaMeta$annotation <- gsub("T4-5a/b","T4",rnaMeta$annotation)
selected_cells_RNA <- row.names(rnaMeta)

cellTypesInGrs_RNA <- unique(rnaMeta$annotation)
exprMat <- as.matrix(GetAssayData(sce, slot = "counts",assay = 'RNA'))
exprMat <- exprMat[,selected_cells_RNA]
dim(exprMat)

exprMat <- t(t(exprMat)/colSums(exprMat))
cells <- selected_cells_RNA

## Keep only TFs 
load("resources/tfs_annotation_dmel_mcv89chip.RData")
allTfs <- tfs_annotation[which(tfs_annotation$tfAnnotationPriority <= 1),]; nrow(allTfs)
allTfs <- intersect(rownames(exprMat), rownames(allTfs)); length(allTfs)

# gplots::venn(list(rownames(rnaMeta), colnames(exprMat)))
exprMatSS <- exprMat[allTfs, cells]
dim(exprMatSS)
# 517 
save(exprMatSS, file=paste0("exprMatSSnorm_tfs_",runName,".RData")) 

## Average by cell type ----
meanExpr <- list()
percExpr <- list()
nCells <- c()
genesDetectedCnt <- list()
genesDetectedPerc <- list()

for(cn in cellTypesInGrs_RNA)
{
  print(cn)
  cn_rna <- cn
  
  cells_rna <- row.names(rnaMeta)[which(rnaMeta$annotation %in% cn_rna)];
  nCells[cn] <- length(cells_rna)

  meanExpr[[cn]] <- rowMeans(exprMatSS[,cells_rna,drop=F])
  percExpr[[cn]] <- rowMeans(exprMatSS[,cells_rna,drop=F]>0)
  genesDetectedPerc[[cn]] <- percExpr[[cn]]

}
meanExprMat <- do.call(cbind, meanExpr)
dim(meanExprMat)

percExprMat <- do.call(cbind, percExpr)
percExprMat <- percExprMat[,colnames(meanExprMat)]

meanExprMat <- meanExprMat[,which(apply(meanExprMat, 2, max, na.rm=T)!= -Inf)] # remove empty columns (S-Shape)

save(meanExprMat, file=paste0("exprMat_",runName,"_meanByCellType.RData")) 
save(percExprMat, file=paste0("percMat_",runName,"_meanByCellType.RData")) 
save(genesDetectedPerc, file=paste0("rna_",runName,"_genesDetected_perc.RData")) 
