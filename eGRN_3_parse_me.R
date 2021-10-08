parse_me <- function(runName,rnaMarkers){
    #### Load and merge: ----
    fileNames <- c(
      "chip_noBg @ auc01"= paste0(runName,"_chipEnrichment_aucMax01.Rds"),
      "chip_noBg @ auc03"= paste0(runName,"_chipEnrichment_aucMax03.Rds"),
      "motifs_noBg @ auc001"= paste0(runName,"_motifEnrichment_noBg_aucMax001.Rds"),
      "motifs_noBg @ auc005"= paste0(runName,"_motifEnrichment_noBg_aucMax005.Rds"),
      "motifs_noBg @ auc01"= paste0(runName,"_motifEnrichment_noBg_aucMax01.Rds"),
      "motifs_noBg @ auc05"= paste0(runName,"_motifEnrichment_noBg_aucMax05.Rds"),
      "motifs_withBg @ auc005"= paste0(runName,"_motifEnrichment_withBackground_aucMax005.Rds"),
      "motifs_withBg @ auc01"= paste0(runName,"_motifEnrichment_withBackground_aucMax01.Rds"),
      "motifs_withBg @ auc05"= paste0(runName,"_motifEnrichment_withBackground_aucMax05.Rds"),
      "motifs_withBg @ auc10"= paste0(runName,"_motifEnrichment_withBackground_aucMax10.Rds")
    )

    enrichmentTable <- NULL
    for(i in seq_along(fileNames))
    {
      fn <- fileNames[i]
      if(file.exists(fn))
      {
        met <- readRDS(paste0(fn))$motifEnrichmentTable
        met$analysis <- names(fileNames)[i]
        colnames(met)[2] <- "motif"
        met <- met[,c("geneSet", "motif", "NES", "TF_annot", "analysis")]
        enrichmentTable <- rbind(enrichmentTable, met)
      }else{
        message(fn,  " does not exist.")
      }
    }
    # setwd("..")

    colnames(enrichmentTable)[1] <- "regionSet"
    colnames(enrichmentTable)[which(colnames(enrichmentTable)=="analysis")] <- "MotifEnrSettings"
    enrichmentTable <- enrichmentTable[order(-enrichmentTable$NES),]
    nrow(enrichmentTable)
    # 3,628,358
    enrichmentTable[1:20,]

    ### Settings
    tmp <- data.frame(do.call(rbind,strsplit(enrichmentTable$MotifEnrSettings, "_| @ ")))
    colnames(tmp) <- c("me_DB", "me_BG", "me_rocThr")
    enrichmentTable <- data.frame(enrichmentTable, tmp)

    enrichmentTable$me_DB <- as.character(enrichmentTable$me_DB)
    enrichmentTable[which(enrichmentTable$me_DB=="chip"),"me_DB"] <- "ChIP"
    enrichmentTable$me_DB <- factor(enrichmentTable$me_DB)

    #### DARset type
    dars <- readRDS(paste0("dars_",runName,"_darsMerged-Formatted.Rds"))
    head(dars)

    dars <- unique(dars[,c("regionSet","cellType", "DAR_contrast", "DAR_within", "DAR_type")]); nrow(dars)

    # DAR set size: 
    load(paste0("dars_",runName,"_darSets.RData")); length(darSets)
    setLen <- lengths(darSets)
    enrichmentTable$nRegions <- setLen[as.character(enrichmentTable$regionSet)]
    head(enrichmentTable)

    library(dplyr)
    enrichmentTable <- enrichmentTable %>%
      left_join(dars)
    head(enrichmentTable)

    enrichmentTable <- data.frame(enrichmentTable)

    #####
    # Annotation to TFs   ----- 

    # By default uses any annotation:
    colnames(enrichmentTable)[which(colnames(enrichmentTable)=="TF_annot")] <- "TF_annot_any"
    TFs <- setNames(strsplit(enrichmentTable$TF_annot_any, "; "), rownames(enrichmentTable))
    enrichmentTable$TFsList_any <- TFs

    ########
    # Split high-low confidence
    library(RcisTarget)
    #motifAnnotations <- readRDS("/staging/leuven/stg_00002/lcb/cbravo/feather/motifAnnot.dmel.mc9nr.withAertsDL.Rds")
    data(motifAnnotations_dmel)
    motifAnnotations <- motifAnnotations_dmel
    # Remove annotations to non-TFs from DB: 
    tfFileName="resources/tfs_annotation_dmel_mcv89chip.RData"
    load(tfFileName)
    tfs <- tfs_annotation[which(tfs_annotation$tfAnnotationPriority <= 1),]; nrow(tfs)
    motifAnnotations <- motifAnnotations[which(as.character(motifAnnotations$TF) %in% rownames(tfs)),]; nrow(motifAnnotations)

    # Add annotations:
    enrichmentTable$TF_annot_1direct <- getMotifAnnotation(as.character(enrichmentTable$motif), motifAnnotations, 
                                                           annotCats = c("directAnnotation"),
                                                           idColumn = "motif", returnFormat="asCharacter",keepAnnotationCategory = FALSE)[as.character(enrichmentTable$motif)]
    enrichmentTable$TF_annot_2ort <- getMotifAnnotation(as.character(enrichmentTable$motif), motifAnnotations, 
                                                        annotCats = c("inferredBy_Orthology","inferredBy_MotifSimilarity_n_Orthology"),
                                                        idColumn = "motif", returnFormat="asCharacter",keepAnnotationCategory = FALSE)[as.character(enrichmentTable$motif)]
    enrichmentTable$TF_annot_3simil <- getMotifAnnotation(as.character(enrichmentTable$motif),motifAnnotations, 
                                                          annotCats = c("inferredBy_MotifSimilarity"),
                                                          idColumn = "motif", returnFormat="asCharacter",keepAnnotationCategory = FALSE)[as.character(enrichmentTable$motif)]

    # use only high conf: 
    enrichmentTable$TFsList_highConf <- getMotifAnnotation(as.character(enrichmentTable$motif), motifAnnotations, 
                                                           annotCats = c("directAnnotation","inferredBy_Orthology","inferredBy_MotifSimilarity_n_Orthology"),
                                                           idColumn = "motif", returnFormat="asCharacter",keepAnnotationCategory = FALSE)[as.character(enrichmentTable$motif)]
    TFs <- setNames(strsplit(enrichmentTable$TFsList_highConf, "; "), rownames(enrichmentTable))
    enrichmentTable$TFsList_highConf <- TFs


    ## Format & save
    enrichmentTable$motif <- factor(enrichmentTable$motif)
    enrichmentTable$cellType <- factor(enrichmentTable$cellType)
    enrichmentTable$DAR_type <- factor(enrichmentTable$DAR_type)
    enrichmentTable <- enrichmentTable[order(-enrichmentTable$NES),]

    enrichmentTable <- RcisTarget::addLogo(enrichmentTable, dbVersion="v9") # v9dl doesnt seem to have enough motifs...
    enrichmentTable[which(enrichmentTable$me_DB=="ChIP"),"logo"] <- ""

    colOrd <- c("cellType", "logo", "NES", "TF_annot_any", 
                "me_DB", "me_BG", "me_rocThr", 
                "DAR_contrast", "DAR_within", "DAR_type", "nRegions", 
                "TF_annot_1direct", "TF_annot_2ort", "TF_annot_3simil", "TFsList_any", "TFsList_highConf",
                "regionSet", 
                "motif", 
                "MotifEnrSettings")

    enrichmentTable <- enrichmentTable[,colOrd, with=F]
    saveRDS(enrichmentTable, file=paste0("motifs_",runName,"_enrichmentTable_1.Rds"))
    
    
    rnaMarkers_subset <- rnaMarkers
    rnaMarkers_subset2 <- rnaMarkers[which(rnaMarkers$cluster %in% c("T4/T5")),]
    rnaMarkers_subset2$cluster <- gsub("T4/T5","T5",rnaMarkers_subset2$cluster)
    rnaMarkers_subset$cluster <- gsub("T4/T5","T4",rnaMarkers_subset$cluster)
    rnaMarkers_subset <- rbind(rnaMarkers_subset,rnaMarkers_subset2)
    saveRDS(rnaMarkers_subset,paste0("RNAmarkers_",runName,".Rds"))

    # load("int/02_cellTypesInGrs_RNA.RData")
    ##### Is the TF expressed in these cells? ----
    load(paste0("rna_",runName,"_genesDetected_perc.RData"))
    genesDetectedPerc <- lapply(genesDetectedPerc, function(x) names(x)[x>=0.10]) 

    ## TF is Marker gene?
    rnaMarkers <- readRDS(paste0("RNAmarkers_",runName,".Rds"))
    rnaMarkers <- rnaMarkers[which(rnaMarkers$avg_logFC>0),]
    rnaMarkers <- rnaMarkers[which(rnaMarkers$pct.ratio > 1.5),]

    ##
    met <- readRDS(paste0("motifs_",runName,"_enrichmentTable_1.Rds"))
    met <- as.data.frame(met)
    # gsWithMatch <- sort(intersect(as.character(unique(met$cellType)), unique(c(names(genesDetectedPerc), rnaMarkers$cluster)))); length(gsWithMatch)
    ctNames <-  as.character(unique(met$cellType))

    met$TFany_expressed <- ""
    met$TFhighConf_expressed <- ""
    met$TFany_marker <- ""
    met$TFhighConf_marker <- ""
    for(ct in ctNames)
    {
      print(ct)
      idx <- which(met$cellType == ct)

      genesDetected <- unique(unname(unlist(genesDetectedPerc[[ct]])))
      met[idx,"TFany_expressed"] <- sapply(met[idx,"TFsList_any"], function(x) paste(unique(x[which(x %in% genesDetected)]), collapse="; "))
      met[idx,"TFhighConf_expressed"] <- sapply(met[idx,"TFsList_highConf"], function(x) paste(unique(x[which(x %in% genesDetected)]), collapse="; "))
      rm(genesDetected)

      rnaCls <- ct
      if(length(rnaCls)>0)
      {
        markerGenes <- unique(as.character(rnaMarkers[which(rnaMarkers$cluster %in% rnaCls),"gene"]))
        met[idx,"TFany_marker"] <- sapply(met[idx,"TFsList_any"], function(x) paste(unique(x[which(x %in% markerGenes)]), collapse="; "))
        met[idx,"TFhighConf_marker"] <- sapply(met[idx,"TFsList_highConf"], function(x) paste(unique(x[which(x %in% markerGenes)]), collapse="; "))
        rm(markerGenes)
      }else{
        print(ct)
      }
    }

    met$TFany_isExpressed <- factor(as.character(met$TFany_expressed!=""))
    met$TFany_isMarker <- factor(as.character(met$TFhighConf_expressed!=""))

    met$TFhighConf_isExpressed <- factor(as.character(met$TFany_expressed!=""))
    met$TFhighConf_isMarker <- factor(as.character(met$TFhighConf_marker!=""))

    colOrd <- c("cellType", "logo", "NES", "TF_annot_any", 
                "me_DB", "me_BG", "me_rocThr", 
                "DAR_contrast", "DAR_within", "DAR_type", "nRegions", 
                "TFany_isExpressed", "TFany_isMarker", "TFany_expressed", "TFany_marker",
                "TFhighConf_isExpressed", "TFhighConf_isMarker", "TFhighConf_expressed", "TFhighConf_marker", 
                "TF_annot_1direct", "TF_annot_2ort", "TF_annot_3simil", "TFsList_any", "TFsList_highConf",
                "regionSet", 
                "motif", 
                "MotifEnrSettings")

    met <- data.frame(met)[,colOrd]
    met <- met[order(met$NES, decreasing=T),]
    saveRDS(met, file=paste0("motifs_",runName,"_enrichmentTable_2_wTFs.Rds"))
}