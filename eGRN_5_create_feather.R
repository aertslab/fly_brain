create_feather <- function(runName){
    files <- list.files(pattern = "signifRegionsList")
    files <- files[grepl(x = files, pattern = runName)]
    files <- files[grepl(x = files, pattern = '.df.RData')]
    signifRegionsList <- list()

    for(file in files){
        dbBg <- gsub(runName,"",file)
        dbBg <- gsub("_signifRegionsList_","",dbBg)
        dbBg <- gsub(".df.RData","",dbBg)
        load(file)
        signifRegionsList.df <- signifRegionsList.df[which(sapply(signifRegionsList.df, class) == "data.frame")]; length(signifRegionsList.df) # remove empty AUCs
        signifRegionsList[[dbBg]] <- signifRegionsList.df
        rm(signifRegionsList.df)

    }
    library(data.table)
    for(ictat in names(signifRegionsList))
    {
      print(ictat)
      for(aucTrh in names(signifRegionsList[[ictat]]))
      {
        print(aucTrh)
        colnames(signifRegionsList[[ictat]][[aucTrh]]) <- c("region", "idEnriched", "regionSet")
        signifRegionsList[[ictat]][[aucTrh]]$me_type <- ictat
        signifRegionsList[[ictat]][[aucTrh]]$rocTrh <- aucTrh    
      }
      signifRegionsList[[ictat]] <- rbindlist(signifRegionsList[[ictat]])
    }
    
    signifRegions <- rbindlist(signifRegionsList)
    rm(signifRegionsList)

    
    signifRegions <- signifRegions[, idEnriched:=as.factor(idEnriched)]
    signifRegions <- signifRegions[, regionSet:=as.factor(regionSet)]
    signifRegions <- signifRegions[, me_type:=as.factor(me_type)]
    signifRegions <- signifRegions[, rocTrh:=as.factor(rocTrh)]
    arrow::write_feather(signifRegions, paste0(runName,"_signifRegions.df.feather"))
    

    library(tibble)
    signifRegions <- arrow::read_feather(file=paste0(runName,"_signifRegions.df.feather"), mmap = TRUE)

    darSetNames <- readRDS(paste0("dars_",runName,"_darsMerged-Formatted.Rds")
    darSetNames <- unique(data.frame(unique(darSetNames[,c("cellType", "DAR_type","regionSet")])))
    darSetNames <- as_tibble(darSetNames)

    signifRegions <- signifRegions %>%
      # filter(rocTrh %in% c("auc01", "auc03", "auc05")) %>% # "auc001", "auc005", "auc10" 
      left_join(darSetNames, by = "regionSet") %>%
      select (c("cellType", "idEnriched", "region", "DAR_type")) %>%  # cell type instead of regionSet (ignores settings & direction...)
      distinct()

    signifRegions <- signifRegions %>% 
      mutate(DAR_type=ifelse(DAR_type %in% "Up (exclusive)", "Up", as.character(DAR_type))) %>% 
      distinct()
                           
    arrow::write_feather(signifRegions, paste0(runName,"_signifRegions.df.simplified.feather")) 
    
    signifRegions <- arrow::read_feather(file=paste0(runName,"_signifRegions.df.feather"), mmap = TRUE)
    ## Filters ----
    # Motif enrichment settings
    library(dplyr)
    signifRegions <- signifRegions %>%
      filter(rocTrh %in% c("auc03", "auc01", "auc05")) %>%      # discards: auc005 auc001 auc10; (auc 03 is only for ChIP?)
      filter(me_type %in% c("noBg", "withBg"))  # discards "ChIP"
    signifRegions <- signifRegions %>%       
      select("region", "idEnriched", "regionSet", "me_type") 
    load("resources/tfs_annotation_dmel_mcv89chip.RData")
    tfs <- tfs_annotation[which(tfs_annotation$tfAnnotationPriority <= 1),]; nrow(tfs)

    # High-confidence TF annotation 
    tfAnnot <- RcisTarget::importAnnotations("/staging/leuven/res_00001/databases/cistarget/motif2tf/motifs-v9-nr.flybase-m0.001-o0.0.tbl"); nrow(tfAnnot)
    tfAnnot <- tfAnnot[apply(tfAnnot[,c("directAnnotation", "inferred_Orthology")], 1, any),]; nrow(tfAnnot) 
    tfAnnot <- tfAnnot[,c("motif", "TF")]
    length(unique(tfAnnot$TF))
    tfAnnot <- tfAnnot[which(as.character(tfAnnot$TF) %in% rownames(tfs)),]; nrow(tfAnnot) 
    length(unique(tfAnnot$TF))

    signifRegions <- signifRegions %>%
      filter(idEnriched %in% tfAnnot$motif) # %>%  distinct() 

    # DARs
    dars <- readRDS("01_darsMerged-Formatted.Rds")
    dars <- unique(dars[,c("regionSet","cellType", "DAR_contrast", "DAR_within", "DAR_type")]); nrow(dars)
    dars <- dars[!dars$DAR_type %in% "Down"]; nrow(dars)
    head(dars)

    signifRegions <- signifRegions %>%
      filter(regionSet %in% dars$regionSet)
    nrow(signifRegions)

    signifRegions <- signifRegions %>%
      left_join(dars)  %>%     
      select("region", "cellType", "idEnriched")     # switch regionSet by cellType, and discards "me_type"
    signifRegions <- signifRegions %>% distinct()
    nrow(signifRegions)
    # 8,744,396 

    # Add annottation
    signifRegions <- signifRegions %>%
      inner_join(tfAnnot, by=c("idEnriched" = "motif")) # inner= keeps only common rows
    nrow(signifRegions)
    # 35,762,191
    signifRegions <- signifRegions %>%
      select(c("region", "TF", "cellType")) %>%
      distinct() 
    nrow(signifRegions)
    # 4,565,288
    length(unique(signifRegions$TF))
    arrow::write_feather(signifRegions, "06_cistromes.feather")

    }