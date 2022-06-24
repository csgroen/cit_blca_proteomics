################################################################################
## Classify
################################################################################
stroggilos_cls <- function(prot_symbol, sa) {
    nps_train <- .nps_train()
    nps_centroids <- nps_train$centroids
    
    #-- Recenter based on NMIBC
    ref_samps <- sa %>% filter(Stage_ == "NMIBC") %>% pull(id)
    ref_values <- rowMeans(prot_symbol[,ref_samps], na.rm = TRUE)
    reref_prot <- lapply(rownames(prot_symbol), function(gene) { prot_symbol[gene,] - ref_values[gene]}) %>%
        bind_rows() %>%
        as.matrix()
    rownames(reref_prot) <- rownames(prot_symbol)
    
    #-- Predict
    nps_pred <- .predict_NPScentroid(reref_prot, nps_centroids)
    pred_tb <- tibble(id = colnames(reref_prot), NPS = nps_pred)
    sa_nps <- sa %>%
        select(id, uPG = CCP_cluster, Stage_) %>%
        left_join(pred_tb) %>%
        mutate(NPS_MIBC = ifelse(Stage_ == "MIBC", "MIBC", NPS))
    table(sa_nps$uPG, sa_nps$NPS)
    # NPS1 NPS2 NPS3
    # A   11    1    1
    # B   10    1    0
    # C    4    2    1
    # D    4   12    0
    # E    2    1   13
    
    return(list(sa_nps = sa_nps, nps_difs = nps_train$nps_diffs))
}
.nps_train <- function() {
    #-- Read protein expression data
    str_prot_exp <- 
        read_excel("data/other_datasets/Stroggilos_IJC/ijc32556-sup-0001-appendixs1.xlsx", sheet = 2, skip = 6) %>%
        select(protein_id = Accession, symbol = `Gene name`, starts_with("Sample"))
    
    str_prot_counts <- str_prot_exp %>%
        select(symbol, starts_with("Sample")) %>%
        mutate(mean_counts = rowMeans(select(., starts_with("Sample")))) %>%
        group_by(symbol) %>%
        slice_max(mean_counts) %>%
        ungroup() %>%
        select(-mean_counts) %>%
        as.data.frame() %>%
        column_to_rownames("symbol") %>%
        as.matrix()
    #---- Scale for training
    str_prot_scaled <- t(scale(t(log2(str_prot_counts+1))))
    
    #-- Read annotation
    str_sa <- read_excel("data/other_datasets/Stroggilos_IJC/ijc32556-sup-0001-tables.xlsx", 
                         sheet = 2, skip = 2) %>%
        select(id = PatientID, NPS = `Consensus Subtype`)
    cls <- str_sa %>% pull(NPS, id)
    
    #-- Top cls markers
    cls_diffs <- findClassDiffs(str_prot_scaled[,names(cls)], cl = cls, pval_cutoff = 1)
    mks <- list(NPS1 = cls_diffs %>% arrange(desc(auc_NPS1)) %>% pull(feat) %>% .[1:10],
                NPS2 = cls_diffs %>% arrange(desc(auc_NPS2)) %>% pull(feat) %>% .[1:10],
                NPS3 = cls_diffs %>% arrange(desc(auc_NPS3)) %>% pull(feat) %>% .[1:10])
    all_mks <- unlist(mks) %>% as.character()
    write.xlsx(cls_diffs, file = "results/tables/Stroggilos_NPS_class_markers.xlsx")
    
    #-- Centroids
    centroids <- sapply(split(names(cls), cls), function(cl) {
        rowMeans(str_prot_scaled[all_mks,cl])
    })
    
    #-- Predict and metrics
    pred <- .predict_NPScentroid(str_prot_scaled[all_mks,names(cls)], centroids)
    
    ncc_cls <- tibble(cls, pred) %>%
        mutate(cls = factor(cls),
               pred = factor(pred))
    # yardstick::conf_mat(ncc_cls, cls, pred)
    # # Truth
    # # Prediction NPS1 NPS2 NPS3
    # # NPS1   15    1    0
    # # NPS2    0   37    0
    # # NPS3    2    4   39
    # yardstick::bal_accuracy(ncc_cls, cls, pred) # 0.942
    # yardstick::f_meas(ncc_cls, cls, pred) # 0.925
    
    return(list(centroids = centroids, nps_diffs = cls_diffs))
}

.predict_NPScentroid <- function(data, centroids) {
    feats <- rownames(centroids)
    fkeep <- feats[feats %in% rownames(data)]
    
    centroid_cors <- cor(centroids[fkeep,], data[fkeep,], use = "pair") %>% t()
    pred <- paste0("NPS", apply(centroid_cors, 1, which.max))
    
    return(pred)
}

################################################################################
## Plot
################################################################################
dir.create("results/suppfig_stroggilos/", showWarnings = FALSE)
plot_stroggilosComp <- function(prot_data, prot_symbol, strog_cls, sa) {
    #-- From ref: 10.1002/ijc.32556
    # stroggilos_prots <- 
    #     list(`Cell cycle` = c("NASP", "RCC2", "CDC37", "YWHAG", "PAICS", "NME2", "GART", 
    #                           "CDC42", "BUB3"),
    #          `Inflammation` = c("STAT1", "STAT3", "SND1", "DHX9", "HMGB1", "HMGB2", "HMGB3", 
    #                             "PTGES3", "RNF213", "TYMP"),
    #          `DNA damage\nresponse` = c("RUVBL2", "PCNA", "PRKDC", "PARP1", "TOP2B", "APEX1"),
    #          `Stromal/EMT` = c("ACTC1", "CNN1", "MFAP4", "ACTA2", "DES", "MYH11", "MYL9",
    #                        "TAGLN", "COL6A3", "COL14A1", 
    #                        "VIM", "COL1A1", "COL1A2", "TFGBI", "CAV1", "NID1", "POSTN", 
    #                    "FLNA"),
    #          `Differentiation` = c("UPK2", "UPK3BL1", "UPK1B", "GPX2", "PDCD4", "SRC", "ADIRF", 
    #                                "FBP1", "FABP4")
    #     )
    
    nps_cls <- c("NPS1", "NPS2", "NPS3")
    stroggilos_prots <- lapply(nps_cls, function(nps) {
        log_var <- paste0("logFC_", nps)
        strog_cls$nps_difs %>%
            slice_max(order_by = !!sym(log_var), n = 30) %>%
            pull(feat)
    })
    names(stroggilos_prots) <- paste0(nps_cls, " markers")
    strog_cls_df <- strog_cls$sa_nps
    
    #-- Get shared proteins
    pts4plot <- lapply(stroggilos_prots, function(pts) {
        pts[pts %in% prot_data$wpAnnot$symbol]
    })
    
    data4plot <- prot_symbol[unlist(pts4plot),] %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("id") %>%
        left_join(select(sa, id, Stage. = StageTa, Subtype, CCP_cluster))
    
    #-- Plot heatmap
    gghm <- data4plot %>%
        left_join(select(strog_cls_df, id, NPS, NPS_MIBC)) %>%
        group_by(CCP_cluster) %>%
        dplyr::rename(uPG = CCP_cluster, Stage = `Stage.`) %>%
        ggheatmap(rowv = pts4plot,
                  hm_colors = viridis::inferno(100),
                  cluster_rows = FALSE,
                  show_dend_col = FALSE,
                  hm_color_values = scales::rescale(c(-5,-1.5,-1,-0.5,-0.2,0,0.2,0.5,1,1.5,4)),
                  show_colnames = FALSE, 
                  colors_title = bquote("Scaled"~log[2](ratio)),
                  dist_method = "spearman",
                  clustering_method = "ward.D2",
                  colv = "id", 
                  fontsize = 8,
                  group_prop = 0.03,
                  dend_prop_col = 0.05,
                  group_lines = TRUE,
                  group_line_color = "white",
                  group_lwd = 0.7,
                  raster = TRUE)
    
    gghm <- gghm %>%
        add_tracks(track_columns = c("Stage", "Subtype", "NPS"),
                   track_colors = list('Stage' = c("Ta"="#fee0d2", "T1"="#fc9272", "MIBC"="#de2d26"), 
                                       'Subtype' = subtype_cols,
                                       'NPS' = c("MIBC" = "grey40",
                                                 "NPS1" = "#acd785",
                                                 "NPS2" = "#a0c7db",
                                                 "NPS3" = "#1e74ae")
                                       # 'NPS_MIBC' = c("MIBC" = "grey30",
                                       #                "NPS1" = "#acd785",
                                       #                "NPS2" = "#a0c7db",
                                       #                "NPS3" = "#1e74ae")
                   ),
                   track_pos = 'top', track_prop = 0.08, fontsize = 8)
    gghm
    ggsave("results/suppfig_stroggilos/fig.pdf", 
           plot = gghm, height = 7, width = 7)
    
    return(list(gghm = gghm, nps_markers = pts4plot))
}


