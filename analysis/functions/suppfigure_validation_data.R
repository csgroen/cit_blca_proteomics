################################################################################
## Stroggilos data: 10.1002/ijc.32556
################################################################################
dir.create("results/suppfig_validation/", showWarnings = FALSE)

get_upg_centroids <- function(prot_symbol, samp_annot, diffProt_uPG, n_feats = 20) {
  #-- Select top markers
  top_mks <- .get_top_markers(diffProt_uPG, n = n_feats, pos_fc = FALSE) %>% unlist()
  
  #-- Get centroids
  cls <- samp_annot %>% pull(uPG, sample) %>% split(names(.), .)
  centroids <- sapply(cls, function(cl) {
    rowMeans(prot_symbol[top_mks,cl], na.rm = TRUE)
  })
  return(centroids)
}

.pp_strog_prot <- function() {
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
  str_prot_log <- log2(str_prot_counts+1)
  return(str_prot_log)
}
stroggilos_cls <- function(prot_symbol, sa) {
  nps_train <- .nps_train()
  nps_centroids <- nps_train$centroids
  nps_centroids %>%
    as.data.frame() %>%
    rownames_to_column("symbol") %>%
    write_csv(file = "results/tables/Stroggilos_NPS_centroids.csv")
  
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
  str_prot <- .pp_strog_prot()
  #---- Scale for training
  str_prot_scaled <- t(scale(t(str_prot)))

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
  ggsave("results/suppfig_validation/stroggilos.pdf", 
         plot = gghm, height = 7, width = 7)
  
  return(list(gghm = gghm, nps_markers = pts4plot))
}




normalize_to_reference <- function(ref_data, new_data) {
  #-- Keep feats
  ref_feats <- rownames(ref_data)
  new_feats <- rownames(new_data)
  keep_feats <- new_feats[new_feats %in% ref_feats]
  
  message(length(keep_feats), " out of ", length(new_feats), " features kept...")
  
  #-- Get gene means and sds
  gene_means <- rowMeans(ref_data[keep_feats,], na.rm = TRUE)
  gene_sds <- rowSds(ref_data[keep_feats,], na.rm = TRUE)
  new_data[keep_feats,]
  
  new_scaled <- t(scale(t(new_data[keep_feats,]), center = TRUE, scale = TRUE))
  
  new_res <- sapply(1:length(gene_means), function(i) {
    (new_scaled[i,] * gene_sds[i]) + gene_means[i]
  }) %>% t()
  rownames(new_res) <- keep_feats
  colnames(new_res) <- colnames(new_data)
  
  return(new_res)
}

################################################################################
## Xu et al data. doi: https://doi.org/10.1186/s13045-022-01291-7
################################################################################
# Mon Aug  8 14:35:09 2022 ------------------------------
.pp_xu_prot <- function() {
  xu_prot_data <- read_excel("data/other_datasets/Xu_JHO/13045_2022_1291_MOESM13_ESM.xlsx", sheet = 8, skip = 2)
  #-- Clean xu_prot_data
  xu_prot_data1 <- xu_prot_data %>%
    mutate(pid = paste0("pid_", 1:n())) %>%
    relocate(pid)
  
  prot_annot <- select(xu_prot_data1, pid, symbol = Symbol)
  xu_prot_data1 <- xu_prot_data1 %>%
    select(-Symbol) %>%
    as.data.frame() %>%
    column_to_rownames("pid")
  
  xu_prot_data <- gexpPreprocess(xu_prot_data1, prot_annot, og_annot = "pid", keep_annot = "symbol")$gexp
}

.pp_xu_sa <- function() {
  xu_samp_annot <- read_excel("data/other_datasets/Xu_JHO/13045_2022_1291_MOESM13_ESM.xlsx", sheet = 2) %>%
    rename_with(.fn = ~ str_replace_all(., " ", "_") %>% str_to_lower())
  xu_mut_annot <- read_excel("data/other_datasets/Xu_JHO/13045_2022_1291_MOESM13_ESM.xlsx", sheet = 3)
  
  # Get FGFR3 annotation
  fgfr3_mut <- xu_mut_annot %>% filter(Hugo_Symbol == "FGFR3") %>% unlist() %>% .[-1]
  fgfr3_tb <- tibble(case_id = names(fgfr3_mut),
                     FGFR3_mutation = fgfr3_mut) %>%
    mutate(FGFR3_mutation = ifelse(!is.na(FGFR3_mutation), "M", "WT")) %>%
    filter(str_detect(case_id, "_T$")) %>%
    mutate(case_id = str_remove(case_id, "_T$"))
  
  xu_samp_annot <- xu_samp_annot %>%
    left_join(fgfr3_tb)
  return(xu_samp_annot)
}

classify_xu_data <- function(upg_centroids, min_cor = 0.2) {
  xu_prot_data <- .pp_xu_prot()

  feats <- rownames(upg_centroids)
  kept_feats <- feats[feats %in% rownames(xu_prot_data)]
  
  #-- Classify
  xu_scaled <- t(scale(t(xu_prot_data)))
  centroid_cors <- cor(xu_scaled[kept_feats,], upg_centroids[kept_feats,], use = "pair")
  
  colnames(centroid_cors) <- paste0("cor_", LETTERS[1:5])
  class_res <- centroid_cors %>% as.data.frame() %>% rownames_to_column("proteome_id")
  
  new_class <- LETTERS[apply(centroid_cors, 1, which.max)]
  max_cors <- apply(centroid_cors, 1, max)
  new_class <- ifelse(max_cors > min_cor, new_class, "unc")
  
  upg_xu <- class_res %>%
    mutate(uPG = new_class,
           max_cor = max_cors) %>%
    relocate(uPG, max_cor, .after = proteome_id) %>%
    tibble()
  
  return(upg_xu)
  
  
}

plots_xu_data <- function(upg_xu, diffProt_uPG) {
  xu_prot_data <- .pp_xu_prot()
  xu_samp_annot <- .pp_xu_sa()
  
  #-- Select top 20 uPG markers
  top_aucs_pos <- .get_top_markers(diffProt_uPG, n = 15, pos_fc = TRUE)
  top_aucs <- top_aucs_pos %>% unlist() %>% unique()
  prots_used <- top_aucs[top_aucs %in% rownames(xu_prot_data)]
  idx <- rowSums(!is.na(xu_prot_data[prots_used,])) > 100
  prots_used <- prots_used[idx]
  
  top_aucs_list <- lapply(top_aucs_pos, function(aucs) {
    aucs[aucs %in% prots_used]
  })
  names(top_aucs_list) <- LETTERS[1:5]
  
  #-- Get annotation
  final_sa <- xu_samp_annot %>%
    mutate(
      months_to_death = as.numeric(months_until_death),
      months_to_last_followup = as.numeric(months_to_last_followup),
      os = ifelse(vital_status == "Alive", 0, 1),
      os_time = case_when(
        !is.na(months_to_death) ~ months_to_death,
        TRUE ~ months_to_last_followup)) %>%
    select(proteome_id, tnm_stage, histologic_grade, proteome_cluster, FGFR3_mutation,
           pfs_status, pfs_months, os, os_time) %>%
    mutate(stage = str_extract(tnm_stage, "T."),
           mibc = case_when(
             stage %in% c("T2", "T3", "T4") ~ "MIBC",
             TRUE ~ stage
           )) %>%
    left_join(upg_xu, by = "proteome_id") %>%
    mutate(pfs_status = ifelse(str_detect(pfs_status, "^0"), 0, 1))
  
  #-- Plot heatmap + class flow
  hm <- .plt_xu_hm(xu_prot_data, final_sa, top_aucs_list)
  plt_cls <- .plt_xu_classFlow(final_sa)
  # plt_surv <- .plt_xu_surv(final_sa)
  # plt_surv2 <- ((plt_surv$plot + labs(y="Progression-free survival")) / 
  #                 (plt_surv$table + labs(y = ""))) + plot_layout(heights = c(3,1))
  
  #-- Stats + survival
  test_res <- final_sa %>%
    select(proteome_id, uPG, FGFR3_mutation, mibc, proteome_cluster, stage) %>%
    doAssociationTests(id_var = "proteome_id", test_var = "uPG")
  
  plt <- (hm | ((plt_cls + guides(fill = "none")) / plot_spacer()) + 
      plot_layout(heights = c(1,1.5))) + plot_layout(widths = c(1.5,1))
  
  ggsave("results/suppfig_validation/xu.pdf", plt, width = 8, height = 6.5)
  
  return(plt)
}

.plt_xu_surv <- function(final_sa) {
  surv_theme <- theme_csg_scatter +
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(),
          axis.line = element_line())
  
  os_plt <- ggsurvplot(survfit(Surv(os_time, os) ~ uPG, 
                               data = final_sa %>% filter(uPG != "unc")), 
                       risk.table = TRUE, pval = TRUE, tables.height = 0.4, 
                       xlim = c(0,60), break.time.by = 12, 
                       ggtheme = surv_theme,
                       fontsize = 3, pval.size = 3, risk.table.title = "")
  
  pfs_plt <- ggsurvplot(survfit(Surv(pfs_months, pfs_status) ~ uPG, 
                               data = final_sa %>% filter(uPG != "unc")), 
                       risk.table = TRUE, pval = TRUE, tables.height = 0.4, 
                       xlim = c(0,60), break.time.by = 12, 
                       ggtheme = surv_theme,
                       fontsize = 3, pval.size = 3, risk.table.title = "")
  
  
  (os_plt$plot + labs(subtitle = "OS")) / (pfs_plt$plot + labs(subtitle = "PFS"))
  return(pfs_plt)
}

.plt_xu_hm <- function(xu_prot_data, sannot, prot_feats) {
  #-- Subset
  prots_used <- unlist(prot_feats)
  xu_prot_data_upg <- xu_prot_data[prots_used,] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("proteome_id")
  
  #-- Join annotation
  xu_prot_data_upg <- sannot %>%
    left_join(xu_prot_data_upg, by = "proteome_id")
  
  #-- Plot
  hm <- 
    xu_prot_data_upg %>%
    rename(`Xu proteome cluster` = proteome_cluster,
           `FGFR3 mutation` = FGFR3_mutation,
           `Stage` = mibc) %>%
    group_by(uPG) %>%
    ggheatmap(colv = "proteome_id",
              rowv = prot_feats,
              raster = TRUE,
              hm_colors = viridis::inferno(100),
              hm_color_values = c(-2.5,-1.5,-1,-0.5,-0.25,0,0.25,0.5,1,1.5,2.5),
              hm_color_limits = c(-2.5,2.5),
              show_colnames = FALSE,
              clustering_method = "ward.D2",
              rows_title = "Protein expression",
              colorbar_dir = "horizontal",
              colors_title = "Scaled factor score",
              scale = TRUE,
              center = TRUE,
              show_dend_row = FALSE,
              show_dend_col = FALSE,
              group_prop = 0.02,
              group_lines = TRUE,
              group_line_color = "white",
              fontsize = 7,
              group_colors = ccp_cols) %>%
    add_tracks(
      track_columns = c('Stage', 'Xu proteome cluster', 'FGFR3 mutation'),
      track_prop = 0.06,
      track_pos = "top", 
      track_colors = list(
        'Stage' = c("Ta"="#fee0d2", "T1"="#fc9272", "MIBC"="#de2d26"),
        'Xu proteome cluster' = 'Set1',
        'FGFR3 mutation' = c("M" = "black", "WT" = "grey80")
      ),
      fontsize = 7
    )
  return(hm)
}

.plt_xu_classFlow <- function(final_sa) {
  cl_comp <- final_sa %>%
    rename(`Xu proteome cluster` = proteome_cluster) %>%
    classFlow("Xu proteome cluster", "uPG", 
              id_var = "proteome_id", alpha = 0.5, 
              label_colors = c("U-I" = "#e41a1c",
                               "U-II" = "#377eb8",
                               "U-III" = "#4daf4a",
                               ccp_cols,
                               "unc" = "grey60"), label_size = 3) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_blank(),
          axis.ticks = element_line(color = "black"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7, face = "bold"))
}

.get_top_markers <- function(diffProt_uPG, 
                             by = c("auc", "lfc")[1],
                             n = 20, 
                             pos_fc = FALSE) {
  if(by == "auc") {
    top_aucs_list <- lapply(LETTERS[1:5], function(upg) {
      auc_var <- paste0("auc_", upg)
      fc_var <- paste0("logFC_", upg)
      if(pos_fc) {
        diffProt_uPG %>%
          filter(!! fc_var > 0) %>%
          slice_max(.data[[auc_var]], n = n) %>%
          pull(symbol)
      } else {
        diffProt_uPG %>%
          slice_max(.data[[auc_var]], n = n) %>%
          pull(symbol)
      }
    })
    return(top_aucs_list)
  } else {
    top_lfc_var <- lapply(LETTERS[1:5], function(upg) {
      auc_var <- paste0("auc_", upg)
      fc_var <- paste0("logFC_", upg)
      if(pos_fc) {
        diffProt_uPG %>%
          filter(!! fc_var > 0) %>%
          slice_max(.data[[fc_var]], n = n) %>%
          pull(symbol)
      } else {
        diffProt_uPG %>%
          slice_max(.data[[fc_var]], n = n) %>%
          pull(symbol)
      }
    })
    return(top_lfc_var)
  }

  
}

################################################################################
## Old
################################################################################

# apoptosis_fgfr3_xu <- function(kegg_paths) {
#   xu_prot <- .pp_xu_prot()
#   xu_samp_annot <- .pp_xu_sa()
#   cl <- pull(xu_samp_annot, FGFR3_mutation, proteome_id)
#   
#   #-- FGFR3
#   difxu_FGFR3 <- findClassDiffs(xu_prot[,names(cl)], cl, type = "non-parametric",
#                                 pval_cutoff = 1, log_transformed = TRUE)
#   lfc_fgfr3 <- pull(difxu_FGFR3, logFC, feat) %>% sort()
#   
#   path_list <- split(kegg_paths$Symbol, kegg_paths$Path_Name)
#   # path_list <- split(react_paths$Symbol, react_paths$Path_Name)
#   gsea_res1 <- fgsea(path_list, lfc_fgfr3)
# 
#   plotEnrichment(path_list[["KEGG_APOPTOSIS"]], lfc_fgfr3)
#   
#   #-- Class E
#   cl <- upg_xu %>% mutate(class_E = ifelse(uPG == "E", "E", "non_E")) %>%
#     pull(class_E, proteome_id)
#   difxu_E <- findClassDiffs(xu_prot[,names(cl)], cl, type = "non-parametric",
#                                 pval_cutoff = 1, log_transformed = TRUE)
#   lfc_E <- pull(difxu_E, logFC, feat) %>% sort()
# 
#   gsea_res2 <- fgsea(path_list, lfc_E)
#   lfc_E
#   
#   pv.out <- pathview::pathview(gene.data = lfc_fgfr3, pathway.id = "04210",
#                                species = "hsa", gene.idtype = "SYMBOL")
#   
# }
# classify_strog <- function(upg_centroids, min_cor = 0.2) {
#   #-- Read protein expression data
#   str_sa <- read_excel("data/other_datasets/Stroggilos_IJC/ijc32556-sup-0001-tables.xlsx", 
#                        sheet = 2, skip = 2) %>%
#     select(id = PatientID, NPS = `Consensus Subtype`, Stage, Grade)
#   str_prot_log <- .pp_strog_prot()
#   
#   #-- "Batch correct" scaled data
#   str_prot_norm <- normalize_to_reference(prot_symbol, str_prot_log)
#   
#   # str_prot_nmibc <- str_prot_log[,str_sa$id]
#   # str_prot_mibc <- str_prot_log[, ! colnames(str_prot_scaled) %in% str_sa$id]
#   # 
#   # nmibc_samps <- samp_annot %>% filter(Stage_ == "NMIBC") %>% pull(sample)
#   # cit_nmibc <- prot_symbol[,nmibc_samps]
#   # 
#   # mibc_samps <- samp_annot %>% filter(Stage_ == "MIBC") %>% pull(sample)
#   # cit_mibc <- prot_symbol[,mibc_samps]
#   # 
#   # str_nmibc_norm <- normalize_to_reference(cit_nmibc, str_prot_nmibc)
#   # str_mibc_norm <- normalize_to_reference(cit_mibc, str_prot_mibc)
#   # 
#   # str_prot_norm <- cbind(str_mibc_norm, str_nmibc_norm)
#   
#   #-- Correlation to centroid
#   kept_feats <- intersect(rownames(str_prot_norm), rownames(upg_centroids))
#   centroid_cors <- cor(str_prot_norm[kept_feats,], upg_centroids[kept_feats,], use = "pair")
#   
#   colnames(centroid_cors) <- paste0("cor_", LETTERS[1:5])
#   class_res <- centroid_cors %>% as.data.frame() %>% rownames_to_column("id")
#   
#   new_class <- LETTERS[apply(centroid_cors, 1, which.max)]
#   max_cors <- apply(centroid_cors, 1, max)
#   new_class <- ifelse(max_cors > min_cor, new_class, "unc")
#   
#   upg_strog <- class_res %>%
#     mutate(uPG = new_class,
#            max_cor = max_cors) %>%
#     relocate(uPG, max_cor, .after = id) %>%
#     tibble()
#   
#   upg_strog <- full_join(upg_strog, str_sa, by = "id") %>%
#     mutate(Stage = ifelse(is.na(Stage), "MIBC", Stage),
#            NPS = ifelse(is.na(NPS), "MIBC", NPS))
#   
# }
# 
# plot_stroggilosHM <- function(upg_strog, upg_centroids) {
#   str_prot <- .pp_strog_prot()
#   pts4plot <- intersect(rownames(upg_centroids), rownames(str_prot))
#   upg_strog <- upg_strog %>%
#     mutate(Stage = ifelse(Stage == "T1 (CIS)", "T1", Stage))
#   
#   #-- Top mks
#   top_aucs_pos <- .get_top_markers(diffProt_uPG, n = 20, pos_fc = TRUE)
#   
#   top_aucs <- top_aucs_pos %>% unlist() %>% unique()
#   prots_used <- top_aucs[top_aucs %in% rownames(str_prot)]
#   idx <- rowSums(!is.na(str_prot[prots_used,])) > 100
#   prots_used <- prots_used[idx]
#   
#   top_aucs_list <- lapply(top_aucs_pos, function(aucs) {
#     aucs[aucs %in% prots_used]
#   })
#   names(top_aucs_list) <- LETTERS[1:5]
#   
#   data4plot <- str_prot[prots_used,] %>%
#     t() %>%
#     as.data.frame() %>%
#     rownames_to_column("id") %>%
#     left_join(upg_strog)
#   
#   
#   #-- Plot heatmap
#   hm <- .plt_strog_hm(data4plot)
#   plt_cf <- .plt_strog_classFlow(upg_strog)
#   
#   plt <- (hm | ((plt_cf + guides(fill = "none")) / plot_spacer()) + 
#             plot_layout(heights = c(1,1.5))) + plot_layout(widths = c(1.5,1))
#   
#   ggsave("results/suppfig_validation/stroggilos_v2.pdf", plt, width = 7, height = 5)
#   
# }
# .plt_strog_hm <- function(data) {
#   gghm <- 
#     data4plot %>%
#     group_by(uPG) %>%
#     ggheatmap(rowv = top_aucs_list,
#               hm_colors = viridis::inferno(100),
#               scale = TRUE, center = TRUE,
#               cluster_rows = FALSE,
#               show_dend_col = FALSE,
#               show_colnames = FALSE, 
#               hm_color_limits = c(-2.5,2.5),
#               hm_color_values = c(-2.5,-1,-0.5,0,0.5,1,2.5),
#               colors_title = bquote("Scaled"~log[2](ratio)),
#               dist_method = "spearman",
#               clustering_method = "ward.D2",
#               colv = "id", 
#               fontsize = 7,
#               group_prop = 0.03,
#               dend_prop_col = 0.05,
#               group_lines = TRUE,
#               group_line_color = "white",
#               group_lwd = 0.7,
#               group_colors = c(ccp_cols, unc = "grey50"),
#               raster = TRUE) %>%
#     add_tracks(track_columns = c("Stage", "NPS"),
#                track_colors = list('Stage' = c("Ta"="#fee0d2", "T1"="#fc9272", "MIBC"="#de2d26"), 
#                                    'NPS' = c("MIBC" = "grey40",
#                                              "NPS1" = "#acd785",
#                                              "NPS2" = "#a0c7db",
#                                              "NPS3" = "#1e74ae")
#                ),
#                track_pos = 'top', track_prop = 0.08, fontsize = 8)
#   return(gghm)
# }
# 
# .plt_strog_classFlow <- function(data) {
#   cl_comp <- data %>%
#     classFlow("NPS", "uPG", 
#               id_var = "id", alpha = 0.5, 
#               label_colors = c("MIBC" = "grey40",
#                                "NPS1" = "#acd785",
#                                "NPS2" = "#a0c7db",
#                                "NPS3" = "#1e74ae",
#                                ccp_cols,
#                                "unc" = "grey60"), label_size = 3) +
#     theme(axis.text = element_text(size = 8),
#           axis.title = element_blank(),
#           axis.ticks = element_line(color = "black"),
#           legend.text = element_text(size = 7),
#           legend.title = element_text(size = 7, face = "bold"))
# }
