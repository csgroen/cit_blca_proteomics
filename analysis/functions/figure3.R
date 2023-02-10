################################################################################
## Differential analysis
################################################################################
protDifferentialAnalysis <- function(prot_data, group_var, p_cutoff = 0.05) {
    wp <- prot_data$wp
    sa <- prot_data$sampAnnot
    pa <- prot_data$wpAnnot

    cl <- structure(pull(sa, !! group_var), names = sa$id)[colnames(wp)]

    prot_clDifs <- findClassDiffs(wp, cl, type = "non-parametric",
                                  pval_cutoff = p_cutoff, log_transformed = TRUE)

    prot_dif_results <- prot_clDifs %>%
        dplyr::rename(protein_id = feat) %>%
        left_join(select(pa, protein_id, symbol, subcell_loc)) %>%
        relocate(symbol, .after = protein_id)

    return(prot_dif_results)
}
mRNADifferentialAnalysis_FGFR3 <- function(mrna_data,
                                           samples = NULL,
                                           data_filter = "all") {
    if(is.null(samples)) {
        samples <- colnames(mrna_data$gexp)
    }

    if(data_filter == "NMIBC") {
        annot <- mrna_data$sampleAnnot %>%
            filter(id %in% samples, mibc == 0)
        design <- annot %>%
            select(id, FGFR3_mutation) %>%
            dummy_cols(select_columns = c("FGFR3_mutation")) %>%
            select(- ends_with("NA"), -FGFR3_mutation) %>%
            as.data.frame() %>%
            filter(across(everything(), ~ ! is.na(.))) %>%
            column_to_rownames("id")

        contrast_mat <- makeContrasts(FGFR3_mutation_M-FGFR3_mutation_WT,
                                      levels = design)
    } else if (data_filter == "MIBC") {
        annot <- mrna_data$sampleAnnot %>%
            filter(id %in% samples, mibc == 1)
        design <- annot %>%
            select(id, FGFR3_mutation) %>%
            dummy_cols(select_columns = c("FGFR3_mutation")) %>%
            select(- ends_with("NA"), -FGFR3_mutation) %>%
            as.data.frame() %>%
            filter(across(everything(), ~ ! is.na(.))) %>%
            column_to_rownames("id")
        contrast_mat <- makeContrasts(FGFR3_mutation_M-FGFR3_mutation_WT,
                                      levels = design)
    } else if (data_filter == "all") {
        annot <- mrna_data$sampleAnnot %>% filter(id %in% samples)
        #-- Make design table + contrasts
        design <- annot %>%
            select(id, FGFR3_mutation, mibc) %>%
            dummy_cols(select_columns = c("FGFR3_mutation")) %>%
            select(- ends_with("NA"), -FGFR3_mutation) %>%
            as.data.frame() %>%
            filter(across(everything(), ~ ! is.na(.))) %>%
            column_to_rownames("id")
        contrast_mat <- makeContrasts(FGFR3_mutation_M-FGFR3_mutation_WT, mibc,
                                      levels = design)
    }

    #-- Fit model
    gexp <- mrna_data$gexp[,rownames(design)]
    fit <- lmFit(gexp, design)
    fit2 <- contrasts.fit(fit,  contrasts = contrast_mat)
    fit2 <- eBayes(fit2)
    degtable <- topTable(fit2, coef="FGFR3_mutation_M - FGFR3_mutation_WT", n=Inf) %>%
        rownames_to_column("symbol")

    return(degtable)
}

################################################################################
## Pathway scores
################################################################################
pathEnrichmentFromFC <- function(diff_res, react_paths,  path_rep = 0.5, id = "symbol") {
    classes <- diff_res %>% select(starts_with("logFC")) %>% colnames() %>% str_remove("^logFC_")
    
    react_ls <- split(react_paths$Symbol, react_paths$Path_ID)
    react_ls <- lapply(react_ls, na.exclude)
    
    #-- Pathway representation
    keep_paths <- sapply(react_ls, function(path) {
        sum(path %in% pull(diff_res, id))/length(path) > path_rep
    })
    message("Keeping ", sum(keep_paths, na.rm = TRUE), " of ", length(react_ls), 
            " pathways with over ", path_rep*100, "% representation...")
    
    react_ls <- react_ls[keep_paths]
    
    #-- Select
    react_pnames <- react_paths %>%
        select(Path_ID, Path_Name) %>%
        filter(Path_ID %in% names(react_ls)) %>%
        distinct()
    
    if(length(classes) == 1) {
        message("Running pathway enrichment for a single class...")
        gr_logFC <- .get_grFC(diff_res, id = id)
        gr_logFC <- gr_logFC[!duplicated(gr_logFC)]
        gsea_scores <- fgsea(react_ls, gr_logFC, nproc = 1, eps = 0)
        # for (i in 1:length(react_ls)) {
        #     message("--- Testing ", names(react_ls)[i])
        #     test <- fgsea(react_ls[i], gr_logFC, nproc = 1)
        # }
        react_path_enrichment <- gsea_scores %>%
            as_tibble() %>%
            left_join(react_pnames, by = c("pathway" = "Path_ID")) %>%
            dplyr::rename(Path_ID = pathway) %>%
            relocate(Path_Name, .after = Path_ID) %>%
            mutate(leadingEdge = map(leadingEdge, paste, collapse = ", ") %>% unlist()) %>%
            arrange(padj)
    } else {
        react_path_enrichment <- lapply(classes, function(gr) {
            message("Running pathway enrichment for ", gr, "...")
            gr_logFC <- .get_grFC(diff_res, gr, id = id)
            gsea_scores <- fgsea(react_ls, gr_logFC) %>%
                as_tibble() %>%
                left_join(react_pnames, by = c("pathway" = "Path_ID")) %>%
                dplyr::rename(Path_ID = pathway) %>%
                relocate(Path_Name, .after = Path_ID) %>%
                mutate(leadingEdge = map(leadingEdge, paste, collapse = ", ") %>% unlist()) %>%
                arrange(padj)
        })
        names(react_path_enrichment) <- classes
    }
    
    return(react_path_enrichment)
}
.get_grFC <- function(data, gr = NULL, id = "symbol") {
    if(is.null(gr)) {
        df_grFC <- data %>%
            select(id = !! id, logFC)
    } else {
        df_grFC <- data %>%
            select(id = !!id, !! paste0("logFC_", gr)) %>%
            group_by(id) %>%
            summarize(across(contains("logFC"), mean), .groups = "drop") %>%
            rename_with(~ str_remove(., paste0("_", gr)))
    }
    gr_logFC <- sort(structure(df_grFC$logFC, names = df_grFC$id), decreasing = TRUE)
    
    return(gr_logFC)
}

join_pFGFR3paths <- function(pathScores_FGFR3, goBPScores_FGFR3, keggScores_FGFR3) {
    all_paths <- pathScores_FGFR3 %>%
        mutate(database = "Reactome") %>%
        bind_rows(mutate(goBPScores_FGFR3, database = "GO:BP")) %>%
        bind_rows(mutate(keggScores_FGFR3, database = "KEGG"))
}



pathScores_fromDEG <- function(degtable, path_table = react_paths) {
    #-- Get react annotation
    react_ls <- split(path_table$Symbol, path_table$Path_ID)
    react_ls <- lapply(react_ls, na.exclude)
    react_pnames <- path_table %>%
        select(Path_ID, Path_Name) %>%
        filter(Path_ID %in% names(react_ls)) %>%
        distinct()
    
    logFC <- structure(degtable$logFC, names = degtable$symbol) %>% sort(decreasing = TRUE)
    
    rgsea <- fgsea(react_ls, logFC) %>%
        as_tibble() %>%
        left_join(react_pnames, by = c("pathway" = "Path_ID")) %>%
        rename(Path_ID = pathway) %>%
        relocate(Path_Name, .after = Path_ID) %>%
        mutate(leadingEdge = map(leadingEdge, paste, collapse = ", ") %>% unlist()) %>%
        arrange(padj)
    
    return(rgsea)
}
join_mFGFR3paths <- function(mRNA_pathScores_FGFR3, mRNA_gobpScores_FGFR3, mRNA_keggScores_FGFR3) {
    all_paths <- mRNA_pathScores_FGFR3 %>%
        mutate(database = "Reactome") %>%
        bind_rows(mutate(mRNA_gobpScores_FGFR3, database = "GO:BP")) %>%
        bind_rows(mutate(mRNA_keggScores_FGFR3, database = "KEGG"))
    
}


################################################################################
## Delta pathway activity
################################################################################
deltaPathwayActivity <- function(protPaths_FGFR3, mrnaPaths_FGFR3,
                                 label = c("R-HSA-611105" = "(RDB) Respiratory electron transport", 
                                           "M23564" = "(GO-BP) ATP synthesis coupled electron transport", 
                                           "R-HSA-1428517" = "The citric acid (TCA) cycle and respiratory electron transport", 
                                           "R-HSA-111465" = "(RDB) Apoptotic cleavage of cellular proteins", 
                                           "R-HSA-75153" = "(RDB) Apoptotic execution phase", 
                                           "KEGG_APOPTOSIS" = "(KEGG) Apoptosis", 
                                           "R-HSA-69190" = "(RDB) DNA strand elongation", 
                                           "M10697" = "(GO-BP) Cell cycle/DNA replication", 
                                           "KEGG_MISMATCH_REPAIR" = "(KEGG) Mismatch repair",
                                           "KEGG_RIBOSOME" = "(KEGG) Ribosome",
                                           "R-HSA-72312" = "(RDB) rRNA processing"
                                 )) {
    mrna_tb <- mrnaPaths_FGFR3 %>%
        select(Path_ID, Path_Name, mrna_padj = padj, mrna_NES = NES)
    
    dnes_tb <- protPaths_FGFR3 %>%
        select(Path_ID, Path_Name, prot_padj = padj, prot_NES = NES, leadingEdge) %>%
        left_join(mrna_tb) %>%
        mutate(deltaNES = prot_NES - mrna_NES,
               abs_dNES = abs(deltaNES),
               dir_dNES = sign(deltaNES)) %>%
        arrange(desc(abs_dNES))
    
    label_annot <- dnes_tb %>%
        right_join(tibble(Path_ID = names(label), label = label)) %>%
        select(Path_ID, Path_Name, label)
    
    dnes4plot <- dnes_tb %>%
        left_join(label_annot)
    
    dnes4plot <- dnes4plot %>%
        mutate(deltaNES = ifelse(prot_padj <= 0.1, deltaNES, NA),
               padj = ifelse(is.na(deltaNES), "ns", "*"))
    
    # Plot
    plt <- ggplot(dnes4plot, aes(x = mrna_NES, y = prot_NES, 
                                 fill = deltaNES, label = label, 
                                 size = padj)) +
        geom_point(aes(alpha = padj), pch = 21) +
        geom_text_repel(min.segment.length = 0, size = 2.5) +
        geom_abline(slope = 1, intercept = 0, lty = "dashed", size = 0.5, color = "black") +
        geom_abline(slope = 1, intercept = 0.5, lty = "dashed", size = 0.5, color = "grey") +
        geom_abline(slope = 1, intercept = -0.5, lty = "dashed", size = 0.5, color = "grey") +
        geom_hline(yintercept = 0, lty = "dotted", size = 0.5, color = "black") +
        geom_vline(xintercept = 0, lty = "dotted", size = 0.5, color = "black") +
        scale_fill_distiller(palette = "PuOr") +
        scale_size_manual(values = c(1.5, 0.2)) +
        scale_alpha_manual(values = c(1, 0.3)) +
        lims(x = c(-3,3), y = c(-3,3)) +
        labs(x = "FGFR3* pathway activity (mRNA)", y = "FGFR3* pathway activity (Protein)") +
        theme_csg_scatter
    
    plt
    ggsave("results/fig3/fig3a.pdf", width = 4.5, height = 3)
    
    return(list(dNES = dnes_tb, label_annot = label_annot))
}

################################################################################
## Enrichment volcanos
################################################################################

deltaPath_volcanos <- function(deltaES) {
    dNES <- deltaES$dNES
    label_annot <- deltaES$label_annot
    
    
    plot_tb <- dNES %>%
        mutate(prot_mlog10 = -log10(prot_padj),
               mrna_mlog10 = -log10(mrna_padj),
               color_p = case_when(
                   prot_mlog10 >= 1 & prot_NES > 0 ~ "Up-regulated",
                   prot_mlog10 > 1 & prot_NES < 0 ~ "Down-regulated",
                   TRUE ~ "ns",
               ),
               color_m = case_when(
                   mrna_mlog10 >= 1 & mrna_NES > 0 ~ "Up-regulated",
                   mrna_mlog10 > 1 & mrna_NES < 0 ~ "Down-regulated",
                   TRUE ~ "ns",
               )) %>%
        left_join(select(mutate(label_annot, show = TRUE), Path_ID, label, show))
    
    mrna_plt <- ggplot(plot_tb, aes(mrna_NES, mrna_mlog10, label = label)) +
        geom_point(size = 0.5, aes(color = color_m)) +
        geom_text_repel(size = 2.5, min.segment.length = 0) +
        geom_hline(yintercept = 1, lty = "dashed") +
        geom_vline(xintercept = 0, lty = "dashed") +
        scale_color_manual(values = c("Up-regulated" = "#ef8a62", "Down-regulated" = "#67a9cf", "ns" = "grey80")) +
        labs(x = "mRNA Pathway Enrichment Score", y = "-log10(adj p)") +
        guides(color = "none") +
        scale_x_continuous(limits = c(-2.5,2.5)) +
        scale_y_continuous(limits = c(0,8)) +
        theme_csg_scatter +
        theme(title = element_text(face = "bold", size = 8))
        
    prot_plt <- ggplot(plot_tb, aes(prot_NES, prot_mlog10, label = label)) +
        geom_point(size = 0.5, aes(color = color_p)) +
        geom_text_repel(size = 2.5, min.segment.length = 0) +
        geom_hline(yintercept = 1, lty = "dashed") +
        geom_vline(xintercept = 0, lty = "dashed") +
        scale_color_manual(values = c("Up-regulated" = "#ef8a62", "Down-regulated" = "#67a9cf", "ns" = "grey80")) +
        labs(x = "Protein Pathway Enrichment Score", y = "-log10(adj p)", color = "Pathway status") +
        scale_x_continuous(limits = c(-2.5,2.5)) +
        scale_y_continuous(limits = c(0,8)) +
        theme_csg_scatter +
        theme(title = element_text(face = "bold", size = 8),
              axis.title.y = element_blank())
        
    nes_volcanos <- (mrna_plt | prot_plt) + 
        plot_layout(guides = "collect") + 
        plot_annotation(tag_levels =  "A")
    
    ggsave("results/fig3/fig3ab.pdf", nes_volcanos, width = 8, height = 3)
    
    return(nes_volcanos)
    
}


# plot_topDEP_FGFR3 <- function(wp_symbol, samp_annot, diffProt_FGFR3, dFC,
#                               n_top = 25) {
#     #-- Get proteins
#     prots_unique <- diffProt_FGFR3 %>%
#         arrange(desc(logFC)) %>%
#         slice(1:n_top) %>%
#         pull(symbol, protein_id)
#     
#     #-- Join data
#     prots4plot <- wp_symbol[prots_unique,] %>%
#         t() %>%
#         as.data.frame() %>%
#         rownames_to_column("id") %>%
#         left_join(select(samp_annot, id, Consensus, UROMOL = NMIBCclass, 
#                          uPG = CCP_cluster, FGFR3_mutation)) %>%
#         dplyr::rename(`FGFR3 mutation` = "FGFR3_mutation")
#     
#     #-- Get scale 
#     quants <- quantile(scale(t(wp_symbol[prots_unique,])), seq(0,1,by=0.01), na.rm = TRUE)
#     
#     #-- Plot
#     hm <- prots4plot %>%
#         group_by(`FGFR3 mutation`) %>%
#         ggheatmap(colv = "id",
#                   scale = TRUE,
#                   hm_colors = viridis::inferno(100),
#                   hm_color_values = scales::rescale(quants[2:100]),
#                   colors_title = "Scaled\nprotein expression",
#                   rowv = rev(as.character(prots_unique)),
#                   cluster_rows = FALSE,
#                   show_dend_row = FALSE,
#                   show_colnames = FALSE, 
#                   group_colors = yn_cols,
#                   group_prop = 0.05,
#                   group_lines = TRUE,
#                   group_line_color = "white",
#                   group_lwd = 1,
#                   fontsize = 10,
#                   raster = TRUE)
#     new_hm <- get_hmPlot(hm) +
#         theme(axis.text.y = element_text(size = 5))
#     hm <- update_hmPlot(hm, new_hm)
#     hm_tracks <- hm %>%
#         add_tracks(track_columns = c("UROMOL", "Consensus", "uPG"),
#                    track_colors = list("uPG" = ccp_cols,
#                                        "UROMOL" = nmibc_cols,
#                                        "Consensus" = consensus2),
#                    track_pos = "top", track_prop = 0.16, fontsize = 10)
#     
#     gene_order <- get_rowLevels(hm_tracks)
#     
#     fc_track <- dFC$dFC %>%
#         filter(gene %in% gene_order) %>%
#         select(gene, mrna_fc, prot_fc) %>%
#         mutate(gene = factor(gene, levels = gene_order)) %>%
#         pivot_longer(cols = mrna_fc:prot_fc, names_to = "data_level", values_to = "fc") %>%
#         mutate(data_level = str_remove(data_level, "_fc")) %>%
#         ggplot(aes(data_level, gene, fill = fc)) +
#         geom_raster() +
#         scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,3.5)) +
#         labs(x = "Data level", y = "", fill = "logFC") +
#         theme_csg_sparse +
#         theme(axis.text.y = element_blank(),
#               axis.ticks.y = element_blank(),
#               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#     
#     hm_panel <- align_to_hm(hm_tracks, fc_track, pos = "right", 
#                             legend_action = "collect", newplt_size_prop = 0.1)
#     
#     ggsave("results/fig3/fig3a.pdf", 
#            plot = hm_panel, width = 6.8, height = 3.3)
#     return(hm_panel)
# }

pathScores_prot_mrna <- function(deltaPathway, 
                                 paths = c("R-HSA-611105" = "(RDB) Respiratory electron transport", 
                                           "M22938" = "(GO-BP) Respiratory electron transport chain", 
                                           "M23564" = "(GO-BP) ATP synthesis coupled electron transport", 
                                           "R-HSA-1428517" = "(RDB) The citric acid (TCA) cycle and respiratory electron transport",
                                           "R-HSA-111465" = "(RDB) Apoptotic cleavage of cellular proteins", 
                                           "R-HSA-75153" = "(RDB) Apoptotic execution phase", 
                                           "KEGG_APOPTOSIS" = "(KEGG) Apoptosis",
                                           "R-HSA-72312" = "(RDB) rRNA processing", 
                                           "KEGG_RIBOSOME" = "(KEGG) Ribosome", 
                                           "R-HSA-156902" = "(RDB) Peptide chain elongation",
                                           "R-HSA-69190" = "(RDB) DNA strand elongation", 
                                           "M10697" = "(GO-BP) Cell cycle/DNA replication", 
                                           "KEGG_DNA_REPLICATION" = "(KEGG) DNA replication", 
                                           "KEGG_MISMATCH_REPAIR" = "(KEGG) KEGG Mismatch repair", 
                                           "R-HSA-5651801" = "(RDB) PCNA-dependent long patch base excision repair")) {
    
    dpath <- tibble(Path_ID = names(paths), label = paths)
    
    path_act <- deltaPathway$dNES %>%
        right_join(dpath) %>%
        mutate(label = factor(label, levels = label))
    
    
    path_act4plot <- path_act %>%
        pivot_longer(cols = c("prot_padj", "mrna_padj"), names_to = "level", values_to = "padj") %>%
        pivot_longer(cols = c("prot_NES", "mrna_NES"), names_to = "level2", values_to = "NES") %>%
        mutate(level = str_remove(level, "_padj"),
               level2 = str_remove(level2, "_NES")) %>%
        filter(level == level2) %>%
        select(-level2) %>%
        mutate(mlogpadj = -log10(padj),
               level = plyr::mapvalues(level, c("mrna", "prot"), c("Transcriptomics", "Proteomics")) %>%
                   factor(levels = c("Transcriptomics", "Proteomics")),
               label = factor(label, levels = rev(paths)),
               signif = ifelse(padj < 0.1, "*", "ns"))
    
    path_plt <- ggplot(path_act4plot, aes(level, label,
                                          color = signif, size = mlogpadj, fill = NES)) +
        geom_point(pch = 21) +
        scale_fill_distiller(type = "div", palette = "RdBu") +
        scale_color_manual(values = c("ns" = "grey60", "*" = "black")) +
        scale_size_continuous(breaks = c(1.3,2,5)) +
        theme_csg_scatter +
        theme(axis.text.x = element_text(angle=45, hjust = 1)) +
        labs(x = "omic", y = "Pathway")
    
    ggsave("results/fig3/fig3b.pdf", path_plt, width = 5.2, height = 4)
    return(path_plt)
}

################################################################################
## Apoptosis pathway plot
################################################################################
apoptosisPathHeatmap <- function(diffProt_FGFR3_unfilt, degtable_fgfr3_all) {
  #-- Apoptosis genes
  apopt_genes <- c("TRAIL" = "TNFSF10", "FAS" = "FASLG", "TNFa" = "TNF",
                   "TRAIL-R1" = "TNFRSF10A", "TRAIL-R2" = "TNFRSF10B", 
                   "TNFalpha" = "TNFRSF1A", "FADD" = "FADD", 
                   "TRADD" = "TRADD", "FLIP" = "CFLAR", "CASP10" = "CASP10", 
                   "CASP8" = "CASP8", "CASP6" = "CASP6", "CASP7" = "CASP7", 
                   "CASP3" = "CASP3", "Bad" = "BAD", "Noxa" = "PMAIP1",
                   "Hrk" = "HRK", "Bcl-2" = "BCL2", "BCL2L1" = "Bcl-XL",  
                   "Mcl-1"= "MCL1", "Bim" = "BCL2L11", "Bid" = "BID",
                   "Puma" = "BBC3", "Bax" = "BAX", "Bak" = "BAK1", "Cyt-c" = "CYCS", 
                   "Apaf-1" = "APAF1", "CASP9" = "CASP9", "SMAC/DIABLO" = "DIABLO",
                   "XIAP" = "XIAP")
  
  pdiff_apopt <- diffProt_FGFR3_unfilt %>%
    filter(symbol %in% apopt_genes) %>%
    select(symbol, prot_logFC = logFC)
  
  mdiff_apopt <- degtable_fgfr3_all %>%
    filter(symbol %in% apopt_genes) %>%
    select(symbol, trans_logFC = logFC)
  
  apopt_fc <- full_join(pdiff_apopt, mdiff_apopt) %>%
    mutate(symbol = factor(symbol, levels = apopt_genes),
           name = plyr::mapvalues(symbol, apopt_genes, names(apopt_genes)))
  
  
  hm <- ggheatmap(apopt_fc, 
                  colv = "name", 
                  rowv = c("trans_logFC", "prot_logFC"),
                  cluster_rows = FALSE, 
                  cluster_cols = FALSE,
                  hm_colors = "RdBu",
                  colors_title = "log2FC\n(FGFR3 mut vs WT)",
                  colorbar_dir = "horizontal",
  )
  
  ggsave(filename = "results/fig3/fig3c.pdf", plot = hm, width = 8, height = 3)
  
  return(hm)
}



