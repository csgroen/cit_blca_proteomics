################################################################################
## Proteomics PCA
################################################################################
proteomicsPCA <- function(prot_data, nfeats = "all") {
    #-- Params
    wp <- prot_data$wp %>% as.data.frame()
    
    if(nfeats != "all" & is.numeric(nfeats)) {
        vars <- apply(wp, 1, var, na.rm = TRUE) %>% sort(decreasing = )
        wp <- wp[names(vars)[1:nfeats],]
    } else if (nfeats == "complete") {
        completeness <- apply(wp, 1, function(v) { sum(is.na(v)) })
        wp_impt <- wp[names(completeness)[completeness == 0], ]
    }
    if (nfeats != "complete") {
        #-- Impute knn
        pp <- preProcess(wp, method = "knnImpute")
        wp_impt <- stats::predict(pp, wp)
    } 
    #-- Do PCA
    pca_res <- prcomp(t(wp_impt))
    
    return(pca_res)
}
plot_proteomicsPCA <- function(prot_data, pca_res, annot_vars, 
                               pcs = c(1,2),
                               gr_colors = NULL,
                               fname = "results/fig2/fig2a.pdf",
                               loadings = NULL,
                               add_ic_loadings = FALSE,
                               width = 4.5,
                               height = 2.5
) {
    annot <- prot_data$sampAnnot[,c("id", annot_vars)]
    if(length(annot_vars) == 2) {
        pc_plot <- pcaVisualization(pca_res, pcs = pcs, annotation = annot, color = annot_vars[1], shape = annot_vars[2],
                                    size = 2)
    } else {
        pc_plot <- pcaVisualization(pca_res, pcs = pcs, annotation = annot, color = annot_vars)
    }
    if(!is.null(gr_colors)) {
        pc_plot <- pc_plot +
            scale_color_manual(values = gr_colors)
    }
    #-- Add loadings
    all_loadings <- pca_res$rotation %>%
        as.data.frame() %>%
        rownames_to_column("protein_id") %>%
        left_join(select(prot_data$wpAnnot, protein_id, symbol)) %>%
        select(symbol, PC1, PC2) %>%
        tibble() 
    if(!is.null(loadings)) {
        #------ Get rotations
        loadings2plot <- filter(all_loadings, symbol %in% loadings)
        
        #------ Add segments
        pc_plot <- pc_plot +
            geom_segment(aes(x = 0, y = 0, xend = PC1*400, yend = PC2*400), 
                         arrow = arrow(length = unit(0.02, "npc")),
                         data = loadings2plot) +
            geom_text(aes(x = PC1*400, y = PC2*400, label = symbol), 
                      data = loadings2plot, size = 2.5, hjust = -0.1, vjust = -0.1)
    }
    #-- Calculate signature loadings
    if(add_ic_loadings) {
        load("data/annotations/biton_ics.RData")
        sigs <- c("IC Smooth muscle" = "M3_SMOOTH_MUSCLE_weight",
                  "IC Cell cycle" = "M7_CELLCYCLE_weight",
                  "IC Immune" = "M8_IMMUNE_weight",
                  "IC Urothelial differentiation" = "M9_UROTHELIALDIFF_weight",
                  "IC Ta pathway" = "M13_BLCAPATHWAYS_weight",
                  "IC Basal-like" = "M6_BASALLIKE_weight")
        
        ic_loads <- lapply(1:length(sigs), function(i) {
            sig <- sigs[[i]]
            loadings_up <- all_loadings %>% filter(symbol %in% biton_ls[[sig]]$top) %>%
                summarize(PC1 = mean(PC1), PC2 = mean(PC2)) %>%
                unlist()
            loadings_down <- all_loadings %>% filter(symbol %in% biton_ls[[sig]]$bottom) %>%
                summarize(PC1 = mean(PC1), PC2 = mean(PC2)) %>%
                unlist()
            loadings_down[is.na(loadings_down)] <- 0
            
            ic_loads <- loadings_up - loadings_down
            tibble(ic = names(sigs)[i],
                   PC1 = ic_loads[1],
                   PC2 = ic_loads[2])
        }) %>% bind_rows()
        pc_plot <- pc_plot +
            geom_segment(aes(x = 0, y = 0, xend = PC1*300, yend = PC2*300), 
                         arrow = arrow(length = unit(0.04, "npc")),
                         color = "grey60", alpha = 0.5, size = 2,
                         data = ic_loads) +
            geom_text(aes(x = PC1*400, y = PC2*400, label = ic), fontface = "bold",
                      data = ic_loads, size = 3, hjust = -0.1, vjust = -0.1)
    }
    ggsave(fname, 
           plot = pc_plot + labs(color = "Consensus/UROMOL subtype", shape = "Stage"), 
           width = width, height = height)
    
    return(pc_plot)
}

export_pcaResults <- function(prot_pca, prot_data, react_paths) {
    wpAnnot <- prot_data$wpAnnot
    prot_loadings <- prot_pca$rotation %>%
        as.data.frame() %>%
        rownames_to_column("protein_id") %>%
        left_join(select(wpAnnot, protein_id, symbol)) %>%
        select(protein_id, symbol, PC1:PC10) %>%
        tibble()
    react_tb <- select(react_paths, pathway_id = Path_ID, pathway = Path_Name) %>% distinct()
    pcs <- paste0("PC", 1:10)
    pc_enrich <- lapply(pcs, function(pc) {
      message("Running enrichment for ", pc)
      path_ls <- split(react_paths$Symbol, react_paths$Path_Name)
      pc_vec <- sort(pull(prot_loadings, {{pc}}, symbol))
      fgseaMultilevel(pathways = path_ls, stats = pc_vec, nproc = 1) %>%
        arrange(padj) %>%
        select(pathway, pval, padj, NES, size, leadingEdge) %>%
        rowwise() %>%
        mutate(leadingEdge = paste0(leadingEdge, collapse = ", ")) %>%
        filter(padj < 0.05) %>%
        tibble() %>%
        left_join(react_tb) %>%
        relocate(pathway_id)
    })
    names(pc_enrich) <- paste0(pcs, "_enrichment")
    pca_tables <- rlist::list.prepend(Proteomics_PCA = prot_loadings,
                        pc_enrich)
    
    # pca_tables <- list(Proteomics_PCA = prot_loadings,
    #                    PC1_enrichment = pc1_enrich,
    #                    PC2_enrichment = pc2_enrich)
    
    openxlsx::write.xlsx(x = pca_tables, file = "results/tables/PCA_proteomics.xlsx",
                         overwrite = TRUE)
    return(pca_tables)
}

################################################################################
## Big heatmap
################################################################################
plot_finalMOFA <- function(final_mofa, mcp_res, wp_symbol, sa) {
    dir.create("results/mofa/finalMofa", showWarnings = FALSE)
    dir.create("results/suppfig_mofa/", showWarnings = FALSE)
    
    mofa_trained <- final_mofa$mofa

    #-- Plots
    mofa_varExp <- .mofa_allVarExp(final_mofa)
    mofa_hm2 <- .mofa_finalHM(final_mofa, mcp_res, sa, wp_symbol)
    mofa_topGenes <- .factor_weights_hm(final_mofa)
    mofa_paths <- .factor_paths(final_mofa)
    mofa_assocs <- .mofa_associations(final_mofa)
    
    row1 <- (mofa_paths | mofa_assocs$others) + plot_layout(widths=c(1.4,1))
    row2 <- wrap_plots(mofa_topGenes, plot_spacer(), ncol = 2, widths = c(5,1))
    
    supp_plot <- (row1 / row2) + plot_layout(heights = c(3.3,1))
    ggsave("results/suppfig_mofa/mofa_suppFig.pdf", 
           supp_plot, width = 10, height = 13)
    
    ggsave("results/suppfig_mofa/mofa_varexp.pdf", mofa_varExp, width = 3.3, height = 2)
    
    sbty <- (mofa_assocs$subtype)
    ggsave("results/fig2/fig2e.pdf", sbty, width = 4.5, height = 2)
    
    return(list(
        mofa_hm = mofa_hm2,
        top_genes = mofa_topGenes,
        paths = mofa_paths,
        assocs = mofa_assocs
    ))
    
}
.mofa_allVarExp <- function(final_mofa) {
    var_exp <- get_variance_explained(final_mofa$mofa)[[2]][[1]]
    rownames(var_exp) <- plyr::mapvalues(rownames(var_exp), final_mofa$facts, final_mofa$fact_names)
    
    var_exp %>%
        as.data.frame() %>%
        rownames_to_column("factor_names") %>%
        pivot_longer(cols = cna:prot, names_to = "data_level", values_to = "var_exp") %>%
        mutate(factor_names = fct_inorder(factor_names)) %>%
        ggplot(aes(factor_names, data_level, fill = var_exp)) +
        geom_tile() +
        scale_fill_distiller(palette = "Blues", direction = 1) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_csg_sparse +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              panel.border = element_rect(color = "black", fill = "transparent")) +
        labs(x = "Factor names", y = "Data level", fill = "Variance\nexplained (%)")
}


.mofa_associations <- function(final_mofa) {
    mofa_trained <- final_mofa$mofa
    factors2show <- paste0("Factor",1:10)
    factors_select <- final_mofa$facts
    factor_names <- structure(plyr::mapvalues(factors2show, 
                                              from = factors_select, 
                                              to = final_mofa$fact_names),
                              names = factors2show)
    
    invert <- structure(rep(1,10), names =  factors2show)
    invert[factors_select] <- final_mofa$invert
    
    #-- Class
    assoc_sty <- associate_factors_with_covariates(mofa_trained, covariates = "Subtype",
                                                   fct_names = factors_select, 
                                                   invert = final_mofa$invert,
                                                   color_limits = c(-2,2),
                                                   size_limits = c(0,8)) +
        scale_y_discrete(labels = rev(factor_names)) +
        theme_mofa
    assoc_upg <- associate_factors_with_covariates(mofa_trained, covariates = "uPG",
                                                   fct_names = factors_select, 
                                                   invert = final_mofa$invert,
                                                   color_limits = c(-2,2),
                                                   size_limits = c(0,8)) +
        scale_y_discrete(labels = rev(factor_names)) +
        theme_mofa
    
    cl_assocs <-  ((assoc_upg + labs(x = "uPG")) | (assoc_sty 
                                                    + labs(x = "Subtype") +
                                                        theme(axis.text.y=element_blank(),
                                                              axis.title.y=element_blank(),
                                                              axis.ticks.y = element_blank()))) + 
        plot_layout(guides = "collect", widths = c(1,2))
    #-- Genomic alterations
    facI_mut <- associate_factors_with_covariates(mofa_trained,
                                                  fct_names = factors2show, 
                                                  invert = invert,
                                                  c("FGFR3_mutation", "TP53_mutation",
                                                    "PIK3CA_mutation", "RAS_mutation"),
                                                  color_limits = c(-1.5,1.5),
                                                  size_limits = c(0,6)) +
        theme_mofa +
        scale_y_discrete(labels = rev(factor_names))
    chr_alts <- c("9q31_1_del", "11p15_3_del", 
                  "1p36_32_del", "13q13_1_amp", 
                  "4q35_2_del", "6p22_3_amp", 
                  "18p11_31_amp", "3p14_2_del", 
                  "11q13_3_amp", "3p25_2_amp", 
                  "8p22_del", "1p12_amp", 
                  "1q23_3_amp")
    chr_labs <- chr_alts %>% 
        str_replace("_([0-9])", ".\\1") %>%
        str_replace_all("_", " ")
    facI_chr <- associate_factors_with_covariates(mofa_trained, 
                                                  covariates =  chr_alts,
                                                  invert = invert, 
                                                  cluster_cols = TRUE,
                                                  color_limits = c(-1.5,1.5),
                                                  size_limits = c(0,6)) +
        theme_mofa +
        scale_x_discrete(labels = structure(chr_labs, names = chr_alts)) +
        scale_y_discrete(labels = rev(factor_names)) +
        labs(x = "CNA")
    
    cna_genes <- c("TP53_loss_bin", "ERBB2_gain_bin",
                   "RB1_loss", "ELF3_loss_bin", "CDKN2A_loss", "PPARG_gain")
    cna_names <- cna_genes %>% str_replace_all("_", " ") %>% str_remove(" bin")
    facI_cnaGenes <- associate_factors_with_covariates(mofa_trained, 
                                                       cna_genes, 
                                                       fct_names = factors2show, 
                                                       invert = invert,
                                                       cluster_cols = TRUE,
                                                       color_limits = c(-1.5,1.5),
                                                       size_limits = c(0,6)) +
        theme_mofa +
        scale_x_discrete(labels = structure(cna_names, names = cna_genes)) +
        scale_y_discrete(labels = factor_names) +
        labs(x = "CNA")
    
    plt_genAlt <- ((((facI_mut + labs(x = "Mutation")) |
                         (facI_cnaGenes + labs(x = "CNA") +
                              theme(axis.text.y = element_blank(),
                                    axis.title.y = element_blank(),
                                    axis.ticks.y = element_blank()))) +
                        plot_layout(widths = c(2,3))) /
                       facI_chr + labs(x = "CNA")) +
        plot_layout(guides = "collect")
    plt_genAlt
    
    #-- trICs 
    trICs <- c("IC Basal-like", "IC Mitochondria", "IC Cell cycle", "IC Urothelial differentiation", 
               "IC Bladder cancer pathways", "IC Neuroendocrine", "IC Smooth muscle", 
               "IC Lymphocytes B&T", "IC Myofibroblasts", "IC Interferon response"
    )
    plt_trIC <- associate_factors_with_covariates(mofa_trained, covariates = trICs,
                                                  invert = invert,
                                                  color_limits = c(-1,1),
                                                  size_limits = c(0,6),
                                                  color_low = "#5e3c99",
                                                  color_high = "#e66101") +
        theme_mofa +
        scale_y_discrete(labels = factor_names) +
        labs(x = "Transcriptomic ICs")
    
    all_associations <- (plt_genAlt / plt_trIC)
    
    return(list(subtype = cl_assocs, others = all_associations))
}

.mofa_finalHM <- function(final_mofa, mcp_res, sa, wp_symbol) {
    mofa_trained <- final_mofa$mofa
    factors2show <- final_mofa$facts
    factor_names <- final_mofa$fact_names
    invert <- final_mofa$invert
    
    #-- Get factor HM
    fct_values <- get_factors(mofa_trained)[[1]]
    fct_mat <- sapply(1:length(factors2show), function(i) {
        fct <- factors2show[i]; inv <- invert[i]
        vals <- fct_values[,fct] * inv
        return(vals)
    })
    colnames(fct_mat) <- factor_names
    pops <- c("T cells", "Cytotoxic lymphocytes", "B lineage", "Neutrophils",
              "Monocytic lineage", "Fibroblasts")
    
    mks <- list("Basal phenotype" = c("KRT6A", "KRT14", "S100A8", "TYMP", "KRT5"),
                "Cell cycle/DNA repair" = c("TOP2A", "TMPO", "FEN1", "UHRF1", "KPNA2"),
                "Lipid metabolism" = c("AKR1C3", "GDPD3", "BCAT2", "ACOX1", "ABCD3"),
                "ECM interactions/Cytoskeleton/Smooth muscle" =  c("FLNC", "TPM2", "TAGLN", 
                                                                   "TPM1", "CNN1"),
                "Ta pathway/Urothelial differentiation" = c("ANXA10", "S100P",
                                                            "IVL", "KRT7", "DHRS2"))
    mks_unlist <- as.character(unlist(mks))
    
    # mks <- c("KRT6A", "KRT14", "SERPINE1", "TOP2A", "MKI67", "AKR1C3", "GDPD3", 
    #          "TAGLN", "DMD", "S100P", "KRT7", "ADIRF")
    
    fct_mat <- rbind(fct_mat, matrix(rep(0, 5), 
                                     nrow = 1, 
                                     dimnames = list("CIT.44", colnames(fct_mat))))
    
    data4hm <- fct_mat %>%
        as.data.frame() %>%
        rownames_to_column("sample") %>%
        left_join(rename(sa, sample = id)) %>%
        left_join(rownames_to_column(mcp_res, "sample")) %>%
        left_join(rownames_to_column(as.data.frame(t(wp_symbol[mks_unlist,])), "sample")) %>%
        rename(uPG = CCP_cluster)
    
    mofa_hm <- .mofa_factor_hm(data4hm, factor_names, pops, mks)
    var_exp <- .variance_exp_plot(mofa_trained, mofa_hm, factors2show, factor_names)
    complete_hm <- align_to_hm(mofa_hm, var_exp, pos = "right", newplt_size_prop = 0.13, legend_action = "collect")
    ggsave("results/fig2/fig2b.pdf", width = 7, height = 6.5)
    
    return(complete_hm)
}

# .factor_weights_hm <- function(final_mofa) {
#     mofa_trained <- final_mofa$mofa
#     factors2show <- final_mofa$facts
#     factor_names <- final_mofa$fact_names
#     invert <- final_mofa$invert
#     mofa_weights <- .plot_mofa_topGenes(mofa_trained, "mofa_mrna10000_cna1500_10f", 
#                                         factors = factors2show, invert = invert, width = 7, height = 7.5)
#     mofa_weights <- lapply(1:length(factor_names), function(i){ 
#         mofa_weights[[i]] + 
#             labs(subtitle = factor_names[i])})
#     mofa_weights[2:5] <- lapply(mofa_weights[2:5], function(plt) { plt + theme(axis.title.y = element_blank())})
#     
#     mofa_top_genes <- wrap_plots(mofa_weights, guides = "collect", nrow = 1) &
#         theme(axis.text.y = element_text(size=5),
#               axis.text.x = element_text(angle=90,hjust=0.5,vjust=0.5))
#     ggsave("results/mofa/finalMofa/factor_topGenes.pdf", width = 8, height = 4)
#     return(mofa_top_genes)
# }

.factor_paths <- function(final_mofa) {
    mofa_trained <- final_mofa$mofa
    factors2show <- paste0("Factor",1:10)
    factors_select <- final_mofa$facts
    factor_names <- plyr::mapvalues(factors2show, from = factors_select, to = final_mofa$fact_names)
    
    invert <- structure(rep(1,10), names =  factors2show)
    invert[factors_select] <- final_mofa$invert
    mofa_paths <- .plot_mofa_factorPathways(mofa_trained, base_name = "finalMofa", 
                                            factors = factors2show, invert = invert, n = 2)
    mofa_paths <- lapply(1:length(factor_names), function(i){ 
        mofa_paths[[i]] + 
            labs(subtitle = factor_names[i])})
    mofa_path_plt <- wrap_plots(mofa_paths, guides = "collect", ncol = 2) &
        theme(plot.subtitle = element_text(size=8,face="bold"),
              axis.text.y=element_text(size=7))
    
    ggsave("results/mofa/finalMofa/factor_paths.pdf", mofa_path_plt, 
           width = 6, height = 10)
    return(mofa_path_plt)
}

sampAnnot_core <- function(mrna, prot, cna_calls, core_samples) {
    #-- Load files
    cna_cyto_file <- "Data/CNA_fromCGH/CN_cytoArm_GISTIC.txt"
    cna_calls_all <- read_tsv("Data/CNA_fromCGH/CN_geneCalls_GISTIC.txt")
    
    sa_mrna <- mrna$sampleAnnot %>% select(-Stage., -Consensus, -NMIBCclass, 
                                           -Subtype, -ERBB2.gain, -FGFR3_mutation,
                                           -PIK3CA_mutation, -RAS_mutation,
                                           -matches("(loss|gain|mut)$"))
    sa_prot <- prot$sampAnnot %>% select(-ICA.smoothMuscle, -matches("(loss|gain)"),
                                         CDKN2A_loss, RB1_loss, PIK3CA_mutation, RAS_mutation,
                                         PPARG_amp)
    
    sa_cyto <- read_tsv(cna_cyto_file) %>%
        filter(str_detect(`Unique Name`, "CN values", negate = TRUE)) %>%
        select(type = `Unique Name`, cytoband = Descriptor, starts_with("CIT")) %>%
        mutate(type = str_remove(type, " .*")) %>%
        pivot_longer(cols = starts_with("CIT"), names_to = "id", values_to = "cna") %>%
        mutate(type = ifelse(type == "Amplification", "amp", "del"),
               cytoband = paste0(cytoband, "_", type),
               cna = factor(cna, levels = c(1,0))) %>%
        select(-type) %>%
        pivot_wider(id_cols = id, names_from = cytoband, values_from = cna)
    
    #-- Genomic instability
    # cna_mat <- cna_calls_all %>%
    #     select(-`Locus ID`, -Cytoband) %>%
    #     as.data.frame() %>%
    #     column_to_rownames("Gene Symbol")
    # 
    # gii_tb <- tibble(id = colnames(cna_mat), 
    #                  Genomic_Instability = apply(cna_mat, 2, function(col) { sum(col != 0)/nrow(cna_mat) }))
    
    genes_per_chr <- cna_calls_all %>%
        mutate(chr = str_sub(Cytoband, 1, 1)) %>%
        pull(`Gene Symbol`, chr)
    genes_per_chr <- split(genes_per_chr, names(genes_per_chr))
    
    cna_mat <- cna_calls_all %>%
        select(-`Locus ID`, -Cytoband) %>%
        as.data.frame() %>%
        column_to_rownames("Gene Symbol")
    
    galt_per_chr <- sapply(genes_per_chr, function(genes) {
        apply(cna_mat[genes,], 2, function(col) { sum(col != 0)/length(genes) })
    })
    
    gii_tb <- tibble(id = colnames(cna_mat), 
                     Genomic_Instability = apply(cna_mat, 2, function(col) { sum(col != 0)/nrow(cna_mat) }),
                     wGII = rowMeans(galt_per_chr))
    
    #-- Select and join
    sa <- left_join(sa_mrna, 
                    select(sa_prot, id, Stage., 
                           uPG = CCP_cluster, Subtype, Consensus, NMIBCclass, FGFR3_mutation,
                           CDKN2A_loss, RB1_loss, PIK3CA_mutation, RAS_mutation, RXRA_mutation,
                           PPARG_mutation, PPARG_amp)) %>%
        left_join(sa_cyto, by = "id") %>%
        left_join(cna_calls, by = c("id" = "sample")) %>%
        left_join(gii_tb, by = "id") %>%
        dplyr::rename(sample = id, stage = Stage.y) %>%
        select(-old, -exclude, -CNV.ID, -inconsensus,-mcl_consensus,-keep,-`OS..Days.`, 
               -survkeep,-NAC_free_Seiler, -response_seiler,-pStage,
               -pred, -ends_with("mut"), -mutation_load, -cisplatin, -dataset, -CT,
               -histology, -CDKN2A_mutation, -CDKN2A_lossMLPA, 
               -CDKN2A_MLPA, -CDKN2A.mlpa, -RB1_MLPA, -RB1_lossMLPA, -FGFR3_amp) %>%
        rename_with(.fn = ~ str_replace_all(., "\\.", "_")) %>%
        mutate(across(matches("[A-Z].*_(gain)$"), ~ ifelse(. > 0, "Gain", "WT"), .names = "{.col}_bin"),
               across(matches("[A-Z].*_(loss)$"), ~ ifelse(. < 0, "Loss", "WT"), .names = "{.col}_bin")) %>%
        filter(sample %in% core_samples) %>%
        select(-CDKN2A_loss_y, -RB1_loss_y) %>%
        rename(CDKN2A_loss = CDKN2A_loss_x, RB1_loss = RB1_loss_x) %>%
        as.data.frame()
    rownames(sa) <- sa$sample
    
    return(sa)
}

.mofa_factor_hm <- function(data4hm, factor_names, pops, mks) {
    #-- Heatmap
    factor_hm <- data4hm %>%
        group_by(uPG) %>%
        select(-Stage) %>%
        dplyr::rename(UROMOL = NMIBCclass, Stage = StageTa) %>%
        # sample = factor(sample, levels = sample_order)) %>%
        ggheatmap(
            colv = "sample",
            rowv = rev(factor_names),
            hm_colors = viridis(100),
            hm_color_limits = c(-2,2),
            group_lines = TRUE,
            group_line_color = "white",
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            show_colnames = FALSE,
            rows_title = "Factor scores",
            colorbar_dir = "horizontal",
            colors_title = "Scaled factor score",
            group_prop = 0.12,
            fontsize = 9) %>%
        add_tracks(track_columns = c("Consensus", "UROMOL", "Stage"), 
                   track_colors = list(Consensus = consensus2,
                                       UROMOL = nmibc_cols,
                                       Stage = c("Ta" = "#fee0d2",
                                                 "T1" = "#fc9272",
                                                 "MIBC" = "#de2d26")),
                   track_pos = "top",
                   fontsize = 9, track_prop = 0.3) %>%
        add_tracks(track_columns = c("TP53_mutation", "FGFR3_mutation", "RAS_mutation",
                                     "CDKN2A_loss", "RB1_loss", "PPARG_amp", "PPARG_mutation",
                                     "RXRA_mutation",
                                     "Genomic_Instability"),
                   track_colors = list(
                       TP53_mutation = yn_cols,
                       FGFR3_mutation = yn_cols,
                       PPARG_mutation = yn_cols,
                       RXRA_mutation = yn_cols,
                       RAS_mutation = c("M" = "black", "WT" = "grey80"),
                       CDKN2A_loss = c("Loss" = "black", "WT" = "grey80"),
                       RB1_loss = c("Loss" = "black", "WT" = "grey80"),
                       PPARG_amp = c("Amp" = "black", "WT" = "grey80"),
                       Genomic_Instability = "OrRd"
                   ),
                   track_prop = 0.45,
                   fontsize = 9) %>%
        add_matrix_track(track_columns = as.character(unlist(mks)),
                         fontsize = 9, 
                         rows_title = "Protein expression",
                         color_limits = c(-4,4),
                         track_colors = inferno(100),
                         pal_dir = -1,
                         track_prop = 0.5) %>%
        add_matrix_track(track_columns = pops,
                         fontsize = 9, 
                         rows_title = "MCPcounter scores",
                         color_limits = c(-2,2),
                         track_colors = magma(100),
                         pal_dir = -1,
                         track_prop = 0.15)
}
.variance_exp_plot <- function(mofa_trained, factor_hm, factors2show, factor_names) {
    var_exp <- mofa_trained@cache $variance_explained$r2_per_factor$group1 %>%
        as.data.frame() %>%
        rownames_to_column("feat") %>%
        filter(feat %in% factors2show) %>%
        mutate(feat = plyr::mapvalues(feat, from = factors2show, to = factor_names)) %>%
        mutate(feat = factor(feat, levels = get_rowLevels(factor_hm))) %>%
        pivot_longer(cols = -feat, names_to = "level", values_to = "variance_exp") %>%
        mutate(variance_exp = variance_exp/100)
    
    var_exp_plt <- ggplot(var_exp, aes(level, feat, fill = variance_exp)) +
        geom_tile() +
        scale_fill_distiller(palette = "Blues", direction = 1, labels = scales::percent) +
        guides(fill = guide_colorbar(direction = "horizontal")) +
        theme_sparse2(base_size = 9) +
        theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
}

.plot_mofa_topGenes <- function(mofa_trained, base_name, factors = NULL, invert = NULL, n = 5,
                                width = 8, height = 5) {
    weights <- get_weights(mofa_trained)
    if(!is.null(factors)) weights <- lapply(weights, function(w) w[,factors])
    if(!is.null(invert)) {
        weights_new <- lapply(weights, function(w) {
            new_w <- sapply(1:ncol(w), function(i) { w[,i]*invert[i]})
            colnames(new_w) <- colnames(w)
            return(new_w)
        } ) 
        weights <- weights_new
    }
    weights_rn <- lapply(names(weights), function(name) {
        mat <- weights[[name]]
        rownames(mat) <- str_remove(rownames(mat), str_glue("_{name}"))
        return(mat)
    })
    names(weights_rn) <- names(weights) 
    plt_fct <- lapply(colnames(weights[[1]]), function(facti) {
        fct_top <- lapply(weights_rn, function(lvl_w) { 
            sorted_genes <- lvl_w[,facti] %>% sort() %>% names()
            c(sorted_genes[1:n], sorted_genes[(length(sorted_genes)-n):length(sorted_genes)])
        })
        fct_top <- Reduce(union, fct_top)
        
        fct_weights <- lapply(names(weights_rn), function(lvl) {
            lvl_w <- weights_rn[[lvl]]
            genes4plot <- intersect(rownames(lvl_w), fct_top)
            tibble(gene = genes4plot, weight = lvl_w[genes4plot,facti], lvl = lvl)
        })
        fct_weights <- do.call(bind_rows, fct_weights)
        order <- fct_weights %>%
            group_by(gene) %>%
            summarize(mean = mean(weight, na.rm = TRUE)) %>%
            arrange(mean) %>%
            pull(gene)
        ggplot(fct_weights, aes(lvl, factor(gene, levels = order), fill = weight)) +
            geom_tile() +
            scale_fill_viridis(option = "inferno", limits = c(-1,1), oob = scales::squish) +
            labs(x = "Data level", y = "Gene symbol", subtitle = facti) +
            theme_csg_sparse +
            theme(axis.text = element_text(size = 8),
                  axis.title = element_text(size=8, face="bold"),
                  plot.subtitle = element_text(size=8, face="bold"))
    })
    
    plt_top_genes <- wrap_plots(plt_fct, nrow = 2, guides = "collect")
    ggsave(str_glue("results/mofa/{base_name}/top_genes.pdf"), width = width, height = height)
    
    return(plt_top_genes)
}

.plot_mofa_QC <- function(mofa_trained, base_name) {
    #-- QC plots
    qc_cor <- plot_factor_cor(mofa_trained)
    qc_plot <- ggcorrplot::ggcorrplot(qc_cor$corr, method = "circle")
    var_plot1 <- plot_variance_explained(mofa_trained) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
    var_plot2 <- plot_variance_explained(mofa_trained, plot_total = T)[[2]]
    
    all_qc_plots <- (qc_plot / (var_plot1 | var_plot2)) + plot_layout(heights = c(1.5,1))
    
    ggsave(str_glue("results/mofa/{base_name}/qcPlots.pdf"),
           all_qc_plots, width = 6, height = 8)
}

.plot_mofa_associations <- function(mofa_trained, base_name) {
    sa <- samples_metadata(mofa_trained)
    #-- Factor sample associations with covariates
    #--------- Subtype
    factI_subtype <- associate_factors_with_covariates(mofa_trained, "Subtype",
                                                       color_limits = c(-2,2),
                                                       size_limits = c(0,12)) +
        theme_mofa
    factI_upg <- associate_factors_with_covariates(mofa_trained, "uPG",
                                                   color_limits = c(-2,2),
                                                   size_limits = c(0,12)) +
        theme_mofa
    sbty <- ((factI_subtype +
                  labs(x = "Subtype")) | 
                 (factI_upg + 
                      labs(x = "uPG") +
                      theme(axis.text.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.ticks.y = element_blank()))) +
        plot_layout(guides = "collect", widths = c(1.2,1))
    ggsave(str_glue("results/mofa/{base_name}/assocSubtype.pdf"),
           sbty, width = 5, height = 4)
    
    #------- Genomic alterations
    facI_mut <- associate_factors_with_covariates(mofa_trained,
                                                  c("FGFR3_mutation", "TP53_mutation"),
                                                  color_limits = c(-2,2),
                                                  size_limits = c(0,8)) +
        theme_mofa
    facI_chr <- associate_factors_with_covariates(mofa_trained, 
                                                  str_subset(colnames(sa), "[0-9].*(amp|del)"),
                                                  cluster_cols = TRUE,
                                                  color_limits = c(-2,2),
                                                  size_limits = c(0,8)) +
        theme_mofa
    facI_cnaGenes <- associate_factors_with_covariates(mofa_trained, 
                                                       str_subset(colnames(sa), "(gain|loss)_bin"),
                                                       cluster_cols = TRUE,
                                                       color_limits = c(-2,2),
                                                       size_limits = c(0,8)) +
        theme_mofa
    
    genAlt <- ((((facI_mut + labs(x = "Mutation")) |
                     (facI_cnaGenes + labs(x = "CNA") +
                          theme(axis.text.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.ticks.y = element_blank()))) +
                    plot_layout(widths = c(1,3))) /
                   facI_chr + labs(x = "CNA")) +
        plot_layout(guides = "collect")
    
    ggsave(str_glue("results/mofa/{base_name}/assocGeneticAlt.pdf"), genAlt, width = 7,
           height = 7)
    
    #------ trIC
    factI_tIC <- associate_factors_with_covariates(mofa_trained, 
                                                   str_subset(colnames(sa), "^IC"),
                                                   color_limits = c(-0.7,0.7),
                                                   cluster_cols = TRUE) +
        theme_mofa
    ggsave(str_glue("results/mofa/{base_name}/assocTranscriptomicICs.pdf"), factI_tIC,
           width = 4, height = 5)
    
    return(list(subtypes = sbty,
                genAlt = genAlt,
                tICs = factI_tIC))
}

.plot_mofa_factorHeatmap <- function(mofa_trained, base_name,
                                     factor_names = "all") {
    #-- Heatmap
    factors <- get_factors(mofa_trained)[[1]]
    if(factor_names != "all") {
        factors <- factors[,factor_names]
    }
    
    sa <- samples_metadata(mofa_trained)
    factor_table <- factors %>%
        as.data.frame() %>%
        rownames_to_column("sample") %>%
        left_join(sa) %>%
        dplyr::rename(Stage = Stage_, UROMOL = NMIBCclass)
    
    factor_hm <- factor_table %>%
        group_by(uPG) %>%
        ggheatmap(colv = "sample",
                  rowv = rev(colnames(factors)),
                  show_colnames = FALSE,
                  cluster_rows = FALSE,
                  hm_color_limits = c(-1.5,1.5),
                  hm_colors = viridis(100),
                  group_lines = TRUE,
                  group_line_color = "white",
                  colorbar_dir = "horizontal",
                  colors_title = "Factor score",
                  raster = TRUE,
                  fontsize = 8,
                  group_prop = 0.06) %>%
        add_tracks(track_columns = c("Consensus", "UROMOL", "Stage"),
                   track_colors = list(
                       Consensus = consensus2,
                       UROMOL = nmibc_cols,
                       Stage = inv_cols
                   ),
                   fontsize = 8,
                   track_pos = "top", 
                   track_prop = 0.2, 
                   legend_action = "collect")
    #----- Variance explained track
    var_exp <- mofa_trained@cache $variance_explained$r2_per_factor$group1 %>%
        as.data.frame() %>%
        rownames_to_column("feat") %>%
        mutate(feat = factor(feat, levels = get_rowLevels(factor_hm))) %>%
        pivot_longer(cols = -feat, names_to = "level", values_to = "variance_exp")
    
    var_exp_plt <- ggplot(var_exp, aes(level, feat, fill = variance_exp)) +
        geom_tile() +
        scale_fill_distiller(palette = "Blues", direction = 1) +
        guides(fill = guide_colorbar(direction = "horizontal")) +
        theme_sparse2(base_size = 8) +
        theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
    
    full_hm <- align_to_hm(factor_hm, var_exp_plt, newplt_size_prop = 0.2,
                           pos = "right", legend_action = "collect")
    
    full_hm
    
    ggsave(str_glue("results/mofa/{base_name}/factorHeatmap.pdf"), 
           plot = full_hm, width = 6, height = 3)
    
    return(full_hm)
}

.plot_mofa_factorPathways <- function(mofa_trained, base_name = NULL, 
                                      factors = NULL, invert = NULL, n = 3) {
    
    #-- Pathway enrichment
    mrna_enrich <- get_enrichment(mofa_trained, level = "mrna", p_cutoff = 1)
    prot_enrich <- get_enrichment(mofa_trained, level = "prot", p_cutoff = 1)
    cna_enrich <- get_enrichment(mofa_trained, level = "cna", p_cutoff = 1)
    
    if(is.null(factors)) {
        factors <- colnames(get_factors(mofa_trained)[[1]])
    }
    if(is.null(invert)) {
        invert <- rep(1, length(factors))
    }
    
    #-- Get top paths 
    top_paths <- lapply(list(mrna_enrich, prot_enrich, cna_enrich), 
                        function(enrich_tb) {
                            paths_up <- enrich_tb %>%
                                group_by(factor) %>%
                                slice_min(mean_diff, n = n) %>%
                                pull(pathway) %>%
                                unique()
                            
                            paths_down <- enrich_tb %>%
                                group_by(factor) %>%
                                slice_max(mean_diff, n = n) %>%
                                pull(pathway) %>%
                                unique()
                            
                            union(paths_up, paths_down)
                        })
    top_paths <- Reduce(union, top_paths)
    
    #-- Filter significant and join
    mrna_enrich_f <- mrna_enrich %>%
        mutate(mlog10 = -log10(padj)) %>%
        select(pathway, factor, mrna_padj = mlog10, mrna_diff = mean_diff) %>%
        filter(pathway %in% top_paths)
    prot_enrich_f <- prot_enrich %>%
        mutate(mlog10 = -log10(padj)) %>%
        select(pathway, factor, prot_padj = mlog10, prot_diff = mean_diff) %>%
        filter(pathway %in% top_paths)
    cna_enrich_f <- cna_enrich %>%
        mutate(mlog10 = -log10(padj)) %>%
        select(pathway, factor, cna_padj = mlog10, cna_diff = mean_diff) %>%
        filter(pathway %in% top_paths)
    
    pathway_enrich <- full_join(mrna_enrich_f, prot_enrich_f, 
                                by = c('factor', 'pathway')) %>% 
        full_join(cna_enrich_f,  by = c('factor', 'pathway'))
    
    #-- Make factor pathway plots
    var_exp <- mofa_trained@cache$variance_explained$r2_per_factor$group1
    slice_vars <- structure(paste0(colnames(var_exp)[apply(var_exp, 1, which.max)], "_padj"),
                            names = rownames(var_exp))[factors]
    
    path_plots <- lapply(1:length(factors), function(i) {
        fact <- factors[i]
        slice_var <- slice_vars[i]
        inv <- invert[i]
        paths2show <- pathway_enrich %>%
            filter(
                factor == fact,
                prot_padj > 2 | mrna_padj > 2 | cna_padj > 2) %>%
            arrange(across(matches(slice_var), dplyr::desc)) %>%
            slice(n = 1:10) %>%
            mutate(
                cna_diff = cna_diff * inv,
                mrna_diff = mrna_diff * inv,
                prot_diff = prot_diff * inv,
                pathway = str_sub(pathway, end = 30),
                pathway = fct_reorder(pathway, prot_diff)) %>%
            pivot_longer(cols = matches("diff|padj"), values_to = "values", names_to = "var") %>%
            mutate(level = str_remove(var, "_diff$|_padj$"),
                   var = str_remove(var, "mrna_|prot_|cna_")) %>%
            pivot_wider(names_from = var, values_from = values, values_fn = mean) %>%
            mutate(mlog10 = -log10(padj)) 
        
        ggplot(paths2show, aes(level, pathway, size = mlog10, fill = diff)) +
            geom_point(pch = 21) +
            scale_fill_distiller(palette = "RdBu", limits = c(-5,5), oob = squish) +
            scale_size_area(limits = c(-1,6), oob = scales::squish) +
            labs(x = "Data level", y = "Pathway", subtitle = fact) +
            theme_mofa
    })
    path_plots <- wrap_plots(path_plots, ncol = 3) +
        plot_layout(guides = "collect")
    
    if(!is.null(base_name)) {
        ggsave(str_glue("results/mofa/{base_name}/pathways.pdf"), width = 10, height = 10)
    }
    
    return(path_plots)
}

gii_B_wilcox <- function(sampAnnot) {
    idx <- sampAnnot$uPG == "B"
    wilcox.test(sampAnnot$Genomic_Instability[idx], sampAnnot$Genomic_Instability[!idx])$p.value
}

################################################################################
## Statistical plots
################################################################################
plot_subtype_upg_assocs <- function(sa) {
    #-- Proportion of subtype per uPG
    sbty <- sa %>%
        select(Subtype, uPG = CCP_cluster)
    
    subtype_n <- sbty %>% group_by(Subtype) %>% count() %>% pull(n, Subtype)
    subtype_n_lvs <- structure(# names = paste0(names(subtype_n), " (n=", subtype_n, ")"),
                               names = names(subtype_n),
                               levels(sa$Subtype))
    sbty1 <- sbty %>%
        mutate(Subtype = fct_recode(Subtype, !!!subtype_n_lvs))
    
    prop <- 
        ggplot(sbty1, aes(Subtype, fill = uPG)) +
        geom_bar() +
        geom_hline(yintercept = 1:16, color = "white", size = 0.4) +
        scale_fill_manual(values = ccp_cols) +
        theme_csg_sparse +
        theme(panel.border = element_rect(fill = "transparent"),
              axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Number of samples") +
        guides(fill = "none")
    
    
    #-- FET subtype
    df4fet <- fastDummies::dummy_cols(sbty)
    colnames(df4fet) <- str_remove(colnames(df4fet), "Subtype_|uPG_")
    subtypes <- c("Ba/Sq", "NE-like", "LumU", "LumNS", "Class 2a", "Class 2b", "Stroma-rich", "LumP", "Class 1")
    upgs <- LETTERS[1:5] 
    
    #-- Signifcance (FET)
    fet_tb <- sapply(subtypes, function(st){
        sapply(upgs, function(upg) {
            fisher.test(makeContingencyTable(df4fet, st, upg), alternative = "greater")$p.value
        })
    }) %>%
        as.data.frame() %>%
        rownames_to_column("uPG") %>%
        pivot_longer(cols=-uPG, names_to = "Subtype", values_to = "FET_pval") %>%
        mutate(
            FET_padj = p.adjust(FET_pval, method = "BH"),
            signif = ifelse(FET_padj < 0.05, "*", "ns"),
            mlog10 = -log10(FET_padj),
            Subtype = factor(Subtype, levels = rev(subtypes)))
    
    fet_plt <- 
        ggplot(fet_tb, aes(uPG, Subtype, size = mlog10, color = signif, fill = mlog10)) +
        geom_point(pch = 21) +
        scale_color_manual(values = c("*" = "black", "ns" = "white")) +
        scale_fill_distiller(palette = "Reds", direction = 1) +
        scale_size(limits = c(-0.5,6)) +
        theme_csg_scatter 
    
    #-- Join
    plt <- (prop / fet_plt) + plot_layout(heights = c(1,1))
    ggsave("results/fig2/fig2c.pdf", height = 4, width = 3.3)
    
    return(plt)
}

plot_path_assocs <- function(exp_table, multiomics_hm) {
    samp_order <- get_colLevels(multiomics_hm$mofa_hm)
    tb4plot <- exp_table %>%
        mutate(Stage = case_when(
            stage == "Ta" ~ "Ta",
            str_detect(stage, "^T1") ~ "T1",
            TRUE ~ "MIBC"
        ),
        Node = ifelse(node == 1, "Yes", "No"),
        Squamous = ifelse(histological_subtype == "squamous differentiation", "Yes", "No"),
        Neuroendocrine = ifelse(histological_subtype == "neuroendocrine", "Yes", "No"),
        sample_id = factor(sample_id, levels = samp_order)) %>%
        select(sample_id, Stage, Grade = grade, 
               Node,
               Squamous, Neuroendocrine,
               Papillary = papillary,
               CIS = cis, 
               `Endophytic growth` = endophytic, 
               class_uPG)
    
    #-- Associations
    test_res <- doAssociationTests(tb4plot, id_var = "sample_id", test_var = "class_uPG")
    
    #-- Plot
    cp_hm <- tb4plot %>%
        select(-class_uPG) %>%
        pivot_longer(cols = -sample_id) %>%
        mutate(name = factor(name, levels = c("Stage", "Grade", "Node",
                                              "CIS", "Papillary", "Endophytic growth",
                                              "Squamous", "Neuroendocrine"))) %>%
        ggplot(aes(sample_id, fct_rev(name), fill = value)) +
        geom_tile() +
        scale_fill_manual(values = c("Ta" = "#fee0d2", "T1" = "#fc9272", "MIBC" = "#de2d26",
                                      "High" = "#e6550d", "Low" = "#fee6ce", yn_cols), 
                          na.value = "white") +
        labs(y = "Clinical/Pathological") +
        theme_csg_sparse +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 7),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank())
    
    ggsave("results/fig2/fig2b_clinPath.pdf", cp_hm, width = 4.5, height = 0.9)
    
    res <- list(clinpath_hmComplement = cp_hm,
               fisher_assocs = test_res)
    return(res)
    
    
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
## uPG pathways
################################################################################
plot_uPG_pathEnrich <- function(pathScores_uPG, pathScores_uPG_KEGG,
                                pathScores_uPG_H, exp_ic_table) {
    #-- List what to plot
    ics2plot <- c("IC Basal-like", "IC Cell cycle", "IC Neuroendocrine", 
                  "IC Urothelial differentiation", "IC Smooth muscle", "IC Ta pathway")
    
    paths2plot <- c("(HM) Interferon gamma" = "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
                    "(KEGG) Ribosome" = "KEGG_RIBOSOME",
                    "(HM) E2F targets" = "HALLMARK_E2F_TARGETS",
                    "(HM) DNA Repair" = "HALLMARK_DNA_REPAIR", 
                    "(HM) Fatty acid metabolism" = "HALLMARK_FATTY_ACID_METABOLISM",
                    "(KEGG) Peroxisome" = "KEGG_PEROXISOME", 
                    "(HM) Myogenesis" = "HALLMARK_MYOGENESIS", 
                    "(KEGG) Focal adhesion" = "KEGG_FOCAL_ADHESION", 
                    "(RDB) Cell-ECM interactions" = "Cell-extracellular matrix interactions", 
                    "(RDB) Respiratory electron transport" = "Respiratory electron transport")
    
    #-- Get IC plot table
    ic_tb <- exp_ic_table %>%
        mutate(feat = ifelse(feat == "IC Bladder cancer pathways", "IC Ta pathway", feat)) %>%
        select(IC = feat, starts_with("logFC"), starts_with("pair")) %>%
        pivot_longer(cols = starts_with("logFC"), names_to = "upg", values_to = "logFC") %>%
        pivot_longer(cols = starts_with("pair"), names_to = "upg2", values_to = "pval") %>%
        mutate(upg = str_remove(upg, "logFC_"),
               upg2 = str_remove(upg2, "pair_pval_")) %>%
        filter(IC %in% ics2plot, upg == upg2) %>%
        select(-upg2) %>%
        mutate(mlogpval = -log10(pval),
               signif = ifelse(pval < 0.05, "*", "ns"),
               IC = factor(IC, levels = rev(ics2plot)))
    
    #-- Get path plot tables and join
    rdb_paths <- .tidy_path_table(pathScores_uPG, paths2plot)
    kegg_paths <- .tidy_path_table(pathScores_uPG_KEGG, paths2plot)
    hallmarks <- .tidy_path_table(pathScores_uPG_H, paths2plot)
    
    paths_tb <- bind_rows(rdb_paths, kegg_paths, hallmarks) %>%
        mutate(Path_Name = plyr::mapvalues(Path_Name, paths2plot, names(paths2plot)),
               Path_Name = factor(Path_Name, levels = rev(names(paths2plot))),
               mlogpadj = -log10(padj),
               signif = ifelse(padj < 0.05, "*", "ns"))
    
    #-- IC plot
    ic_plt <- ggplot(ic_tb, aes(upg, IC, fill = logFC, color = signif, size = mlogpval)) +
        geom_point(pch = 21) +
        scale_fill_distiller(palette = "RdBu") +
        scale_color_manual(values = c("*" = "black", "ns" = "grey60")) +
        scale_size_area(limits = c(0, 8), max_size = 4, oob = scales::squish) +
        theme_csg_scatter +
        labs(y = "Transcriptomic IC") +
        guides(size = "none") +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank())
    
    #-- Path plot
    path_plt <- ggplot(paths_tb, aes(uPG, Path_Name, fill = NES, color = signif, size = mlogpadj)) +
        geom_point(pch = 21) +
        scale_fill_distiller(palette = "PuOr") +
        scale_color_manual(values = c("*" = "black", "ns" = "grey60")) +
        scale_size_area(limits = c(0, 8), max_size = 4, oob = scales::squish, breaks = c(0, 1.33, 8)) +
        theme_csg_scatter +
        labs(x = "uPG", y = "Pathway")
    
    panel_plt <- ( ic_plt / path_plt ) + plot_layout(guides = "collect", heights = c(length(ics2plot),
                                                                        length(paths2plot)))
    
    ggsave("results/fig2/figb_2.pdf", panel_plt, width = 4, height = 3.5)
    
}

.tidy_path_table <- function(path_table, paths2plot) {
    lapply(names(path_table), function(upg) {
        path_table[[upg]] %>%
            select(Path_Name, padj, NES) %>%
            mutate(uPG = upg)
    }) %>%
        bind_rows() %>%
        filter(Path_Name %in% paths2plot)
}

subtype_ccp_flow <- function(samp_annot_63) {
  plt_flow <- samp_annot_63 %>%
    rename(uPG = CCP_cluster) %>%
    mutate(Subtype = factor(Subtype, levels = c("Ba/Sq", "LumU", "NE-like", "LumNS", 
                                        "Class 2a", "Class 2b", "Stroma-rich",
                                        "LumP", "Class 1"))) %>%
    classFlow("Subtype", "uPG", 
              id_var = "id", alpha = 0.5, 
              label_colors = c(subtype_cols, ccp_cols),
              label_size = 2) +
      guides(fill = "none") +
      theme(axis.text = element_text(size = 8),
            axis.title = element_text(size = 8, face = "bold"),
            axis.title.y = element_blank(),
            axis.ticks = element_line(color = "black"), 
            axis.line.y = element_blank()) +
    labs(x = "Class type")
  
  ggsave("results/fig2/fig2_sankey.pdf", plot = plt_flow, width = 2.5, height = 2.5)
  
}

get_subset_data <- function(prot_data, stage = "NMIBC") {
    ids <- prot_data$sampAnnot$id[prot_data$sampAnnot$Stage. == stage]
    new_prot_data <- prot_data
    new_prot_data$wp <- new_prot_data$wp[,ids]
    new_prot_data$pep <- new_prot_data$pep[,ids]
    new_prot_data$sampAnnot <- new_prot_data$sampAnnot[ids,]
    return(new_prot_data)
}

