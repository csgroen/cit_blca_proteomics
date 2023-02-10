library(patchwork)
library(ggrepel)

plt_upg_volcanos <- function(diffProt_uPG, react_paths, hallmark_paths, kegg_paths) {
    prot_highlights <- list(
        "A" = c("KRT6A", "KRT14", "S100A8", "KRT5", "KRT16"),
        "B" = c("TOP2A", "TMPO", "FEN1", "MCM5", "MCM6", "KNPA2", "MCM2", "PRIM1", "STAG2"),
        "C" = c("AKR1C3", "GDPD3", "OSBP", "ACAA2"),
        "D" = c("FLNC", "TPM2", "TAGLN", "TPM1", "CNN1"),
        "E" = c("ANXA10", "S100P", "IVL", "KRT7", "DHIRS2", "AGR2")
    )
    #-- Get path highlights
    all_path_tars <- bind_rows(select(react_paths, Path_Name, Symbol), 
                               select(hallmark_paths, Path_Name, Symbol),
                               select(kegg_paths, Path_Name, Symbol))
    path_highlights <- list(
        "A" = c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                "KEGG_RIBOSOME"),
        "B" = c("HALLMARK_E2F_TARGETS", "HALLMARK_DNA_REPAIR"),
        "C" = c("HALLMARK_FATTY_ACID_METABOLISM", "KEGG_PEROXISOME"),
        "D" = c("HALLMARK_MYOGENESIS", "KEGG_FOCAL_ADHESION", "Cell-extracellular matrix interactions"),
        "E" = c("Respiratory electron transport")
    )
    
    path_tars_highlights <- lapply(LETTERS[1:5], function(upg) {
        paths <- path_highlights[[upg]]
        path_symbols <- all_path_tars %>% filter(Path_Name %in% paths) %>% pull(Symbol) %>% unique()
    })
    names(path_tars_highlights) <- LETTERS[1:5]
    all_highlights <- lapply(LETTERS[1:5], function(upg) {
        union(path_tars_highlights[[upg]], prot_highlights[[upg]])
    })
    names(all_highlights) <- LETTERS[1:5]

    #-- Plot
    upg_volcanoplots <- lapply(LETTERS[1:5], .upg_volcano, diffProt_uPG, all_highlights)
    upg_volcanoplots[[2]] <- upg_volcanoplots[[2]] + lims(y = c(0, 4.5))
    upg_volcanoplots[[3]] <- upg_volcanoplots[[3]] + lims(y = c(0, 2))
    all_volcanos <- wrap_plots(upg_volcanoplots, ncol = 2)
    all_volcanos
    ggsave("results/suppfig_volcanos/a.pdf", all_volcanos, width = 7, height = 8)
    
}

.upg_volcano <- function(subtype, diffProt_uPG, all_highlights) {
    highlight <- all_highlights[[subtype]]
    pval_var <- paste0("pair_pval_", subtype)
    logfc_var <- paste0("logFC_", subtype)
    diffprot <- diffProt_uPG %>%
        select(symbol, pval = !! pval_var, log2fc = !! logfc_var) %>%
        mutate(mlog10pval = -log10(pval)) %>%
        mutate(
            type = case_when(
                log2fc > 0 & pval < 0.05 ~ "up",
                log2fc < 0 & pval < 0.05 ~ "down",
                TRUE ~ "neither"),
            type2 = case_when(
                log2fc > 0 & pval < 0.05 ~ "up",
                log2fc < 0 & pval < 0.05 ~ "down",
                TRUE ~ NA_character_),
            abs_logfc = abs(log2fc))
    
    to_label_fc <- diffprot %>% 
        group_by(type) %>%
        group_by(type2) %>%
        slice_max(abs_logfc, n = 8) %>%
        pull(symbol)
    
    to_label_pval <- diffprot %>% 
        group_by(type) %>%
        ungroup() %>%
        slice_max(mlog10pval, n = 10) %>%
        slice(1:10) %>%
        pull(symbol)
    
    to_label_auto <- union(to_label_fc, to_label_pval)[1:15]
    to_label <- union(to_label_auto, to_label_pval)
    
    diffprot %>%
        mutate(label = ifelse(symbol %in% to_label, symbol, NA),
               highlight = ifelse(label %in% highlight, "highlight", "no")) %>%
        ggplot(aes(log2fc, mlog10pval, label = label)) +
        geom_point(aes(color = type), size = 0.8) +
        geom_text_repel(aes(color = highlight), size = 2, max.overlaps = Inf) + 
        geom_hline(yintercept = 1.33, lty = "dashed") +
        geom_vline(xintercept = 0, lty = "dashed") +
        scale_color_manual(values = c("up" = "#ef8a62", "down" = "#67a9cf", "neither" = "grey80",
                                      "highlight" = "#ca0020", "no" = "black")) +
        # scale_color_manual(values = c("highlight" = "#ca0020", "no" = "black")) +
        labs(x = "log2(FC)", y = "-log10(p-value)", title = paste0("uPG-", subtype)) +
        guides(color = "none", fill = "none") +
        theme_csg_scatter +
        theme(title = element_text(face = "bold", size = 8))
}

################################################################################
## P-value heatmap MCPcounter
################################################################################
plt_mcpCounter_pvals <- function(samp_annot_63, mcp_res) {
    mcp_res_tb <- mcp_res %>%
        as.data.frame() %>%
        rownames_to_column("id") %>%
        rowwise() %>%
        ungroup()
    
    tab4assoc <- mcp_res_tb %>%
        left_join(select(samp_annot_63, id, uPG = CCP_cluster))
    
    #-- Tests
    pval_mat_higher <- sapply(LETTERS[1:5], function(cls) {
        sapply(colnames(mcp_res), function(var) {
            #-------- Get in/out class
            cls_values <- tab4assoc %>% pull(uPG, id)
            in_cls <- names(cls_values)[cls_values == cls]
            out_cls <- names(cls_values)[cls_values != cls]
            
            #-------- Get values
            var_values <- tab4assoc %>% pull(!!{{var}}, id)
            #-------- Test
            wilcox.test(var_values[in_cls], var_values[out_cls], alternative = "greater")$p.value/2
        })
    })
    pval_mat_less <- sapply(LETTERS[1:5], function(cls) {
        sapply(colnames(mcp_res), function(var) {
            #-------- Get in/out class
            cls_values <- tab4assoc %>% pull(uPG, id)
            in_cls <- names(cls_values)[cls_values == cls]
            out_cls <- names(cls_values)[cls_values != cls]
            
            #-------- Get values
            var_values <- tab4assoc %>% pull(!!{{var}}, id)
            #-------- Test
            wilcox.test(var_values[in_cls], var_values[out_cls], alternative = "less")$p.value/2
        })
    })
    idx <- pval_mat_less < pval_mat_higher
    
    mlogpval_mat <- -log10(pval_mat_higher)
    mlogpval_mat[idx] <- log10(pval_mat_less)[idx]
    mcp_pvals <- ggheatmap(mlogpval_mat,
              cluster_rows = TRUE,
              cluster_cols = FALSE,
              hm_color_values = c(-6,-2,-1,0,1,2,6),
              hm_color_breaks = c(-6,-2,0,2,6),
              hm_colors = "RdBu",
              colors_title = "-log10(p-value)")

    
    # mcp_box <- mcp_res_tb %>%
    #     pivot_longer(cols = -id, names_to = "Population", values_to = "MCPcounter Score") %>%
    #     left_join(select(samp_annot_63, id, uPG = CCP_cluster)) %>%
    #     ggplot(aes(uPG, `MCPcounter Score`, fill = uPG)) +
    #     facet_wrap(~ Population, scales = "free_y") +
    #     geom_violin(alpha = 0.2) +
    #     geom_boxplot(width = 0.2) +
    #     stat_compare_means(size = 2.5) +
    #     theme_csg_scatter
    # ggsave("results/suppfig_volcanos/b.pdf", mcp_box, width = 7.5, height = 4)
    ggsave("results/suppfig_volcanos/b.pdf", mcp_pvals, width = 4.5, height = 2.5)
} 

################################################################################
## uPG - IC differences
################################################################################
upg_ic_table <- function(mrna_data, samp_annot_63) {
    #-- Get IC table
    ic_table <- mrna_data$sampleAnnot %>%
        select(id, starts_with("IC")) %>%
        right_join(select(samp_annot_63, id, uPG = CCP_cluster)) %>%
        select(-`IC Smooth muscle`) %>%
        rename(`IC Smooth muscle` = ICA.smoothMuscle)
    
    ic4hm <- ic_table %>%
        group_by(uPG) %>%
        summarize(across(starts_with("IC"), mean, na.rm = TRUE))

    #-- Differential test
    ic_mat <- ic_table %>% select(-uPG) %>% as.data.frame() %>% column_to_rownames("id") %>% as.matrix() %>% t()
    cls <- structure(ic_table$uPG, names = ic_table$id)
    
    diff_ics <- findClassDiffs(ic_mat, cl = cls, type = "non-parametric")
    return(diff_ics)
}

plt_uPG_PDL1 <- function(prot_data, unfilt_prot, mrna_data, samp_annot_63) {
    #-- Check if PDL1 or PD1 proteins were identified
    prot_data$wpAnnot %>%
        tibble() %>%
        filter(symbol %in% c("CD274", "PDCD1"))
    
    unfilt_prot$wpAnnot %>%
        tibble() %>%
        filter(symbol %in% c("CD274", "PDCD1"))
    
    #-- PDL1 prot quantification + PDL1 expression
    pdl1_quant <- ifelse(!is.na(unfilt_prot$wp["Q9NZQ7",]), "PD-L1+", "PD-L1-")
    pdl1_exp <- mrna_data$gexp["CD274",]
    pdl1_exp <- tibble(PDL1_exp = pdl1_exp, id = names(pdl1_exp))
    
    sa <- samp_annot_63 %>%
        mutate(PDL1_quant = pdl1_quant) %>%
        left_join(pdl1_exp)
    
    #-- Get summaries + stats
    table(sa$PDL1_quant, sa$CCP_cluster)
    sa <- sa %>%
        mutate(CCP_A = ifelse(CCP_cluster == "A", "A", "Other"))
    
    fisher.test(makeContingencyTable(sa, "CCP_A", "PDL1_quant"))$p.value
    
    pdl1_a_d <- sa %>%
        filter(Consensus == "Ba/Sq", CCP_cluster %in% c("A", "D")) %>%
        ggplot(aes(CCP_cluster, PDL1_exp, fill = CCP_cluster)) +
        geom_violin(alpha = 0.4) +
        geom_boxplot(width = 0.3) +
        geom_quasirandom() +
        stat_compare_means(size = 3) +
        scale_fill_manual(values = ccp_cols) +
        guides(fill = "none") +
        labs(x = "uPG", y = "PD-L1 gene expression") +
        theme_csg_scatter
    
    ggsave("results/suppfig_volcanos/extra_baAD.pdf", pdl1_a_d, width = 2, height = 1.5)
    
    
    #-- Plot
    plt_pdl1_prot <- sa %>%
        group_by(PDL1_quant, CCP_cluster) %>%
        summarize(n = n()) %>%
        ggplot(aes(fct_rev(CCP_cluster), n, fill = PDL1_quant)) +
        geom_col(color = "black", size = 0.3, width = 0.5) +
        annotate("text", x = 5, y = 16, label = "***", angle = 90) +
        scale_fill_manual(values = c("PD-L1+" = "#de2d26", "PD-L1-" = "#fee0d2")) +
        scale_x_discrete(position = "top") +
        labs(y = "Number of samples", fill = "PD-L1 protein\n quantification") +
        coord_flip() +
        theme_csg_hist +
        theme(axis.line.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.major.x = element_line(color = "grey90", linetype = "dotted"),
              panel.grid.major.y = element_line(color = "grey90", linetype = "solid", size = 0.5))
    plt_pdl1_prot
    plt_pdl1_exp <- ggplot(sa, aes(fct_rev(CCP_cluster), PDL1_exp, fill = CCP_cluster)) +
        geom_violin(alpha = 0.3) +
        geom_boxplot(width = 0.2) +
        geom_quasirandom(size = 0.5) +
        stat_compare_means(comparisons = list(c("A", "D")),
                           label = "p.signif") +
        stat_compare_means(label.y = 5.4, label.x = 0.5, size = 3) +
        labs(x = "uPG", y = "PD-L1 (CD274) gene expression") +
        guides(fill = "none") +
        coord_flip() +
        theme_csg_scatter
    plt_pdl1_exp
    
    #-- Join plots
    plt_pdl1 <- (plt_pdl1_exp | plt_pdl1_prot) +
        plot_layout(widths = c(1,0.3), guides = "collect")
    
    ggsave("results/suppfig_volcanos/c.pdf", plt_pdl1, width = 6, height = 2.3)
    return(plt_pdl1)
}


################################################################################
## Boxplots MOFA
################################################################################
# plt_mofa_boxplots <- function(samp_annot_63, final_mofa) {
#     mofa_trained <- final_mofa$mofa
#     factors2show <- final_mofa$facts
#     factor_names <- final_mofa$fact_names
#     invert <- final_mofa$invert
#     
#     #-- Get factor HM
#     fct_values <- MOFA2::get_factors(mofa_trained)[[1]]
#     fct_mat <- sapply(1:length(factors2show), function(i) {
#         fct <- factors2show[i]; inv <- invert[i]
#         vals <- fct_values[,fct] * inv
#         return(vals)
#     })
#     colnames(fct_mat) <- factor_names
#     
#     fct_box <- fct_mat %>%
#         as.data.frame() %>%
#         rownames_to_column("id") %>%
#         pivot_longer(cols = -id, names_to = "Factor", values_to = "MOFA Score") %>%
#         left_join(select(samp_annot_63, id, uPG = CCP_cluster)) %>%
#         ggplot(aes(uPG, `MOFA Score`, fill = uPG)) +
#         facet_wrap(~ Factor, scales = "free_y",  ncol = 5) +
#         geom_violin(alpha = 0.2) +
#         geom_boxplot(width = 0.2) +
#         stat_compare_means(size = 2) +
#         theme_csg_scatter
#     fct_box
#     ggsave("results/suppfig_volcanos/c.pdf", fct_box, width = 8, height = 1.5)
# } 

