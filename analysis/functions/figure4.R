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
    
    ggsave(filename = "results/fig4/fig4b.pdf", plot = hm, width = 8, height = 3)
    
    return(hm)
}
################################################################################
## siFGFR3
################################################################################
source("aux_functions/pathway_enrich_funs.R")
source("aux_functions/pathenrich_plots2.R")
dir.create("results/tables", showWarnings = FALSE)

siFGFR3_plot <- function(
    paths2plot = c("(H) GLYCOLYSIS", "(H) OXIDATIVE PHOSPHORYLATION", "(RDB) APOPTOTIC EXECUTION PHASE",
                   "(RDB) RESPIRATORY ELECTRON TRANSPORT", "(RDB) APOPTOSIS"),
    # paths2plot =  c("(RDB) INNATE IMMUNE SYSTEM", "(RDB) INTERFERON SIGNALING", 
    #                                      "(RDB) CELL CYCLE", "(RDB) FOXO MEDIATED TRANSCRIPTION", 
    #                                      "(H) MYC TARGETS V2", "(H) GLYCOLYSIS", 
    #                                      "(H) OXIDATIVE PHOSPHORYLATION", "(H) MTORC1 SIGNALING"),
    pcutoff = 0.05, width = 10, height = 5) {
    
    #-- Read data
    clines <- readxl::excel_sheets("data/experiments/siFGFR3/200416-MGHU3-RT112-UMUC14-siF3-Log2FC-for-gseaPreRanked.xlsx")
    sifgfr3_fc <- lapply(clines, function(i) {
        read_excel("data/experiments/siFGFR3/200416-MGHU3-RT112-UMUC14-siF3-Log2FC-for-gseaPreRanked.xlsx", 
                   sheet = i) })
    clines <- str_remove(clines, " siFGFR3")
    names(sifgfr3_fc) <- clines
    
    #-- Conform to format for fgsea (named, sorted vector of lfc)
    sifgfr3_fc <- lapply(sifgfr3_fc, function(sitable) {
        sitable <- sitable %>%
            mutate(log2FC = as.numeric(log2FC),
                   abs_lfc = abs(log2FC)) %>%
            group_by(hgnc_symbol) %>%
            top_n(1, abs_lfc)
        lfc_vec <- structure(sitable$log2FC, names = sitable$hgnc_symbol)
        sort(lfc_vec, decreasing = TRUE)
    })
    
    ################################################################################
    ## Run enrichment
    ################################################################################
    #-- Reactome collection
    reactome_paths <- get_pathway_list("C2", "CP:REACTOME", "^REACTOME_")
    reactome_enrich <- lapply(sifgfr3_fc, run_fgsea, reactome_paths$list, reactome_paths$annot)
    #-------------------------------------------------------------------------------
    # Cross Reactome
    path_ids <- sapply(reactome_enrich, function(tb) {
        tb %>%
            filter(padj < 0.05) %>%
            pull(pathway_id)
    }) %>% reduce(intersect)
    
    intersect_negenrich <- lapply(names(reactome_enrich), function(cl) {
        tb <- reactome_enrich[[cl]]
        tb <- tb %>%
            filter(pathway_id %in% path_ids) %>%
            select(pathway_id, pathway_name, padj, NES)
        colnames(tb)[3:4] <- paste0(cl, "_", c("padj", "NES"))
        return(tb)
    }) %>% reduce(full_join)
    
    #-------------------------------------------------------------------------------
    #-- Hallmark collection
    hallmarks <- get_pathway_list("H", "", "HALLMARK_")
    hallmark_enrich <- lapply(sifgfr3_fc, run_fgsea, hallmarks$list, hallmarks$annot)
    
    #-- Cross halmarks
    path_ids <- sapply(hallmark_enrich, function(tb) {
        tb %>%
            filter(padj < 0.05) %>%
            pull(pathway_id)
    }) %>% reduce(intersect)
    
    intersect_negenrich <- lapply(names(hallmark_enrich), function(cl) {
        tb <- hallmark_enrich[[cl]]
        tb <- tb %>%
            filter(pathway_id %in% path_ids) %>%
            select(pathway_name, padj, NES)
        colnames(tb)[2:3] <- paste0(cl, "_", c("padj", "NES"))
        return(tb)
    }) %>% reduce(full_join)
    
    ################################################################################
    ## Select results
    ################################################################################
    reactome_enrich <- lapply(reactome_enrich, function(tb) {
        tb %>% 
            mutate(pathway = paste0("(RDB)_", pathway_name))
    })
    hallmark_enrich <- lapply(hallmark_enrich, function(tb) {
        tb %>% 
            mutate(pathway = paste0("(H)_", pathway_name))
    })
    
    #-- Join
    joined_enrich <- lapply(names(reactome_enrich), function(cl) {
        bind_rows(reactome_enrich[[cl]], hallmark_enrich[[cl]])
    })
    names(joined_enrich) <- names(reactome_enrich)
    
    enrich_plot <- .FGFR3enrich_plot(joined_enrich, paths2plot, clines)
    
    ggsave(filename = "results/fig4/fig4a.pdf", 
           plot = enrich_plot, width = 5, height = 4)
    
    return(list(plot = enrich_plot, all_enrich = joined_enrich))
}

.FGFR3enrich_plot <- function(joined_enrich, paths2plot = paths2plot, clines) {
    venn <- as_grob(.venn_pathEnrich(joined_enrich, direction = "down", 
                                     pcutoff = 0.05, colors = c("#4393c3", "#92c5de", "#d1e5f0")))
    path_scores <- lapply(clines, function(cl) {
        joined_enrich[[cl]] %>%
            select(pathway, NES, padj) %>%
            rename_with(~ paste0(., "_", cl), NES:padj)
    }) %>% reduce(inner_join, by = "pathway")
    
    downreg4plot <- path_scores %>%
        mutate(pathway = .fix_path_names(pathway),
               NES_median = matrixStats::rowMedians(as.matrix(select(path_scores, starts_with("NES")))),
               NES_ci = rowSds(as.matrix(select(path_scores, starts_with("NES"))))*2,
               NES_high = NES_median + NES_ci,
               NES_low = NES_median - NES_ci,
               padj_median = rowMedians(as.matrix(select(path_scores, starts_with("padj")))),
               log10padj = -log10(padj_median)) %>%
        filter(pathway %in% paths2plot)
    
    dot_plot <- 
        ggplot(downreg4plot, aes(NES_median, fct_reorder(pathway, NES_median, .desc = TRUE))) +
        geom_errorbar(aes(xmin = NES_low, xmax = NES_high), width = .3) +
        geom_point(aes(size = log10padj, fill = NES_median), color = "black", pch = 21) +
        geom_vline(xintercept = 0, lty = "dashed") +
        scale_fill_gradient(low = "#084594", high = "#c6dbef", 
                            limits = c(-2,0), breaks = c(-2,-1,0)) +
        labs(x = "Normalized enrichment score\n(NES)", y = "Pathways", color = "Median NES", 
             size = expression(-log[10](p[adj]))) +
        scale_radius(limits = c(-1,5), range = c(2,6),
                     breaks = c(0,1.3,2,5)) +
        theme_csg_scatter +
        labs(fill = "Median NES") +
        theme(panel.border = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_line(linetype = "dotted", color = "grey30"),
              legend.position = "bottom")
    
    grid <- plot_grid(venn, dot_plot, ncol = 1, rel_heights = c(1,1.5))
    return(grid)
}

export_siFGFR3_table <- function(siFGFR3) {
    all_enrich <- siFGFR3$all_enrich
    sifgfr3_signifEnrich <- lapply(names(all_enrich), function(cl) { 
        mutate(all_enrich[[cl]], cell_line = cl)}) %>%
        bind_rows() %>%
        group_by(pathway) %>%
        filter(all(padj < 0.05)) %>%
        pivot_wider(id_cols = pathway,
                    names_from = cell_line,
                    values_from = c(NES, padj))
    
    openxlsx::write.xlsx(sifgfr3_signifEnrich, 
                         file = "results/tables/siFGFR3_pathways.xlsx", 
                         overwrite = TRUE)
    return(sifgfr3_signifEnrich)
}
