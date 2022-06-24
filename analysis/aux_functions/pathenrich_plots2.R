################################################################################
## Pathway enrich plots
################################################################################
#-- Venn diagram
.venn_pathEnrich <- function(path_enrich_list, direction = "up", 
                            colors = NULL, pcutoff = pcutoff, ...) {
    if(direction == "up") {
        reg_paths <- lapply(path_enrich_list, function(enrich_tab) {
            enrich_tab %>%
                filter(padj < pcutoff, NES > 0) %>%
                pull(pathway)
        }) 
    } else {
        reg_paths <- lapply(path_enrich_list, function(enrich_tab) {
            enrich_tab %>%
                filter(padj < pcutoff, NES < 0) %>%
                pull(pathway)
        })
    }
    euler_df <- data.frame(pathway = reduce(reg_paths, union))
    euler_df <- euler_df %>%
        mutate(map(reg_paths, function(paths) { euler_df$pathway %in% paths }) %>%
                   as.data.frame())
    plot(venn(select(euler_df, -pathway)), 
         quantities = TRUE, fills = colors, ... = ...)
    # return(list(plt, euler_df[apply(euler_df[,-1], 1, sum) == ncol(euler_df)-1,1]))
}

.dotplot_pathEnrich <- function(enrich_results, clines, direction = 1, ntop = 30, paths2plot = NULL) {
    if(direction > 0) {
        upreg_scores <- lapply(clines, function(cl) {
            enrich_results[[cl]] %>%
                filter(NES > 0) %>%
                select(pathway, NES, padj) %>%
                rename_with(~ paste0(., "_", cl), NES:padj)
        }) %>% reduce(inner_join, by = "pathway")
        
        upreg4plot <- upreg_scores %>%
            mutate(pathway = .fix_path_names(pathway),
                   NES_median = matrixStats::rowMedians(as.matrix(select(upreg_scores, starts_with("NES")))),
                   NES_ci = rowSds(as.matrix(select(upreg_scores, starts_with("NES"))))*2,
                   NES_high = NES_median + NES_ci,
                   NES_low = NES_median - NES_ci,
                   padj_median = rowMedians(as.matrix(select(upreg_scores, starts_with("padj")))),
                   log10padj = -log10(padj_median))
        if(is.null(ntop)) {
            upreg4plot <- upreg4plot %>% filter(pathway %in% paths2plot)
            
        } else {
            upreg4plot <- upreg4plot %>% top_n(ntop, log10padj)
        }
        
        plot <- ggplot(upreg4plot, aes(NES_median, fct_reorder(pathway, NES_median))) +
            geom_errorbar(aes(xmin = NES_low, xmax = NES_high), width = .3) +
            geom_point(aes(size = log10padj, color = NES_median)) +
            scale_color_gradient(low = "#fcbba1", high = "#cb181d") +
            labs(x = "Normalized enrichment score\n(NES)", y = "Pathways", color = "Median NES", size = expression(-log[10](p[adj]))) +
            lims(x = c(0, max(upreg4plot$NES_high))) +
            theme_light() +
            theme(panel.border = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.grid.major.x = element_line(linetype = "dotted", color = "grey30"),
                  axis.text = element_text(color = "black"),
                  legend.position = "bottom")
    } else {
        downreg_scores <- lapply(clines, function(cl) {
            enrich_results[[cl]] %>%
                filter(NES < 0) %>%
                select(pathway, NES, padj) %>%
                rename_with(~ paste0(., "_", cl), NES:padj)
        }) %>% reduce(inner_join, by = "pathway")
        
        downreg4plot <- downreg_scores %>%
            mutate(pathway = .fix_path_names(pathway),
                   NES_median = rowMedians(as.matrix(select(downreg_scores, starts_with("NES")))),
                   NES_ci = rowSds(as.matrix(select(downreg_scores, starts_with("NES"))))*2,
                   NES_high = NES_median + NES_ci,
                   NES_low = NES_median - NES_ci,
                   padj_median = rowMedians(as.matrix(select(downreg_scores, starts_with("padj")))),
                   log10padj = -log10(padj_median))
        
        if(is.null(ntop)) {
            downreg4plot <- downreg4plot %>% filter(pathway %in% paths2plot)
            
        } else {
            downreg4plot <- downreg4plot %>% top_n(ntop, log10padj)
        }
        
        plot <- ggplot(downreg4plot, aes(NES_median, fct_reorder(pathway, NES_median, .desc = TRUE))) +
            geom_errorbar(aes(xmin = NES_low, xmax = NES_high), width = .3) +
            geom_point(aes(size = log10padj, color = NES_median)) +
            scale_color_gradient(low = "#084594", high = "#c6dbef") +
            labs(x = "Normalized enrichment score\n(NES)", y = "Pathways", color = "Median NES", size = expression(-log[10](p[adj]))) +
            lims(x = c(min(downreg4plot$NES_low), 0)) +
            theme_light() +
            theme(panel.border = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.grid.major.x = element_line(linetype = "dotted", color = "grey30"),
                  axis.text = element_text(color = "black"),
                  legend.position = "bottom")
    }
    return(plot)
}

.fix_path_names <- function(pathway_names) {
    #-- Get gene names for auto-replacement
    genes_upper <- names(as.list(org.Hs.eg.db::org.Hs.egSYMBOL2EG)) %>% unique()
    genes_title <- paste0("\b", str_to_title(genes_upper), "\b")
    #-- Make patterns and replacements
    patterns <- c("Dna", "Rna", "Mrna", "Trna", "Rrna", "Gtp", "Atp", "G2m", "G1s",
                  genes_title)
    replacements <- c("DNA", "RNA", "mRNA", "tRNA", "rRNA", "GTP", "ATP", "G2M", "G1S",
                      genes_upper)
    #-- Fix names
    pathway_names <- str_replace_all(pathway_names, "_", " ") %>% str_trim()
    # str_to_title() %>% str_trim()
    # pathway_names <- reduce2(patterns, replacements, .init = pathway_names, str_replace_all)
    # 
    return(pathway_names)
}

siFGFR3_plot_enrichmentPanel <- function(enrich_res, paths2plot = NULL, clines = NA, updown_rel_heights = c(1,1), 
                                 venn_dot_rel_widths = c(1,3), ntop = 15) {
    #-- Get Venns
    venn1 <- as_grob(.venn_pathEnrich(enrich_res, direction = "up",
                                     pcutoff = 0.05, colors = c("#d6604d", "#f4a582", "#fddbc7")))
    venn2 <- as_grob(.venn_pathEnrich(enrich_res, direction = "down", 
                                     pcutoff = 0.05, colors = c("#4393c3", "#92c5de", "#d1e5f0")))
    #-- Get dotplots
    if(is.null(ntop)) {
        dot1 <- .dotplot_pathEnrich(enrich_res, clines, paths2plot = paths2plot, ntop = NULL)
        dot2 <- .dotplot_pathEnrich(enrich_res, clines, direction = -1,  paths2plot = paths2plot, ntop = NULL)
    } else {
        dot1 <- .dotplot_pathEnrich(enrich_res, clines, ntop = ntop)
        dot2 <- .dotplot_pathEnrich(enrich_res, clines, direction = -1,  ntop = ntop)
    }
    
    padjs <- c(dot1$data$log10padj, dot2$data$log10padj)
    dot1 <- dot1 + scale_size_continuous(breaks = pretty(padjs, 6)[1:6], 
                                         limits = c(min(padjs), max(padjs))) +
        guides(size = FALSE)
    dot2 <- dot2 + scale_size_continuous(breaks = pretty(padjs, 6)[1:6], 
                                         limits = c(min(padjs), max(padjs)))
    
    #-- Venn panel
    left <- plot_grid(venn1, venn2, NULL, scale = 0.7, ncol = 1, rel_heights = c(2,2,1),
                      labels = c("Up-regulated", "Down-regulated"), label_size = 10)
    #-- Dotplot panel
    paths1 <- plot_grid(dot1 + theme(axis.title.x = element_blank()),
                        dot2, ncol = 1, 
                        rel_heights =  updown_rel_heights, align = "v")
    
    #-- Join venn
    panel1 <- plot_grid(left, paths1, nrow = 1, rel_widths = venn_dot_rel_widths)
    return (panel1)
}
