plot_physChemHist <- function(fname = "results/suppfig_bioinfoqc/c.pdf", save = TRUE,
                              width = 8, height = 2.5) {
    load("data/annotations/Proteomics_SubCell_PhysChem_info.RData")
    
    physchem4plot <- prot_info %>%
        mutate(filt = ifelse(Filtered, "Filtered", "Unfiltered"),
               nobs = ifelse(Filtered, 2, 1),
               logMW = log2(MW)) %>%
        uncount(nobs) %>%
        mutate(filt = ifelse(duplicated(PID), "Unfiltered", filt)) %>%
        select(PID, Scaled_Solubility, Avg_pI, logMW, filt) %>%
        rename(`Scaled solubility` = Scaled_Solubility, `Average pI` = Avg_pI, `log2(Molecular weight)` = logMW) %>%
        pivot_longer(cols = `Scaled solubility`:`log2(Molecular weight)`)
    
    
    physchem_plot <- ggplot(physchem4plot, aes(value, color = fct_rev(filt), fill = fct_rev(filt))) +
        facet_wrap(~ name, scales = "free_x") +
        geom_histogram(alpha = 0.5, bins = 30, position = "identity") +
        theme_csg_hist +
        scale_color_manual(values = c("Filtered" = "#045a8d", "Unfiltered" = "grey20")) +
        scale_fill_manual(values = c("Filtered" = "#74a9cf", "Unfiltered" = "grey70")) +
        labs(x = "", y = "Protein count", color = "Filter status", fill = "Filter status") +
        theme_csg_sparse +
        theme(panel.border = element_rect(fill = NA, size = 0.8))
    
    if(save) ggsave(filename = fname, plot = physchem_plot, width = width, height = height)
    return(physchem_plot)
}
# Supp Fig 1B
prot_function_dynamicRange <- function(unfilt_prot, 
                                       fname = "results/suppfig_bioinfoqc/d.pdf") {
    #-- Get "big terms" for molecular function
    library(msigdbr)
    go_cc <- msigdbr(category = "C5", subcategory = "CC")
    go_mf <- msigdbr(category = "C5", subcategory = "MF")
    
    go_mf$gs_name %>% unique() %>% view()
    
    bigterms1 <- go_cc %>% 
        mutate(big_term = case_when(
            gs_name == "GOCC_RIBOSOME" ~ "Ribosomal",
            gs_name == "GOCC_MITOCHONDRION" ~ "Mitochondrial",
            str_detect(gs_name, "CYTOSKELETON") ~ "Cytoskeleton",
            TRUE ~ NA_character_
        )) %>%
        filter(!is.na(big_term)) %>%
        relocate(big_term) %>%
        select(-gs_cat, -gs_subcat) %>%
        arrange(big_term)
    
    bigterms2 <- go_mf %>% 
        mutate(big_term = case_when(
            str_detect(gs_name, "TRANSCRIPTION_FACTOR") ~ "TF",
            str_detect(gs_name, "TYROSINE_KINASE") ~ "Tyrosine Kinase",
            str_detect(gs_name, "PHOSPHATASE") ~ "Phosphatase",
            str_detect(gs_name, "RECEPTOR") ~ "Receptor",
            TRUE ~ NA_character_
        )) %>%
        filter(!is.na(big_term)) %>%
        relocate(big_term) %>%
        select(-gs_cat, -gs_subcat) %>%
        arrange(big_term)
    
    bigterms <- bind_rows(bigterms1, bigterms2)
    
    #-- Query unfiltered proteome
    wp_unf <- unfilt_prot$wp
    
    bigcc_protinfo <- bigterms %>%
        left_join(select(unfilt_prot$wpAnnot, gene_symbol = symbol, protein_id)) %>%
        left_join(tibble(protein_id = rownames(wp_unf),
                         samp_prop = rowSums(wp_unf > 0, na.rm = TRUE))) %>%
        mutate(samp_prop = ifelse(is.na(samp_prop), 0, samp_prop))
    
    #-- Get proportions of term / samples quantified
    all_terms <- bigcc_protinfo$big_term %>% unique()
    dr4plot <- lapply(all_terms, function(term) {
        term_sp <- bigcc_protinfo %>%
            filter(big_term == term) %>%
            arrange(desc(samp_prop)) %>%
            pull(samp_prop)
        
        tibble(term = term, 
               quantile = 1 - seq(0,1,0.01),
               q_sp = quantile(term_sp, seq(0,1,0.01)))
    }) %>% 
        bind_rows() %>%
        mutate(term = factor(term, levels = c("Receptor", "TF", "Tyrosine Kinase",
                                              "Cytoskeleton", "Phosphatase",
                                              "Mitochondrial",  "Ribosomal"))) %>%
        filter(!is.na(term))
    #-- Join into plot
    dr_plot <- ggplot(dr4plot, aes(quantile, term, fill = q_sp)) +
        geom_raster() +
        scale_fill_viridis() +
        scale_x_continuous(labels = scales::percent) +
        labs(y = "Protein function", x = "Fraction of gene set quantified",
             fill = "Quantified in\nnumber of samples") +
        theme_csg_sparse +
        theme(legend.position="bottom")
    
    ggsave(dr_plot, file = fname,
           width = 3, height = 2.5)
    
    return(dr_plot)
}

#-- Supp Figure 1C
plot_subCell_counts <- function(fname = "results/suppfig_bioinfoqc/e.pdf",
                                save = TRUE, width = 9, height = 5) {
    load("data/annotations/Proteomics_SubCell_PhysChem_info.RData")
    #-- locs to plot
    locs <- c("Cytoplasm", "Nucleus", "Not annotated", "Membrane", "Mitochondrion",
              "Cytoplasm, cytosol", "Cytoplasm, cytoskeleton",
              "Nucleus, nucleolus", "Mitochondrion inner membrane", "Mitochondrion matrix",
              "Secreted", "Golgi apparatus membrane", "Mitochondrion outer membrane", "Nucleus, nucleoplasm",
              "Lysosome", 'Cell junction', "Cell projection")
    
    #-- Loc annots
    loc_annot_unf <- .getLocAnnot(prot_info)
    loc_annot_filt <- .getLocAnnot(prot_info %>% filter(Filtered))
    
    #-- Plots
    filt_plot <- .plot_subcellLoc(loc_annot_filt) +
        labs(title = "Filtered, 3,781 proteins", x = "") +
        guides(fill = FALSE)
    unfilt_plot <- .plot_subcellLoc(loc_annot_unf) +
        labs(title = "Unfiltered, 15,714 proteins")
    
    #-- Panel
    leg <- get_legend(unfilt_plot)
    unfilt_plot <- unfilt_plot + guides(fill = FALSE)
    
    subcell_plot <- plot_grid(plot_grid(unfilt_plot, filt_plot, scale = 0.9), leg, rel_heights = c(6,1), nrow = 2)
    
    #-- Return
    if(save)  ggsave(filename = fname, plot = subcell_plot, width = width, height = height)
    return(subcell_plot)
}
.getLocAnnot <- function(prot_info) {
    prot_info %>%
        group_by(Cell_Loc2) %>%
        count() %>%
        mutate(loc_overall = case_when(
            Cell_Loc2 %in% c("Cell junction", "Cell projection", "Membrane") ~ "Membrane",
            str_detect(Cell_Loc2, "Cytoplasm") ~ "Cytoplasm",
            Cell_Loc2 %in% c("Endoplasmic reticulum", "Golgi apparatus membrane",
                             "Lysosome") ~ "Organelle",
            str_detect(Cell_Loc2, "Nucleus") ~ "Nucleus",
            str_detect(Cell_Loc2, "Mitochondrion") ~ "Mitochondrion",
            Cell_Loc2 == "Not annotated" ~ NA_character_,
            TRUE ~ Cell_Loc2)) %>%
        rename(loc = Cell_Loc2)
}

.plot_subcellLoc <- function(loc_annot) {
    ggplot(loc_annot, aes(fct_reorder(loc, n), n)) +
        geom_col(aes(fill = loc_overall)) +
        rcartocolor::scale_fill_carto_d(palette = "Pastel") +
        coord_flip() +
        theme_csg_hist +
        theme(panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_line(linetype = "dotted"),
              legend.position = "bottom") +
        labs(x = "Protein subcellular location", y = "Number of proteins", fill = "")
}
#-- Supp Figure 1D
plot_protFilter <- function(unfilt_prot, min_pep = 3, cutoff = 0.66, x_by = 0.033,
                            save = TRUE, fname = "results/suppfig_bioinfoqc/a.pdf",
                            width = 4, height = 3, bar_color = "#360057") {
    nsamps <- ncol(unfilt_prot$pep)
    samp_perc <- seq(0, 1, by = x_by)
    
    #-- Get representation pep/quant
    prep <- sapply(1:nrow(unfilt_prot$wp), function(i) {
        sum(!is.na(unfilt_prot$wp[i,]) & unfilt_prot$pep[i,] >= min_pep)/nsamps
    })
    
    nprot_df <- tibble(sample_rep = samp_perc * 100,
                       nprots = sapply(samp_perc, function(sp) sum(prep > sp))
    )
    
    #-- Make cutoff legend
    cutoff_perc <- cutoff * 100
    cutoff_text <- paste0("Cutoff:\n", cutoff_perc, "% of samples\n",
                          filter(nprot_df, sample_rep == cutoff_perc)$nprots, " proteins")
    
    #-- Plot
    plot_pfilter <- ggplot(nprot_df, aes(sample_rep, nprots)) +
        geom_col(fill = bar_color, color = "white") +
        theme_csg_hist +
        labs(x = "% of samples with protein quantification",
             y = "Number of proteins\n(>= 3 peptides)") +
        scale_y_continuous(breaks = seq(0, 7500, length.out = 6)) +
        scale_x_continuous(expand = c(0,0)) +
        geom_segment(x = cutoff_perc, xend = cutoff_perc, y = 0, yend = 7500,
                     lty = "dashed", color = "black", size = 0.4) +
        annotate("text", x = cutoff_perc + 1, y = quantile(nprot_df$nprots, 0.8),
                 label = cutoff_text, hjust = 0, size = 3.5)
    
    if(save) {
        ggsave(filename = fname, plot = plot_pfilter, width = width, height = height)
    }
    
    return(list(pfilter_df = nprot_df, pfilter_plot = plot_pfilter))
}
#-- Supp Figure 1E
protein_overlap <- function(prot_data, prot_data_bystage,
                            fname = "results/suppfig_bioinfoqc/b.pdf") {
    mibc_data <- prot_data_bystage$mibc
    nmibc_data <- prot_data_bystage$nmibc
    ptotal <- tibble(id = prot_data$wp %>% rownames(), Total = 1)
    pmibc <- tibble(id = mibc_data$wp %>% rownames(), MIBC = 1)
    pnmibc <- tibble(id = nmibc_data$wp %>% rownames(), NMIBC = 1)
    
    #-- Intercept
    all_p <- full_join(ptotal, pmibc) %>% full_join(pnmibc) %>%
        mutate(across(Total:NMIBC, ~ ifelse(is.na(.), 0, .))) %>%
        as.data.frame() %>%
        column_to_rownames("id")
    
    #-- Plots
    pdf(fname, width = 4, height = 4)
    print(plot(eulerr::venn(all_p),
               fills = c("grey30", "#fb8072", "#80b1d3")))
    dev.off()
    
    return(all_p)
}
