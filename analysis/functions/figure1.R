################################################################################
## Protein/mRNA correlations
################################################################################
# Calculate correlations
calculate_protmrna_cor <- function(prot_data, mrna_data) {
    #-- Get data for calculation
    common_ids <- intersect(prot_data$wpAnnot$symbol, rownames(mrna_data$gexp))
    common_samps <- intersect(prot_data$sampAnnot$id, mrna_data$sampleAnnot$id)
    
    prot_ids <- prot_data$wpAnnot %>% dplyr::filter(symbol %in% common_ids)
    prot_ids <- structure(prot_ids$protein_id, names = prot_ids$symbol)
    
    wp_filt <- prot_data$wp[prot_ids,common_samps]
    rownames(wp_filt) <- names(prot_ids)
    gexp_filt <- mrna_data$gexp[names(prot_ids),common_samps]
    
    #-- Calculate spearman correlation for each gene (complete obs)
    prot_mrna_cor <- tibble(Gene = names(prot_ids),
                            pm_rho = sapply(1:nrow(wp_filt), function(i) {
                                cor(wp_filt[i,], gexp_filt[i,], method = "spearman", use = "complete")
                            }),
                            pm_cor = sapply(1:nrow(wp_filt), function(i) {
                                cor(wp_filt[i,], gexp_filt[i,], method = "pearson", use = "complete")
                            }),
                            )
    #-- Add p-values
    nsamp <- ncol(mrna_data$gexp)
    prot_mrna_cor <- prot_mrna_cor %>%
        mutate(cor_pval = .cor_pval(pm_rho, nsamp),
               cor_padj = p.adjust(cor_pval, "BH"))
    
    return(prot_mrna_cor)
}
.cor_pval <- function(r,n) {
    t <- r*sqrt((n-2)/(1-r^2))
    p <- 1 - pt(t, n - 1)
    return(p)
}
#-------------------------------------------------------------------------------
# Pathway correlation distributions
path_protmrna_cor <- function(prot_mrna_cor, path_table) {
    #-- Correlation by pathway (is the distribution of values different)
    paths <- path_table$Path_ID %>% unique()
    
    #-- Calculate distributions
    pathway_distributions <- lapply(paths, function (path) {
        genesInPath <- path_table %>% filter(Path_ID == path) %>% pull(Symbol)
        path_name <- path_table %>% filter(Path_ID == path) %>% pull(Path_Name) %>% unique()
        genesPresent <- genesInPath[genesInPath %in% prot_mrna_cor$Gene]
        pathReprs <- length(genesPresent)/length(genesInPath)
        
        pathRhos <- prot_mrna_cor %>%
            filter(Gene %in% genesPresent) %>%
            pull(pm_rho)
        if(pathReprs > 0) KS_pValue <- ks.test(pathRhos, prot_mrna_cor$pm_rho)$p.value
        else KS_pValue <- 1
        
        return(data.frame("Path_ID" = path, "Pathway_Name" = path_name, 
                          "Pathway_Representation" = pathReprs, "KS_pValue" = KS_pValue, 
                          "Median_rho" = median(pathRhos),
                          "Pathway_Size" = length(genesInPath)))
    }) %>% bind_rows()
    
    #-- P-adjust and order
    pathway_rhodist <- pathway_distributions %>%
        tibble() %>%
        mutate(KS_pAdj = p.adjust(KS_pValue, "BH")) %>%
        arrange(KS_pValue, desc(Pathway_Representation))
    
    #-- Add genes
    pathway_rhodist <- pathway_rhodist %>%
        mutate(Higher = Median_rho >= 0.47) %>%
        filter(Pathway_Representation > 0)
    
    path_ls <- split(path_table$Symbol, path_table$Path_Name)
    pathway_rhodist$Genes <- sapply(1:nrow(pathway_rhodist), function(i) {
        path <- pathway_rhodist$Pathway_Name[i]
        path_dir <- pathway_rhodist$Higher[i]
        path_genes <- path_ls[[path]]
        path_cors <- prot_mrna_cor %>% filter(Gene %in% path_genes) %>% pull(pm_rho, Gene)
        genes <- sort(path_cors, decreasing = path_dir) %>% names()
        paste0(genes, collapse = ", ")
    })
    pathway_rhodist <- select(pathway_rhodist, -Higher)
    
    return(pathway_rhodist)
}
#-------------------------------------------------------------------------------
#-- Plot Figure 1C
plot_protmrna_cor <- function(prot_mrna_cor, pathway_rhodist, path_table, 
                              paths4plot = NULL, rel_heights = c(1,1.5,2.2),
                              rel_widths = c(6,0.5,3)) {
    
    #-- Get pathway tables ready for plotting
    paths4plot_names <- structure(pathway_rhodist$Pathway_Name, 
                                  names = pathway_rhodist$Path_ID)[paths4plot]
    names(paths4plot_names) <- NULL
    mean_rho <- mean(prot_mrna_cor$pm_rho)
    
    path_select <- path_table %>%
        filter(Path_ID %in% paths4plot) %>%
        left_join(pathway_rhodist, by = "Path_ID") %>%
        mutate(mean2ref = case_when(
            (Median_rho > mean_rho & KS_pAdj < 0.05) ~ "higher",
            (Median_rho < mean_rho & KS_pAdj < 0.05) ~ "lower",
            TRUE ~ "ns"),
            Path_Name = factor(Path_Name, levels = paths4plot_names)) %>%
        select(path_id = Path_ID, path_name = Path_Name, Gene = Symbol, path_rep = Pathway_Representation, mean2ref)
    
    pm_cors <- prot_mrna_cor %>%
        mutate(path = "Overall") %>%
        left_join(path_select)
    
    #-- Get count tables
    path_count <- path_select %>%
        group_by(path_name) %>%
        dplyr::count() %>%
        mutate(n_txt = paste0("(n = ", n, ")")) %>%
        left_join(select(pathway_rhodist, path_name = Pathway_Name, path_rep = Pathway_Representation)) %>%
        ungroup() %>%
        mutate(path_name = factor(path_name, levels = paths4plot_names))
    
    #-- Get KS tables
    path_ks <- pathway_rhodist %>%
        filter(Path_ID %in% paths4plot) %>%
        mutate(
            Pathway_Name = factor(Pathway_Name, levels = paths4plot_names),
            signif = case_when(
                KS_pAdj < 0.001 ~ "***",
                KS_pAdj < 0.01 ~ "**",
                KS_pAdj < 0.05 ~ "*",
                TRUE ~  "ns"),
            x = 0) %>%
        select(path_name = Pathway_Name, x, signif) 
    
    #-- Overall distribution ridge
    ov_plot <- ggplot(pm_cors, aes(x = pm_rho)) +
        geom_histogram(binwidth = 0.05, alpha = 0.6) +
        geom_density(aes(y = 0.05*..count..), size = 0.7) +
        geom_vline(xintercept = mean_rho, lty = "dashed") +
        lims(x = c(-0.5, 1)) +
        labs(y = "Number of proteins") +
        theme_csg_sparse +
        theme(panel.grid.major.x = element_line(color = "grey60", linetype = "dotted"),
              axis.text.y = element_text(size=8),
              axis.title.y = element_text(size=8, face= "bold"),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              axis.line.y = element_line()) +
        coord_cartesian(clip = "off")
    
    #-- Pathway plots
    path_ridges <- .plot_protcor_ridges(pm_cors, y = "path_name", ylab = "Pathway", mean_rho = mean_rho) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    path_rep <- .plot_protcor_pathrep(path_count) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    path_stars <- .plot_protmrna_signif(path_ks, y = "path_name") + theme(axis.title.x = element_blank())
    
    #-- Get signatures
    biton_ic_tb <- .get_bitonICAs()
    biton_ics <- lapply(1:length(biton_ic_tb), function(i)
        structure(biton_ic_tb[,i], names = rownames(biton_ic_tb)))
    names(biton_ics) <- colnames(biton_ic_tb)
    
    biton_ls <- .get_topICgenes(biton_ics, nsd = 3)
    save(biton_ls, file = "data/annotations/biton_ics.RData")
    
    #------------ Added IC signatures from other sources
    m10_11p15 <- read_tsv("data/Biton_ICA_signatures/M10_11p15_5_top.csv", col_names = c("gene", "score"))
    m18_neuro <- read_tsv("data/Biton_ICA_signatures/M18_NEUROENDOCRINE_top.tsv", col_names = c("gene", "score"))
    
    biton_ls[["M10_11p15.5"]][["top"]] <- m10_11p15 %>% filter(score > 2.5) %>% pull(gene)
    biton_ls[["M10_11p15.5"]][["bottom"]] <- m10_11p15 %>% filter(score < -2.5) %>% pull(gene)
    biton_ls[["M18_NEUROENDOCRINE"]][["top"]] <- m18_neuro %>% filter(score > 2.5) %>% pull(gene)
    biton_ls[["M18_NEUROENDOCRINE"]][["bottom"]] <- m18_neuro %>% filter(score < -2.5) %>% pull(gene)
    
    #-- Make list
    biton_ls <- lapply(biton_ls, unlist)
    names(biton_ls) <- c("IC Myofibroblasts", "IC BLCA pathways", "IC Stress", 
                         "IC Smooth Muscle", "IC Mitochondrial Translation", "IC Interferon",
                         "IC Basal-like", "IC Cell cycle", "IC Immune", "IC Urothelial differentiation",
                         "IC 11p15.5", "IC Neuroendocrine")
    # load("Data/base47.RData")
    # sigs <- rlist::list.append(biton_ls, BASE47 = base47, 
    #                            "FGFR3 co-expressed genes" = c("FGFR3", "TP63", "IRS1", "WNT7B", "INKA1", "CAPNS2", "SEMA4B", 
    #                                                           "DUOXA1", "C16orf74", "ZNF385A", "SMAD3", "SLC2A9", "NSG1", "CLCA4", 
    #                                                           "SYTL1", "PLCH2", "SSH3", "PTPN13", "DUOX1", "TMPRSS4"))
    sigs <- biton_ls
    sig_rhodist <- sig_protmrna_cor(prot_mrna_cor, sigs)
    order <- sig_rhodist %>% arrange() %>% pull(Signature) %>% rev()
    
    #-- Sig rhodist
    unlist_sigs <- lapply(sigs, unlist) %>% unlist()
    sig_cors <- tibble(sig_name = rep(names(sigs), sapply(sigs, length)), Gene = unlist_sigs) %>%
        mutate(weight_dir = ifelse(str_detect(names(unlist_sigs), "top"), "up", "down")) %>%
        left_join(prot_mrna_cor) %>%
        filter(!is.na(pm_rho)) %>%
        left_join(sig_rhodist, by = c("sig_name" = "Signature")) %>%
        mutate(mean2ref = case_when(
            (Median_rho > mean_rho & KS_pAdj < 0.05) ~ "higher",
            (Median_rho < mean_rho & KS_pAdj < 0.05) ~ "lower",
            TRUE ~ "ns"),
            sig_name = factor(sig_name, levels = order)) %>%
        select(sig_name, Gene, pm_rho, path_rep = Pathway_Representation, mean2ref, weight_dir)
    
    save(sig_cors, file = "data/annotations/biton_cors.RData")
    
    #-- Signature counts / signif
    sig_count <- sig_cors %>%
        group_by(sig_name) %>%
        count() %>%
        mutate(n_txt = paste0("(n = ", n, ")")) %>%
        left_join(select(sig_rhodist, sig_name = Signature, path_rep = Pathway_Representation)) %>%
        ungroup() %>%
        mutate(path_name = factor(sig_name, levels = order))
    
    sig_ks <- sig_rhodist %>%
        mutate(
            path_name = factor(Signature, levels = order),
            signif = case_when(
                KS_pAdj < 0.001 ~ "***",
                KS_pAdj < 0.01 ~ "**",
                KS_pAdj < 0.05 ~ "*",
                TRUE ~  "ns"),
            x = 0) %>%
        select(path_name, x, signif) 
    
    #-- Signature plots
    sig_ridges <- .plot_protcor_ridges(sig_cors, y = "sig_name", ylab = "Signature", mean_rho = mean_rho)
    sig_stars <- .plot_protmrna_signif(sig_ks, y = "path_name") 
    sig_rep <- .plot_protcor_pathrep(sig_count)
    
    #-- Make grid
    pmcor_plot <- ((ov_plot / path_ridges / sig_ridges) + plot_layout(heights = rel_heights) |
                       ((plot_spacer() / path_stars / sig_stars) + plot_layout(heights = rel_heights)) | 
                       (plot_spacer() / path_rep / sig_rep) + plot_layout(heights = rel_heights)) +
        plot_layout(widths = rel_widths)
    
    ggsave(filename = "results/fig1/fig1b.pdf", plot = pmcor_plot, 
           width = 7, height = 4.5)
    
    return(pmcor_plot)
}
.plot_protcor_ridges <- function(pm_cors, x = "pm_rho", y, ylab = "", mean_rho) {
    y <- sym(y)
    x <- sym(x)
    pm_cors %>%
        filter(!is.na(!! y)) %>%
        ggplot(aes(!!x, y = !! y, fill = mean2ref)) +
        geom_density_ridges(alpha = 0.5) +
        geom_vline(xintercept = mean_rho, lty = "dashed") +
        scale_x_continuous(limits = c(-0.5,1)) +
        scale_fill_manual(values = c("higher" = "#ca0020", "ns" = "grey60", "lower" = "#0571b0")) +
        theme_csg_sparse +
        theme(panel.grid.major.x = element_line(color = "grey60", linetype = "dotted")) +
        guides(fill = "none", color = "none") +
        labs(x = expression("mRNA/Protein correlation"~italic(rho)~"distribution"), y = ylab) +
        theme(axis.title = element_text(size=8, face="bold")) +
        coord_cartesian(clip = "off")
}

.plot_protcor_pathrep <- function(var_count, xlab = "% of pathway represented") {
    var_count %>%
        ggplot(aes(path_name, path_rep)) +
        geom_bar(stat = "identity") +
        labs(y = xlab) +
        scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        scale_x_discrete(labels = var_count$n_txt) +
        coord_flip() +
        theme_csg_sparse +
        theme(axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              axis.line.x = element_blank(),
              panel.grid.major.x = element_line(color = "grey60", linetype = "dotted"),
              legend.position = "bottom")
}

.plot_protcor_counts <- function(var_count, x, fill = NULL, fill_lab = "",
                                 xlab = "Number of proteins") {
    x <- sym(x)
    if(!is.null(fill)) {
        count_plt <- var_count %>%
            ggplot(aes(!! x, n, fill = !! sym(fill))) +
            geom_bar(stat = "identity") +
            scale_fill_gradient2(low = "#f7f7f7", high = "#252525") +
            labs(y = xlab, fill = fill_lab)
    } else {
        count_plt <- var_count %>%
            ggplot(aes(!! x, n)) +
            geom_bar(stat = "identity", size = 0.5) +
            labs(y = xlab)
    }
    
    count_plt +
        scale_y_continuous(breaks = c(0,200,400,600,800), limits = c(0,800)) +
        scale_x_discrete(labels = var_count$n_txt) +
        coord_flip() +
        theme_csg_sparse +
        theme(axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              axis.line.x = element_blank(),
              panel.grid.major.x = element_line(color = "grey60", linetype = "dotted"),
              legend.position = "bottom")
}

.plot_protmrna_signif <- function(ks_ptb, y) {
    ks_ptb %>%
        ggplot(aes(x, !! sym(y), label = signif)) +
        geom_text(size = 4) +
        theme_void() +
        theme(axis.title.x = element_text(size=8, face="bold")) +
        labs(x = expression("KS adj."~italic("p")), y = "")
}
.get_bitonICAs <- function() {
    signatures <- list.files("data/Biton_ICA_signatures/") %>% str_subset(".rnk") %>% str_remove(".rnk$")
    metagenes <- lapply(signatures, function(sig) {
        read_tsv(str_glue("data/Biton_ICA_signatures/{sig}.rnk"), 
                 col_names = c("gene", str_glue("{sig}_weight")))
    })
    metagenes <- Reduce(full_join, metagenes)
    
    metagene_mat <- metagenes %>%
        as.data.frame() %>%
        column_to_rownames("gene")
    return(metagene_mat)
}
.get_topICgenes <- function(ic_ls, weight = 4, nsd = NA) {
    if(!is.na(nsd)) {
        ic_nsds <- sapply(ic_ls, sd, na.rm = TRUE) * nsd
        ic_means <- sapply(ic_ls, mean, na.rm = TRUE)
        ic_thresh1 <- ic_means + ic_nsds
        ic_thresh2 <- ic_means - ic_nsds
    } else {
        ic_thresh1 <- weight
        ic_thresh2 <- -weight
    }
    top_genes <-  lapply(1:length(ic_ls), function(i) {
        icv <- ic_ls[[i]]
        idx <- icv > ic_thresh1[i]
        idx <- ifelse(is.na(idx), FALSE, idx)
        genes <- names(icv)[idx] })
    bottom_genes <- lapply(1:length(ic_ls), function(i) {
        icv <- ic_ls[[i]]
        idx <- icv < ic_thresh2[i]
        idx <- ifelse(is.na(idx), FALSE, idx)
        genes <- names(icv)[idx] })
    top_ls <- lapply(1:length(ic_ls), function(i) {
        list(top = top_genes[[i]], bottom = bottom_genes[[i]])
    })
    names(top_ls) <- names(ic_ls)
    return(top_ls)
}
sig_protmrna_cor <- function(prot_mrna_cor, sigs) {
    # Correlation by signature (is the distribution of values different)
    
    #-- Calculate distributions
    sig_distributions <- lapply(1:length(sigs), function (i) {
        sig <- sigs[[i]]
        genesPresent <- sig[sig %in% prot_mrna_cor$Gene]
        sigReprs <- length(genesPresent)/length(sig)
        
        sigRhos <- prot_mrna_cor %>%
            filter(Gene %in% genesPresent) %>%
            pull(pm_rho)
        if(sigReprs > 0) KS_pValue <- ks.test(sigRhos, prot_mrna_cor$pm_rho)$p.value
        else KS_pValue <- 1
        
        return(data.frame("Signature" = names(sigs)[i], 
                          "Pathway_Representation" = sigReprs, "KS_pValue" = KS_pValue, 
                          "Median_rho" = median(sigRhos),
                          "Pathway_Size" = length(sig)))
    }) %>% bind_rows()
    
    #-- P-adjust and order
    sig_rhodist <- sig_distributions %>%
        tibble() %>%
        mutate(KS_pAdj = p.adjust(KS_pValue, "BH")) %>%
        arrange(KS_pValue, desc(Pathway_Representation))
    
    return(sig_rhodist)
}
################################################################################
## CNA correlations
################################################################################
cna_cors <- function(cna_gr, mrna_data, prot_data, sampAnnot) {
    load("data/annotations/Ensembl99_geneAnnot.RData")
    consensus_genes <- intersect(intersect(cna_gr$gene_name,
                                           mrna_data$geneAnnot$symbol), 
                                 prot_data$wpAnnot$symbol)
    consensus_genes <- prot_data$wpAnnot %>% 
        filter(symbol %in% consensus_genes) %>% 
        pull(symbol, protein_id)
    
    mrna <- mrna_data$gexp[consensus_genes,sampAnnot$sample]
    prot <- prot_data$wp[names(consensus_genes),sampAnnot$sample]
    rownames(prot) <- consensus_genes
    cna <- cna_gr %>% 
        filter(gene_name %in% consensus_genes) %>% 
        select(gene_name, starts_with("CIT")) %>%
        as.data.frame() %>%
        column_to_rownames("gene_name")
    cna <- cna[consensus_genes,sampAnnot$sample]
    cna <- as.matrix(cna)
    
    prot_cors <- sapply(consensus_genes, function(gene) cor(cna[gene,], prot[gene,], 
                                                            method = "spearman", use = "complete"))
    mrna_cors <- sapply(consensus_genes, function(gene) cor(cna[gene,], mrna[gene,], 
                                                            method = "spearman", use = "complete"))
    
    ks.test(prot_cors, mrna_cors, alternative = "greater")
    
    cna_cors <- tibble(gene = consensus_genes,
                       prot = prot_cors,
                       prot_padj = p.adjust(.cor_pval(prot, 62), method = "BH"),
                       mrna  = mrna_cors,
                       mrna_padj = p.adjust(.cor_pval(mrna, 62), method = "BH")) %>%
        left_join(select(gene_ensembl, seqnames, start, end, gene = gene_name)) %>%
        filter(!is.na(start)) %>%
        mutate(seqnames = factor(seqnames, c(1:22, "X", "Y"))) %>%
        filter(!is.na(seqnames)) %>%
        group_by(gene) %>%
        slice(1) %>%
        ungroup() %>%
        arrange(seqnames, start) %>%
        mutate(pos = 1:n())
    
    return(cna_cors)
}

compare_cna_cors <- function(cna_cor_res, cna_gr) {
    #-- Compare correlations
    cna4plot <- cna_cor_res %>%
        left_join(select(cna_gr, gene = gene_name, cytoband)) %>%
        mutate(which_high = case_when(
                   mrna > prot+0.2 ~ "Higher mRNA correlation",
                   prot > mrna+0.2 ~ "Higher protein correlation",
                   TRUE ~ "Neither"))
    
    #-- Plot
    scatter <- cna4plot %>%
        ggplot(aes(mrna, prot, color = which_high)) +
        geom_point(size = 0.3) +
        geom_abline(slope = 1, intercept = 0, size = 0.3) +
        lims(x = c(-0.4, 0.75), y = c(-0.4, 0.75)) +
        scale_color_manual(values = c("Higher mRNA correlation" = "#998ec3",
                                      "Higher protein correlation" = "#f1a340",
                                      "Neither" = "grey60")) +
        labs(x = "CNA/mRNA correlations", y = "CNA/Protein correlations", color = "") +
        guides(color = "none") +
        annotate("text", x = -0.3, y = 0.5, label = "Higher protein\ncorrelations", size = 3,
                 hjust = 0) +
        annotate("text", x = 0.25, y = -0.25, label = "Higher mRNA\ncorrelations", size = 3,
                 hjust = 0) +
        theme_csg_scatter
    #-- Top histograms
    histo1 <- ggplot(cna4plot, aes(mrna)) +
        geom_density(fill = "white") +
        geom_vline(xintercept = 0.28, lty = "dashed", color = "black") +
        lims(x = c(-0.4, 0.75)) +
        annotate("text", x = 0.3, y = 0.5, label = "mean CNA/mRNA\ncor = 0.28", hjust = 0, size = 3) +
        theme_void()
    histo2 <- ggplot(cna4plot, aes(prot)) +
        geom_density(fill = "white") +
        geom_vline(xintercept = 0.19, lty = "dashed", color = "black") +
        lims(x = c(-0.4, 0.75)) +
        coord_flip() +
        annotate("text", x = 0.3, y = 0, label = "mean CNA/Protein\ncor = 0.19", size = 3, hjust = 0) +
        theme_void()
    #-- Join plots
    left <- (histo1 / scatter) + plot_layout(heights = c(0.3,1))
    right <-  (plot_spacer() / histo2) + plot_layout(heights = c(0.3,1))
    
    fig1c <- (left | right) + 
        plot_layout(widths = c(1,0.3), guides = "collect")
    
    ggsave("results/fig1/fig1c.pdf", fig1c, width = 3.3, height = 2.5)
    
    return(fig1c)
}

cna_cors_genomePlots <- function(cna_cor_res, cna_calls, core_samples) {
    #-- Get chr limits
    line_pos <- cna_cor_res %>%
        group_by(seqnames) %>%
        summarize(n = n()) %>%
        mutate(end_pos = cumsum(n),
               mid_pos = end_pos - n/2)
    #-- Get significant correlations (>0.3)
    cna_plt3_signif <- cna_cor_res %>%
        pivot_longer(prot:mrna, names_to = "omic", values_to = "rho") %>%
        mutate(signif = ifelse(rho >= 0.3, "rho ≥ 0.3", "rho < 0.3"),
               omic = factor(omic, levels = c("prot", "mrna"))) %>%
        ggplot(aes(pos, omic)) +
        geom_raster(aes(fill = signif)) +
        scale_fill_manual(values = c("rho ≥ 0.3" = "grey50", "rho < 0.3" = "white")) +
        geom_vline(xintercept = line_pos$end_pos, color = "black", size = 0.3) +
        scale_x_continuous(breaks = line_pos$mid_pos, labels = line_pos$seqnames) +
        scale_y_discrete(labels = c("CNA/Protein", "CNA/mRNA")) +
        labs(x = "", y = "Correlation", fill = "Significance") +
        theme_csg_sparse +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank())
    #-- Get frequency of alterations in genes from CNA calls (GISTIC) 
    cna_call_mat <- cna_calls %>%
        mutate(symbol = `Gene Symbol`) %>%
        select(symbol, starts_with("CIT")) %>%
        as.data.frame() %>%
        column_to_rownames("symbol") %>%
        as.matrix() %>%
        .[,core_samples]
    alt_freq <- tibble(
        gene = rownames(cna_call_mat),
        gain = rowSums(cna_call_mat == 1)/ncol(cna_call_mat),
        amp = rowSums(cna_call_mat == 2)/ncol(cna_call_mat),
        loss = rowSums(cna_call_mat == -1)/ncol(cna_call_mat),
        dd = rowSums(cna_call_mat == -2)/ncol(cna_call_mat))
    
    cna_alt_tb <- left_join(cna_cor_res, alt_freq) %>%
        rowwise() %>%
        mutate(alt_freq = gain + amp + loss + dd) %>%
        ungroup()
    
    #-- Plot gains/losses
    gain_loss <- 
        cna_alt_tb %>%
        pivot_longer(gain:dd, names_to = "type", values_to = "freq") %>%
        mutate(freq = ifelse(type %in% c("loss", "dd"), -freq, freq),
               event_type = ifelse(type %in% c("loss", "dd"), "loss", "amp")) %>%
        ggplot(aes(pos, freq, fill = type)) +
        facet_wrap(~ event_type, nrow = 2, scales = "free_y") +
        geom_area() +
        scale_fill_manual(values = c("gain" = "#ef8a62", "loss" = "#67a9cf",
                                     "amp" = "#ca0020", "dd" = "#0571b0")) +
        geom_vline(xintercept = line_pos$end_pos, color = "black", size = 0.3) +
        scale_x_continuous(breaks = line_pos$mid_pos, labels = line_pos$seqnames) +
        scale_y_continuous(labels = function(x) scales::percent(abs(x)),
                           breaks = c(-0.4,-0.2,0,0.2,0.4)) +
        theme_csg_sparse +
        theme(strip.text = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
        labs(x = "Chromosomic position", y = "Frequency")
    
    #-- Join plots
    fig1e <- (cna_plt3_signif / gain_loss) + plot_layout(heights = c(0.5,1))
    
    ggsave("results/fig1/fig1e.pdf", fig1e, width = 7, height = 2.5)
    
    return(list(fig = fig1e, cna_alt_tb = cna_alt_tb))
}

cna_freqAlt_compare <- function(cna_genome_plots) {
    cna_alt_tb <- cna_genome_plots$cna_alt_tb
    cna_alt_tb2 <- cna_alt_tb %>%
        mutate(
            alt_cat = case_when(
                alt_freq < 0.3 ~ paste0("<30%\n(n=", sum(cna_alt_tb$alt_freq < 0.3), ")"),
                alt_freq < 0.5 ~ paste0("30-50%\n(n=", sum(cna_alt_tb$alt_freq < 0.5 &
                                                               cna_alt_tb$alt_freq > 0.3), ")"),
                alt_freq >= 0.5 ~ paste0(">50%\n(n=", sum(cna_alt_tb$alt_freq >= 0.5), ")")),
            alt_cat = factor(alt_cat, levels = c("<30%\n(n=503)", "30-50%\n(n=2294)", ">50%\n(n=383)"), ordered = TRUE)) %>%
        pivot_longer(cols = prot:mrna, names_to = "omic", values_to = "cor")
    fig1d <- ggplot(cna_alt_tb2, aes(alt_cat, cor, fill = alt_cat)) +
        facet_wrap(~ omic, ncol = 1) +
        geom_violin(alpha = 0.3) +
        geom_boxplot(width = 0.2) +
        stat_compare_means(size = 3, label.y.npc = "bottom") +
        scale_fill_brewer(palette = "OrRd") +
        labs(x = "CNA frequency", y = "CNA correlation", fill = "CNA frequency") +
        theme_csg_scatter
    
    ggsave(fig1d, file = "results/fig1/fig1d.pdf", width = 3.5, height = 3)
    
    return(fig1d)
    
}

median_pep_quant <- function(unfilt_prot) {
    pep <- unfilt_prot$pep
    prop_3pep <- rownames(pep)[]
    idx <- apply(pep, 1, function(prot) prot > 3)
    prot_3pep <- rownames(pep)[idx] %>% na.exclude()
    mean(rowSums(pep[prop_3pep,], na.rm = TRUE))
    median(rowSums(pep[prop_3pep,], na.rm = TRUE))

}

plot_pp_pm_complex_cors <- function(prot_data, prot_mrna_cor) {
  # CORUM ------------------
  complex_data <- read_tsv("data/annotations/CORUM_humanComplex_20220725es.txt")
  
  complex_long <- complex_data %>%
    select(complex_name = ComplexName, subunits = `subunits(UniProt IDs)`) %>%
    mutate(subunits = str_split(subunits, ";")) %>%
    unnest(cols = subunits)
  complex_list <- split(complex_long$subunits, complex_long$complex_name)
  
  ## Get interactions ----------
  prots <- prot_data$wpAnnot$protein_id
  p_complex <- lapply(complex_list, function(complex_prots) prots %in% complex_prots)
  
  all_interactions <- lapply(complex_list, function(complex_prots) {
    if(length(complex_prots) < 2) {
      return(tibble())
    } else {
      interactions <- combn(complex_prots, 2) %>%
        t() %>%
        as.data.frame() %>%
        tibble()
    }
  }) %>%
    bind_rows()
  
  all_interactions <- all_interactions %>%
    mutate(pair = map2_chr(V1,V2, function(x,y) { paste(sort(c(x,y)), collapse = "~")}))
  
  # Calculate prot/prot cors -----------
  cached_file <- "cached_results/prot_prot_cor.RData"
  if(file.exists(cached_file)) {
    load(cached_file)
  } else {
    prot_prot_cors <- cor(t(prot_data$wp), t(prot_data$wp), use = "pair", method = "spearman")
    prot_prot_cors <- prot_prot_cors %>%
      replace(upper.tri(.), NA) %>%
      reshape2::melt(na.rm = TRUE) %>%
      tibble() %>%
      rename(prot1 = Var1, prot2 = Var2, prot_cor = value)
    
    save(prot_prot_cors_tb, file = "cached_results/prot_prot_cor.RData")
  }
  
  # Join tables (janky but works) ---------
  pp_cors1 <- prot_prot_cors %>%
    inner_join(select(all_interactions, prot1 = V1, V2), by = "prot1") %>%
    filter(prot2 == V2) %>%
    distinct()
  pp_cors2 <- prot_prot_cors %>%
    inner_join(select(all_interactions, prot1 = V2, V1), by = "prot1") %>%
    filter(prot2 == V1) %>%
    distinct()
  pp_cors3 <- prot_prot_cors %>%
    inner_join(select(all_interactions, prot2 = V1, V2), by = "prot2") %>%
    filter(prot1 == V2) %>%
    distinct()
  pp_cors4 <- prot_prot_cors %>%
    inner_join(select(all_interactions, prot2 = V2, V1), by = "prot2") %>%
    filter(prot1 == V1) %>%
    distinct()
  
  # Get complex protein/protein cors --------
  complex_cors <- bind_rows(pp_cors1, pp_cors2, pp_cors3, pp_cors4) %>%
    select(-V1, -V2, -prot_cor) %>%
    mutate(status = "in same complex")
  
  pp_cors_complex <- prot_prot_cors %>%
    left_join(complex_cors, by = c("prot1", "prot2")) %>%
    filter(prot1 != prot2) %>%
    mutate(status = ifelse(is.na(status), "not in same complex", status) %>%
             factor(levels = c("in same complex", "not in same complex"))) 
  
  mean_pp_cors <- pp_cors_complex %>%
    group_by(status) %>%
    summarize(mean_rho = mean(prot_cor) %>% signif(3))
  
  # Get proteins that are in complex prot/mRNA cors -------
  prots_in_complex <- unique(unique(complex_cors$prot1), unique(complex_cors$prot2))
  length(prots_in_complex)
  
  wp_annot <- prot_data$wpAnnot %>%
    select(protein_id, Gene = symbol)
  
  pm_cors <- prot_mrna_cor %>%
    left_join(wp_annot) %>%
    mutate(status = ifelse(protein_id %in% prots_in_complex, "in a complex", "not in a complex"))
  
  mean_pm_cors <- pm_cors %>%
    group_by(status) %>%
    summarize(mean_rho = mean(pm_rho) %>% signif(3))
  
  #-- KS tests
  pm_cor_ls <- pm_cors %>% pull(pm_cor, status) %>% split(., names(.))
  pp_cor_ls <- pp_cors_complex %>% pull(prot_cor, status) %>% split(., names(.))
  
  ks.test(pm_cor_ls$`in a complex`, pm_cor_ls$`not in a complex`)
  ks.test(pp_cor_ls$`in same complex`, pm_cor_ls$`not in a complex`)
  
  ks_label <- "KS p-value < 2.2e-16"
  
  #-- Plot
  plt_pp_cor_complex <- ggplot(pp_cors_complex, aes(x = prot_cor, y = status, fill = status)) +
    stat_density_ridges(alpha = 0.5, quantile_lines = TRUE, quantiles = 2) +
    geom_text(data = mean_pp_cors, aes(y = status, x = mean_rho, label = mean_rho), 
              hjust = 0, nudge_x = 0.05, nudge_y = 0.3, size = 3) +
    annotate("text", label = ks_label, size = 3, x = Inf, y = Inf, hjust = 1.5) +
    guides(fill = "none") +
    labs(x =  expression("Protein/protein correlation"~italic(rho)~"distribution"), y = "Protein type") +
    scale_x_continuous(position = "top", limits = c(-0.5,1)) +
    coord_cartesian(clip = "off") +
    theme_csg_sparse +
    theme(axis.title = element_text(size=8, face="bold"),
          panel.grid.major.x = element_line(color = "grey60", linetype = "dotted"))
  
  plt_pm_cor_complex <- ggplot(pm_cors, aes(x = pm_rho, y = status, fill = status)) +
    stat_density_ridges(alpha = 0.5, quantile_lines = TRUE, quantiles = 2) +
    geom_text(data = mean_pm_cors, aes(y = status, x = mean_rho, label = mean_rho), 
              hjust = 0, nudge_x = 0.05, nudge_y = 0.3, size = 3) +
    annotate("text", label = ks_label, size = 3, x = Inf, y = -Inf, hjust = 1.5) +
    guides(fill = "none") +
    labs(x =  expression("mRNA/protein correlation"~italic(rho)~"distribution"), y = "Protein type") +
    lims(x = c(-0.5,1)) +
    coord_cartesian(clip = "off") +
    theme_csg_sparse +
    theme(axis.title = element_text(size=8, face="bold"),
          panel.grid.major.x = element_line(color = "grey60", linetype = "dotted"))
  
  plt_pp_pm_densities <- plt_pp_cor_complex / plt_pm_cor_complex
  
  ggsave("results/suppfig_bioinfoqc/f.pdf", plt_pp_pm_densities, width = 5, height = 3)
  
  return(list(pp_cors = pp_cors_complex, pm_cors = pm_cors, plt = plt_pm_cor_complex))
}

cna_cor_example <- function(gene, cna_gr, cna_calls, mrna_data, prot_symbol, core_samples) {
  gene_title <- cna_gr %>%
    filter(gene_name == gene) %>%
    mutate(name = paste(c(gene_name, cytoband), collapse = " | ")) %>%
    pull(name)
  
  cna_values <- cna_gr %>%
    filter(gene_name == gene) %>%
    select(starts_with('CIT')) %>%
    unlist() %>%
    .[core_samples]
  cna_call <- cna_calls %>%
    filter(`Gene Symbol` == !! gene) %>%
    select(starts_with('CIT')) %>%
    unlist() %>%
    .[core_samples]
  prot_values <- prot_symbol[gene,core_samples]
  mrna_values <- mrna_data$gexp[gene,core_samples]
  
  tibble(
    id = core_samples,
    cna = cna_values,
    cna_call = cna_call,
    mrna = mrna_values,
    prot = prot_values) %>%
    pivot_longer(cols = mrna:prot, names_to = "level", values_to = "exp") %>%
    # mutate(cna_call = case_when(
    #   cna_call == -2 ~ "del",
    #   cna_call == -1 ~ "loss",
    #   cna_call == 0 ~ "no cnv",
    #   cna_call == 1 ~ "gain"
    # ) %>% factor(levels = c("del","loss", "no cnv", "gain"))) %>%
    mutate(cna_call = factor(cna_call, levels = -2:2)) %>%
    ggplot(aes(cna, exp)) +
    facet_wrap(~ level, scales = "free") +
    geom_point(aes(fill = cna_call), pch = 21) +
    geom_smooth(method = "lm", se = FALSE, lty = "dotted", color = "black", size = 0.5) +
    scale_fill_manual(values = c("-2" = "#0571b0", "-1" = "#92c5de",
                                 "0" = "#f7f7f7", "1" = "#f4a582",
                                 "2" = "#ca0020")) +
    stat_cor(aes(label = ..r.label..), size = 3, method = "spearman") +
    labs(title = gene_title, x = "CNA smoothed value", y = "Expression", fill = "CNA call") +
    theme_csg_scatter
}


