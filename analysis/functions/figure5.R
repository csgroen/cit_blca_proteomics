## TRAIL mechanism in FGFR3* BLCA ----------------------
dir.create("results/fig5/", showWarnings = FALSE)
################################################################################
## Check levels of FAS/TRAIL/TNF
################################################################################
plot_siFGFR3_apoptPath <- function() {
    #-- Read data
    clines <- readxl::excel_sheets("data/experiments/siFGFR3/200925-UMUC14-MGHU3-RT112-Pool_Transcriptomic_results-LIMMA-siFGFR3vsLipo.xlsx")
    sifgfr3_fc <- lapply(clines, function(i) {
        read_excel("data/experiments/siFGFR3/200925-UMUC14-MGHU3-RT112-Pool_Transcriptomic_results-LIMMA-siFGFR3vsLipo.xlsx", 
                   sheet = i) %>%
            mutate(logFC = as.numeric(logFC),
                   adj.P.Val = as.numeric(adj.P.Val)) %>%
            select(symbol, logFC, padj = adj.P.Val)})
    clines <- str_remove(clines, "-LIMMA-siFGFR3vsLipo")
    names(sifgfr3_fc) <- clines
    
    #-- Get genes
    # receptors <- c("TRAIL-R1" = "TNFRSF10A", "TRAIL-R2" = "TNFRSF10B", "FAS" = "FAS", "TNF-R1" = "TNFRSF1A")
    genes <- c("TRAIL-R1" = "TNFRSF10A", "c-FLIP" = "CFLAR", "TRAIL-R2" = "TNFRSF10B")
    
    #-- Plot
    data4plot <- lapply(clines, function(cl) { 
        fc <- sifgfr3_fc[[cl]]
        fc %>%
            filter(symbol %in% genes) %>%
            mutate(cell_line = cl,
                   gene_name = plyr::mapvalues(symbol, genes, names(genes)))
    }) %>%
        bind_rows() %>%
        arrange(gene_name, logFC) %>%
        group_by(symbol) %>%
        mutate(agg_pval = .stouffer_adjust(padj)) %>%
        ungroup()
    
    plt <- ggplot(data4plot, aes(gene_name, logFC)) +
        geom_violin(alpha = 0.5, fill = "grey80") +
        stat_summary(fun = median, geom = "crossbar", width = 0.8, size = 0.2) +
        geom_quasirandom(aes(fill = cell_line, shape = cell_line), size = 2) +
        geom_hline(yintercept = 0, lty = "dashed") +
        scale_y_continuous(limits = c(-1.8,1)) +
        scale_fill_manual(values = c("MGHU3" = "#1170aa", "RT112" = "#fc7d0b", "UMUC14" = "#5fa2ce")) +
        scale_shape_manual(values = c("MGHU3" = 24, "RT112" = 21, "UMUC14" = 24)) +
        theme_csg_scatter +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
        labs(x = "Gene", y = "log2FC after siFGFR3", fill = "Cell line")
    
    ggsave("results/fig5/fig5b.pdf", plt, width = 3.3, height = 2.2)
    return(plt)
}

.stouffer_adjust <- function(pvals, weights = NULL) {
    if(is.null(weights)) weights <- rep(1, length(pvals))
    zp <- (qnorm(pvals, lower.tail = FALSE) %*% weights) / sqrt(sum(weights^2))
    norm_p <- pnorm(zp, lower.tail = FALSE)[1]
    return(norm_p)
}

## siFGFR3 rescue cell viability --------------------
plot_rescue_siFGFR3 <- function(pcols = c("grey70", "grey30")) {
    #-- Read data
    umuc14_exp <- read_excel("data/experiments/TRAIL/201125-Pooled results-Rescue TRAIL with FGFR3 KD-ALL_Stat-and-SourceData.xlsx",
                             sheet = 2, skip = 2)
    colnames(umuc14_exp) <- c("TRAIL", 
                              paste0("siFGFR3n1_", 1:4),
                              paste0("siFGFR3n3_", 1:4))
    mghu3_exp <- read_excel("data/experiments/TRAIL/201125-Pooled results-Rescue TRAIL with FGFR3 KD-ALL_Stat-and-SourceData.xlsx",
                            sheet = 3, skip = 2)
    colnames(mghu3_exp) <- c("TRAIL", 
                             paste0("siFGFR3n1_", 1:3),
                             paste0("siFGFR3n3_", 1:3))
    
    #-- Join and reshape
    .pp_siRescue <- function(dat, cl) {
        dat %>%
            pivot_longer(cols = starts_with("siFGFR3"), 
                         names_to = "siFGFR3", 
                         values_to = "Relative_viability") %>%
            mutate(replicate = str_remove(siFGFR3, "siFGFR3n._") %>% as.numeric(),
                   siFGFR3 = str_replace(siFGFR3, "siFGFR3n", "siFGFR3#") %>% str_remove("_.*$"),
                   TRAIL = str_remove(TRAIL, "TRAIL "),
                   cell_line = cl)
    }
    exp_data <- bind_rows(.pp_siRescue(umuc14_exp, "UM-UC-14"), 
                          .pp_siRescue(mghu3_exp, "MGH-U3")) %>%
        mutate(cell_line = factor(cell_line, levels = c("UM-UC-14", "MGH-U3")),
               TRAIL = str_remove(TRAIL, " ng/ml"),
               siFGFR3 = case_when(
                   siFGFR3 == "siFGFR3#1" ~ "siFGFR3#I",
                   siFGFR3 == "siFGFR3#3" ~ "siFGFR3#II"
               ))
    
    #-- Stats
    exp_data_mean <- exp_data %>%
        group_by(TRAIL, siFGFR3, cell_line) %>%
        summarize(
            pval = t.test(Relative_viability, mu = 0, var.equal = TRUE)$p.value,
            Relative_viability = max(Relative_viability),
        ) %>%
        mutate(signif = ifelse(pval < 0.05, "*", "ns"))
    
    si_rescue_plot <-
        ggplot(exp_data, aes(TRAIL, Relative_viability, fill = siFGFR3)) +
        facet_wrap(~ cell_line, scales = "free_x") +
        stat_summary(geom = "crossbar", fun = "mean", size = 0.2, 
                     width = 0.8, position = position_dodge(width=1)) +
        geom_quasirandom(size = 3, shape = 21, dodge.width=1) +
        geom_text(aes(label = signif, vjust = -0.1), position=position_dodge(width=1),
                  data = exp_data_mean, size = 5) +
        geom_hline(yintercept = 0, lty = "dashed") +
        scale_fill_manual(values = pcols) +
        scale_y_continuous(limits = c(-2, 6), breaks = c(-1, 0, 2, 4, 6)) +
        labs(y = "Relative cell viability\nvs siCTL (log\u2082FC)", x = "TRAIL (ng/ml)", 
             fill = "") +
        theme_csg_scatter +
        theme(legend.position = 'bottom',
              panel.grid.major.y = element_line(linetype = "dashed", color = "grey50"))
    
    ggsave("results/fig5/fig5a.pdf", 
           plot = si_rescue_plot, width = 3.8, height = 2.7)
    
    return(si_rescue_plot)
}
.pp_siRescue <- function(dat, cl) {
    dat %>%
        pivot_longer(cols = starts_with("siFGFR3"), 
                     names_to = "siFGFR3", 
                     values_to = "Relative_viability") %>%
        mutate(replicate = str_remove(siFGFR3, "siFGFR3n._") %>% as.numeric(),
               siFGFR3 = str_replace(siFGFR3, "siFGFR3n", "siFGFR3#") %>% str_remove("_.*$"),
               TRAIL = str_remove(TRAIL, "TRAIL "),
               cell_line = cl)
}