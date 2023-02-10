################################################################################
## TRAIL plots
################################################################################
se <- function(x) sqrt(var(x)/length(x))

plot_trailSensitivity <- function(cell_lines = c("MGH-U3", "UM-UC-14", "RT112", "RT4", "UM-UC-9"),
                                  # cell_lines = c("MGH-U3", "UM-UC-14", "RT112", "RT4", "UM-UC-9", "SCaBER"),
                                  # cl_colors =  ggthemes::tableau_color_pal("Color Blind")(6)[c(1,5,2,6,3,4)],
                                  fig_height = 4,
                                  cl_colors =  ggthemes::tableau_color_pal("Color Blind")(6)[c(1,5,2,6,3)]
                                  ) {
    names(cl_colors) <- cell_lines
    #-- Read raw data
    exp_data <- read_excel("data/experiments/TRAIL/Source data - 200720-Results_compilation-200629c - 200706a - 200715a.xlsx",
                           skip = 1)
    #-- Rename to pivot
    cl_cols <- lapply(c("RT112", "MGH-U3", "UM-UC-14", "RT4", "SCaBER", "UM-UC-9"), 
                      function(cl) { paste(cl, sep = "_", 1:3)}) %>% unlist()
    colnames(exp_data) <- c("TRAIL_concentration", cl_cols)
    #-- Pivot, add annotations
    exp_pivot <- exp_data %>%
        pivot_longer(cols = cl_cols, 
                     names_to = "replicate", values_to = "cell_viability") %>%
        mutate(cell_line = str_remove(replicate, "_.*")) %>%
        filter(cell_line %in% cell_lines) %>%
        group_by(TRAIL_concentration, cell_line) %>%
        summarize(mean_cell_viability = mean(cell_viability), 
                  cell_viability_se = se(cell_viability)) %>%
        mutate(FGFR3_alteration = case_when(
            cell_line %in% c("UM-UC-14", "MGH-U3") ~ "Mutation",
            cell_line %in% c("RT112", "RT4") ~ "Fusion",
            TRUE ~ "WT"),
            # cell_line = factor(cell_lines, levels = cell_lines),
            # TRAIL_concentration = ifelse(TRAIL_concentration <= 0, 0, TRAIL_concentration),
            cell_viability_high = mean_cell_viability + cell_viability_se,
            cell_viability_low = mean_cell_viability - cell_viability_se)
    
    cell_viability_plot <-
        ggplot(exp_pivot, aes(TRAIL_concentration, mean_cell_viability, shape = FGFR3_alteration)) +
        geom_line(aes(group = cell_line), color = "grey50") +
        geom_errorbar(aes(ymin = cell_viability_low, ymax = cell_viability_high), width = 0.1) +
        geom_point(aes(fill = cell_line), size = 3) +
        scale_x_continuous(trans=pseudo_log_trans(base = 10), breaks = c(0, 10, 100, 1000)) +
        scale_fill_manual(values = cl_colors, guide = guide_legend(override.aes = list(pch = 22))) +
        scale_shape_manual(values = c(21, 24, 22)) +
        labs(x = "TRAIL (ng/ml)", y = "Relative cell viability (%)", 
             shape = "FGFR3 status", fill = "Cell line") +
        theme_csg_scatter +
        theme(legend.position = "bottom", legend.box="vertical",
              panel.grid.major.y = element_line(linetype = "dashed", color = "grey50"))
    
    ggsave("results/fig4/fig4a.pdf", 
           plot = cell_viability_plot, width = 3.5, height = fig_height)
    return(cell_viability_plot)
}
dir.create("results/suppfig_apop/", showWarnings = FALSE)

plot_rescue_erdafitinib <- function(pcols = brewer.pal(3, "Reds")) {
    #-- Read data
    umuc14_data <- read_excel("data/experiments/TRAIL/200831b-200907b-200914a- UMUC14-MGHU3-treated-ERDAandTRAIL-CellTiterGlo48h_SourceData.xlsx",
                              skip = 2)
    colnames(umuc14_data) <- c("Erdafitinib",
                               paste0("Control_", 1:3),
                               paste0("TRAIL_30_", 1:3))
    
    mghu3_data <- read_excel("data/experiments/TRAIL/200831b-200907b-200914a- UMUC14-MGHU3-treated-ERDAandTRAIL-CellTiterGlo48h_SourceData.xlsx",
                             skip = 2, sheet = 2)
    colnames(mghu3_data) <- c("Erdafitinib",
                              paste0("Control_", 1:3),
                              paste0("TRAIL_300_", 1:2))
    
    umuc14_data <- umuc14_data %>%
        select(Erdafitinib, starts_with("Control"), starts_with("TRAIL_30")) %>%
        pivot_longer(cols = matches("Control|TRAIL"), 
                     names_to = "TRAIL", 
                     values_to = "Relative_viability") %>%
        mutate(TRAIL = ifelse(str_detect(TRAIL, "TRAIL"), "30 ng/ml", "0 ng/ml"),
               cell_line = "UM-UC-14")
    
    mghu3_data <- mghu3_data %>%
        select(Erdafitinib, starts_with("Control"), starts_with("TRAIL_300")) %>%
        pivot_longer(cols = matches("Control|TRAIL"), 
                     names_to = "TRAIL", 
                     values_to = "Relative_viability") %>%
        mutate(TRAIL = ifelse(str_detect(TRAIL, "TRAIL"), "300 ng/ml", "0 ng/ml"),
               cell_line = "MGH-U3")
    
    exp_data <- bind_rows(umuc14_data, mghu3_data) %>%
        mutate(cell_line = factor(cell_line, levels = c("UM-UC-14", "MGH-U3"))) 
    
    exp_data2plot <- exp_data %>%
        group_by(cell_line, Erdafitinib, TRAIL) %>%
        summarize(mean_viability = mean(Relative_viability),
                  low_viability = mean(Relative_viability) - se(Relative_viability),
                  high_viability = mean(Relative_viability) + se(Relative_viability)) %>%
        filter(Erdafitinib > 0) %>%
      mutate(Erdafitinib = Erdafitinib * 100)
    
    plt <- ggplot(exp_data2plot, aes(Erdafitinib, mean_viability, fill = TRAIL)) +
        facet_wrap(~ cell_line) +
        geom_line(color = "grey20") +
        geom_errorbar(aes(ymin = low_viability, ymax = high_viability), width = 0.1) +
        geom_point(size = 3, pch = 21) +
        scale_x_log10() +
        scale_fill_manual(values = pcols) +
        labs(x = "Erdafitinib (nM)", y = "Cell viability (% relative to untreated)") +
        theme_csg_scatter +
        theme(panel.grid.major.y = element_line(linetype = "dashed", color = "grey50"))
    
    ggsave("results/fig4/fig4b.pdf", 
           plot = plt, width = 5, height = 2.5)
    return(plt)
}

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
    
    ggsave("results/fig4/fig4d.pdf", plt, width = 3.3, height = 2.2)
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
  
  ggsave("results/fig5/fig5c.pdf", 
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
