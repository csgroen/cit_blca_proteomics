################################################################################
## TRAIL plots
################################################################################
se <- function(x) sqrt(var(x)/length(x))

plot_trailSensitivity <- function(cell_lines = c("MGH-U3", "UM-UC-14", "RT112", "RT4", "UM-UC-9", "SCaBER"),
                                  cl_colors =  ggthemes::tableau_color_pal("Color Blind")(6)[c(1,5,2,6,3,4)]) {
    names(cl_colors) <- cell_lines
    #-- Read raw data
    exp_data <- read_excel("data/experiments/TRAIL/Source data - 200720-Results_compilation-200629c - 200706a - 200715a.xlsx",
                           skip = 1)
    #-- Rename to pivot
    colnames(exp_data) <- c("TRAIL_concentration", 
                            paste("RT112", sep = "_", 1:3),
                            paste("MGH-U3", sep = "_", 1:3),
                            paste("UM-UC-14", sep = "_", 1:3),
                            paste("RT4", sep = "_", 1:3),
                            paste("SCaBER", sep = "_", 1:3),
                            paste("UM-UC-9", sep = "_", 1:3))
    #-- Pivot, add annotations
    exp_pivot <- exp_data %>%
        pivot_longer(cols = RT112_1:`UM-UC-9_3`, names_to = "replicate", values_to = "cell_viability") %>%
        mutate(cell_line = str_remove(replicate, "_.*")) %>%
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
        scale_shape_manual(values = c(21, 24, 25)) +
        labs(x = "TRAIL (ng/ml)", y = "Relative cell viability (%)", 
             shape = "FGFR3 status", fill = "Cell line") +
        theme_csg_scatter +
        theme(legend.position = "bottom", legend.box="vertical",
              panel.grid.major.y = element_line(linetype = "dashed", color = "grey50"))
    
    ggsave("results/fig5/fig5a.pdf", 
           plot = cell_viability_plot, width = 3.5, height = 4.5)
    return(cell_viability_plot)
}

plot_trailBirina <- function(cell_lines = c("UM-UC-14", "MGH-U3", "RT112"),
                             rel_widths = c(1,1,0.5)) {
    p1 <- .plot_trailBirina_lines(cell_lines)
    p2 <- .plot_trailBirina_synergy(cell_lines)
    leg1 <- get_legend(p1)
    leg2 <- get_legend(p2)
    p1 <- p1 + guides(fill = FALSE)
    p2 <- p2 + guides(fill = FALSE)
    leg <- plot_grid(leg1, leg2, ncol = 1)
    grid_plt <- plot_grid(p1, p2, leg, nrow = 1, align = "h", rel_widths = rel_widths)
    
    ggsave("results/fig5/fig5d.pdf", 
           plot = grid_plt, width = 6, height = 6)
    
    return(grid_plt)
    
}

.plot_trailBirina_lines <- function(cell_lines =  c("UM-UC-14", "MGH-U3", "RT112"),
                                    pal = "Reds") {
    exp_data <- read_excel("data/experiments/TRAIL/UMUC14_MGHU3-RT112-ZIP-data-source.xlsx", sheet = 2)
    colnames(exp_data) <- c("cell_line", "drug1", "drug2", "concentration_trail", "concentration_birina", "cell_viability", "synergy")
    
    exp_data <- exp_data %>%
        mutate(cell_viability = 100 - cell_viability,
               cell_line = factor(cell_line, levels = cell_lines))
    
    lplot <- exp_data %>%
        filter(concentration_birina > 0) %>%
        mutate(concentration_trail = format(concentration_trail, digits = 3)) %>%
        ggplot(aes(concentration_birina, cell_viability)) +
        facet_wrap(~ cell_line, ncol = 1) +
        geom_line(aes(group = concentration_trail), color = "grey50") +
        geom_point(aes(fill = concentration_trail), pch = 21, size = 2.5) +
        scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
        scale_fill_brewer(palette = pal) +
        # scale_fill_viridis_d(option = pal) +
        labs(x = "Birinapant (nM)", fill = "TRAIL (ng/ml)", y = "Relative cell viability (%)") +
        theme_csg_scatter +
        theme(panel.grid.major.y = element_line(linetype = "dashed", color = "grey50"))
    
    return(lplot)
}

.plot_trailBirina_synergy <- function(cell_lines =  c("UM-UC-14", "MGH-U3", "RT112"),
                                      pal = "plasma", contour = TRUE) {
    #-- Read data
    exp_data <- read_excel("data/experiments/TRAIL/UMUC14_MGHU3-RT112-ZIP-data-source.xlsx", sheet = 2)
    colnames(exp_data) <- c("cell_line", "drug1", "drug2", "concentration_trail", "concentration_birina", "cell_viability", "synergy")
    ZIP_scores = c("UM-UC-14" = 12.988, "MGH-U3" = 19.031, "RT112" = 13.921)
    
    exp_data <- exp_data %>%
        mutate(cell_viability = 100 - cell_viability,
               cell_line = factor(cell_line, levels = cell_lines))
    
    #-- Plot
    synergy2d <- .plot_synergy2d(exp_data, cell_lines, 
                                 labels = paste0("ZIP synergy score: ", ZIP_scores[cell_lines]), 
                                 pal = pal, contour = TRUE)
    
    return(synergy2d)
}

.plot_synergy2d <- function(exp_data, cell_lines, labels = NULL, pal = "plasma", contour = FALSE) {
    test <- 1
    # Get surfaces
    surfs <- lapply(cell_lines, function(cl) {
        #-- Get data
        cl_data <- exp_data %>%
            filter(cell_line == cl, concentration_birina > 0, concentration_trail > 0) %>%
            mutate(across(starts_with("concentration"), log10)) %>%
            as.data.frame()
        
        #-- Fit surface
        surf_df <- .fit3dsurface(cl_data, "concentration_birina", 
                                 "concentration_trail", "synergy") %>%
            mutate(cell_line = cl)
        return(surf_df)
    }) %>% bind_rows() %>%
        mutate(cell_line = factor(cell_line, levels = cell_lines))
    
    #-- Make labeller
    if(is.null(labels))
        facet_names <- as_labeller(structure(cell_lines, names = cell_lines))
    else
        facet_names <- as_labeller(structure(labels, names = cell_lines))
    
    #-- Plot 2D
    synergy2d <- 
        ggplot(surfs, aes(concentration_birina, concentration_trail, fill = synergy)) + 
        facet_wrap(~ cell_line, ncol = 1, labeller = facet_names) +
        geom_raster() + 
        scale_fill_viridis_c(option = pal) +
        scale_x_continuous(expand = c(0,0),
                           labels = math_format(10^.x),
                           name = "Birinapant (nM)") +
        scale_y_continuous(expand = c(0,0), 
                           breaks = c(0, 0.568636235841012, 
                                      1.04575749056067, 1.52287874528034, 2, 
                                      2.47712125471966, 2.95424250943932), 
                           labels = c("0","3.7","11.1","33.3","100","300","900"),
                           name = "TRAIL (ng/mL)", limits = c(log10(3.7), log10(900))) +
        labs(y = "TRAIL (ng/mL)", fill = "Synergy score") +
        theme_csg_sparse +
        theme(panel.grid = element_blank())
    
    if(contour) {
        synergy2d <- 
            synergy2d +
            geom_contour(aes(z = synergy), color = "white")
    }
    
    return(synergy2d)
}
dir.create("results/suppfig_apop/", showWarnings = FALSE)
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
    
    ggsave("results/suppfig_apop/a.pdf", 
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
        filter(Erdafitinib > 0)
    
    plt <- ggplot(exp_data2plot, aes(Erdafitinib, mean_viability, fill = TRAIL)) +
        facet_wrap(~ cell_line) +
        geom_line(color = "grey20") +
        geom_errorbar(aes(ymin = low_viability, ymax = high_viability), width = 0.1) +
        geom_point(size = 3, pch = 21) +
        scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
        scale_fill_manual(values = pcols) +
        labs(x = "Erdafitinib (μM)", y = "Relative cell viability (%)") +
        theme_csg_scatter +
        theme(panel.grid.major.y = element_line(linetype = "dashed", color = "grey50"))
    
    ggsave("results/fig5/fig5b.pdf", 
           plot = plt, width = 5, height = 2.5)
    return(plt)
}

################################################################################
## Check levels of FAS/TRAIL/TNF
################################################################################
plot_siFGFR3_apoptPath <- function() {
    #-- Read data
    clines <- readxl::excel_sheets("data/experiments/siFGFR3/200416-MGHU3-RT112-UMUC14-siF3-Log2FC-for-gseaPreRanked.xlsx")
    sifgfr3_fc <- lapply(clines, function(i) {
        read_excel("data/experiments/siFGFR3/200416-MGHU3-RT112-UMUC14-siF3-Log2FC-for-gseaPreRanked.xlsx", 
                   sheet = i) })
    clines <- str_remove(clines, " siFGFR3")
    names(sifgfr3_fc) <- clines
    
    #-- Get genes
    # receptors <- c("TRAIL-R1" = "TNFRSF10A", "TRAIL-R2" = "TNFRSF10B", "FAS" = "FAS", "TNF-R1" = "TNFRSF1A")
    genes <- c("TRAIL-R1" = "TNFRSF10A", "c-FLIP" = "CFLAR")
    
    #-- Plot
    plt <- lapply(clines, function(cl) { 
        fc <- sifgfr3_fc[[cl]]
        fc %>%
            filter(hgnc_symbol %in% genes) %>%
            mutate(cell_line = cl)
    }) %>%
        bind_rows() %>%
        mutate(log2FC = as.numeric(log2FC)) %>%
        arrange(hgnc_symbol, log2FC) %>%
        ggplot(aes(hgnc_symbol, log2FC)) +
        geom_violin(alpha = 0.5, fill = "grey80") +
        stat_summary(fun = median, geom = "crossbar", width = 0.8, size = 0.2) +
        geom_quasirandom(aes(fill = cell_line), size = 2, pch = 21) +
        geom_hline(yintercept = 0, lty = "dashed") +
        scale_y_continuous(limits = c(-1.8,1)) +
        scale_fill_manual(values = c("MGHU3" = "#1170aa", "RT112" = "#fc7d0b", "UMUC14" = "#5fa2ce")) +
        theme_csg_scatter +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
        labs(x = "Gene", y = "log2FC after siFGFR3", fill = "Cell line")
    
    ggsave("results/fig5/fig5c.pdf", plt, width = 3, height = 2.2)
    return(plt)
}