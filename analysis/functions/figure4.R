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
                  cell_viability_se = se(cell_viability),
                  cell_ci95_high = .ci95(cell_viability, bound = "upper"),
                  cell_ci95_low = .ci95(cell_viability, bound = "lower")) %>%
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
        geom_errorbar(aes(ymin = cell_ci95_low, ymax = cell_ci95_high), width = 0.1) +
        geom_point(aes(fill = cell_line), size = 3) +
        scale_x_continuous(trans=pseudo_log_trans(base = 10), breaks = c(0, 10, 100, 1000)) +
        scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
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
                high_viability = mean(Relative_viability) + se(Relative_viability),
                low_viability_95ci = .ci95(Relative_viability, "lower"),
                high_viability_95ci = .ci95(Relative_viability, "upper")
                ) %>%
        filter(Erdafitinib > 0) %>%
      mutate(Erdafitinib = Erdafitinib * 100)
    
    plt <- ggplot(exp_data2plot, aes(Erdafitinib, mean_viability, fill = TRAIL)) +
        facet_wrap(~ cell_line) +
        geom_line(color = "grey20") +
        geom_errorbar(aes(ymin = low_viability_95ci, ymax = high_viability_95ci), width = 0.1) +
        geom_point(size = 3, pch = 21) +
        scale_x_log10() +
        scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
        scale_fill_manual(values = pcols) +
        labs(x = "Erdafitinib (nM)", y = "Cell viability (% relative to untreated)") +
        theme_csg_scatter +
        theme(panel.grid.major.y = element_line(linetype = "dashed", color = "grey50"))
    
    ggsave("results/fig4/fig4b.pdf", 
           plot = plt, width = 5, height = 2.5)
    return(plt)
}



.ci95 <- function(vec, bound) {
    t_score <- qt(p=0.05/2, df = length(vec)-1, lower.tail = FALSE)
    margin_error <- t_score * se(vec)
    if(bound == "lower")
    {
        return(mean(vec) - margin_error)
    } else {
        return(mean(vec) + margin_error)
    }
}
