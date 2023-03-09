# Figure 6: Birinapant + TRAIL combo ---------------------------
dir.create("results/fig6/", showWarnings = FALSE, recursive = TRUE)

## Cell viability + synergy -------------------------
plot_trailBirina_synergy <- function(cell_lines = c("UM-UC-14", "MGH-U3", "RT112"),
                             rel_widths = c(1,1,0.5)) {
    p1 <- .plot_trailBirina_lines(cell_lines)
    p2 <- .plot_trailBirina_synergy(cell_lines)
    leg1 <- get_legend(p1)
    leg2 <- get_legend(p2)
    p1 <- p1 + guides(fill = FALSE)
    p2 <- p2 + guides(fill = FALSE)
    leg <- plot_grid(leg1, leg2, ncol = 1)
    grid_plt <- plot_grid(p1, p2, leg, nrow = 1, align = "h", rel_widths = rel_widths)
    
    ggsave("results/fig6/fig6b.pdf", 
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
    
    #-- Synergy plot
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



## Cell viability combo --------------------
plot_trailBirina_caspActivity  <- function(cell_lines = c("MGH-U3", "RT112", "RT4", "UM-UC-14")) {
    #-- Get caspase glow data
    casp_glo <- read_excel("data/experiments/TRAIL/200907a-200831a-200824-200818-Combi CaspActivity and Viab-Stat-and-SourceData.xlsx",
               sheet = 3, skip = 2)
    
    colnames(casp_glo) <- c("cell_line",
                            paste0("birina_", 1:4),
                            paste0("trail_10_", 1:4),
                            paste0("trail_30_", 1:4),
                            paste0("trail_100_", 1:4),
                            paste0("birina_trail_10_", 1:4),
                            paste0("birina_trail_30_", 1:4),
                            paste0("birina_trail_100_", 1:4))
    casp_glo <- casp_glo %>% filter(cell_line %in% cell_lines)
    
    casp_glo_tidy <- casp_glo %>%
        pivot_longer(cols = -cell_line, names_to = "treatment", values_to = "casp_act") %>%
        mutate(trail = case_when(
                   str_detect(treatment, "trail_10_") ~ 10,
                   str_detect(treatment, "trail_30_") ~ 30,
                   str_detect(treatment, "trail_100_") ~ 100,
                   TRUE ~ 0
               ),
               birina = ifelse(str_detect(treatment, "birina"), 100, 0),
               cell_line = factor(cell_line)) %>%
        relocate(trail, birina, .before = casp_act) %>%
        select(-treatment) %>%
        bind_rows(data.frame(cell_line = rep(casp_glo$cell_line, each = 4),
                             trail = 0,
                             birina = 0,
                             casp_act = 1)) %>%
        mutate(birina = factor(birina),
               trail = factor(trail))
    
    #-- Summarize for plot
    casp_glo_summ <- casp_glo_tidy %>%
        group_by(cell_line, trail, birina) %>%
        summarize(mean_casp_act = mean(casp_act, na.rm = TRUE), 
                  low_casp_act = mean(casp_act, na.rm = TRUE) - sd(casp_act, na.rm = TRUE),
                  high_casp_act = mean(casp_act, na.rm = TRUE) + sd(casp_act, na.rm = TRUE)) %>%
        ungroup()
    
    # T-test results
    ttests <- casp_glo_tidy %>%
        group_by(cell_line, trail) %>%
        summarize(ttest = list(t.test(casp_act ~ birina))) %>%
        mutate(ttest = map(ttest, tidy)) %>%
        unnest(cols = c(ttest)) %>%
        select(cell_line, trail, p.value) %>%
        mutate(psignif = case_when(
            p.value < 0.001 ~ "***",
            p.value < 0.01 ~ "**",
            p.value < 0.05 ~ "*",
            TRUE ~ ""
        ))
        
    #-- Plot
    casp_plt <- ggplot(casp_glo_summ, aes(trail, mean_casp_act, fill = birina)) +
        facet_wrap(~ cell_line, nrow = 1) +
        geom_errorbar(aes(ymin = low_casp_act, ymax = high_casp_act, group = birina), 
                      width = 0.3, position = position_dodge(width = 0.7)) +
        geom_col(alpha = 0.4, color = "black", position = "dodge", width = 0.7) +
        geom_point(aes(x = factor(trail), y = casp_act, fill = factor(birina)), pch = 21, size = 1.5, 
                   position = position_jitterdodge(), data = casp_glo_tidy) +
        stat_compare_means(aes(x = trail, y = casp_act), method = "t.test", label = "p.signif",
                           label.y = 12, size = 4, data = casp_glo_tidy) +
        scale_y_continuous(limits = c(0, 14), expand = c(0,0), breaks = c(1,4,7,10)) +
        scale_fill_manual(values = c("#d7cefe", "#3f0595")) +
        labs(x = "TRAIL (ng/mL)", y = "Relative Caspase 3/7 activity",
             fill = "Birinapant (nM)") +
        theme_csg_scatter +
        theme(panel.grid.major.y = element_line(linetype = "dashed", color = "grey50"))

    
    # ggsave("results/fig6/fig6a1.pdf", casp_plt, width = 7, height = 2)
    return(casp_plt)
    
}

plot_trailBirina_cellViability  <- function(cell_lines = c("MGH-U3", "RT112", "RT4", "UM-UC-14")) {
  #-- Get caspase glow data
  cv_glo <- read_excel("data/experiments/TRAIL/200907a-200831a-200824-200818-Combi CaspActivity and Viab-Stat-and-SourceData.xlsx",
                         sheet = 1, skip = 2)
  
  colnames(cv_glo) <- c("cell_line",
                          paste0("birina_", 1:4),
                          paste0("trail_10_", 1:4),
                          paste0("trail_30_", 1:4),
                          paste0("trail_100_", 1:4),
                          paste0("birina_trail_10_", 1:4),
                          paste0("birina_trail_30_", 1:4),
                          paste0("birina_trail_100_", 1:4))
  cv_glo <- cv_glo %>% filter(cell_line %in% cell_lines)
  
  cv_glo_tidy <- cv_glo %>%
    pivot_longer(cols = -cell_line, names_to = "treatment", values_to = "cell_viability") %>%
    mutate(trail = case_when(
      str_detect(treatment, "trail_10_") ~ 10,
      str_detect(treatment, "trail_30_") ~ 30,
      str_detect(treatment, "trail_100_") ~ 100,
      TRUE ~ 0
    ),
    birina = ifelse(str_detect(treatment, "birina"), 100, 0),
    cell_line = factor(cell_line)) %>%
    relocate(trail, birina, .before = cell_viability) %>%
    select(-treatment) %>%
    bind_rows(data.frame(cell_line = rep(cv_glo$cell_line, each = 4),
                         trail = 0,
                         birina = 0,
                         cell_viability = 100)) %>%
    mutate(birina = factor(birina),
           trail = factor(trail))
  
  #-- Summarize for plot
  cv_glo_summ <- cv_glo_tidy %>%
    group_by(cell_line, trail, birina) %>%
    summarize(mean_cell_viability = mean(cell_viability, na.rm = TRUE), 
              low_cell_viability = mean(cell_viability, na.rm = TRUE) - sd(cell_viability, na.rm = TRUE),
              high_cell_viability = mean(cell_viability, na.rm = TRUE) + sd(cell_viability, na.rm = TRUE)) %>%
    ungroup()
  
  # Welch results
  ttests <- cv_glo_tidy %>%
    group_by(cell_line, trail) %>%
    summarize(ttest = list(t.test(cell_viability ~ birina))) %>%
    mutate(ttest = map(ttest, tidy)) %>%
    unnest(cols = c(ttest)) %>%
    select(cell_line, trail, p.value) %>%
    mutate(psignif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ))
  
  #-- Plot
  cv_plt <- 
    ggplot(cv_glo_summ, aes(trail, mean_cell_viability, fill = birina)) +
    facet_wrap(~ cell_line, nrow = 1) +
    geom_errorbar(aes(ymin = low_cell_viability, ymax = high_cell_viability, group = birina), 
                  width = 0.3, position = position_dodge(width = 0.7)) +
    geom_col(alpha = 0.4, color = "black", position = "dodge", width = 0.7) +
    geom_point(aes(x = factor(trail), y = cell_viability, fill = factor(birina)), pch = 21, size = 1.5, 
               position = position_jitterdodge(), data = cv_glo_tidy) +
    stat_compare_means(aes(x = trail, y = cell_viability), method = "t.test", label = "p.signif",
                       label.y = 110, size = 4, data = cv_glo_tidy) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 120), breaks = c(0,25,50,75,100)) +
    scale_fill_manual(values = c("#d7cefe", "#3f0595")) +
    labs(x = "TRAIL (ng/mL)", y = "Relative cell viability (%)",
         fill = "Birinapant (nM)") +
    theme_csg_scatter +
    theme(panel.grid.major.y = element_line(linetype = "dashed", color = "grey50"))
  # ggsave("results/fig6/fig6a2.pdf", cv_plt, width = 7.1, height = 2)
  return(cv_plt)
}

plot_trailBirina <- function(cell_lines = c("MGH-U3", "RT112", "RT4", "UM-UC-14")) {
  casp_plt <- plot_trailBirina_caspActivity(cell_lines)
  cv_plt <- plot_trailBirina_cellViability(cell_lines)
  
  
  plt_trail_birina <- (casp_plt + theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank())) / 
    (cv_plt + theme(strip.background = element_blank(),
                    strip.text = element_blank())) +
    plot_layout(guides = "collect")
  
  ggsave("results/fig6/fig6a.pdf", plt_trail_birina, width = 7, height = 4)
  
}
