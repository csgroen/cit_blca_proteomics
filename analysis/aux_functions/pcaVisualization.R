pcaVisualization <- function(pca_res, annotation, pcs = c(1, 2), 
                             color, shape = NA, id_name = "id", size = 2) {
    #-- Summary statistics (for variance)
    prop_var <- summary(pca_res)$importance[2,]
    prop_var_plt <- paste0(" (", format(prop_var * 100, digits = 2, trim = TRUE), "%)")
    
    #-- Plot table
    plot_df <- pca_res$x[,pcs] %>% as.data.frame() %>% rownames_to_column(id_name)
    annotation <- as.data.frame(annotation)
    #-- Make plot
    if (is.na(shape)) {
        plot_df <- left_join(plot_df, annotation[,c(id_name, color)], by = id_name)
        colnames(plot_df) <- c("id", "PCa", "PCb", "color")
        pca_plot <- ggplot(plot_df) +
            geom_point(aes(PCa, PCb, color = color), size = size)
    } else if (color != shape) {
        plot_df <- left_join(plot_df, annotation[,c(id_name, color, shape)], by = id_name)
        colnames(plot_df) <- c("id", "PCa", "PCb", "color", "shape")
        
        pca_plot <- ggplot(plot_df) +
            geom_point(aes(PCa, PCb, color = color, shape = shape), size = size)
    } else {
        plot_df <- left_join(plot_df, annotation[,c(id_name, color)], by = id_name)
        colnames(plot_df) <- c("id", "PCa", "PCb", "color")
        
        pca_plot <- ggplot(plot_df) +
            geom_point(aes(PCa, PCb, color = color, shape = color), size = size)
    }
    #-- Tweak styling
    pca_plot <- pca_plot +
        geom_hline(yintercept = 0, lty = "dashed", color = "grey50") +
        geom_vline(xintercept = 0, lty = "dashed", color = "grey50") +
        labs(x = paste0("PC ", pcs[1], prop_var_plt[pcs[1]]),
             y = paste0("PC ", pcs[2], prop_var_plt[pcs[2]]),
             color = color, shape = shape) +
        theme_csg_scatter +
        theme(panel.grid.minor = element_blank())
    return(pca_plot)
}
