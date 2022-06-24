require(ggalluvial)
#' Sankey diagram / parallel plot using ggalluvial
#' 
#' Takes a tidy `sampAnnot` with at least three columns: `classv1`, `classv2`
#' and `id_var`, where `classv1` is a categorical variable used as the first
#' stratum and `classv2` is another categorical variable used as the second stratum.
#' 
#' @param sampAnnot a tidy table with at least 3 columns, described below.
#' @param classv1 a character with the name of the column containing the
#' categorical variable used as the first stratum. Origin point of the "lode".
#' @param classv2 a character with the name of the column containing the
#' categorical variable used as the second stratum. Destination point of the "lode".
#' @param id_var a character variable with the name of the column containing unique
#' ids for each observation
#' @param group_var optional character variable to use for facetting.
#' @param flip a boolean. If `flip = TRUE`, the strata will be arranged vertically.
#' @param alpha controls lode transparecy. See: {ggalluvial::geom_lode}
#' @param label_colors an optional named vector of colors for the strata and lodes
#' @param label_type a character, one of: 'label', 'text' or 'none'. If 'label',
#' `geom_label` is used for strata names, if 'text', `geom_text` is used. If 'none',
#' strata are left unlabelled.
#' @param label_size a number, the size of the text label of the strata
#' @param label_fill a color, used for the background of the label if `label_type == "label"`
#' @param save a boolean, whether to save or not.
#' @param fpath a character with the path of the output directory
#' @param width width of the output pdf
#' @param height height of the output pdf
#' 
require(ggalluvial)
classFlow <- function(sampAnnot, classv1, classv2, id_var = "id", 
                      group_var = NULL,
                      flip = FALSE,
                      alpha = 0.1, 
                      label_colors = NULL,
                      label_type = c("label", "text", "none")[1], 
                      label_size = 3.5, 
                      label_fill = "white",
                      text_color = "black", 
                      save = FALSE, fpath = "./",
                      width = 6, height = 5) {
    if(is.null(label_colors)) {
        lv1 <- unique(sampAnnot[,classv1]) %>% unlist() %>% sort() %>% as.character()
        lv2 <- unique(sampAnnot[,classv2]) %>% unlist() %>% as.character()
        label_colors <- c(structure(ggsci::pal_npg()(length(lv1)), names = lv1),
                          structure(ggsci::pal_jco()(length(lv2)), names = lv2))
    }
    
    if(is.null(group_var)) {
        cls_v1 <- pull(sampAnnot, !!classv1)
        cls_v2 <- pull(sampAnnot, !!classv2)
        if(is.factor(cls_v1)) {
            cls_v1_lvs <- levels(cls_v1)
        } else {
            cls_v1_lvs <- unique(cls_v1)
        }
        if(is.factor(cls_v2)) {
            cls_v2_lvs <- levels(cls_v2)
        } else {
            cls_v2_lvs <- unique(cls_v2)
        }
        sa4plot <-  sampAnnot %>%
            dplyr::select(all_of(id_var), all_of(classv1), all_of(classv2)) %>%
            dplyr::rename(id = !! id_var) %>%
            pivot_longer(cols = c(classv1, classv2),
                         names_to = "Class_Type", values_to = "Class_Label") 
        if(flip) {
            sa4plot <- sa4plot %>%
                mutate(Class_Type = factor(Class_Type, levels = c(classv2, classv1)),
                       Class_Label = factor(Class_Label, levels = c(rev(cls_v1_lvs), rev(cls_v2_lvs))))
        } else {
            sa4plot <- sa4plot %>%
                mutate(Class_Type = factor(Class_Type, levels = c(classv1, classv2)),
                       Class_Label = factor(Class_Label, levels = c(cls_v1_lvs, cls_v2_lvs)))
        }
        
        #-- ggalluvial (long lode)
        aluv_plot <- ggplot(sa4plot, aes(x = Class_Type, stratum = Class_Label, 
                                         alluvium = id, label = Class_Label)) +
            geom_flow(aes(fill = Class_Label), stat = "alluvium", 
                      lode.guidance = "frontback", aes.flow = ifelse(flip, "backward", "forward"),
                      alpha = alpha) +
            geom_stratum(aes(fill = Class_Label)) +
            scale_y_continuous(expand = c(0,0)) +
            scale_color_manual(values = label_colors) +
            scale_fill_manual(values = label_colors) +
            scale_x_discrete(expand = c(0.1,0)) +
            theme_light() +
            theme(axis.text = element_text(size = 10, color = "black"),
                  axis.line.y = element_line(color = "black"),
                  axis.ticks.y = element_line(color = "black"),
                  legend.position = "bottom",
                  panel.grid = element_blank(),
                  panel.border = element_blank())
        
    } else {
        sa4plot <-  sampAnnot %>%
            dplyr::select(all_of(id_var), all_of(group_var), all_of(classv1), all_of(classv2)) %>%
            dplyr::rename(id = !! id_var, group = !! group_var) %>%
            pivot_longer(cols = c(classv1, classv2),
                         names_to = "Class_Type", values_to = "Class_Label") %>%
            mutate(Class_Type = case_when(
                flip ~ factor(Class_Type, levels = c(classv2, classv1)),
                TRUE ~ factor(Class_Type, levels = c(classv1, classv2))))
        
        #-- ggalluvial (long lode)
        aluv_plot <- ggplot(sa4plot, aes(x = Class_Type, stratum = Class_Label, alluvium = id, label = Class_Label)) +
            facet_wrap(~ group, scales = "free_y") +
            geom_flow(aes(fill = Class_Label, color = Class_Label), 
                      stat = "alluvium", lode.guidance = "frontback",
                      aes.flow = ifelse(flip, "backward", "forward"),
                      alpha = alpha) +
            geom_stratum(aes(fill = Class_Label)) +
            scale_x_discrete(expand = c(0.1,0.1)) +
            scale_y_continuous(expand = c(0,0)) +
            scale_color_manual(values = label_colors) +
            scale_fill_manual(values = label_colors) +
            theme_light() +
            theme(axis.text = element_text(size = 10, color = "black"),
                  axis.line.y = element_line(color = "black"),
                  axis.ticks.y = element_line(color = "black"),
                  legend.position = "bottom",
                  panel.grid = element_blank(),
                  panel.border = element_blank(), 
                  strip.text = element_text(color = "black"))
    }
    if(label_type == "label") {
        aluv_plot <- aluv_plot + geom_label(stat = "stratum", size = label_size, 
                                            fill = label_fill, color = text_color)
    } else if (label_type == "text") {
        aluv_plot <- aluv_plot + geom_text(stat = "stratum", size = label_size, color = text_color)
    }
    if(flip) {
        aluv_plot <- aluv_plot + coord_flip()
    }
    
    if(save) {
        if(is.null(group_var)) fname <- paste0(fpath, "ClassCompare_", classv1, "_", classv2, ".pdf")
        else fname <- paste0(fpath, "ClassCompare_", classv1, "_", classv2, "_by_", group_var, ".pdf")
        ggsave(filename = fname, aluv_plot, width = width, height = height)
    }
    return(aluv_plot)
    
}