# Supp Figure 4
################################################################################
## MOFA-related functions
################################################################################
runMOFA <- function(mrna_data, prot_data, cna, core_samples,
                    n_factors, n_transcripts, n_cna,
                    cache = TRUE) {
    dir.create("results/mofa", showWarnings = FALSE)
    base_name <- str_glue("mofa_mrna{n_transcripts}_cna{n_cna}_{n_factors}f")
    dir.create(str_glue("results/mofa/{base_name}/"), showWarnings = FALSE)
    fname <- paste0("cached_results/mofa/", base_name, "_trained.RData")
    euler_plot <- str_glue("results/mofa/{base_name}/eulerFeatures.pdf")
    
    #-- Feature selection
    gene_cv <- apply(mrna_data,1,cv)
    idx_mrna <- order(gene_cv, decreasing = TRUE)[1:n_transcripts]
    
    cna_cv <- apply(cna,1,cv)
    idx_cna <- order(cna_cv, decreasing = TRUE)[1:n_cna]
    idx_cna <- rownames(cna)[idx_cna]
    
    #-- Get data
    cna_mat <- as.matrix(cna[idx_cna,core_samples])
    mrna_mat <- as.matrix(mrna_data[idx_mrna,core_samples])
    prot_mat <- as.matrix(prot_data[,core_samples])
    
    #-- Plot data overview
    df_tr <- tibble(feature = rownames(mrna_mat), transcriptomics = 1)
    df_pr <- tibble(feature = rownames(prot_data), proteomics = 1)
    df_cna <- tibble(feature = rownames(cna_mat), cna = 1)
    
    all_feats <- df_tr %>%
        full_join(df_pr) %>%
        full_join(df_cna) %>%
        mutate(across(-feature, ~ ifelse(is.na(.), 0, .))) %>%
        filter(!is.na(feature)) %>%
        as.data.frame() %>%
        column_to_rownames("feature")
    
    feature_plot <- euler(all_feats) %>% plot(quantities = TRUE,
                                              fills = c("#8da0cb", "#66c2a5", "#fc8d62"))
    pdf(euler_plot, width = 4, height = 4)
    Sys.sleep(1)
    feature_plot
    dev.off()
    
    if(file.exists(fname) & cache) {
        load(fname)
    } else {
        #-- Create object
        rownames(cna_mat) <- paste0(rownames(cna_mat), "_cna")
        rownames(mrna_mat) <- paste0(rownames(mrna_mat), "_mrna")
        rownames(prot_mat) <- paste0(rownames(prot_mat), "_prot")
        
        data_list <- list(
            cna = cna_mat,
            mrna = mrna_mat,
            prot = prot_mat)
        
        mofa_obj <- create_mofa(data_list)
        
        #-- Model options
        data_opts <- get_default_data_options(mofa_obj)
        model_opts <- get_default_model_options(mofa_obj)
        model_opts$num_factors <- n_factors
        train_opts <- get_default_training_options(mofa_obj)
        train_opts$convergence_mode <- "slow"
        
        #-- Train MOFA
        mofa_obj <- prepare_mofa(mofa_obj,
                                 data_options = data_opts,
                                 model_options = model_opts,
                                 training_options = train_opts)
        mofa_trained <- run_mofa(mofa_obj, use_basilisk = TRUE)
        save(mofa_trained, file = fname)
    }
    return(list(mofa_obj = mofa_trained,
                base_name = base_name))
}

conformMOFAres <- function(mofa_res) {
    obj <- seq(1,length(mofa_res),by=2)
    name <- seq(2,length(mofa_res),by=2)
    mofa_res2 <- lapply(1:(length(mofa_res)/2), 
                        function(i) { list(mofa_obj = mofa_res[[obj[i]]], 
                                           base_name = mofa_res[[name[i]]])})
    return(mofa_res2)
    
}

characterizeMOFA <- function(mofa_res_list, sampAnnot = sampAnnot) {
    mofa_interps <- lapply(mofa_res_list, function(mofa_res) {
        mofa_trained <- mofa_res$mofa_obj
        base_name <- mofa_res$base_name
        
        samples_metadata(mofa_trained) <- sampAnnot
        
        mofa_qc <- .plot_mofa_QC(mofa_trained, base_name)
        mofa_assocs <- .plot_mofa_associations(mofa_trained, base_name)
        mofa_hm <- .plot_mofa_factorHeatmap(mofa_trained, base_name)
        
        mofa_paths_data <- .plot_mofa_factorPathways(mofa_trained, base_name)
        
        mofa_top_genes <- .plot_mofa_topGenes(mofa_trained, base_name)
        
        mofa_chr <- list(qc = mofa_qc,
                         associations = mofa_assocs,
                         factor_heatmap = mofa_hm,
                         factor_pathways = mofa_paths_data,
                         factor_top_genes = mofa_top_genes)
        return(mofa_chr)

    })
    names(mofa_interps) <- sapply(mofa_interps, "[[", "base_name")
    return(mofa_interps)
    
}

finalizeMOFA <- function(mofa_select,
                         sampAnnot, 
                         factors2show, 
                         factor_names, 
                         invert) {
    fname <- str_glue("cached_results/mofa/{mofa_select}_trained.RData")
    load(fname)
    
    samples_metadata(mofa_trained) <- sampAnnot
    
    list(mofa = mofa_trained,
         facts = factors2show,
         fact_names = factor_names,
         invert = invert)
}

.mofa_allVarExp <- function(final_mofa) {
    var_exp <- MOFA2::get_variance_explained(final_mofa$mofa)[[2]][[1]]
    rownames(var_exp) <- plyr::mapvalues(rownames(var_exp), final_mofa$facts, final_mofa$fact_names)
    
    var_exp %>%
        as.data.frame() %>%
        rownames_to_column("factor_names") %>%
        pivot_longer(cols = cna:prot, names_to = "data_level", values_to = "var_exp") %>%
        mutate(factor_names = fct_inorder(factor_names)) %>%
        ggplot(aes(factor_names, data_level, fill = var_exp)) +
        geom_tile() +
        scale_fill_distiller(palette = "Blues", direction = 1) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_csg_sparse +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              panel.border = element_rect(color = "black", fill = "transparent")) +
        labs(x = "Factor names", y = "Data level", fill = "Variance\nexplained (%)")
}


.mofa_associations <- function(final_mofa) {
    mofa_trained <- final_mofa$mofa
    factors2show <- paste0("Factor",1:10)
    factors_select <- final_mofa$facts
    factor_names <- structure(plyr::mapvalues(factors2show, 
                                              from = factors_select, 
                                              to = final_mofa$fact_names),
                              names = factors2show)
    
    invert <- structure(rep(1,10), names =  factors2show)
    invert[factors_select] <- final_mofa$invert
    
    #-- Class
    assoc_sty <- associate_factors_with_covariates(mofa_trained, covariates = "Subtype",
                                                   fct_names = factors_select, 
                                                   invert = final_mofa$invert,
                                                   color_limits = c(-2,2),
                                                   size_limits = c(0,8)) +
        scale_y_discrete(labels = rev(factor_names)) +
        theme_mofa
    assoc_upg <- associate_factors_with_covariates(mofa_trained, covariates = "uPG",
                                                   fct_names = factors_select, 
                                                   invert = final_mofa$invert,
                                                   color_limits = c(-2,2),
                                                   size_limits = c(0,8)) +
        scale_y_discrete(labels = rev(factor_names)) +
        theme_mofa
    
    cl_assocs <-  ((assoc_upg + labs(x = "uPG")) | (assoc_sty 
                                                    + labs(x = "Subtype") +
                                                        theme(axis.text.y=element_blank(),
                                                              axis.title.y=element_blank(),
                                                              axis.ticks.y = element_blank()))) + 
        plot_layout(guides = "collect", widths = c(1,2))
    #-- Genomic alterations
    facI_mut <- associate_factors_with_covariates(mofa_trained,
                                                  fct_names = factors2show, 
                                                  invert = invert,
                                                  c("FGFR3_mutation", "TP53_mutation",
                                                    "PIK3CA_mutation", "RAS_mutation"),
                                                  color_limits = c(-1.5,1.5),
                                                  size_limits = c(0,6)) +
        theme_mofa +
        scale_y_discrete(labels = rev(factor_names))
    chr_alts <- c("9q31_1_del", "11p15_3_del", 
                  "1p36_32_del", "13q13_1_amp", 
                  "4q35_2_del", "6p22_3_amp", 
                  "18p11_31_amp", "3p14_2_del", 
                  "11q13_3_amp", "3p25_2_amp", 
                  "8p22_del", "1p12_amp", 
                  "1q23_3_amp")
    chr_labs <- chr_alts %>% 
        str_replace("_([0-9])", ".\\1") %>%
        str_replace_all("_", " ")
    facI_chr <- associate_factors_with_covariates(mofa_trained, 
                                                  covariates =  chr_alts,
                                                  invert = invert, 
                                                  cluster_cols = TRUE,
                                                  color_limits = c(-1.5,1.5),
                                                  size_limits = c(0,6)) +
        theme_mofa +
        scale_x_discrete(labels = structure(chr_labs, names = chr_alts)) +
        scale_y_discrete(labels = rev(factor_names)) +
        labs(x = "CNA")
    
    cna_genes <- c("TP53_loss_bin", "ERBB2_gain_bin",
                   "RB1_loss", "ELF3_loss_bin", "CDKN2A_loss", "PPARG_gain")
    cna_names <- cna_genes %>% str_replace_all("_", " ") %>% str_remove(" bin")
    facI_cnaGenes <- associate_factors_with_covariates(mofa_trained, 
                                                       cna_genes, 
                                                       fct_names = factors2show, 
                                                       invert = invert,
                                                       cluster_cols = TRUE,
                                                       color_limits = c(-1.5,1.5),
                                                       size_limits = c(0,6)) +
        theme_mofa +
        scale_x_discrete(labels = structure(cna_names, names = cna_genes)) +
        scale_y_discrete(labels = factor_names) +
        labs(x = "CNA")
    
    plt_genAlt <- ((((facI_mut + labs(x = "Mutation")) |
                         (facI_cnaGenes + labs(x = "CNA") +
                              theme(axis.text.y = element_blank(),
                                    axis.title.y = element_blank(),
                                    axis.ticks.y = element_blank()))) +
                        plot_layout(widths = c(2,3))) /
                       facI_chr + labs(x = "CNA")) +
        plot_layout(guides = "collect")
    plt_genAlt
    
    #-- trICs 
    trICs <- c("IC Basal-like", "IC Mitochondria", "IC Cell cycle", "IC Urothelial differentiation", 
               "IC Bladder cancer pathways", "IC Neuroendocrine", "IC Smooth muscle", 
               "IC Lymphocytes B&T", "IC Myofibroblasts", "IC Interferon response"
    )
    plt_trIC <- associate_factors_with_covariates(mofa_trained, covariates = trICs,
                                                  invert = invert,
                                                  color_limits = c(-1,1),
                                                  size_limits = c(0,6),
                                                  color_low = "#5e3c99",
                                                  color_high = "#e66101") +
        theme_mofa +
        scale_y_discrete(labels = factor_names) +
        labs(x = "Transcriptomic ICs")
    
    all_associations <- (plt_genAlt / plt_trIC)
    
    return(list(subtype = cl_assocs, others = all_associations))
}

.mofa_finalHM <- function(final_mofa, mcp_res, sa, wp_symbol) {
    mofa_trained <- final_mofa$mofa
    factors2show <- final_mofa$facts
    factor_names <- final_mofa$fact_names
    invert <- final_mofa$invert
    
    #-- Get factor HM
    fct_values <- MOFA2::get_factors(mofa_trained)[[1]]
    fct_mat <- sapply(1:length(factors2show), function(i) {
        fct <- factors2show[i]; inv <- invert[i]
        vals <- fct_values[,fct] * inv
        return(vals)
    })
    colnames(fct_mat) <- factor_names
    pops <- c("T cells", "Cytotoxic lymphocytes", "B lineage", "Fibroblasts")
    
    mks <- list("Basal phenotype" = c("KRT6A", "KRT14", "S100A8", "TYMP", "KRT5"),
                "Cell cycle/DNA repair" = c("TOP2A", "TMPO", "FEN1", "UHRF1", "KPNA2"),
                "Lipid metabolism" = c("AKR1C3", "GDPD3", "ADIRF", "ACOX1", "ABCD3"),
                "ECM interactions/Cytoskeleton/Smooth muscle" =  c("FLNC", "TPM2", "TAGLN", 
                                                                   "TPM1", "CNN1"),
                "Ta pathway/Urothelial differentiation" = c("ANXA10", "S100P",
                                                            "IVL", "KRT7", "DHRS2"))
    mks_unlist <- as.character(unlist(mks))
    
    # mks <- c("KRT6A", "KRT14", "SERPINE1", "TOP2A", "MKI67", "AKR1C3", "GDPD3", 
    #          "TAGLN", "DMD", "S100P", "KRT7", "ADIRF")
    
    data4hm <- fct_mat %>%
        as.data.frame() %>%
        rownames_to_column("sample") %>%
        left_join(sa) %>%
        left_join(rownames_to_column(mcp_res, "sample")) %>%
        left_join(rownames_to_column(as.data.frame(t(wp_symbol[mks_unlist,])), "sample"))
    
    mofa_hm <- .mofa_factor_hm(data4hm, factor_names, pops, mks)
    var_exp <- .variance_exp_plot(mofa_trained, mofa_hm, factors2show, factor_names)
    complete_hm <- align_to_hm(mofa_hm, var_exp, pos = "right", newplt_size_prop = 0.13, legend_action = "collect")
    ggsave("results/mofa/finalMofa/factor_hm_complete2.pdf", width = 7, height = 6.3)
    
    return(complete_hm)
}

.factor_weights_hm <- function(final_mofa) {
    mofa_trained <- final_mofa$mofa
    annot_fcts <- final_mofa$facts
    factor_names <- final_mofa$fact_names
    factors <- structure(paste0("Factor", 1:10), names = paste0("Factor", 1:10))
    idx <- factors %in% annot_fcts
    names(factors)[idx] <- factor_names
    invert <- rep(1, 10)
    invert[idx] <- final_mofa$invert
    mofa_weights <- .plot_mofa_topGenes(mofa_trained, "mofa_mrna10000_cna1500_10f", 
                                        factors = factors, invert = invert)
    mofa_weights <- lapply(1:length(factors), function(i){ 
        mofa_weights[[i]] + 
            labs(subtitle = names(factors)[i])})
    mofa_weights[c(2:5, 7:10)] <- lapply(mofa_weights[c(2:5, 7:10)], function(plt) { plt + theme(axis.title.y = element_blank())})
    
    mofa_top_genes <- wrap_plots(mofa_weights, guides = "collect", nrow = 2) &
        theme(axis.text.y = element_text(size=5),
              axis.text.x = element_text(angle=90,hjust=0.5,vjust=0.5))
    mofa_top_genes
    ggsave("results/mofa/finalMofa/factor_topGenes.pdf", mofa_top_genes, width = 8, height = 8.5)
    return(mofa_top_genes)
}

.factor_paths <- function(final_mofa) {
    mofa_trained <- final_mofa$mofa
    factors2show <- paste0("Factor", 1:10)
    factors_select <- final_mofa$facts
    factor_names <- plyr::mapvalues(factors2show, from = factors_select, to = final_mofa$fact_names)
    
    invert <- structure(rep(1,10), names =  factors2show)
    invert[factors_select] <- final_mofa$invert
    mofa_paths <- .plot_mofa_factorPathways(mofa_trained, base_name = "finalMofa", 
                                            factors = factors2show, invert = invert, n = 2)
    mofa_paths <- lapply(1:length(factor_names), function(i){ 
        mofa_paths[[i]] + 
            labs(subtitle = factor_names[i])})
    mofa_path_plt <- wrap_plots(mofa_paths, guides = "collect", ncol = 2) &
        theme(plot.subtitle = element_text(size=8,face="bold"),
              axis.text.y=element_text(size=7))
    
    ggsave("results/mofa/finalMofa/factor_paths.pdf", mofa_path_plt, 
           width = 6, height = 10)
    return(mofa_path_plt)
}
.mofa_factor_hm <- function(data4hm, factor_names, pops, mks) {
    #-- Heatmap
    factor_hm <- data4hm %>%
        group_by(uPG) %>%
        dplyr::rename(`Stage` = Stage_, UROMOL = NMIBCclass) %>%
        # sample = factor(sample, levels = sample_order)) %>%
        ggheatmap(
            colv = "sample",
            rowv = rev(factor_names),
            hm_colors = viridis(100),
            hm_color_limits = c(-2,2),
            group_lines = TRUE,
            group_line_color = "white",
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            show_colnames = FALSE,
            rows_title = "Factor scores",
            colorbar_dir = "horizontal",
            colors_title = "Scaled factor score",
            group_prop = 0.12,
            fontsize = 9) %>%
        add_tracks(track_columns = c("Consensus", "UROMOL", "Stage"), 
                   track_colors = list(Consensus = consensus2,
                                       UROMOL = nmibc_cols,
                                       `Stage` = inv_cols),
                   track_pos = "top",
                   fontsize = 9, track_prop = 0.3) %>%
        add_tracks(track_columns = c("TP53_mutation", "FGFR3_mutation", "RAS_mutation",
                                     "CDKN2A_loss", "RB1_loss", "PPARG_gain", "PPARG_mutation",
                                     "RXRA_mutation",
                                     "Genomic_Instability"),
                   track_colors = list(
                       TP53_mutation = yn_cols,
                       FGFR3_mutation = yn_cols,
                       PPARG_mutation = yn_cols,
                       RXRA_mutation = yn_cols,
                       RAS_mutation = c("M" = "black", "WT" = "grey80"),
                       CDKN2A_loss = c("Loss" = "black", "WT" = "grey80"),
                       RB1_loss = c("Loss" = "black", "WT" = "grey80"),
                       PPARG_gain = c("Gain" = "black", "WT" = "grey80"),
                       Genomic_Instability = "OrRd"
                   ),
                   track_prop = 0.45,
                   fontsize = 9) %>%
        add_matrix_track(track_columns = as.character(unlist(mks)),
                         fontsize = 9, 
                         rows_title = "Protein expression",
                         color_limits = c(-4,4),
                         track_colors = inferno(100),
                         pal_dir = -1,
                         track_prop = 0.5) %>%
        add_matrix_track(track_columns = pops,
                         fontsize = 9, 
                         rows_title = "MCPcounter scores",
                         color_limits = c(-2,2),
                         track_colors = magma(100),
                         pal_dir = -1,
                         track_prop = 0.1)
}
.variance_exp_plot <- function(mofa_trained, factor_hm, factors2show, factor_names) {
    var_exp <- mofa_trained@cache $variance_explained$r2_per_factor$group1 %>%
        as.data.frame() %>%
        rownames_to_column("feat") %>%
        filter(feat %in% factors2show) %>%
        mutate(feat = plyr::mapvalues(feat, from = factors2show, to = factor_names)) %>%
        mutate(feat = factor(feat, levels = get_rowLevels(factor_hm))) %>%
        pivot_longer(cols = -feat, names_to = "level", values_to = "variance_exp") %>%
        mutate(variance_exp = variance_exp/100)
    
    var_exp_plt <- ggplot(var_exp, aes(level, feat, fill = variance_exp)) +
        geom_tile() +
        scale_fill_distiller(palette = "Blues", direction = 1, labels = scales::percent) +
        guides(fill = guide_colorbar(direction = "horizontal")) +
        theme_sparse2(base_size = 9) +
        theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
}

.plot_mofa_topGenes <- function(mofa_trained, base_name, factors = NULL, invert = NULL, n = 10,
                                width = 9, height = 7) {
    weights <- MOFA2::get_weights(mofa_trained)
    if(!is.null(factors)) weights <- lapply(weights, function(w) w[,factors])
    if(!is.null(invert)) {
        weights_new <- lapply(weights, function(w) {
            new_w <- sapply(1:ncol(w), function(i) { w[,i]*invert[i]})
            colnames(new_w) <- colnames(w)
            return(new_w)
        } ) 
        weights <- weights_new
    }
    weights_rn <- lapply(names(weights), function(name) {
        mat <- weights[[name]]
        rownames(mat) <- str_remove(rownames(mat), str_glue("_{name}"))
        return(mat)
    })
    names(weights_rn) <- names(weights) 
    plt_fct <- lapply(colnames(weights[[1]]), function(facti) {
        fct_top <- lapply(weights_rn, function(lvl_w) { 
            sorted_genes <- lvl_w[,facti] %>% sort() %>% names()
            c(sorted_genes[1:n], sorted_genes[(length(sorted_genes)-n):length(sorted_genes)])
        })
        fct_top <- Reduce(union, fct_top)
        
        fct_weights <- lapply(names(weights_rn), function(lvl) {
            lvl_w <- weights_rn[[lvl]]
            genes4plot <- intersect(rownames(lvl_w), fct_top)
            tibble(gene = genes4plot, weight = lvl_w[genes4plot,facti], lvl = lvl)
        })
        fct_weights <- do.call(bind_rows, fct_weights)
        order <- fct_weights %>%
            group_by(gene) %>%
            summarize(mean = mean(weight, na.rm = TRUE)) %>%
            arrange(mean) %>%
            pull(gene)
        ggplot(fct_weights, aes(lvl, factor(gene, levels = order), fill = weight)) +
            geom_tile() +
            scale_fill_viridis(option = "inferno", limits = c(-1,1), oob = scales::squish) +
            labs(x = "Data level", y = "Gene symbol", subtitle = facti) +
            theme_csg_sparse +
            theme(axis.text = element_text(size = 7),
                  axis.title = element_text(size=8, face="bold"),
                  plot.subtitle = element_text(size=8, face="bold"))
    })
    
    plt_top_genes <- wrap_plots(plt_fct, nrow = 2, guides = "collect")
    ggsave(str_glue("results/mofa/{base_name}/top_genes.pdf"), width = width, height = height)
    
    return(plt_top_genes)
}

.plot_mofa_QC <- function(mofa_trained, base_name) {
    #-- QC plots
    qc_cor <- plot_factor_cor(mofa_trained)
    qc_plot <- ggcorrplot::ggcorrplot(qc_cor$corr, method = "circle")
    var_plot1 <- plot_variance_explained(mofa_trained) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
    var_plot2 <- plot_variance_explained(mofa_trained, plot_total = T)[[2]]
    
    all_qc_plots <- (qc_plot / (var_plot1 | var_plot2)) + plot_layout(heights = c(1.5,1))
    
    ggsave(str_glue("results/mofa/{base_name}/qcPlots.pdf"),
           all_qc_plots, width = 6, height = 8)
}

.plot_mofa_associations <- function(mofa_trained, base_name) {
    sa <- samples_metadata(mofa_trained)
    #-- Factor sample associations with covariates
    #--------- Subtype
    factI_subtype <- associate_factors_with_covariates(mofa_trained, "Subtype",
                                                       color_limits = c(-2,2),
                                                       size_limits = c(0,12)) +
        theme_mofa
    factI_upg <- associate_factors_with_covariates(mofa_trained, "uPG",
                                                   color_limits = c(-2,2),
                                                   size_limits = c(0,12)) +
        theme_mofa
    sbty <- ((factI_subtype +
                  labs(x = "Subtype")) | 
                 (factI_upg + 
                      labs(x = "uPG") +
                      theme(axis.text.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.ticks.y = element_blank()))) +
        plot_layout(guides = "collect", widths = c(1.2,1))
    ggsave(str_glue("results/mofa/{base_name}/assocSubtype.pdf"),
           sbty, width = 5, height = 4)
    
    #------- Genomic alterations
    facI_mut <- associate_factors_with_covariates(mofa_trained,
                                                  c("FGFR3_mutation", "TP53_mutation"),
                                                  color_limits = c(-2,2),
                                                  size_limits = c(0,8)) +
        theme_mofa
    facI_chr <- associate_factors_with_covariates(mofa_trained, 
                                                  str_subset(colnames(sa), "[0-9].*(amp|del)"),
                                                  cluster_cols = TRUE,
                                                  color_limits = c(-2,2),
                                                  size_limits = c(0,8)) +
        theme_mofa
    facI_cnaGenes <- associate_factors_with_covariates(mofa_trained, 
                                                       str_subset(colnames(sa), "(gain|loss)_bin"),
                                                       cluster_cols = TRUE,
                                                       color_limits = c(-2,2),
                                                       size_limits = c(0,8)) +
        theme_mofa
    
    genAlt <- ((((facI_mut + labs(x = "Mutation")) |
                     (facI_cnaGenes + labs(x = "CNA") +
                          theme(axis.text.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.ticks.y = element_blank()))) +
                    plot_layout(widths = c(1,3))) /
                   facI_chr + labs(x = "CNA")) +
        plot_layout(guides = "collect")
    
    ggsave(str_glue("results/mofa/{base_name}/assocGeneticAlt.pdf"), genAlt, width = 7,
           height = 7)
    
    #------ trIC
    factI_tIC <- associate_factors_with_covariates(mofa_trained, 
                                                   str_subset(colnames(sa), "^IC"),
                                                   color_limits = c(-0.7,0.7),
                                                   cluster_cols = TRUE) +
        theme_mofa
    ggsave(str_glue("results/mofa/{base_name}/assocTranscriptomicICs.pdf"), factI_tIC,
           width = 4, height = 5)
    
    return(list(subtypes = sbty,
                genAlt = genAlt,
                tICs = factI_tIC))
}

.plot_mofa_factorHeatmap <- function(mofa_trained, base_name,
                                     factor_names = "all") {
    #-- Heatmap
    factors <- MOFA2::get_factors(mofa_trained)[[1]]
    if(factor_names != "all") {
        factors <- factors[,factor_names]
    }
    
    sa <- samples_metadata(mofa_trained)
    factor_table <- factors %>%
        as.data.frame() %>%
        rownames_to_column("sample") %>%
        left_join(sa) %>%
        dplyr::rename(UROMOL = NMIBCclass)
    
    factor_hm <- factor_table %>%
        group_by(uPG) %>%
        ggheatmap(colv = "sample",
                  rowv = rev(colnames(factors)),
                  show_colnames = FALSE,
                  cluster_rows = FALSE,
                  hm_color_limits = c(-1.5,1.5),
                  hm_colors = viridis(100),
                  group_lines = TRUE,
                  group_line_color = "white",
                  colorbar_dir = "horizontal",
                  colors_title = "Factor score",
                  raster = TRUE,
                  fontsize = 8,
                  group_prop = 0.06) %>%
        add_tracks(track_columns = c("Consensus", "UROMOL", "Stage"),
                   track_colors = list(
                       Consensus = consensus2,
                       UROMOL = nmibc_cols,
                       Stage = inv_cols
                   ),
                   fontsize = 8,
                   track_pos = "top", 
                   track_prop = 0.2, 
                   legend_action = "collect")
    #----- Variance explained track
    var_exp <- mofa_trained@cache $variance_explained$r2_per_factor$group1 %>%
        as.data.frame() %>%
        rownames_to_column("feat") %>%
        mutate(feat = factor(feat, levels = get_rowLevels(factor_hm))) %>%
        pivot_longer(cols = -feat, names_to = "level", values_to = "variance_exp")
    
    var_exp_plt <- ggplot(var_exp, aes(level, feat, fill = variance_exp)) +
        geom_tile() +
        scale_fill_distiller(palette = "Blues", direction = 1) +
        guides(fill = guide_colorbar(direction = "horizontal")) +
        theme_sparse2(base_size = 8) +
        theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
    
    full_hm <- align_to_hm(factor_hm, var_exp_plt, newplt_size_prop = 0.2,
                           pos = "right", legend_action = "collect")
    
    full_hm
    
    ggsave(str_glue("results/mofa/{base_name}/factorHeatmap.pdf"), 
           plot = full_hm, width = 6, height = 3)
    
    return(full_hm)
}

.plot_mofa_factorPathways <- function(mofa_trained, base_name = NULL, 
                                      factors = NULL, invert = NULL, n = 3) {
    
    #-- Pathway enrichment
    mrna_enrich <- get_enrichment(mofa_trained, level = "mrna", p_cutoff = 1)
    prot_enrich <- get_enrichment(mofa_trained, level = "prot", p_cutoff = 1)
    cna_enrich <- get_enrichment(mofa_trained, level = "cna", p_cutoff = 1)
    
    if(is.null(factors)) {
        factors <- colnames(get_factors(mofa_trained)[[1]])
    }
    if(is.null(invert)) {
        invert <- rep(1, length(factors))
    }
    
    #-- Get top paths 
    top_paths <- lapply(list(mrna_enrich, prot_enrich, cna_enrich), 
                        function(enrich_tb) {
                            paths_up <- enrich_tb %>%
                                group_by(factor) %>%
                                slice_min(mean_diff, n = n) %>%
                                pull(pathway) %>%
                                unique()
                            
                            paths_down <- enrich_tb %>%
                                group_by(factor) %>%
                                slice_max(mean_diff, n = n) %>%
                                pull(pathway) %>%
                                unique()
                            
                            union(paths_up, paths_down)
                        })
    top_paths <- Reduce(union, top_paths)
    
    #-- Filter significant and join
    mrna_enrich_f <- mrna_enrich %>%
        mutate(mlog10 = -log10(padj)) %>%
        select(pathway, factor, mrna_padj = mlog10, mrna_diff = mean_diff) %>%
        filter(pathway %in% top_paths)
    prot_enrich_f <- prot_enrich %>%
        mutate(mlog10 = -log10(padj)) %>%
        select(pathway, factor, prot_padj = mlog10, prot_diff = mean_diff) %>%
        filter(pathway %in% top_paths)
    cna_enrich_f <- cna_enrich %>%
        mutate(mlog10 = -log10(padj)) %>%
        select(pathway, factor, cna_padj = mlog10, cna_diff = mean_diff) %>%
        filter(pathway %in% top_paths)
    
    pathway_enrich <- full_join(mrna_enrich_f, prot_enrich_f, 
                                by = c('factor', 'pathway')) %>% 
        full_join(cna_enrich_f,  by = c('factor', 'pathway'))
    
    #-- Make factor pathway plots
    var_exp <- mofa_trained@cache$variance_explained$r2_per_factor$group1
    slice_vars <- structure(paste0(colnames(var_exp)[apply(var_exp, 1, which.max)], "_padj"),
                            names = rownames(var_exp))[factors]
    
    path_plots <- lapply(1:length(factors), function(i) {
        fact <- factors[i]
        slice_var <- slice_vars[i]
        inv <- invert[i]
        paths2show <- pathway_enrich %>%
            filter(
                factor == fact,
                prot_padj > 2 | mrna_padj > 2 | cna_padj > 2) %>%
            arrange(across(matches(slice_var), dplyr::desc)) %>%
            slice(n = 1:10) %>%
            mutate(
                cna_diff = cna_diff * inv,
                mrna_diff = mrna_diff * inv,
                prot_diff = prot_diff * inv,
                pathway = str_sub(pathway, end = 30),
                pathway = fct_reorder(pathway, prot_diff)) %>%
            pivot_longer(cols = matches("diff|padj"), values_to = "values", names_to = "var") %>%
            mutate(level = str_remove(var, "_diff$|_padj$"),
                   var = str_remove(var, "mrna_|prot_|cna_")) %>%
            pivot_wider(names_from = var, values_from = values, values_fn = mean) %>%
            mutate(mlog10 = -log10(padj)) 
        
        ggplot(paths2show, aes(level, pathway, size = mlog10, fill = diff)) +
            geom_point(pch = 21) +
            scale_fill_distiller(palette = "RdBu", limits = c(-5,5), oob = squish) +
            scale_size_area(limits = c(-1,6), oob = scales::squish) +
            labs(x = "Data level", y = "Pathway", subtitle = fact) +
            theme_mofa
    })
    # path_data <- lapply(path_plots, function(plt) { plt$data }) %>% bind_rows()
    path_plots <- wrap_plots(path_plots, ncol = 3) +
        plot_layout(guides = "collect")
    
    if(!is.null(base_name)) {
        ggsave(str_glue("results/mofa/{base_name}/pathways.pdf"), width = 8, height = 10)
    }
    
    return(path_plots)
}

mcpCounter <- function(mrna_data, all_samps) {
    mcp_res <- MCPcounter.estimate(mrna_data$gexp[,all_samps], featuresType = "HUGO_symbols") %>%
        t() %>%
        scale() %>%
        as.data.frame()
    return(mcp_res)
}

export_mofa <- function(final_mofa, interpret_mofa) {
    
    all_mofa_weights <- get_weights(final_mofa$mofa)
    fct_rpl <- structure(final_mofa$facts, names = final_mofa$fact_names)
    fct_names <- structure(paste0("Factor", 1:10), names = paste0("Factor", 1:10))
    fct_invert <- structure(rep(1, 10), names = fct_names)
    fct_invert[fct_rpl] <- final_mofa$invert
    fct_names[fct_rpl] <- names(fct_rpl)
    all_mofa_weights2 <- lapply(all_mofa_weights, function(mofa_weights) {
        rownames(mofa_weights) <- str_remove(rownames(mofa_weights), "_.*")
        mofa_weights <- sapply(1:10, function(i) {
            mofa_weights[,i] * fct_invert[i]
        })
        colnames(mofa_weights) <- fct_names
        return(mofa_weights)
    })
    names(all_mofa_weights2) <- c("CNA", "Transcriptomics", "Proteomics")
    mofa_weights <- lapply(all_mofa_weights2, function(tb) {
        tb %>%
            as.data.frame() %>%
            rownames_to_column("symbol") %>%
            mutate(symbol = str_remove_all(symbol, "_.*")) %>%
            tibble()
    })
    openxlsx::write.xlsx(mofa_weights, "results/tables/MOFA_allWeights.xlsx")
    
}

export_mofaAssociations <- function(interpret_mofa) {
    ic_table <- interpret_mofa$`NULL`$associations$tICs$data
    gen_alt <- interpret_mofa$`NULL`$associations$genAlt$data
    
    #-- IC
    ic_table <- ic_table %>%
        mutate(
            mlogpval = ifelse(is.infinite(mlogpval), NA, mlogpval),
            padj = 10^(-mlogpval)) %>%
        select(factor = feat, IC = var, Spearmans_rho = value, padj)
    #-- CNA
    gen_alt <- gen_alt %>%
        mutate(
            mlogpval = ifelse(is.infinite(mlogpval), NA, mlogpval),
            padj = 10^(-mlogpval)) %>%
        select(factor = feat, CNA = var, log2FC = value, padj)
    
    #-- Paths
    path_plts <- interpret_mofa$`NULL`$factor_pathways$patches$plots
    path_data <- lapply(path_plts, function(plt) plt$data) %>% bind_rows()
    paths <- path_data %>%
        mutate(source = "Reactome",
               factor = factor(factor, levels = paste0("Factor", 1:10)),
               level = factor(level, levels = c("cna", "mrna", "prot")),
               padj = 10^-padj) %>%
        arrange(factor, pathway, level) %>%
        select(-mlog10) %>%
        rename(differential_ES = diff)
        
    openxlsx::write.xlsx(list(Pathway_enrichment = paths, 
                              IC_correlation = ic_table, 
                              CNA_association = gen_alt), "results/tables/MOFA_pathway_IC_CNA_associations.xlsx")
    
}