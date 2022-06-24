associate_factors_with_covariates <- function(mofa_obj, covariates, 
                                              fct_names = NULL,
                                              invert = NULL,
                                              id_var = "sample",
                                              cluster_rows = FALSE, cluster_cols = FALSE,
                                              color_high = "#ca0020", color_low = "#0571b0",
                                              palette = NULL,
                                              color_limits = NULL,
                                              size_limits = NULL) {
    #-- Get basic info
    factors <- get_factors(mofa_obj)[[1]] %>% t()
    if(is.null(fct_names)) {
        fct_names <- rownames(factors)
    }
    factors <- factors[fct_names,]
    
    if(!is.null(invert)) {
        factors <- sapply(1:nrow(factors), function(i) factors[i,] * invert[i]) %>%
            t()
        rownames(factors) <- fct_names
    }
    
    covar_tb <- mofa_obj@samples_metadata %>%
        rownames_to_column("id") %>%
        select(id, {{covariates}})
    
    #-- Get variable type and do tests (compare for factors, cor for numeric)
    var_type <- sapply(covariates, function(var) {
        class(pull(covar_tb, {{var}})) %in% c("character", "factor")
    })
    var_levels <- sapply(covariates, function(var) { 
        pull(covar_tb, {{var}}) %>% unique() %>% na.exclude() %>% length()
    })
    
    all_vtype <- ifelse(all(var_type) & all(var_levels == 2), "binary", "other")
    all_vtype <- ifelse(all(var_type) & all(var_levels > 2) & length(covariates) == 1, "categorical", all_vtype)
    all_vtype <- ifelse(!any(var_type), "numeric", all_vtype)
    
    #-- Run tests
    if(all_vtype == "binary") {
        message("Plotting binary variables...")
        factor_diffs <- lapply(covariates, function(var) {
            cl_var <- covar_tb %>% pull({{var}}, id)
            
            .findClassDiffs(factors, cl_var, pval_cutoff = 1, type = "non-parametric") %>%
                mutate(feat = factor(feat, levels = fct_names),
                       var = var,
                       mlogpval = -log10(padj)) %>%
                arrange(feat) %>%
                select(feat, var, mlogpval, value = logFC)
        }) %>% bind_rows()
        
    } else if (all_vtype == "numeric") {
        message("Plotting numerical variables...")
        
        factor_diffs <- lapply(covariates, function(var) {
            cl_var <- covar_tb %>% pull({{var}}, id)
            tibble(
                feat = factor(fct_names, levels = fct_names),
                value = cor(t(factors[,names(cl_var)]), cl_var, method = "spearman", use = "complete")[,1]) %>%
                mutate(padj = .cor_pval(abs(value), length(cl_var)),
                       mlogpval = -log10(padj),
                       var = var) %>%
                select(feat, var, mlogpval, value)
        }) %>%
            bind_rows()
        
    } else if (all_vtype == "categorical") {
        cl_var <- covar_tb %>% pull({{covariates}}, id)
        factor_diffs <- .findClassDiffs(factors, cl_var, pval_cutoff = 1, type = "non-parametric") %>%
            select(feat, starts_with("logFC"), starts_with("pair_pval")) %>%
            pivot_longer(cols = starts_with("logFC"), names_to = "var", values_to = "value") %>%
            pivot_longer(cols = starts_with("pair_pval"), names_to = "var2", values_to = "pval") %>%
            mutate(var = str_remove(var, "logFC_"),
                   var2 = str_remove(var2, "pair_pval_")) %>%
            filter(var == var2) %>%
            mutate(mlogpval = -log10(pval)) %>%
            select(feat, var, mlogpval, value)
    } else {
        stop("Can't mix binary, categorical and continuous variables for this visualization")
    }
    fdf <- factor_diffs %>%
        filter(!is.na(feat)) %>%
        pivot_wider(id_cols = feat, names_from = var, values_from = value) %>%
        as.data.frame() %>%
        column_to_rownames("feat") %>%
        as.matrix() %>%
        .[fct_names,]
    
    
    #-- Cluster for visualization
    if(cluster_rows) {
        factor_order <- fct_names[hclust(dist(fdf), method = "ward.D2")$order] %>% rev()
    } else {
        factor_order <- fct_names %>% rev()
    }
    if(cluster_cols) {
        covar_order <- colnames(fdf)[hclust(dist(t(fdf)), method = "ward.D2")$order]
    } else {
        if(all_vtype == "categorical") {
            covar_order <- levels(covar_tb[,covariates])
        } else {
            covar_order <- covariates
            
        }
    }
    
    fill_name <- ifelse(all_vtype == "categorical" | all_vtype == "binary", "log2FC", "Spearman's rho")
    #-- Plot
    plt <- factor_diffs %>%
        mutate(feat = factor(feat, levels = factor_order),
               var = factor(var, levels = covar_order)) %>%
        ggplot(aes(var, feat, fill = value, size = mlogpval)) +
        geom_point(color = "black", pch = 21) +
        labs(x = "Covariates", y = "Factor", fill = fill_name, size = "-log10(adj.\np-value)") +
        scale_size_area(limits = size_limits, oob = scales::squish) +
        theme_scatter() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    if(is.null(palette)) {
        plt <- plt + 
            scale_fill_gradient2(low = color_low, 
                                 high = color_high,
                                 oob = scales::squish, 
                                 limits = color_limits)
    } else if (palette %in% c("viridis", "magma", "plasma", "inferno", "cividis")) {
        plt <- plt + scale_fill_viridis(option = palette,
                                        oob = scales::squish, 
                                        limits = color_limits)
    } else {
        plt <- plt + scale_fill_distiller(palette = palette,
                                          oob = scales::squish, 
                                          limits = color_limits)
    }
    
    return(plt)
}

#' data = data.frame or matrix of biological data
#' cl = class vector, samples as names
#'
# AUC = U/(n0 * n1) where U is the Mann-Whitney U stat, n0 and n1 are nobs (target, non-target, resp.)
# using moderated t-test 

.findClassDiffs <- function(data, cl, type = c("parametric", "non-parametric")[1],
                           pval_cutoff = 0.05, p_adj = "fdr", log_transformed = TRUE,
                           add_paired = TRUE) {
    #-- Checks of parameters
    if(!all(colnames(data) %in% names(cl))) stop("`cl` must be named with the sample ids in `data`.")
    if(!is.factor(cl)) cl <- as.factor(cl)
    ngr <- length(levels(cl))
    if(ngr < 2) stop("`cl` must have at least 2 levels for testing.")
    data <- data[,names(cl)]
    
    #-- Define test type
    test_type <- case_when(
        type == "parametric" & ngr > 2 ~ "anova",
        type == "non-parametric" & ngr > 2 ~ "kw",
        type == "parametric" & ngr == 2 ~ "ttest",
        TRUE ~ "wilcox"
    )
    test_name <- switch(test_type,
                        "anova" = "ANOVA test",
                        "kw" = "Kruskal-Wallis test",
                        "ttest" = "Moderated t-test",
                        "wilcox" = "Mann-Whitney U-test")
    
    #-- Run test
    message(paste0("-- Running the ", test_name, "..."))
    pvals <- .runStatTests(data, cl, test_type)
    
    #-- Make results table
    results <- tibble(feat = rownames(data),
                      pval = pvals,
                      padj = p.adjust(pvals, method = p_adj)) %>%
        arrange(padj, pval)
    
    #-- Add more information: Class AUC + logFC
    grs <- split(names(cl), cl)
    grs <- grs[levels(cl)]
    grs <- grs[sapply(grs, length) > 0]
    ngr <- length(grs)
    message("-- Adding class AUCs...")
    class_aucs <- .classAUCs(data, grs)
    message("-- Adding logFC...")
    log_fcs <- .classLogFC(data, grs, ngr, log_transformed)
    if(add_paired & ngr > 2) {
        message("-- Adding paired tests...")
        paired_pvals <- sapply(grs, function(in_gr) {
            out_gr <- setdiff(names(cl), in_gr)
            cl_gr <- structure(c(rep("A", length(in_gr)),rep("B", length(out_gr))),
                               names = c(in_gr, out_gr)) %>% as.factor()
            tt <- ifelse(type == "parametric", "ttest", "wilcox")
            in_pvals <- .runStatTests(data, cl = cl_gr, test_type = tt)
            p.adjust(in_pvals, "BH")
        })
        colnames(paired_pvals) <- paste0("pair_pval_", names(grs))
        class_ppvals <- paired_pvals %>%
            as.data.frame() %>%
            rownames_to_column("feat")
    }
    #-- Join and return
    results_table <- results %>%
        left_join(class_aucs, by = "feat") %>%
        left_join(log_fcs, by = "feat") %>%
        filter(padj < pval_cutoff)
        # rename_with(~ paste(test_type, c("pval", "padj"), sep = "_"), 2:3)
    
    if(add_paired & ngr > 2) results_table <- left_join(results_table, class_ppvals, by = "feat")
    
    return(results_table)
}

# Compute test p-value
.runStatTests <- function(data, cl, test_type) {
    if(test_type == "ttest") {
        # pvals <- sapply(rownames(data), function(feat) {
        #     fsplit <- split(data[feat,], cl)
        #     if(sum(!is.na(fsplit[[1]])) > 1 & sum(!is.na(fsplit[[2]])) > 1)
        #         t.test(fsplit[[1]], fsplit[[2]])$p.value
        #     else
        #         NA
        # })
        #-- Moderated p-value
        gr1 <- split(names(cl), cl)[[1]]; gr2 <- split(names(cl), cl)[[2]]
        n1 <- length(gr1); n2 <- length(gr2)
        dx <- rowMeans(data[,gr1], na.rm = TRUE)-rowMeans(data[,gr2], na.rm = TRUE)
        stderr <- sqrt( (rowVars(data[,gr1], na.rm = TRUE)*(n1-1) + rowVars(data[,gr2], na.rm = TRUE)*(n2-1)) / (n1+n2-2) * ( 1/n1 + 1/n2 ))
        mod.stderr <- (stderr + median(stderr, na.rm = TRUE)) / 2
        t.mod = dx / mod.stderr
        p.mod = 2*pt( -abs(t.mod), n1+n2-2 )
        names(p.mod) <- rownames(data)
        return(p.mod)
    } else if(test_type == "wilcox") {
        pvals <- sapply(rownames(data), function(feat) {
            fsplit <- split(data[feat,names(cl)], cl)
            if(sum(!is.na(fsplit[[1]])) > 1 & sum(!is.na(fsplit[[2]])) > 1)
                wilcox.test(fsplit[[1]], fsplit[[2]])$p.value
            else
                NA
        })
    } else if(test_type == "kw") {
        pvals <- sapply(rownames(data), function(feat) {
            splt <- split(data[feat, names(cl)], cl)
            nobs_gr <- sapply(splt, function(sp) { sum(!is.na(sp)) })
            if(sum(nobs_gr > 0) <= 1) return(NA)
            kruskal.test(data[feat,], cl)$p.value
        })
    } else if (test_type == "anova") {
        pvals <- sapply(rownames(data), function(feat) {
            res <- aov(data[feat,] ~ cl)
            summary(res)[[1]][["Pr(>F)"]][1]
        })
    }
}

# Compute AUC for each class (or AUC for variable)
.classAUCs <- function(data, grs) {
    #-- Class > 2 groups
    if(length(grs) > 2) {
        class_aucs <- sapply(grs, function(in_gr) {
            out_gr <- setdiff(colnames(data), in_gr)
            aucs <- sapply(rownames(data), function(feat) {
                val_in <- data[feat, in_gr]
                val_out <- data[feat, out_gr]
                if(sum(!is.na(val_in)) > 1 & sum(!is.na(val_out))) {
                    u_val <- wilcox.test(val_in, val_out)$statistic
                    names(u_val) <- NULL
                    return(u_val / (length(in_gr) * length(out_gr)))
                } else {
                    return(NA)
                }
            })
        })
        colnames(class_aucs) <- paste0("auc_", colnames(class_aucs))
        class_aucs <- class_aucs %>%
            as.data.frame() %>%
            rownames_to_column("feat") %>%
            mutate(across(starts_with("auc"), function(aucs) { 
                ifelse(aucs < 0.5, 1 - aucs, aucs)
            })) %>%
            as_tibble()
    }
    else {
        us <- sapply(rownames(data), function(feat) {
            fsplit <- list(data[feat, grs[[1]]], data[feat, grs[[2]]])
            if(sum(!is.na(fsplit[[1]])) > 1 & sum(!is.na(fsplit[[2]])) > 1)
                wilcox.test(fsplit[[1]], fsplit[[2]])$statistic
            else
                NA
        })
        aucs <- us/(length(grs[[1]]) * length(grs[[2]]))
        class_aucs <- tibble(feat = rownames(data), auc = ifelse(aucs < 0.5, 1 - aucs, aucs))
    }
    
    return(class_aucs)
}

# Compute logFC against class v others (or exp v control)
.classLogFC <- function(data, grs, ngr, log_transformed) {
    if(length(grs) > 2) {
        log_fcs <- sapply(grs, function(in_gr) {
            out_gr <- setdiff(colnames(data), in_gr)
            vals1 <- rowMeans(data[,in_gr], na.rm = TRUE)
            vals2 <- rowMeans(data[,out_gr], na.rm = TRUE)
            if(log_transformed) return(rowMeans(data[,in_gr], na.rm = TRUE) - rowMeans(data[,out_gr], na.rm = TRUE))
            else stop("Not supported yet.")
        }) %>%
            as.data.frame() %>%
            rownames_to_column("feat") %>%
            rename_with(~ paste0("logFC_", .), 2:(ngr+1)) %>%
            as_tibble()
    } else {
        if(log_transformed) fc <- rowMeans(data[,grs[[1]]], na.rm = TRUE) - rowMeans(data[,grs[[2]]], na.rm = TRUE)
        else stop("Not supported yet.")
        
        log_fcs <- tibble(feat = rownames(data), logFC = fc)
    }
    return(log_fcs)
}

get_enrichment <- function(mofa_trained, level = "mrna", 
                           p_cutoff = 0.05, top_n = NULL) {
    load("data/annotations/Reactome_pathways.RData")
    reactome_mat <- react_paths %>%
        select(Path_Name, Symbol) %>%
        filter(!is.na(Symbol)) %>%
        mutate(Symbol = paste0(Symbol, "_", level),
               in_path = 1) %>%
        pivot_wider(id_cols = Path_Name, names_from = Symbol, values_from = in_path, values_fn = sum) %>%
        mutate(across(ends_with(level), ~ ifelse(is.na(.), 0, 1))) %>%
        as.data.frame() %>%
        column_to_rownames("Path_Name") %>%
        as.matrix()
    
    gseaP <- run_enrichment(mofa_trained, view = level,
                            feature.sets = reactome_mat, 
                            sign = "positive")
    gseaN <- run_enrichment(mofa_trained, view = level,
                            feature.sets = reactome_mat, 
                            sign = "negative")
    
    factors <- colnames(gseaP$pval.adj)
    res_tb <- lapply(factors, function(fact) {
        resP <- .get_enrichResTable(gseaP, fact, p_cutoff)
        resN <- .get_enrichResTable(gseaN, fact, p_cutoff) %>%
            mutate(mean_diff = -mean_diff)
        gsea_res <- bind_rows(resP, resN)
    } ) %>% bind_rows() %>%
        mutate(factor = factor(factor, levels = factors))
    
    if(!is.null(top_n)) {
        top <- res_tb %>%
            group_by(factor) %>%
            slice_max(mean_diff, n = top_n/2)
        bottom <- res_tb %>%
            group_by(factor) %>%
            slice_min(mean_diff, n = top_n/2)
        
        res_tb <- bind_rows(top, bottom) %>% 
            ungroup(factor) %>%
            arrange(factor, desc(mean_diff))
    }
    return(res_tb)
    
}

.get_enrichResTable <- function(res_table, fact, p_cutoff) {
    resP <- tibble(pathway = rownames(res_table$pval.adj), 
                   padj = res_table$pval.adj[,fact],
                   mean_diff = res_table$set.statistics[,fact]) %>%
        filter(padj < p_cutoff) %>%
        mutate(factor = fact) %>%
        relocate(factor) %>%
        arrange(padj)
}
.cor_pval <- function(r,n) {
    t <- r*sqrt((n-2)/(1-r^2))
    p <- 1 - pt(t, n - 1)
    return(p)
}
