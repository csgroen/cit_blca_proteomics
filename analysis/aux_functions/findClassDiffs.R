#' data = data.frame or matrix of biological data
#' cl = class vector, samples as names
#'
# AUC = U/(n0 * n1) where U is the Mann-Whitney U stat, n0 and n1 are nobs (target, non-target, resp.)
# using moderated t-test 

findClassDiffs <- function(data, cl, type = c("parametric", "non-parametric")[1],
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
        filter(padj < pval_cutoff) %>%
        rename_with(~ paste(test_type, c("pval", "padj"), sep = "_"), 2:3)
    
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
