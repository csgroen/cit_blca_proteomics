annotatedHeatmap <- function(data, sa, group_var = NULL, id_var = NULL,
                             annot_vars = NULL,
                             impute = TRUE, kimpt = 3,
                             quantile_breaks = TRUE, bks_lims = c(1,101),
                             hc_dist_method = "euclidean", hc_dist_p = 2, hc_hclust_method = "complete",
                             hc_cor_use = "everything", hc_merge_height = NA, 
                             feature_cluster = NULL,
                             ...) {
    #-- Check to impute
    if(impute & any(is.na(data))) {
        impt <- impute::impute.knn(data, k = kimpt)$data
    } else if (!any(is.na(data))) {
        impt <- data
    } else {
        stop("There are `NA`s that need to be imputed or dealt with previous to function call.")
    }
    #-- Quantile breaks?
    if(quantile_breaks) {
        bks <- quantile(impt, seq(0,1,0.01))[bks_lims[1]:bks_lims[2]]
    } else {
        bks <- seq(min(impt), max(impt), length.out = 101)[bks_lims[1]:bks_lims[2]]
    }
    #-- Harmonize sa
    if (is.null(id_var)) {
        sa <- rownames_to_column(sa, "id")
    } else if (id_var != "id") {
        sa <- rename(id = !! id_var)
    }
    #-- Make annotation
    if(!is.null(annot_vars)) {
        annot_col <- sa %>%
            select(id, !! annot_vars) %>%
            as.data.frame() %>%
            column_to_rownames("id")
    } else {
        annot_col <- NULL
    }
    #-- Semi-supervised or not?
    #-------- Normal
    if (is.null(group_var)) {
        #-- Calculate dist
        s_hc <- hierarch_cluster(impt, dist_method = hc_dist_method,
                                    dist_p = hc_dist_p, cor_use = hc_cor_use,
                                    hclust_method = hc_hclust_method)
        if(is.null(feature_cluster)) {
            f_hc <- hierarch_cluster(t(impt), dist_method = hc_dist_method,
                                     dist_p = hc_dist_p, cor_use = hc_cor_use,
                                     hclust_method = hc_hclust_method)
        } else {
            f_hc <- feature_cluster
        }
        
        data4plot <- data
    #------ Semi-supervised
    } else {
        s_groups <- split(sa$id, sa[,group_var])
        s_hcss <- hclust_semisupervised(t(impt), groups = s_groups, 
                                        dist_method = hc_dist_method, 
                                        dist_p = hc_dist_p, 
                                        cor_use = hc_cor_use,
                                        hclust_method = hc_hclust_method, 
                                        merge_height = hc_merge_height)
        if(is.null(feature_cluster)) {
            f_hc <- hierarch_cluster(s_hcss$data, dist_method = hc_dist_method,
                                     dist_p = hc_dist_p, cor_use = hc_cor_use,
                                     hclust_method = hc_hclust_method)
        } else {
            f_hc <- feature_cluster
        }

        
        data4plot <- data[colnames(s_hcss$data),rownames(s_hcss$data)]
        s_hc <- s_hcss$hclust
    }
    
    #-- Make heatmap
    if(is.null(annot_col)) {
        hm <- pheatmap(data4plot,
                       breaks = bks,
                       border_color = NA, 
                       cluster_rows = f_hc, 
                       cluster_cols = s_hc,
                       ...)
    } else {
        hm <- pheatmap(data4plot,
                       breaks = bks,
                       border_color = NA, 
                       cluster_rows = f_hc, 
                       cluster_cols = s_hc,
                       annotation_col = annot_col,
                       ...)
    }
    return(hm)
}

hierarch_cluster <- function(data, dist_method, dist_p, cor_use, hclust_method) {
    if(dist_method %in% c("pearson", "spearman", "kendall")) {
        disto <- as.dist(1 - cor(data, use = cor_use, method = dist_method))
    } else {
        disto <- dist(t(data), method = dist_method)
    }
    return(hclust(disto, method = hclust_method))
}

#' Semi-supervised hierarchical clustering
#' 
#' Semi-supervised hierarchical clustering by chosen groups with hclust.
#'
#' @param data a data.frame to be clustered by rows
#' @param groups a list of vectors. If we unlist(groups), all elements must be
#'   present in the rownames of data. Each vector in the list will be treated as
#'   a separate group for the hierarchical clustering, and rejoined in order at
#'   the end.
#' @param dist_method a distance computation method. Must be one of "euclidean", 
#' "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman"
#' @param dist_p the power of the Minkowski distance, if chosen dist_method is "minkowski"
#' @param hclust_method an agglomeration method. Should be a method supported by
#' hclust, one of:  "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), 
#' "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param cor_use If using correlation as distance, chooses the method for computing
#' covariances in the presence of missing values. See ? cor.
#' 
#' @return hclust_semisupervised returns a list. The first element of the list
#' is the data, reordered so that the merged hclust object will work. The second
#' element is the result of the semi-supervised hierarchical clustering.

hclust_semisupervised <- function(data, groups, dist_method = "euclidean",
                                  dist_p = 2, hclust_method = "complete", cor_use = "everything",
                                  merge_height = NA) {
    #-- Group checks
    if (!all(unlist(groups) %in% rownames(data))) {
        stop("vectors in `groups` must contain rownames of `data`")
    }
    if(anyDuplicated(unlist(groups))){
        stop("`groups` can't have elements duplicated within or in different 
             groups")
    }
    
    #-- dist_method checks
    alldists <- c("euclidean", "maximum", "manhattan", "canberra", "binary", 
                  "minkowski", "pearson", "spearman", "kendall")
    if(!(dist_method %in% alldists)) {
        stop("`dist_method` must be one of: `euclidean`, `maximum`, `manhattan`, 
             `canberra`, `binary`, `minkowski`, `pearson`, `spearman` or `kendall`.")
    }
    
    #-- Use check
    if(dist_method %in% c("pearson", "spearman", "kendall")) {
        alluse <- c("everything", "all.obs", "complete.obs", "na.or.complete", 
                    "pairwise.complete.obs")
        if(!(cor_use %in% alluse)) {
            stop('`cor_use` must be one of "everything", "all.obs", "complete.obs", 
                 "na.or.complete", or "pairwise.complete.obs".')
        }
    }
    
    #-- Get groups with 1 member
    g_size <- sapply(groups, length)
    if (any(g_size == 1)) {
        # s_groups <- groups[g_size == 1]
        groups <- groups[g_size != 1]
        
    }
    
    #-- Make distance matrices
    if (dist_method %in% c("pearson", "spearman", "kendall")) {
        
        distlist <- lapply(groups, function(group) {
            as.dist(1 - cor(t(data[group,]), method = dist_method, use = cor_use))
        })
    } else {
        distlist <- lapply(groups, function (group) {
            dist(data[group,], method = dist_method, p = dist_p)
        })
    }
    #-- Use hclust
    hclist <- lapply(distlist, hclust, method = hclust_method)
    
    
    hc <- .merge_hclust(hclist, height = merge_height)
    
    #-- Join groups with one element
    if(exists("s_groups")) {
        # s_groups <- unlist(s_groups)
        # if (s_groups > 1) s_hc <- hclust(dist(data[s_groups,]))
        # else {
        #     s_hc <- list(merge = matrix(c(1), ncol = 1), height = max(hc$height + sd(hc$height)),
        #          order = 1, labels = s_groups, call = NA, method = NA,
        #          dist.method = NA)
        #     class(s_hc) <- "hclust"
        # }
        # hc <- .merge_hclust(list(hc, s_hc), height = merge_height)
        # data_reordered <- data[c(unlist(groups), s_groups),]
    } else {
        data_reordered <- data[unlist(groups),]
    }
    
    data_reordered <- data[unlist(groups),]
    
    return(list(data = data_reordered, 
                hclust = hc))
}

.merge_hclust <- function(hclist, height) {
    #-- Check
    if(!is.list(hclist)) {
        stop("`hclist` must be a list.")
    }
    if(!all(sapply(hclist, class) == "hclust")){
        stop("All objects in `hclist` must be `hclust-class`")
    }
    
    #-- Merge
    d <- as.dendrogram(hclist[[1]])
    for (i in 2:length(hclist)) {
        if(is.na(height)) { 
            d <- merge(d, as.dendrogram(hclist[[i]]), adjust = "add.max")
        } else {
            d <- merge(d, as.dendrogram(hclist[[i]]), adjust = "add.max", height = height)
        }
    }
    as.hclust(d)
}