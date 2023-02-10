clusterProteins <- function(prot_data, kmin, kmax, methods = c("pam"),
                            distance = "spearman", impute = FALSE, impute_method = "knnImpute",
                            res_dir = "results/protein_clustering/", make_plots = TRUE,
                            annot_col = prot_data$sampAnnot %>%
                                dplyr::select(id,
                                              Consensus, 
                                              Stage.,
                                              FGFR3_mutation, 
                                              TP53_mutation) %>% column_to_rownames("id"),
                            annot_colors = list(Consensus = consensus2,
                                                Stage. = inv_cols,
                                                FGFR3_mutation = yn_cols,
                                                TP53_mutation = yn_cols),
                            n_jobs = 6) {
    dir.create(res_dir, showWarnings = FALSE)
    #-- Impute for clustering
    if(impute) {
        pp <- preProcess(prot_data$wp, method = impute_method, k = 3)
        impt_data <- predict(pp, newdata = prot_data$wp) %>% t()
        use <- "everything"
    } else {
        impt_data <- prot_data$wp %>% t()
        use <- "pairwise.complete.obs"
    }

    #-- Get clusters
    if(file.exists("cached_results/protein_clustering.RData")) {
        message("Using cached results...")
        load("cached_results/protein_clustering.RData")
    } else {
        #-- Calculate distances
        diss <- as.dist(1 - cor(impt_data, method = distance, use = use))
        clustermq::register_dopar_cmq(n_jobs=n_jobs)
        hclust_met <- methods %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", 
                                     "median", "centroid")
        
        if(any(hclust_met)) {
            prot_clust <- hclust(diss, method = methods[hclust_met][1])
        }
        all_clust_sols <- lapply(methods, function(cmethod) {
            if(cmethod == "pam") {
                pam_grs <- foreach(k=kmin:kmax, .combine = cbind, .export = c("diss")) %dopar%
                    (pam(diss, k = k)$clustering)
                colnames(pam_grs) <- paste0("pam_k=", kmin:kmax)
                return(pam_grs)
            } else {
                hc_grs <- sapply(kmin:kmax, function (k) {
                    cutree(prot_clust, k = k)
                })
                colnames(hc_grs) <- paste0(cmethod, "_k=", kmin:kmax)
                return(hc_grs)
            }
        }) %>% purrr::reduce(cbind)
        # save(prot_clust, all_clust_sols, file = "cached_results/protein_clustering.RData")
        save(all_clust_sols, file = "cached_results/protein_clustering.RData")
        
    }

    #-- Cluster rows
    sdist <- as.dist(1 - cor(t(impt_data), method = distance, use = "complete"))
    sample_clust <- hclust(sdist, method = "ward.D2")
    
    #-- Visual inspection
    if (make_plots) {
        if(!dir.exists(res_dir)) dir.create(res_dir)
        sapply(colnames(all_clust_sols), function(sol_name) {
            clust_sol <- all_clust_sols[,sol_name]
            annotatedHeatmap(t(prot_data$wp), 
                             sa = data.frame(id = names(clust_sol),
                                             clust = clust_sol), 
                             id_var = "id",
                             group_var = "clust",
                             silent = TRUE,
                             filename = filePath(res_dir, paste0("protein_clustering_", sol_name, ".pdf")),
                             width = 6,
                             height = 6,
                             cutree_col = max(clust_sol),
                             show_rownames = FALSE)
        }) 
    }
    return(list(protein_clusts = all_clust_sols, sample_clust = sample_clust))
}
enrich_protClusters <- function(prot_clusters, prot_data2, path_table,
                                type = "chosen_solution", solution = NULL) {
    if(type == "all_solutions") {
        prot_cluster_tb <- tibble(protein_id = rownames(prot_clusters),
                                  pclust = prot_clusters[,solution]) %>%
            left_join(prot_data2$wpAnnot) %>%
            select(protein_id, symbol, pclust)
        #-- ORA
        k <- max(prot_clusters)
        pclust_enrich <- lapply(1:k, function(i) {
            prots <- prot_cluster_tb %>% filter(pclust == i) %>%
                pull(symbol)
            pathwayORA(prots, path_table) %>%
                filter(padj < 0.05)
        })
        names(pclust_enrich) <- paste0("cluster", 1:k)
        
        pclust_df <- lapply(names(pclust_enrich), function(pc_name) {
            pclust_enrich[[pc_name]] %>% mutate(cluster = pc_name,
                                                solution = solution)
        }) %>% bind_rows()
        return(pclust_df)
    } else {
        prot_cluster_tb <- tibble(protein_id = names(prot_clusters),
                                  pclust = prot_clusters) %>%
            left_join(prot_data2$wpAnnot) %>%
            select(protein_id, symbol, pclust)
        #-- ORA
        k <- max(prot_clusters)
        pclust_enrich <- lapply(1:k, function(i) {
            prots <- prot_cluster_tb %>% filter(pclust == i) %>%
                pull(symbol)
            pathwayORA(prots, path_table) %>%
                filter(padj < 0.05)
        })
        names(pclust_enrich) <- paste0("cluster", 1:k)
        return(pclust_enrich)
    }
    
}

plot_proteomeFull <- function(prot_data2, prot_clusters, chosen_solution = "pam_k=8") {
    prot_groups <- prot_clusters
    prot_groups <- split(names(prot_groups), prot_groups)
    var <- prot_data2$sampAnnot %>% select(id, uPG = CCP_cluster)
    
    wp <- prot_data2$wp
    wp4plot <- wp %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("id") %>%
        left_join(var) %>%
        tibble() %>%
        group_by(uPG)
    
    hm <- ggheatmap(wp4plot,
                    colv = "id",
                    hm_colors = viridis::inferno(100),
                    hm_color_limits = c(-1.5,1.5),
                    hm_color_values = scales::rescale(c(-1.5,-0.75,-0.5,-0.25,-0.1,0,0.1,0.25,0.5,0.75,1.5)),
                    colors_title = "Protein expression",
                    rowv = prot_groups,
                    show_dend_row = FALSE,
                    show_rownames = FALSE,
                    show_colnames = FALSE,
                    raster = TRUE,
                    group_label = TRUE,
                    group_track = TRUE,
                    group_lines = TRUE,
                    group_line_color = "white",
                    group_prop = 0.05)
    
    ggsave("results/suppfig_ccp/d.pdf", 
           hm,
           width = 6,
           height = 4)
    
    return(hm)
    
}

proteinPathNetwork <- function(prot_clusters, prot_data, path_table,
                           ng = 8, solution = "pam_k=8", gr_cols = ggsci::pal_npg()(8),
                           top_n_community = 2) {
    wpAnnot <- prot_data$wpAnnot
    #-- Calculate additional parameters
    mads <- apply(prot_data$wp, 1, mad, na.rm = TRUE)
    mads <- tibble(protein_id = names(mads), mad = mads)
    path_n <- path_table %>%
        filter(Symbol %in% wpAnnot$symbol) %>%
        group_by(Symbol) %>%
        select(Path_ID, symbol = Symbol) %>%
        count() %>%
        ungroup() %>%
        rename(path_n = n)
    
    #-- Get groups
    pgrs <- prot_clusters$protein_clusts[,solution]
    pgr_ls <- split(names(pgrs), pgrs)
    
    # #-- Make group functional networks
    # prot_groups_network <- map2(pgr_ls, gr_cols, .group_network, wpAnnot, mads, path_n,
    #                             path_table, top_n_community)
    # if(!dir.exists("results/protein_clustering/")) dir.create("results/protein_clustering/")
    # save(prot_groups_network, file = "results/protein_clustering/protein_group_networks.RData")
    
    #-- Overall network
    set.seed(0)
    total_network <- .overall_network(wpAnnot, path_table, pgrs)
    
    #-- Save plots and tables
    ggsave(plot = total_network$total_netplot, 
           filename = "results/suppfig_clusterprots/c.pdf",
           width = 6, height = 3.5)
    #-- Return
    return(total_network)
    
}

.overall_network <- function(wpAnnot, path_table, pgrs) {
    #-- Get protein group network
    wpAnnot[names(pgrs),"Prot_Group"] <- pgrs 
    
    #-- Get network paths for all prots
    total_paths <- wpAnnot %>%
        right_join(select(path_table, Path_ID, symbol = Symbol)) %>%
        select(symbol1 = symbol, Path_ID, group1 = Prot_Group)
    
    #-- Make group adjancency table
    total_adj <- total_paths %>%
        left_join(rename(total_paths, symbol2 = symbol1, group2 = group1), by = "Path_ID") %>%
        filter(symbol1 != symbol2, 
               group1 != group2) %>%
        group_by(group1, group2) %>%
        summarize(weight = n()) %>%
        rowwise() %>%
        mutate(gpair = paste(sort(c(group1, group2)), collapse = "~")) %>%
        ungroup() %>%
        filter(!duplicated(gpair)) %>%
        select(-gpair) %>%
        arrange(group1)
    
    #-- Make network and annotate
    tnet <- network(select(total_adj, group1, group2))
    tnet %e% "weight" <- total_adj$weight
    set.seed(0)
    tndf <- ggnetwork(tnet, weights = "weight")
    
    #-- Network plot
    tnet_plot <- ggplot(tndf, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(aes(size = weight, alpha = weight), color = "grey70", curvature = 0.2) +
        geom_nodes(size = 10, color = "grey50") +
        geom_text(aes(label = vertex.names)) +
        scale_size_continuous(range = c(0.3, 5)) +
        labs(size = "Number of connections", alpha = "Number of connections") +
        theme_blank()
    
    return(list(total_network = tnet, total_netplot = tnet_plot))
}

upg_pathConnections <- function(pathScores_uPG_KEGG) {
    #-- Get pathway connections
    upgs <- LETTERS[1:5]
    upg_combos <- combn(upgs, 2) %>%  t()
    upg_combo_info <- apply(upg_combos, 1, .get_shared_paths, pathScores_uPG_KEGG) %>% 
        bind_rows() %>%
        filter(nshared > 0)
    
    #-- Add manual annotation
    edge_annot <- c("↑ Ribosomal, Cell Cycle\n↓ Oxidative Phosphorylation, Peroxisome",
                    "↑ ECM interactions",
                    "↓ ECM",
                    "↓ Smooth muscle",
                    "↑ Fatty acid metabolism\n↓ Spliceosome, Ribosome",
                    "↓ DNA replication",
                    "↓ Cell cycle, DNA repair")
    
    upg_combo_info <- upg_combo_info %>%
        mutate(edge_annot = edge_annot)
    
    net <- network(select(upg_combo_info, upg1, upg2, nshared, edge_annot))
    
    #-- Plot network
    set.seed(2); upg_net <- ggnetwork(net) %>%
        ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(aes(size = nshared), color = "grey60", alpha = 0.5) +
        geom_nodes(aes(color = vertex.names), size = 10) +
        geom_nodetext(aes(label = vertex.names)) +
        geom_edgetext(aes(label = edge_annot), size = 3, alpha = 0) +
        guides(color = "none") +
        labs(size = "Number of shared\nKEGG pathway activities") +
        theme_void() +
        theme(legend.title = element_text(size = 8, face = "bold"),
              legend.text = element_text(size = 8))
    
    ggsave("results/suppfig_clusterprots/b.pdf", upg_net, device=cairo_pdf, width = 5, height = 3)
    
    return(upg_net)
}

.get_shared_paths <- function(combo, pathScores) {
    #-- Get paths
    c1_paths <- pathScores[[combo[1]]] %>% filter(padj < 0.05)
    c2_paths <- pathScores[[combo[2]]] %>% filter(padj < 0.05)
    #-- Shared up
    c1_up <- c1_paths %>% filter(NES > 0) %>% pull(Path_ID)
    c2_up <- c2_paths %>% filter(NES > 0) %>% pull(Path_ID)
    shared_up <- intersect(c1_up, c2_up)
    #-- Shared down
    c1_down <- c1_paths %>% filter(NES < 0) %>% pull(Path_ID)
    c2_down <- c2_paths %>% filter(NES < 0) %>% pull(Path_ID)
    shared_down <- intersect(c1_down, c2_down)
    
    nshared <- length(shared_up) + length(shared_down)
    
    tibble(upg1 = combo[1], upg2 = combo[2], nshared = nshared, 
           shared_up = paste(shared_up, collapse = ", "),
           shared_down = paste(shared_down, collapse = ", "))
}

protClustChoice <- function(enrich_all) {
    path_unique_props <- lapply(enrich_all, function(solution) {
        enriched_paths <- split(solution$Path_ID, solution$cluster)
        all_cls <- unique(solution$cluster)
        unique_prop <- sapply(all_cls, function(i) {
            other_cls <- setdiff(all_cls, i)
            cl_paths <- enriched_paths[[i]]
            other_paths <- purrr::reduce(enriched_paths[other_cls], union)
            n_unique <- length(setdiff(cl_paths, other_paths))
            n_unique/length(cl_paths)
        })
        tibble(cl = all_cls, path_unique_prop = unique_prop, solution = solution$solution[1])
    }) %>% bind_rows() %>%
        mutate(k = str_remove(solution, "pam_k=") %>% factor(levels = as.character(5:15)),
               k_num = as.numeric(k))
    prop_meds <- path_unique_props %>%
        group_by(k) %>%
        summarize(path_unique_prop = median(path_unique_prop))
    
    ggplot(path_unique_props, aes(k, path_unique_prop)) +
        geom_violin(aes(fill = k)) +
        geom_quasirandom() +
        ggforce::geom_bspline0(aes(group=1), data = prop_meds) +
        stat_summary(geom = "crossbar", fun = "median", size = 0.3) +
        scale_fill_viridis_d() +
        theme_csg_scatter +
        scale_y_continuous(labels = scales::percent) +
        labs(y = "Proportion of uniquely enriched pathways") +
        guides(fill = "none")

}

export_protAnnot <- function(unfilt_prot, prot_data, prot_clusters) {
    filtered <- rownames(prot_data$wp)
    prot_cls <- prot_clusters$protein_clusts[,"pam_k=8"]
    prot_cls_tb <- tibble(protein_id = names(prot_cls), protein_cluster = prot_cls) %>%
        mutate(protein_cluster = factor(protein_cluster))
    
    n3p <- tibble(nsamples_3peptides = rowSums(unfilt_prot$pep > 3, na.rm = TRUE),
           long_id = rownames(unfilt_prot$pep))
    
    wp_annot <- unfilt_prot$wpAnnot
    wp_annot %>%
        tibble() %>%
        rename(long_id = ID) %>%
        left_join(n3p) %>%
        mutate(in_filtered = protein_id %in% filtered) %>%
        left_join(prot_cls_tb) %>%
        write_csv("results/tables/CITproteomics_proteinAnnotation.csv")
}
getClusters <- function(prot_clusters, version = "pam_k=8", reorder = c(5,2,1,6,4,3,7,8)) {
    cls <- prot_clusters$protein_clusts[,version]
    remapped <- plyr::mapvalues(cls, 1:8, reorder)
    return(remapped)
}

plot_silhouetteWidth <- function(multiomics_hm, prot_data) {
  samp_order <- get_colLevels(multiomics_hm$mofa_hm)
  # samp_dist <- dist(t(prot_data$wp[,samp_order]))
  samp_dist <- as.dist(1 - cor(prot_data$wp[,samp_order], method = "spearman", use = "complete"))
  ccp_cls <- pull(prot_data$sampAnnot, CCP_cluster, id)
  samp_dist_mat <- as.matrix(samp_dist)
  
  sil_width <- silhouette(as.numeric(ccp_cls[samp_order]), samp_dist) %>%
    as.data.frame() %>%
    mutate(id = samp_order)
  
  sil_width4plot <- 
    sil_width %>%
    select(id, sil_width) %>%
    left_join(select(prot_data$sampAnnot, id, CCP_cluster), by = "id") %>%
    group_by(CCP_cluster) %>%
    mutate(id = fct_reorder(id, sil_width, .desc = TRUE)) %>%
    arrange(id)
  
  cl_summary <- sil_width4plot %>%
    ungroup() %>%
    mutate(n = 1:n()) %>%
    group_by(CCP_cluster) %>%
    summarize(max_n = min(n) - 0.5,
              mean_n = mean(n))
  cl_lines <- filter(cl_summary, CCP_cluster != "A")
  plt_silhouette <- ggplot(sil_width4plot) +
    geom_col(aes(id, sil_width, fill = CCP_cluster), width = 1) +
    geom_segment(aes(x = max_n, xend = max_n, y = -Inf, yend = Inf, group = CCP_cluster), 
                 lty = "dotted",
                 data = filter(cl_lines, CCP_cluster != "A")) +
    geom_hline(aes(yintercept = 0), lty = "dotted") +
    geom_text(aes(x = mean_n, y = 0.45, label = CCP_cluster), data = cl_summary,
              size = 3, fontface = "bold") +
    scale_y_continuous(limits = c(-.1,0.5), breaks = c(-0.1,0,0.1,0.25,0.5)) +
    labs(x = "", y = "Silhouette width") +
    guides(fill = "none") +
    theme_csg_sparse +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(color = "black", fill = NA))
  plt_silhouette
  
  ggsave("results/suppfig_ccp/fig_silhouette.pdf", plt_silhouette, width = 2.5, height = 2)
}


