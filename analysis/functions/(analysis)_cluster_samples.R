# library(ConsensusClusterPlus)
run_sample_CCP <- function(prot_data, k, linkage, cluster_alg = "hc", nfeats = Inf,
                           dist = "spearman", verbose = TRUE, 
                           out_dir = "cached_results/sample_clustering/", 
                           ...) {
    current_dir <- getwd()
    if(!dir.exists(out_dir)) dir.create(out_dir)
    setwd(out_dir)
    #-- Create dir name
    dir_name <- case_when(
        is.infinite(nfeats) & cluster_alg == "hc" ~ paste("CCP_k", k, linkage, dist, sep = "_"),
        is.infinite(nfeats) & cluster_alg != "hc" ~ paste("CCP_k", k, cluster_alg, dist, sep = "_"),
        is.finite(nfeats) & cluster_alg == "hc" ~ paste("CCP_k", k, linkage, dist, nfeats, "feats", sep = "_"),
        TRUE ~ paste("CCP_k", k, cluster_alg, dist, nfeats, "feats", sep = "_")
    )
    #-- Select features
    if(is.infinite(nfeats)) {
        wp <- prot_data$wp
    } else {
        vars <- apply(prot_data$wp, 1, mad, na.rm = TRUE) %>% sort(decreasing = TRUE)
        topn <- names(vars)[1:nfeats]
        wp <- prot_data$wp[topn,]
    }
    #-- Run CCP
    if(cluster_alg == "hc") {
        if(!dir.exists(dir_name)) {
            message("Running CCP on ", dir_name)
            ccp <- ConsensusClusterPlus(d = prot_data$wp, maxK = k, 
                                        reps = 1000, pItem = 0.8, pFeature = 0.95,
                                        innerLinkage = linkage, finalLinkage = linkage, distance = dist,
                                        corUse = "pair", writeTable = TRUE, plot = "pdf", title = dir_name, 
                                        ...)
        } 
    } else {
        if(!dir.exists(dir_name)) {
            ccp <- ConsensusClusterPlus(d = prot_data$wp, maxK = k, 
                                        reps = 1000, pItem = 0.8, pFeature = 0.95,
                                        innerLinkage = linkage, finalLinkage = linkage, distance = dist,
                                        corUse = "pair", writeTable = TRUE, plot = "pdf", title = dir_name,
                                        clusterAlg = cluster_alg,
                                        ...)
        }
    }
    setwd(current_dir)
    return(out_dir)
}

source("aux_functions/doAssociationTests.R")
upg_assoc_tests <- function(sa, mcp_res_tb) {
    mcp_df <- mcp_res_tb %>%
        as.data.frame() %>%
        rownames_to_column("sample")
    
    sa <- sa %>%
        select(sample, uPG, matches("mutation$|gain_bin$|loss_bin$|Regulon$|^IC "), 
               CDKN2A_loss, RB1_loss, PPARG_gain, Genomic_Instability, Lund) %>%
        mutate(across(matches("mutation"), ~ ifelse(. == "Yes", "Mut", "WT"))) %>%
        left_join(mcp_df, by = "sample")
    assoc_tests <- doAssociationTests(sa, id_var = "sample", test_var = "uPG")
    return(assoc_tests)
}

consensus_nmibc_mibc_clusters_compared <- function(samp_annot_63) {
    nmibc_cls <- read_csv("extras/nmibc_ccp_res/CCP_k_8_pam_spearman/CCP_k_8_pam_spearman.k=4.consensusClass.csv", 
             col_names = c("id", "new_CCP")) %>%
        mutate(new_CCP = paste0("nmibc_", new_CCP),
               stage = "NMIBC")
    
    nmibc_k3 <-  read_csv("extras/nmibc_ccp_res/CCP_k_8_pam_spearman/CCP_k_8_pam_spearman.k=3.consensusClass.csv", 
                          col_names = c("id", "new_CCPk3")) %>%
        mutate(new_CCPk3 = paste0("nmibc_", new_CCPk3),
               stage = "NMIBC")
    
    left_join(nmibc_cls, nmibc_k3) %>%
        group_by(new_CCP, new_CCPk3) %>%
        summarize(n = n())
    
    mibc_cls <- read_csv("extras/mibc_ccp_res/CCP_k_8_pam_spearman/CCP_k_8_pam_spearman.k=4.consensusClass.csv", 
                          col_names = c("id", "new_CCP")) %>%
        mutate(new_CCP = paste0("mibc_", new_CCP),
               stage = "MIBC")
    
    mibc_pal <- structure(c("#FB6A4A", "#E25477", "#B05592", "#745993"), names = paste0("mibc_", 1:4))
    nmibc_pal <- structure(c("#67A9CF", "#3BC1D9", "#2CD6CE", "#67E8B2"), names = paste0("nmibc_", 1:4))
    
    
    new_cls <- bind_rows(nmibc_cls, mibc_cls)  %>%
        left_join(select(samp_annot_63, id, uPG = CCP_cluster, Consensus, NMIBCclass)) %>%
        mutate(new_CCP = factor(new_CCP, levels = c("mibc_1", "mibc_2", "nmibc_3", "mibc_3",
                                                    "nmibc_1", "mibc_4", "nmibc_2", "nmibc_4")))
    
    plt_nmibc_pr <- new_cls %>%
        filter(stage == "NMIBC") %>% 
        classFlow("uPG", "new_CCP", 
              label_colors = c(ccp_cols, mibc_pal, nmibc_pal),
              alpha = 0.5, label_size = 2.5) +
        guides(fill = "none") +
        labs(x = "Class type") +
        theme_csg_sparse
    
    plt_mibc_pr <- new_cls %>%
        filter(stage == "MIBC") %>% 
        classFlow("uPG", "new_CCP", 
                  label_colors = c(ccp_cols, mibc_pal, nmibc_pal),
                  alpha = 0.5, label_size = 2.5) +
        guides(fill = "none") +
        labs(x = "Class type") +
        theme_csg_sparse
    
    plt_nmibc_tr <- new_cls %>%
        filter(stage == "NMIBC") %>%
        mutate(NMIBCclass = ifelse(NMIBCclass == "Class 1", "Class 1/3", NMIBCclass) %>%
                   factor(levels = c("Class 2a", "Class 2b", "Class 1/3")),
               new_CCP = factor(new_CCP, levels = c("nmibc_3", "nmibc_1", "nmibc_2", "nmibc_4"))) %>%
        classFlow("NMIBCclass", "new_CCP", alpha = 0.5, label_size = 2.5,
                  label_colors = c("Class 1/3" = "#80bf7c", "Class 2a" = "#790d16", "Class 2b" = "#b57824",
                                   nmibc_pal)) +
        guides(fill = "none") +
        labs(x = "Class type") +
        theme_csg_sparse
    
    plt_mibc_tr <- new_cls %>%
        filter(stage == "MIBC") %>%
        mutate(new_CCP = factor(new_CCP, levels = c("mibc_1", "mibc_2", "mibc_3", "mibc_4")),
               Consensus = factor(Consensus, levels = c("Ba/Sq", "LumU", "NE-like", 
                                                        "LumNS", "Stroma-rich", "LumP"))) %>%
        classFlow("Consensus", "new_CCP", alpha = 0.5, label_size = 2.5,
                  label_colors = c(consensus_pretty,
                                  mibc_pal)) +
        guides(fill = "none") +
        labs(x = "Class type") +
        theme_csg_sparse
    
    plt_left <- ((plt_mibc_pr + labs(subtitle = "MIBC") + theme(axis.title.x = element_blank())) / 
                     plt_nmibc_pr + labs(subtitle = "NMIBC")) +
        plot_layout(heights = c(1.5,1))
    plt_right <- ((plt_mibc_tr + theme(axis.title.x = element_blank())) / plt_nmibc_tr) +
        plot_layout(heights = c(1.5,1))
    
    plt_left | plt_right
    

}

