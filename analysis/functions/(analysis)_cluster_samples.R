run_sample_CCP <- function(prot_data, k, linkage, cluster_alg = "hc", nfeats = Inf,
                           dist = "spearman", verbose = TRUE, 
                           out_dir = "cached_results/sample_clustering/", 
                           ...) {
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
        setwd(out_dir)
        if(!dir.exists(dir_name)) {
            message("Running CCP on ", dir_name)
            ccp <- ConsensusClusterPlus(d = prot_data$wp, maxK = k, 
                                        reps = 1000, pItem = 0.8, pFeature = 0.95,
                                        innerLinkage = linkage, finalLinkage = linkage, distance = dist,
                                        corUse = "pair", writeTable = TRUE, plot = "pdf", title = dir_name, 
                                        ...)
        } 
    } else {
        setwd("cached_results/sample_clustering/")
        if(!dir.exists(dir_name)) {
            ccp <- ConsensusClusterPlus(d = prot_data$wp, maxK = k, 
                                        reps = 1000, pItem = 0.8, pFeature = 0.95,
                                        innerLinkage = linkage, finalLinkage = linkage, distance = dist,
                                        corUse = "pair", writeTable = TRUE, plot = "pdf", title = dir_name,
                                        clusterAlg = cluster_alg,
                                        ...)
        }
    }
    setwd("../../")
    return(NULL)
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
