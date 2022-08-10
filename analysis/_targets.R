# install.packages(c("targets"))
# see: https://books.ropensci.org/targets/
library(targets)
setwd("~/Desktop/CITproteomics")

R.utils::sourceDirectory("functions", modifiedOnly = FALSE)

list(
    # Load data ----------
    tar_target(mrna_data, load_file("data/CIT_mrna_data.RData")),
    tar_target(prot_data, load_file("data/CIT_proteomics_data.RData")),
    tar_target(unfilt_prot, load_file("data/CIT_unfiltered_proteomics.RData")),
    tar_target(lowfilt_prot, load_file("data/CIT_lowfilt_proteomics.RData")),
    tar_target(prot_symbol, load_file("data/CIT_proteomics_symbol.RData")),
    tar_target(prot_data_bystage, load_file("data/CIT_proteomics_dataByStage.RData")),
    tar_target(cna_calls, read_tsv("data/GISTIC/CN_geneCalls_GISTIC.txt")),
    tar_target(cna_gr, load_file("data/CIT_GenomicRanges_CNAsmooth.RData")),
    tar_target(cna_mat, load_file("data/CIT_CNAmatrix.RData")),
    tar_target(react_paths, load_file("data/annotations/Reactome_pathways.RData")),
    tar_target(gobp_paths, load_file("data/annotations/GOBP_pathways.RData")),
    tar_target(hallmark_paths, load_file("data/annotations/Hallmark_pathways.RData")),
    tar_target(kegg_paths, load_file("data/annotations/KEGG_pathways.RData")),
    
    tar_target(samp_annot, load_file("data/CIT_sampleAnnot1.RData")),
    tar_target(samp_annot_63, load_file("data/CIT_sampleAnnot2.RData")),
    tar_target(core_samples, load_file("data/CIT_coreSamples.RData")),
    
    tar_target(all_samples, colnames(prot_data$wp)),
    tar_target(mcp_res, mcpCounter(mrna_data, all_samples)),
    
    # Supp Figure: QC --------------
    tar_target(nprot_filter, plot_protFilter(unfilt_prot)), #A
    tar_target(prot_ovrl, protein_overlap(prot_data, prot_data_bystage)), #B
    tar_target(physchem_plot, plot_physChemHist()), #C
    tar_target(dynamic_range_prot, prot_function_dynamicRange(unfilt_prot)), #D
    tar_target(subcellCounts_plot, plot_subCell_counts()), #E
    tar_target(complex_cors, plot_pp_pm_complex_cors(prot_data, prot_mrna_cor)), #F

    
    # Figure 1: Data description ------------------
    ## B) Protein/mRNA correlation ------------------------------
    tar_target(prot_mrna_cor, calculate_protmrna_cor(prot_data, mrna_data)),
    tar_target(prot_mrna_pathdist, 
               path_protmrna_cor(prot_mrna_cor, path_table = react_paths)),
    tar_target(pmcor_plots, plot_protmrna_cor(
        prot_mrna_cor, prot_mrna_pathdist, react_paths,
        paths4plot = c("R-HSA-72766", "R-HSA-72312", "R-HSA-1428517",
                       "R-HSA-72172", "R-HSA-8953854", "R-HSA-169911",
                       "R-HSA-8978868", "R-HSA-909733", "R-HSA-6805567"))),

    ## C-E) CNA/mRNA | CNA/Prot correlations --------------
    tar_target(cna_cor_res, cna_cors(cna_gr, mrna_data, prot_data, samp_annot)),
    tar_target(cna_cors_compare, compare_cna_cors(cna_cor_res, cna_gr)), # C
    tar_target(cna_cors_genome, cna_cors_genomePlots(cna_cor_res, cna_calls, core_samples)), #E
    tar_target(cna_freq, cna_freqAlt_compare(cna_cors_genome)), #D
    
    
    ## Sample clustering: ConsensusClusterPlus -----------------
    tar_target(ccp_linkage, c("ward.D2", "complete")),
    tar_target(ccp_cluster_alg, c("hc", "pam")),
    tar_target(ccp_nfeats, c(Inf, 1000)),
    tar_target(ccp_runs,
             run_sample_CCP(prot_data, k = 8, ccp_linkage, ccp_cluster_alg, ccp_nfeats),
             pattern = cross(ccp_linkage, ccp_cluster_alg, ccp_nfeats)),
    tar_target(upg_assocs, upg_assoc_tests(samp_annot, mcp_res)),

    # Supp Figure: Protein clustering --------------------
    tar_target(prot_clusters,
               clusterProteins(prot_data, kmin = 5, kmax = 15,
                               res_dir = "results/protein_clustering")),
    tar_target(pc_solutions, colnames(prot_clusters$protein_clusts)),
    tar_target(enrich_all_RDB,
               enrich_protClusters(prot_clusters$protein_clusts, prot_data, react_paths, solution = pc_solutions,
                                   type = "all_solutions"),
               pattern = map(pc_solutions), iteration = "list"),
    tar_target(prot_clust_choice, protClustChoice(enrich_all_RDB)),
    tar_target(prot_clusters_final, getClusters(prot_clusters, version = "pam_k=8")),
    tar_target(enrich_prot_clusters_RDB,
               enrich_protClusters(prot_clusters_final, prot_data, react_paths)),
    tar_target(enrich_prot_clusters_HM,
               enrich_protClusters(prot_clusters_final, prot_data, hallmark_paths)),
    tar_target(enrich_prot_clusters_KEGG,
               enrich_protClusters(prot_clusters_final, prot_data, kegg_paths)),
    tar_target(hm_fullProteome, plot_proteomeFull(prot_data, prot_clusters_final)), # D
    
    # MultiOmics Factor Analysis --------------------------------
    ## MOFA params -------------------------
    tar_target(n_factors, c(10,15)),
    tar_target(n_cna, c(1500, 5000)),
    tar_target(n_transcripts, c(5e3,1e4)),
    ## Run MOFA ----------------------------
    tar_target(mofa_res,
               runMOFA(mrna_data$gexp, prot_symbol, cna_mat,
                       core_samples = core_samples, n_transcripts = n_transcripts,
                       n_cna = n_cna, n_factors = n_factors),
                pattern = cross(n_transcripts, n_cna, n_factors)),
    tar_target(mofa_res_list,
               conformMOFAres(mofa_res)),
    tar_target(interpret_mofa,
               characterizeMOFA(mofa_res_list, samp_annot)),
    ## Add annotations --------------------------
    tar_target(final_mofa,
               finalizeMOFA(mofa_select = "mofa_mrna10000_cna1500_10f",
                            sampAnnot = samp_annot,
                            factors2show = paste0("Factor", c(1,2,3,5,6)),
                            factor_names = c("Ta pathway",
                                             "DNA repair",
                                             "Stroma",
                                             "Lipid metabolism",
                                             "Immune system"),
                            invert = c(1,1,-1,1,-1)
                            )),
    tar_target(exp_mofa, export_mofa(final_mofa, interpret_mofa)),
    tar_target(exp_mofa2, export_mofaAssociations(interpret_mofa)),

    # Figure 2 -------------------
    ## A) Proteomics PCA -----------------------------
    tar_target(prot_pca, proteomicsPCA(prot_data)),
    tar_target(plt_prot_pca, plot_proteomicsPCA(prot_data,
                                                prot_pca, 
                                                annot_vars = c("Subtype", "Stage."), 
                                                gr_colors = subtype_cols)), #A
    tar_target(plt_prot_pca_loadings, plot_proteomicsPCA(prot_data,
                                                prot_pca, 
                                                annot_vars = c("Subtype", "Stage."), 
                                                loadings = c("S100A8", "TNC", "CD14",
                                                             "S100P", "ADIRF", "KRT7",
                                                             "KRT6A", "FLNC", "TAGLN",
                                                             "TOP2A", "CDK1"),
                                                add_ic_loadings = TRUE,
                                                fname = "results/suppfig_pca/a.pdf",
                                                gr_colors = subtype_cols,
                                                width = 6.5, height = 4)), # Supp 2
    ## B/E) uPG heatmap / MOFA/subtype -----------------
    tar_target(multiomics_hm, plot_finalMOFA(final_mofa, mcp_res, prot_symbol, samp_annot_63)),
    ## C) Association uPG/Subtype ---------------------
    tar_target(plt_upg_sbty, plot_subtype_upg_assocs(samp_annot_63)),
    ## D) uPG pathway associations ----------------------
    tar_target(plt_upg_assocs, plot_uPG_pathEnrich(pathScores_uPG, pathScores_uPG_KEGG,
                                                   pathScores_uPG_H, exp_ic_table)),
    
    # Figure 3 ------------------------
    ## A) FGFR3 differential protein expression ----------------------
    ### Proteomics -----------------
    tar_target(diffProt_FGFR3,
               protDifferentialAnalysis(prot_data, "FGFR3_mutation", p_cutoff = 1)),
    tar_target(pathScores_FGFR3, 
               pathEnrichmentFromFC(diffProt_FGFR3, react_paths, path_rep = 0.3)),
    tar_target(goBPScores_FGFR3, 
               pathEnrichmentFromFC(diffProt_FGFR3, gobp_paths,  path_rep = 0.3)),
    tar_target(keggScores_FGFR3, 
               pathEnrichmentFromFC(diffProt_FGFR3, kegg_paths,  path_rep = 0.3)),
    tar_target(hallmarkScores_FGFR3, 
               pathEnrichmentFromFC(diffProt_FGFR3, hallmark_paths,  path_rep = 0.3)),
    tar_target(protPaths_FGFR3,
               join_pFGFR3paths(pathScores_FGFR3, goBPScores_FGFR3, keggScores_FGFR3)),
    tar_target(exp_fgfr3_kegg, openxlsx::write.xlsx(keggScores_FGFR3, "results/tables/pathEnrichment_FGFR3_KEGG_0.3.xlsx")),
    tar_target(exp_fgfr3_h,  openxlsx::write.xlsx(hallmarkScores_FGFR3, "results/tables/pathEnrichment_FGFR3_Hallmarks_0.3.xlsx")),
    
    ### Transcriptomics ----------------
    tar_target(degtable_fgfr3, mRNADifferentialAnalysis_FGFR3(mrna_data, samples = all_samples)),
    tar_target(degtable_fgfr3_filt, degtable_fgfr3 %>% filter(symbol %in% rownames(prot_symbol))),
    tar_target(mRNA_pathScores_FGFR3, pathScores_fromDEG(degtable_fgfr3_filt, react_paths)),
    tar_target(mRNA_gobpScores_FGFR3, pathScores_fromDEG(degtable_fgfr3_filt, gobp_paths)),
    tar_target(mRNA_keggScores_FGFR3, pathScores_fromDEG(degtable_fgfr3_filt, kegg_paths)),

    tar_target(mrnaPaths_FGFR3, join_mFGFR3paths(mRNA_pathScores_FGFR3,
                                                 mRNA_gobpScores_FGFR3,
                                                 mRNA_keggScores_FGFR3)),
    tar_target(all_paths,
               bind_rows(react_paths,
                         mutate(bind_rows(gobp_paths, kegg_paths),
                                Entrez = as.character(Entrez))) %>%
                   bind_rows(kegg_paths)),
    ## A) Differential pathways (NES) --------------------
    tar_target(plt_FGFR3_pm_paths, pathScores_prot_mrna(deltaPathway_FGFR3)),
    
    ## B) Delta pathway (dNES) --------------
    tar_target(deltaPathway_FGFR3, deltaPathwayActivity(protPaths_FGFR3, mrnaPaths_FGFR3)),
    tar_target(exp_fgfr3, export_fgfr3_enrichment(diffProt_FGFR3, degtable_fgfr3, mrnaPaths_FGFR3, 
                                                  protPaths_FGFR3, deltaPathway_FGFR3)),
    tar_target(exp_dpath_fgfr3, openxlsx::write.xlsx(deltaPathway_FGFR3$dNES, file = "results/tables/SuppTable_deltaPathway_FGFR3.xlsx")),
    
    ## C) Apoptosis pathway heatmap -----------------
    tar_target(diffProt_FGFR3_lowfilt,
               protDifferentialAnalysis(lowfilt_prot, "FGFR3_mutation", p_cutoff = 1)),
    tar_target(apop_hm,
               apoptosisPathHeatmap(diffProt_FGFR3_lowfilt, degtable_fgfr3)),
    
    # Supp Figure: Validation --------------------
    tar_target(upg_centroids,
               get_upg_centroids(prot_symbol, samp_annot, diffProt_uPG)),
    tar_target(strog_cls,
               stroggilos_cls(prot_symbol, samp_annot_63)),
    tar_target(strog_prots,
               plot_stroggilosComp(prot_data, prot_symbol, strog_cls, samp_annot_63)),
    tar_target(upg_xu,
               classify_xu_data(upg_centroids)),
    tar_target(plt_xu,
               plots_xu_data(upg_xu, diffProt_uPG)),
    
    # Figure 4 ---------------------
    tar_target(trail_sensitivity, plot_trailSensitivity()), # A
    tar_target(fgfr3_erda, plot_rescue_erdafitinib()), # B
    tar_target(fgfr3_rescue, plot_rescue_siFGFR3()), # C
    tar_target(apop_genes_sifgfr3, plot_siFGFR3_apoptPath()), # D
    
    tar_target(trail_birina, plot_trailBirina()), #A
    tar_target(trail_birina_synergy, plot_trailBirina_synergy()), # B
    
    # Extras -------------------
    ## uPG differential ----------------
    tar_target(diffProt_uPG, protDifferentialAnalysis(prot_data, "CCP_cluster", p_cutoff = 1)),
    tar_target(pathScores_uPG, pathEnrichmentFromFC(diffProt_uPG, react_paths, path_rep = 0.3)),
    tar_target(pathScores_uPG_KEGG, pathEnrichmentFromFC(diffProt_uPG, kegg_paths, path_rep = 0)),    
    tar_target(pathScores_uPG_H, pathEnrichmentFromFC(diffProt_uPG, hallmark_paths, path_rep = 0)),
    tar_target(exp_ic_table, upg_ic_table(mrna_data, samp_annot_63)),
    
    tar_target(exp_upg, export_upg_enrichment(diffProt_uPG, pathScores_uPG, pathScores_uPG_H,
                                     pathScores_uPG_KEGG, exp_ic_table)),
    
    ## Supp Fig: MCP + volcanos ------------
    tar_target(upg_volcanos, plt_upg_volcanos(diffProt_uPG, react_paths, hallmark_paths, kegg_paths)),
    tar_target(mcp_pvals, plt_mcpCounter_pvals(samp_annot_63, mcp_res)),

    ## Subtype differential ---------------
    tar_target(diffProt_Subtype,
               protDifferentialAnalysis(prot_data, "Subtype", p_cutoff = 1)),
    tar_target(pathScores_Subtype, pathEnrichmentFromFC(diffProt_Subtype, react_paths, path_rep = 0.3)),
    tar_target(exp_diff_subtype, openxlsx::write.xlsx(diffProt_Subtype, 
                                                      file = "results/tables/diffProteomics_Subtype.xlsx.xlsx")),
    
    # Tables -------------------
    ## Data ---------------------
    tar_target(exp_prot, exportProtTables(prot_data, unfilt_prot, 
                                          mcp_res, prot_clusters_final,
                                          multiomics_hm)),
    ## QC + Cors ---------------
    tar_target(ic_cors, get_ic_cors()),
    tar_target(exp_pm, write.xlsx(
        list(
            Genes = rename(prot_mrna_cor, symbol = Gene, spearman_cor = pm_rho),
            Pathways = prot_mrna_pathdist,
            ICs = ic_cors
            ), 
                                file = "results/tables/protmRNACorrelation.xlsx",
                                overwrite = TRUE)),
    tar_target(exp_cna_cors, write.xlsx(
        rename(cna_cor_res, cna_prot_cor = prot, cna_mrna_cor = mrna) %>%
            select(-pos), 
        file = "results/tables/cnaCorrelations.xlsx",
        overwrite = TRUE)),
    ## PCA + Differential Protein Expression ---------------
    tar_target(exp_pca, export_pcaResults(prot_pca, prot_data, react_paths)),
    tar_target(exp_dp_upg, write.xlsx(diffProt_uPG, 
                                      "results/tables/fig2b_diffProteomics_uPG.xlsx")),
    tar_target(exp_pclust, export_enrichprotClust(enrich_prot_clusters_RDB, enrich_prot_clusters_HM,
                                                  enrich_prot_clusters_KEGG))
    # tar_target(exp_sifgfr3, export_siFGFR3_table(siFGFR3))
)
