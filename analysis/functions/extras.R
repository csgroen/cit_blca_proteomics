exportProtTables <- function(prot_data, unfilt_prot, mcp_res, prot_clusters,
                             multiomics_hm) {
    
    # Sample metadata
    #-- Prepare additional data
    mcp_data <- mcp_res %>%
        rownames_to_column("sample_id") %>%
        rename_with(.cols = -sample_id, 
                    .fn = ~ str_replace_all(paste0("MCPcounter_", .), " ", "_"))
    
    mofa_data <- multiomics_hm$mofa_hm$data %>%
        select(sample_id = observations, `Ta pathway`:`Immune system`) %>%
        rename_with(.cols = -sample_id, 
                    .fn = ~ str_replace_all(paste0("MOFAscore_", .), " |/|\n", "_"))
    
    biton_annot <- read_excel("data/1-s2.0-S2211124714009048-mmc4.xlsx")
    sample_table <- prot_data$sampAnnot %>%
        select(-age, -Sex, -sex) %>%
        left_join(select(biton_annot, id = Sample, age, sex = gender,
                         histological_type = histologicalType), by = "id") %>%
        select(
            sample_id = id,
            stage = Stage,
            stage_inv = Stage.,
            grade = Grade,
            node, 
            age,
            sex,
            histological_type,
            class_Consensus = Consensus,
            class_UROMOL = NMIBCclass,
            class_Subtype = Subtype,
            # class_Lund = Lund,
            class_uPG = CCP_cluster,
            TP53_mutation,
            FGFR3_mutation,
            ends_with("mutation"),
            TP53_mutation_type,
            FGFR3_mutation_type,
            CDKN2A_loss,
            RB1_loss,
            PPARG_gistic_amp = PPARG_amp,
            starts_with("IC ")) %>%
        mutate(
        across(matches("mutation$|loss$|amp$"), function(x) {
            case_when(x == "Loss" ~ "Yes",
                      x == "Amp" ~ "Yes",
                      x == "M" ~ "Yes",
                      x == "WT" ~ "No",
                      x == "Yes" ~ "Yes",
                      x == "No" ~ "No",
                      is.na(x) ~ "NA")
        }),
        class_UROMOL = ifelse(class_UROMOL == "Class 1", "Class 1/3", class_UROMOL)) %>%
        left_join(mcp_data) %>%
        left_join(mofa_data) %>%
        rename_with(~ str_replace(., "_mutation$", "") %>% paste0("mutation_", .),
                    .cols = matches("mutation$")) %>%
        rename_with(~ str_replace(., "_mutation_type$", "") %>% paste0("mutation_type_", .),
                    .cols = matches("_mutation_type$")) %>%
        rename_with(~ paste0("cna_", .),
                    .cols = matches("loss$|amp$")) %>%
        rename_with(~ str_replace_all(., " ", "_")) %>%
        rename_with(~ str_replace_all(., "IC_", "ic_")) %>%
        rename_with(~ str_replace_all(., "MOFAscore", "mofaScore")) %>%
        rename_with(~ str_replace_all(., "MCPcounter", "mcpCounter")) %>%
        tibble()
    
    openxlsx::write.xlsx(sample_table, "results/tables/CITproteomics_samplesMetadata.xlsx")
    
    # Protein expression
    write.csv(prot_data$wp, "results/tables/CITproteomics_proteinExpression.csv")
    write.csv(unfilt_prot$wp, "results/tables/CITproteomics_unfilteredProteinExpression.csv")

    # Peptide counts
    wp_annot <- unfilt_prot$wpAnnot
    pep <- unfilt_prot$pep[wp_annot$ID,]
    rownames(pep) <- wp_annot$protein_id
    write.csv(unfilt_prot$pep, "results/tables/CITproteomics_peptideCounts.csv")
    
    # Protein annotation
    #-------- Added annotation
    median_peps <- matrixStats::rowMedians(pep, na.rm = TRUE)
    expr_prop <- apply(unfilt_prot$wp[wp_annot$protein_id,], 1, 
                       function(prot_peps) sum(!is.na(prot_peps)))/63
    
    
    wp_annot_tb <- wp_annot %>%
        tibble() %>%
        select(-ID) %>%
        mutate(passed_filter = protein_id %in% rownames(prot_data$wp)) %>%
        left_join(tibble(
            protein_id = names(prot_clusters),
            protein_cluster = prot_clusters)) %>%
        mutate(median_npeptides_quantification = median_peps,
               prop_samples_expr = expr_prop)
    write_csv(wp_annot_tb, "results/tables/CITproteomics_proteinAnnotation.csv")
    
    return(TRUE)
}

exportTableFromList <- function(list, file) {
    names(list) <- str_replace_all(names(list), "\\/|\\-", ".")
    list <- lapply(list, as.data.frame)
    openxlsx::write.xlsx(list, file = file, overwrite = TRUE, na.string = "NA")
}

export_enrichprotClust <- function(enrich_prot_clusters_RDB, enrich_prot_clusters_HM,
                                   enrich_prot_clusters_KEGG) {
    joined_tables <- lapply(1:length(enrich_prot_clusters_HM), function(i) {
        enrich_prot_clusters_HM[[i]] %>%
            mutate(Database = "Hallmarks") %>%
            bind_rows(mutate(enrich_prot_clusters_KEGG[[i]], Database = "KEGG")) %>%
            bind_rows(mutate(enrich_prot_clusters_RDB[[i]], Database = "Reactome")) %>%
            tibble()
    })
    names(joined_tables) <- names(enrich_prot_clusters_HM)
    
    openxlsx::write.xlsx(joined_tables, "results/tables/proteinClusters_pathways.xlsx")
    return(joined_tables)

    
}
export_upg_enrichment <- function(diffProt_uPG, pathScores_uPG, pathScores_uPG_H,
                                  pathScores_uPG_KEGG, exp_ic_table) {
    
    react <- .join_pathscores(pathScores_uPG)
    hallmark <- .join_pathscores(pathScores_uPG_H) %>% select(-Path_ID)
    kegg <- .join_pathscores(pathScores_uPG_KEGG) %>% select(-Path_ID)
    
    all_tables <- list(
        `Differential proteins` = diffProt_uPG,
        Reactome = react,
        Hallmark = hallmark,
        KEGG = kegg,
        `Transcriptomic IC` = exp_ic_table
    )
    
    openxlsx::write.xlsx(all_tables, "results/tables/diffProteomics_uPG_paths.xlsx")
    return(all_tables)
}

export_fgfr3_enrichment <- function(diffProt_FGFR3, degtable_fgfr3, mrnaPaths_FGFR3, 
                                    protPaths_FGFR3, deltaPathway_FGFR3) {
    tables <- list(`Differential proteins` = diffProt_FGFR3,
                   `Differential genes` = degtable_fgfr3, 
                   `Enrichment scores (mRNA)` = mrnaPaths_FGFR3,
                   `Enrichment scores (proteomics)` = protPaths_FGFR3,
                   `Delta ES` = deltaPathway_FGFR3$dNES
    )
    openxlsx::write.xlsx(tables, file = "results/tables/diffProteomics_FGFR3_paths.xlsx")
    
    return(tables)
}

get_ic_cors <- function() {
    load("data/annotations/biton_cors.RData")
    sig_cors %>%
        select(IC = sig_name, gene_symbol = Gene, `Spearman_cor` = pm_rho, direction = weight_dir)
}

.join_pathscores <- function(tb) {
    lapply(LETTERS[1:5], function(upg) { 
        tb[[upg]] %>%
            mutate(uPG = upg) %>%
            relocate(uPG)
    }) %>% bind_rows()
}

