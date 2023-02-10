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
    path_annot <- read_excel("data/CIT_pathUpdates.xlsx", na = "NA")
    # path_annot2 <- .get_path_annot()
    sample_table <- prot_data$sampAnnot %>%
        select(-age, -Sex, -sex) %>%
        left_join(select(biton_annot, id = Sample, age, sex = gender), by = "id") %>%
        left_join(path_annot, by = c(id = "sample_id")) %>%
        select(
            sample_id = id,
            stage = Stage,
            stage_inv = Stage.,
            grade = Grade,
            node, 
            age,
            sex,
            histological_subtype,
            papillary:endophytic,
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
            histological_subtype = case_when(
                histological_subtype == "classic" ~ "conventional",
                histological_subtype == "urothelial with squamous component" ~ "urothelial with squamous differentiation",
                TRUE ~ histological_subtype),
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
        # left_join(path_annot2, by = "sample_id") %>%
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
        # relocate(papillary:stroma, .after = histological_subtype)
    
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
    
    return(sample_table)
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

# export_ClinSummary <- function(samp_annot_63) {
#     samp_annot_63 %>%
#         mutate(Stage = str_sub(Stage, 1, 2),
#                Histologi)
# }


plt_apopt_receptors <- function(degtable_fgfr3, apop_hm) {
    apop_genes <- c("TNFRSF10A", "TNFRSF10B", "FASLG", "TNFRSF1A")
    tb <- degtable_fgfr3 %>%
        tibble() %>%
        filter(symbol %in% apop_genes) %>%
        mutate(mlogpval = -log10(P.Value)) %>%
        left_join(apop_hm$data, by = "symbol") %>%
        mutate(
            observations = as.character(observations),
            observations = ifelse(observations == "TNFalpha", "TNF-R1", observations) %>%
                   factor(levels = c("TNF-R1", "FAS", "TRAIL-R2", "TRAIL-R1")))
    
    ggplot(tb, aes(logFC, observations, fill = mlogpval)) +
        geom_vline(xintercept = 0, lty = "dashed") +
        geom_segment(aes(xend = 0, yend = observations), size = 0.5, lty = "dotted") +
        geom_point(pch = 21, color = "black", size = 3) +
        scale_fill_distiller(palette = "Reds", direction = 1, 
                             breaks = c(0,0.3,0.7,1,1.3,1.7),
                             labels = c(1,0.5,0.2,0.1,0.05,0.02),
                             limits = c(0,1.7)) +
        scale_x_continuous(limits = c(-0.1,0.6)) +
        labs(x = "FGFR3 mutated tumors vs WT - mRNA log2FC",
             y = "Apoptosis membrane receptors", 
             fill = "p-value") +
        theme_csg_scatter
    
    # Get cell line ------------
    ccle_meta <- depmap::metadata_21Q1()
    ccle_blad <- ccle_meta %>%
        filter(lineage == "urinary_tract")
    
    ccle_tpm <- depmap::TPM_21Q1() %>% filter(depmap_id %in% ccle_blad$depmap_id)
    ccle_mut <- depmap::mutationCalls_21Q1() %>% filter(depmap_id %in% ccle_blad$depmap_id)
    fgfr3_mut <- ccle_mut %>%
        filter(gene_name == "FGFR3", var_class != "Silent") %>%
        mutate(FGFR3_mutation = "M") %>%
        select(depmap_id, FGFR3_mutation)
    ccle_blad <- ccle_blad %>%
        select(depmap_id, cell_line_name) %>%
        left_join(fgfr3_mut) %>%
        mutate(FGFR3_mutation = ifelse(is.na(FGFR3_mutation), "WT", FGFR3_mutation))
    
    ccle_fgfr3_apop_exp <- ccle_tpm %>%
        filter(gene_name %in% apop_genes) %>%
        left_join(ccle_blad) %>%
        select(cell_line, gene_name, FGFR3_mutation, rna_expression)
    
    ggplot(ccle_fgfr3_apop_exp, aes(rna_expression, gene_name, fill = FGFR3_mutation)) +
        geom_boxplot() +
        stat_compare_means()
        
    

}

upg_check_morphology <- function(exp_table) {
    table2test <- exp_table %>%
        select(sample_id, class_uPG, stage, grade,  histological_subtype:endophytic) %>%
        mutate(
        squamous = case_when(
            str_detect(histological_subtype, "squamous") ~ "squamous",
            !is.na(histological_subtype) ~ "other",
            TRUE ~ NA_character_
        ),
        neuroendocrine = case_when(
            str_detect(histological_subtype, "neuroendocrine") ~ "neuroendocrine",
            !is.na(histological_subtype) ~ "other",
            TRUE ~ NA_character_)
        )
    
    res1 <- doAssociationTests(table2test, "sample_id", "class_uPG")
    
    
    tb_bin <- table2test %>%
        select(class_uPG, stage, grade, squamous, neuroendocrine) %>%
        fastDummies::dummy_cols(select_columns = "class_uPG") %>%
        mutate(Ta = ifelse(stage == "Ta", "Ta", "other"),
               T1 = ifelse(str_detect(stage, "T1"), "T1", "other"),
               T2 = ifelse(str_detect(stage, "T2"), "T2", "other"),
               `T3/T4` = ifelse(str_detect(stage, "T3|4"), "T3/T4", "other"))
    
    upgs <- paste0("class_uPG_", LETTERS[1:5])
    vars <- c("grade", "squamous", "neuroendocrine", "Ta", "T1", "T2", "T3/T4")
    
    combs <- expand_grid(var1 = upgs, var2 = vars) %>% as.data.frame() %>%
        rowwise() %>%
        mutate(pval = fisher.test(makeContingencyTable(tb_bin, var1, var2), alternative = "greater")$p.value, 
               padj = p.adjust(pval),
               mlog10p = -log10(padj))
    
    plt <- ggplot(combs, aes(var1, var2, fill = mlog10p, size = mlog10p)) +
        geom_point()
    
    return(list(plt = plt, assoc_tests = res1))
        
    
}

.get_path_annot <- function() {
    path_annot <- read_excel("extras/CIT_Pathology091017.xls", sheet = 2) %>%
        .[-c(1,2),] %>%
        mutate(sample_id = paste0("CIT.", str_remove_all(`nom initial ->`, "[A-Z]"))) %>%
        select(sample_id, papillary = `Papilles en surface`, cis = `CIS en surface`, pattern = Pattern,
               necrosis = `nécrose`, inflamation_LP = `inflammation LP`,
               inflammation_PN = `inflammation PN`, stroma)
    path_annot2 <-
        path_annot %>%
        mutate(
        across(c(papillary, cis, necrosis), ~ case_when(
            str_to_lower(.) == "non" ~ "no",
            str_to_lower(.) == "oui" ~ "yes",
            TRUE ~ .)),
        across(inflamation_LP:stroma, as.ordered),
        pattern = case_when(
            str_detect(pattern, "endophytique") ~ "endophytic",
            str_detect(pattern, "travées") ~ "trabecular",
            str_detect(pattern, "larges massifs") ~ "large and irregular islets",
            str_detect(pattern, "plages") ~ "discohesive sheets",
            TRUE ~ NA_character_
        )) %>%
        filter(sample_id %in% prot_data$sampAnnot$id)
}
