################################################################################
## Pathway enrichment functions
################################################################################
library(msigdbr)
library(fgsea)
get_pathway_list <- function(category, subcategory, name_fix = "-") {
    path_table <- msigdbr(category = category, subcategory = subcategory) %>%
        mutate(pathway_name = str_remove(gs_name, name_fix))
    
    pathlist <- split(path_table$gene_symbol, path_table$pathway_name)
    path_annot <- path_table %>%
        select(pathway_id = gs_exact_source, pathway_name) %>%
        distinct()
    return(list(list = pathlist, annot = path_annot))
}

run_fgsea <- function(lfc_vec, pathlist, path_annot) {
    fgsea(pathlist, lfc_vec, eps = 0) %>%
        as_tibble() %>%
        arrange(padj) %>%
        mutate(leadingEdge = map_chr(leadingEdge, ~ paste(., collapse = ", "))) %>%
        left_join(path_annot, by = c("pathway" = "pathway_name")) %>%
        relocate(pathway_id, pathway_name = pathway)
}

export_pathTable <- function(enrich_results, file) {
    overlap <- lapply(clines, function(cl) {
        enrich_results[[cl]] %>%
            filter(padj < 0.05) %>%
            select(pathway, NES, padj) %>%
            rename_with(~ paste0(., "_", cl), NES:padj)
    }) %>% reduce(inner_join, by = "pathway")
    
    enrich_results[["Intersection_significant"]] <- overlap
    
    write.xlsx(enrich_results, file)
}
