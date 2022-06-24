source("aux_functions/makeContingencyTable.R")
pathwayORA <- function(genes, pathway_table, id = "Symbol", 
                       pcutoff = 0.05) {
    #-- Get pathway members and names
    path_ls <- split(pull(pathway_table, !! id), pathway_table$Path_ID)
    path_names <- pathway_table %>% select(Path_ID, Path_Name) %>% distinct() %>% pull(Path_Name, Path_ID)
    
    #-- Universe size
    univ <- pathway_table %>% pull(!! id) %>% unique() %>% length()
    
    enrich_res <- lapply(names(path_ls), function(id) {
        path <- path_ls[[id]]
        ginpath <- sum(genes %in% path)
        gopath <- length(genes) - ginpath
        opath <- length(path) - ginpath
        rest <- univ - ginpath - gopath - opath
        ctg <- matrix(c(ginpath, opath, gopath, rest), nrow = 2)
        pfish <- fisher.test(ctg, alternative = "greater")$p.value
        
        data.frame(Path_ID = id, Path_Name = path_names[id], pval = pfish, 
                   GeneRatio = paste0(ginpath, '/', length(genes)),
                   BgRatio = paste0(opath, '/', univ - length(path)),
                   Genes = paste(genes[genes %in% path], collapse = ", "))
    }) %>% bind_rows()
    
    enrich_res %>%
        mutate(padj = p.adjust(pval, "BH")) %>%
        arrange(padj) %>%
        filter(padj < pcutoff)
    
} 