results_folders <- function() {
    folders_to_create <- c("results",
                           paste0("results/", 
                                  c(paste0("fig", 1:5),
                                    paste0("suppfig_", c("mofa", "ccp",
                                                         "bioinfoqc", "apop", "stroggilos",
                                                         "volcanos", "pca")),
                                    "mofa",
                                    "protein_clustering",
                                    "tables")))
    sapply(folders_to_create, dir.create, showWarnings = FALSE)
}