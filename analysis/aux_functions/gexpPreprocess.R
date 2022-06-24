#' Pre-process a gene expression matrix to have different gene annotation
#' 
#' @param gexp A gene expression matrix, with samples in columns and genes in rows
#' @param gene_annotation A data.frame of gene annotation, refering to the rows of
#' `gexp`. It must cotidyselectntain at least two columns: one for the `keep_annot` and one
#' for the `og_annot`
#' @param keep_annot The name of the column in `gene_annotation` that contains
#' the gene names that you would like to put as rownames of your `gexp`, i.e. the
#' kept annotation
#' @param og_annot The name of the column in `gene_annotation` that contains the
#' gene names currently being used by `gexp`, i.e. the original annotation.
#' @param keep_stat a character with a function name for the statistical measure 
#' used for keeping a gene when `keep_annot` doesn't map uniquely to `og_annot`.
#' 
#' @example
# set.seed(0)
# gexp <- matrix(rnorm(1000000), nrow = 10000,
#                dimnames = list(paste0("ID", 1:10000), paste0("Sample_", 1:100)))
# #-- Make fake annotation
# gene_annotation <- data.frame(Original_ID = paste0("ID", 1:10000),
#                               New_ID = paste0("NewID", sample(1:9000, 10000, replace = TRUE)),
#                               stringsAsFactors = FALSE)
# new_data <- gexpPreprocess(gexp, gene_annotation, keep_annot = "New_ID", og_annot = "Original_ID")
# head(new_data$gexp)
# head(new_data$gene_annotation)

gexpPreprocess <- function(gexp, gene_annotation, keep_annot = "Symbol",
                           og_annot = "ProbeID", keep_stat = "cv") {
    #-- Get important annotation
    gannot <- gene_annotation %>%
        dplyr::select(!! og_annot, !! keep_annot)
    samples <- colnames(gexp)
    
    #-- Calculate the keep_stat
    stat <- apply(gexp, 1, eval(parse(text=keep_stat)))
    gexp_int <- cbind(gexp, stat)
    
    #-- Filter IDs by the new annotation based on max stat per group
    gexp_int <- gexp_int %>%
        as.data.frame() %>%
        rownames_to_column(og_annot) %>%
        left_join(gannot, by = og_annot) %>%
        dplyr::filter(!is.na(!!sym(keep_annot))) %>%
        group_by(!!sym(keep_annot)) %>%
        top_n(1, stat) %>%
        dplyr::slice(1) %>%
        ungroup()
    
    #-- Filter annotation
    gannot <- gene_annotation %>%
        dplyr::filter(eval(parse(text=og_annot)) %in% pull(gexp_int, !!og_annot))
    
    #-- Re-assemble gexp
    gexp_res <- gexp_int %>%
        dplyr::select(-!!og_annot, -stat) %>%
        as.data.frame() %>%
        column_to_rownames(keep_annot) %>%
        as.matrix()
    
    return(list(gexp = gexp_res, gene_annotation = gannot))
}

# Coefficient of variation function
cv <- function (v) { sd(v, na.rm = TRUE)/mean(v, na.rm = TRUE) }
