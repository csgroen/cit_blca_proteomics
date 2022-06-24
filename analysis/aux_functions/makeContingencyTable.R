makeContingencyTable <- function(data, var1, var2) {
    #-- Count and pivot
    ctg_table <- data %>%
        group_by(!! sym(var1), !! sym(var2)) %>%
        dplyr::count() %>%
        pivot_wider(names_from = var1,
                    values_from = n) %>%
        ungroup() %>%
        mutate(across(matches(var2), ~ replace_na(as.character(.), "NA"))) %>%
        as.data.frame() %>%
        column_to_rownames(var2) %>%
        as.matrix()
    names(dimnames(ctg_table)) <- c(var2, var1)
    #-- Correct NAs to 0
    ctg_table[is.na(ctg_table)] <- 0
    
    return(ctg_table)
    
}