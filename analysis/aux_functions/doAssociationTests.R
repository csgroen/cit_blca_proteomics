library(tidyselect)
doAssociationTests <- function(data, id_var = "id", test_var = "class") {
    data <- as.data.frame(data)
    #-- Get variable types
    to_test <- setdiff(colnames(data), c(id_var, test_var))
    vtype <- sapply(to_test, function(col) {
        if(all(data[,col] %in% c("0","1", NA))) return("binary")
        else if(is.numeric(data[,col])) return("numeric")
        else return("categorical") 
    })
    #-- Get test_varn and id_varn
    # test_var <- colnames(data)[test_var]
    
    #-- Get pvalues
    all_pvals <- lapply(1:length(vtype), function(i) {
        vt <- vtype[i]
        var <- names(vtype)[i]
        
        if(vt %in% c("binary", "categorical")) {
            #-- Make FET table
            fet_tb <- data %>%
                dplyr::select(all_of(test_var), all_of(var)) %>%
                group_by(!! sym(test_var), !! sym(var)) %>%
                dplyr::summarize(n=n()) %>%
                pivot_wider(id_cols = test_var, names_from = var, values_from = n) %>%
                column_to_rownames(test_var)
            fet_tb[is.na(fet_tb)] <- 0
            
            gen_pval <- fisher.test(fet_tb, simulate.p.value = TRUE)$p.value
            
            #-- Group v non-group
            class_pvals <- sapply(1:nrow(fet_tb), function(j) {
                gr <- fet_tb[j,]
                ngr <- colSums(fet_tb[-j,])
                if(vt == "binary") fisher.test(rbind(gr,ngr))$p.value
                else fisher.test(rbind(gr,ngr), simulate.p.value = TRUE)$p.value
            })
            names(class_pvals) <- rownames(fet_tb)
            
            return(c(General = gen_pval, class_pvals))
            
        } else {
            #-- Prepare KW vectors
            if (all(is.na(data[,var]))) {
                gen_pv <- NA
            } else {
                gen_pv <- kruskal.test(data[,var], data[,test_var])$p.value
                dunn_res <- DescTools::DunnTest(data[,var], as.factor(data[,test_var]))[[1]] %>%
                    as.data.frame() %>%
                    rownames_to_column("comparison") %>%
                    arrange(comparison)
                gen_pv <- list(kw = gen_pv, dunn = dunn_res)
            }
            
            return(gen_pv)
        }
    })
    names(all_pvals) <- to_test
    stat_tests <- plyr::mapvalues(vtype, from = c("binary", "categorical", "numeric"),
                                  to = c("Fisher's test", "Fisher's test", "Kruskal-Wallis"))
    
    
    return(list(pvalues = all_pvals, test_type = stat_tests))
    
}