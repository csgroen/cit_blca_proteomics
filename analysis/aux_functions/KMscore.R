KMscore <- function(annot, sig_name, time, event, qbks = 4, qnames = NULL, 
                    palette = NULL, xlim = c(0,24), 
                    save = TRUE, fpath = "Results/", bk_time_by = 6,
                    width = 4, height = 6, rev_pal = TRUE) {
    
    #-- Make or check qnames
    if(is.null(qnames)) {
        if(length(sig_name) > 1) { 
            stop("Must provide break names for multiple variable evaluations.")
        }
        qnames <- paste0("Q", 1:qbks)
    } else {
        # if(length(qnames) != qbks*2) {
        #     stop("Number of `qnames` must be the same as the `qbks`.")
        # }
    }
    
    #-- Make palette if none was provided
    if(is.null(palette)) {
        if(qbks == 2) {
            palette <- rev(RColorBrewer::brewer.pal(3, "RdBu"))[c(1,3)]
        } else {
            palette <- rev(RColorBrewer::brewer.pal(qbks*length(sig_name), "RdBu"))
            
        }
    }
    
    #-- Select variables from table
    sel_annot <- select(annot, id, !! time, !! event, !!! sig_name) 
    
    if(length(sig_name) == 1) {
        #-- Rename and discretize
        annot4plot <- sel_annot %>%
            rename(sig = !! sig_name, time = !! time, event = !! event) %>%
            mutate(sig_qt = arules::discretize(sig, breaks = qbks, labels = qnames))
        
        #-- Fit survival model
        sf <- survfit(Surv(time, event) ~ sig_qt, data = annot4plot) 
        names(sf$strata) <- qnames
        
        #-- Make survival plot
        sig_km <- ggsurvplot(sf, data = annot4plot, pval = T, 
                             palette = palette, risk.table = TRUE, xlim = xlim,
                             break.time.by = bk_time_by, fontsize = 4, tables.height = 0.3)
        
        if(save) {
            cowplot::plot_grid(sig_km$plot, sig_km$table, nrow = 2, rel_heights = c(0.7,0.3)) +
                ggsave(paste0(fpath, "KMplot_", sig_name, ".pdf"), width = width, height = height)
        }
        
        return(list(survfit = sf, plot = sig_km))
        
    } else {
        #-- Rename and discretize
        annot4plot <- sel_annot %>%
            rename(sig1 = sig_name[1], sig2 = sig_name[2], time = !! time, event = !! event) %>%
            mutate(sig1_qt = arules::discretize(sig1, breaks = qbks, labels = qnames[1:qbks]),
                   sig2_qt = arules::discretize(sig2, breaks = qbks, labels = qnames[(qbks+1):(qbks*2)]))
        
        #-- Fit survival model
        sf <- surv_fit(Surv(time, event) ~ sig1_qt + sig2_qt, data = annot4plot)
        
        ptxts <- c()
        for(qname in qnames[1:2]) {
            sub_annot <- filter(annot4plot, sig1_qt == qname)
            surv_res <- surv_fit(Surv(time, event) ~ sig2_qt, data =  sub_annot)
            ptxts[qname] <- surv_pvalue(surv_res)$pval.txt
        }
        ptxts <- paste0(names(ptxts), ", ", ptxts)
        
        names(sf$strata) <- names(sf$strata) %>% 
            str_remove("sig1_qt=") %>% str_remove("sig2_qt=") %>% str_replace(",", " |")
        
        sigs_km <- ggsurvplot(sf, data = annot4plot, pval = FALSE, 
                              palette = palette, risk.table = TRUE, xlim = xlim,
                              break.time.by = 6, fontsize = 4, tables.height = 0.3)
        
        sigs_km$plot <- sigs_km$plot +
            annotate("text", x = 0, y = 0.25, label = ptxts[1], hjust = 0, size = 4.5) +
            annotate("text", x = 0, y = 0.18, label = ptxts[2], hjust = 0, size = 4.5)
        
        if(save) {
            cowplot::plot_grid(sigs_km$plot, sigs_km$table, nrow = 2, rel_heights = c(0.7,0.3), align = "v") +
                ggsave(paste0(fpath, "KMplot_", sig_name[1], "_", sig_name[2], ".pdf"), width = width, height = height)
        }
        return(list(survfit = sf, plot = sigs_km))
        
    }
}