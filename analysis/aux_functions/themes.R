theme_csg_hist <- theme_minimal() +
    theme(panel.grid = element_blank(), 
          panel.grid.major.y = element_line(color = "grey30", linetype = "dotted"),
          axis.line.x = element_line(),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 8),
          axis.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8, face = "bold"))
theme_csg_scatter <- theme_light() +
    theme(panel.border = element_rect(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(size = 8, color = "black"), 
          axis.title = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "transparent", color = "black"),
          strip.text = element_text(color = "black"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8, face = "bold"))
theme_csg_sparse <- theme_minimal() +
    theme(panel.grid = element_blank(), 
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(size = 8, color = "black"),
          axis.title = element_text(size = 8, face="bold"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8, face = "bold"))
theme_mofa <- theme_csg_scatter +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))