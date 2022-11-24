# ==================================================================

#                   Assemble and export figure 2

# ==================================================================





# ==================================================================
# 1. Source plot files
# ==================================================================
 
 # To combine as 1 figure
 library(ggpubr)
 
 # Script paths
 fold_bac <- file.path(getwd(), "code", "figures", "plots-bac")
 fold_euk <- file.path(getwd(), "code", "figures", "plots-euk")
 
 # These files have to be run in this order
 #1
 source(file.path(fold_bac, "ordination-plots-env-bac.R"))
 #2
 source(file.path(fold_euk, "ordination-plots-env-euk.R"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Combine the plots as one figure
# ==================================================================

 # Combine into one figure
 fig2 <- ggarrange(
     plot1, NULL, plot2,
     plot3, NULL, plot4,
     nrow = 2, ncol = 3,
     widths = c(1, 0.1, 1, 1, 0.1, 1),
     common.legend = TRUE,
     legend = "top",
     labels = c("(A)", "", "(B)", "(C)", "", "(D)"),
     hjust = c(-0.1, -0.1),
     vjust = c(1.5, 1.5)
 )

 # Add labels to the figure
 fig2 <- annotate_figure(
   fig2,
   left = text_grob(
     "Bacteria", hjust = -2.2, vjust = -0.5,
     face = "bold", size = "12", rot = 90)
 )

 fig2 <- annotate_figure(
   fig2,
   left = text_grob(
     "Eukaryotes", hjust = 2.4, vjust = 1.1,
     face = "bold", size = "12", rot = 90)
 )

 # Export the figure
 path <- file.path(getwd(), "outputs")
 ggexport(
    fig2,
    filename = file.path(path, "figures", "figure2.png"),
    width = 2500,
    height = 2200,
    res = 300
 )

# ==================================================================
# ==================================================================