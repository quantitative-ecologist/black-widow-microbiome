# ==================================================================

#                   Assemble and export figure 1

# ==================================================================





# ==================================================================
# 1. Source plot files
# ==================================================================
 
 # To combine as 1 figure
 library(ggpubr)
 
 # Script paths
 path1 <- file.path(getwd(), "code", "plots-bac")
 path2 <- file.path(getwd(), "code", "plots-euk")
 
 # These files have to be run in this order
 #1
 source(file.path(path1, "glm-plots-bac.R"))
 #2
 source(file.path(path2, "glm-plots-euk.R"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Combine the plots as one figure
# ==================================================================

 # Combine into one figure
 fig1 <- ggarrange(
     plot1, NULL, plot2,
     plot4, NULL, plot5,
     nrow = 2, ncol = 3,
     widths = c(1, 0.1, 1, 1, 0.1, 1),
     common.legend = TRUE,
     legend = "top",
     labels = c("(A)", "", "(B)", "(C)", "", "(D)"),
     hjust = c(-0.1, -0.1),
     vjust = c(1.5, 1.5)
 )

 # Add labels to the figure
 fig1 <- annotate_figure(
   fig1,
   left = text_grob(
     "Bacteria", hjust = -1.8, vjust = -0.5,
     face = "bold", size = "12", rot = 90)
 )

 fig1 <- annotate_figure(
   fig1,
   left = text_grob(
     "Eukaryotes", hjust = 2.4, vjust = 1.1,
     face = "bold", size = "12", rot = 90)
 )

 # Export the figure
 ggexport(
    fig1,
    filename = file.path(path, "figures", "figure1.png"),
    width = 2000,
    height = 1800,
    res = 300
 )
 
# ==================================================================
# ==================================================================