# ==================================================================

#                       Plot the copula models

# ==================================================================





# ==================================================================
# 1. Import libraries and data
# ==================================================================
 

# Import libraries -------------------------------------------------
 
 library(ggplot2)
 library(corrplot)
 library(ecoCopula)



# Import data ------------------------------------------------------
 
 # Folder
 path <- file.path(getwd(), "env-folder", 
                   "env-data", 'env-data-clean')
 

 # Data files
 # Community matrices
 comm_bw <- readRDS(file.path(path, "env-bac-comm-bw.rds"))
 comm_w <- readRDS(file.path(path, "env-bac-comm-w.rds"))
 
 # Metadata
 meta_bw <- readRDS(file.path(path, "env-bac-metadata-bw.rds"))
 meta_w <- readRDS(file.path(path, "env-bac-metadata-w.rds"))

 # Taxa
 taxa_bw <- readRDS(file.path(path, "env-bac-taxa-bw.rds"))
 taxa_w <- readRDS(file.path(path, "env-bac-taxa-w.rds"))


 # Delete rare communities for web samples
 dim(comm_w)
 comm_w <- comm_w[, colSums(comm_w) > 1]
 dim(comm_w)
 
 # Load the models


# ==================================================================
# ==================================================================





# ==================================================================
# 2. Plot the model for spiders
# ==================================================================


# Biplot spiders model ---------------------------------------------
 
 plot(bw_lv, biplot = TRUE)
 
 # Prepare plotting parameters
 alpha <- 2
 site_res <- data.frame(bw_lv$scores, meta_bw)
 sp_res <- data.frame(bw_lv$loadings,
                      ASV = colnames(comm_bw))
 fam <- data.frame(ASV = rownames(taxa_bw),
                   family = taxa_bw[,6])
 sp_res <- merge(sp_res, fam, by = "ASV")

# Draw the plot
biplot1 <- ggplot() +
    geom_point(data = site_res,
       aes(x = Factor1,
           y = Factor2,
           color = sample_env),
       size = 3) +
    geom_text(data = sp_res,
       aes(x = Factor1*alpha,
           y = Factor2*alpha,
           label = family),
       size = 2.5) +
    scale_x_continuous(breaks = seq(-2, 2, 0.5),
                       limits = c(-2.5, 2.5)) +
    scale_y_continuous(breaks = seq(-2, 2, 0.5),
                       limits = c(-2.5, 2.5)) +
    scale_fill_manual(values = c("#E69F00",
                                  "#666666")) +
    scale_color_manual(name = "Environment :",
                       values = c("#E69F00",
                                  "#666666")) +
    xlab("\nFactor 1") + ylab("Factor 2\n") +                              
    theme_bw() + 
    theme(legend.position = "top",
          panel.grid = element_blank())



# Correlation matrix spiders model ---------------------------------

 # Prepare color gradient
 COL2(diverging = c("RdBu", "BrBG", "PiYG",
                    "PRGn", "PuOr", "RdYlBu"),
      n = 200)
 
 # Extract correlation matrix
 cormat <- bw_lv$sigma
 
 # Diagonal matrix with hclust no squares
 corrplot(cormat,
          type = "lower",
          diag = FALSE,
          method = "circle",
          order = "hclust",
          tl.srt = 45,
          tl.col = "black",
          tl.cex = 0.7,
          col = COL2("RdBu", 10),
          cl.ratio = 0.1)

 # Matrix with hclust and squares
 #corrplot(cormat,
 #         #type = "lower",
 #         #diag = FALSE,
 #         #diag = TRUE,
 #         method = "circle",
 #         order = "hclust",
 #         addrect = 3,
 #         #tl.srt = 45,
 #         tl.col = "black",
 #         tl.cex = 0.7,
 #         col = COL2("RdBu", 10),
 #         cl.ratio = 0.1)

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Plot the model for webs
# ==================================================================



# ==================================================================
# ==================================================================