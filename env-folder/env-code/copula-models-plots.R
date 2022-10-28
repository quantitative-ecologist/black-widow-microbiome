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
 library(RColorBrewer)

 

# Setup folder paths -----------------------------------------------
 
 # Setup paths to import data
 path1 <- file.path(getwd(), "env-folder", 
                   "env-data", 'env-data-clean')
 path2 <- file.path(getwd(), "env-folder", 
                   "env-outputs")
 
 # Setup a path to export the plots
 path3 <- file.path(getwd(), "env-folder", "env-outputs")



# Import data ------------------------------------------------------
 
 # Data files
 # Community matrices
 comm_bw <- readRDS(file.path(path1, "env-bac-comm-bw.rds"))
 comm_w <- readRDS(file.path(path1, "env-bac-comm-w.rds"))
 
 # Metadata
 meta_bw <- readRDS(file.path(path1, "env-bac-metadata-bw.rds"))
 meta_w <- readRDS(file.path(path1, "env-bac-metadata-w.rds"))

 # Taxa
 taxa_bw <- readRDS(file.path(path1, "env-bac-taxa-bw.rds"))
 taxa_w <- readRDS(file.path(path1, "env-bac-taxa-w.rds"))
 
 # Load the models
 lvm_bw <- readRDS(file.path(path2, "lvm_bw.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Plot the model for spiders
# ==================================================================


# Prepare some custom plotting options -----------------------------
 
 custom_theme <- theme(
     axis.text.y   = element_text(size = 13, color = "black"),
     axis.text.x   = element_text(size = 13, color = "black"),
     axis.title.y  = element_text(size = 14),
     axis.title.x  = element_text(size = 14),
     axis.ticks.length = unit(.15, "cm"),
     panel.background = element_blank(),
     panel.grid.major = element_blank(), 
     panel.grid.minor = element_blank(),
     axis.line = element_line(colour = "black"),
     panel.border = element_rect(colour = "black",
                                 fill = NA,
                                 size = 0.5),
     legend.key = element_rect(fill = "transparent"),
     legend.title = element_text(size = 14),
     legend.text = element_text(size = 14),
     legend.position = "top"
 )



# Prepare plot data ------------------------------------------------

 alpha <- 2
 site_res <- data.frame(lvm_bw$scores, meta_bw)
 site_res$sample_env <- as.factor(site_res$sample_env)
 sp_res <- data.frame(lvm_bw$loadings,
                      ASV = colnames(comm_bw))
 fam <- data.frame(ASV = rownames(taxa_bw),
                   family = taxa_bw[,5])
 sp_res <- merge(sp_res, fam, by = "ASV")



# Calculate ellipses for each environment --------------------------
 
 # Mean Factor axis for each environment
 lvm_mean <- aggregate(site_res[,1:2], list(site_res$sample_env), mean)
 

 # Compute the function to calculate ellipses
 veganCovEllipse <-
  function (cov,
            center = c(0, 0),
            scale = 1,
            npoints = 100)
  {
    theta <- (0:npoints) * 2 * pi / npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }


 # Create a data frame of the ellipses results for ggplot
 df_ell <- data.frame()
  for (g in levels(site_res$sample_env)) {
    df_ell <- rbind(df_ell,
                    cbind(as.data.frame(with(
                      site_res[site_res$sample_env == g, ],
                      veganCovEllipse(cov.wt(
                        cbind(Factor1, Factor2),
                        wt = rep(1 / length(Factor1), length(Factor2))
                      )$cov,
                      center = c(mean(Factor1), mean(Factor2)))
                    )),
                    sample_env = g))
  }



# Compute the biplot and export ------------------------------------

 # Compute the plot
 plot1 <- ggplot() +
     geom_point(
        data = site_res,
        aes(x = Factor1,
            y = Factor2,
            color = sample_env,
            shape = sample_env),
        size = 3
     ) +
     geom_path(
        data = df_ell,
        aes(x = Factor1, y = Factor2, color = sample_env),
        size = 1,
        linetype = "dashed",
        show.legend = FALSE
     ) +
     scale_x_continuous(
        breaks = seq(-2, 2, 1),
        limits = c(-2.5, 2.5)
     ) +
     scale_y_continuous(
        breaks = seq(-2, 2, 1),
        limits = c(-2.5, 2.5)
     ) +
     scale_fill_manual(values = c("#E69F00",
                                  "#666666")) +
     scale_color_manual(values = c("#E69F00",
                                   "#666666")) +
     labs(color = "Environment :",
          shape = "Environment :") +
     xlab("\nFactor 1") + ylab("Factor 2\n") +
     custom_theme
 
 # Save the plot
 ggsave(plot1, file = file.path(path3, "lvm-biplot-bw.png"))



# Correlation matrix spiders model ---------------------------------

 # Prepare color gradient
 COL2(diverging = c("RdBu", "BrBG", "PiYG",
                    "PRGn", "PuOr", "RdYlBu"),
      n = 200)
 
 # Extract correlation matrix
 cormat <- lvm_bw$sigma
 
 # Produce and export the correlation plot
 png(file.path(path3, "corrplot_bw.png"),
     res = 300,
     width = 1500,
     height = 1500)
 
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

 dev.off()

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Plot the model for webs
# ==================================================================



# ==================================================================
# ==================================================================