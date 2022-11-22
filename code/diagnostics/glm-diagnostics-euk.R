# ==================================================================

#                   Plot the ASV richness plots

# ==================================================================




# ==================================================================
# 1. Import data, libraries, and model fits
# ==================================================================

# Libraries --------------------------------------------------------
 
 library(brms)
 library(ggplot2)
 library(ggpubr)



# Data -------------------------------------------------------------
 
 # Folder paths
 path1 <- file.path(getwd(), "data", "data-clean-env")

 # Metadata
 meta_bw <- readRDS(file.path(path1, "metadata-env-euk-bw.rds"))
 meta_w <- readRDS(file.path(path1, "metadata-env-euk-w.rds"))



# Model fits -------------------------------------------------------

 # Folder path
 path2 <- file.path(getwd(), "outputs")
 
 # Models
 fit1 <- readRDS(file.path(path2, "fits-euk", "glm-env-euk-bw.rds"))
 fit2 <- readRDS(file.path(path2, "fits-euk", "glm-env-euk-w.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Plot the model diagnostics
# ==================================================================
 

# Trace plots and parameter distributions --------------------------

 # Inspect if the chains have converged
 plot(fit1)
 plot(fit2)



# Posterior predictive checks --------------------------------------
 
 # Check distributions
 pp1 <- pp_check(fit1, ndraws = 500)
 pp2 <- pp_check(fit2, ndraws = 500)
  
 # Effects plot (mean variance plot)
 stat1 <- pp_check(fit1, type = "stat_2d")
 stat2 <- pp_check(fit2, type = "stat_2d")
 
 
 # Predicted means
 mean1 <- pp_check(fit1, type = "stat", stat = "mean")
 mean2 <- pp_check(fit2, type = "stat", stat = "mean")
 
 
 # Error scatter
 e_scat1 <- pp_check(fit1, type = "error_scatter_avg")
 e_scat2 <- pp_check(fit2, type = "error_scatter_avg")
 

# Model assumptions ------------------------------------------------
 
 # Homogeneity of variance
 boxplot(residuals(fit1) ~ meta_bw$sample_env)
 boxplot(residuals(fit2) ~ meta_w$sample_env)
 
 # Normal distribution of residuals
 hist(residuals(fit1))
 hist(residuals(fit2))

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Save the plots into figures
# ==================================================================


# Arrange figures --------------------------------------------------
 
 # Model 1
 stat_fig1 <- ggarrange(pp1, stat1,
                        mean1, e_scat1,
                        ncol = 2, nrow = 2)
 
 # Model 2
 stat_fig2 <- ggarrange(pp2, stat2,
                        mean2, e_scat2,
                        ncol = 2, nrow = 2)
 


 # Export the figures ----------------------------------------------

 # Folder path
 path3 <- file.path(getwd(), "outputs", "diagnostics")

 ggexport(stat_fig1,
          filename = file.path(path3, "glm-diag-env-euk-bw.png"),
          width = 3000, height = 2500, res = 300)
 
 ggexport(stat_fig2,
          filename = file.path(path3, "glm-diag-env-euk-w.png"),
          width = 3000, height = 2500, res = 300)

# ==================================================================
# ==================================================================