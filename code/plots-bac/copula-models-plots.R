# ==================================================================

#                       Plot the copula models

# ==================================================================





# ==================================================================
# 1. Import libraries and data
# ==================================================================
 

# Import libraries -------------------------------------------------
 
 library(ggplot2)
 library(viridis)
 library(ecoCopula)



# Setup folder paths -----------------------------------------------
 
 # Setup paths to import data
 path1 <- file.path(getwd(), "data")
 path2 <- file.path(getwd(), "outputs", "fits")



# Import data ------------------------------------------------------
 
 # Data files
 # Community matrices
 comm1 <- readRDS(
  file.path(
    path1,
    "env-data-clean",
    "env-bac-comm-bw.rds"
  )
 )
 
 comm2 <- readRDS(
  file.path(
    path1,
    "diet-data-clean",
    "diet-bac-comm-bw.rds"
  )
 )


 # Metadata
 meta1 <- readRDS(
  file.path(
    path1,
    "env-data-clean",
    "env-bac-metadata-bw.rds"
  )
 )

 meta2 <- readRDS(
  file.path(
    path1,
    "diet-data-clean",
    "diet-bac-metadata-bw.rds"
  )
 )


 # Taxa
 taxa1 <- readRDS(
  file.path(
    path1,
    "env-data-clean",
    "env-bac-taxa-bw.rds"
  )
 )
 
 taxa2 <- readRDS(
  file.path(
    path1,
    "diet-data-clean",
    "diet-bac-taxa-bw.rds"
  )
 )


 # Load the models
 lvm1 <- readRDS(
  file.path(
    path2,
    "env-bac-lvm-bw.rds"
  )
 )

 lvm2 <- readRDS(
  file.path(
    path2,
    "diet-bac-lvm2-bw.rds"
  )
 )

# ==================================================================
# ==================================================================






# ==================================================================
# 2. Prepare a custom theme
# ==================================================================
 
 custom_theme <- theme(
     axis.text.y   = element_text(size = 13,
                                  color = "black"),
     axis.text.x   = element_text(size = 13,
                                  color = "black"),
     axis.title.y  = element_text(size = 14),
     axis.title.x  = element_text(size = 14),
     axis.ticks.length = unit(.15, "cm"),
     panel.background = element_blank(),
     panel.grid.major = element_blank(), 
     panel.grid.minor = element_blank(),
     axis.line = element_line(colour = "black"),
     panel.border = element_rect(colour = "black",
                                 fill = NA,
                                 linewidth = 0.5),
     legend.key = element_rect(fill = "transparent"),
     legend.title = element_text(size = 14),
     legend.text = element_text(size = 14),
     legend.position = "top"
 )

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Plot the model for env spiders
# ==================================================================


# Prepare plot data ------------------------------------------------

 alpha1 <- 0.7 # ou 2 comme dans l'exemple?
 site_res1 <- data.frame(lvm1$scores, meta1)
 site_res1$sample_env <- as.factor(site_res1$sample_env)
 sp_res1 <- data.frame(lvm1$loadings,
                      ASV = colnames(comm1))
 fam1 <- data.frame(ASV = rownames(taxa1),
                   family = taxa1[,5])
 sp_res1 <- merge(sp_res1, fam1, by = "ASV")



# Calculate ellipses for each environment --------------------------
 
 # Mean Factor axis for each environment
 lvm_mean1 <- aggregate(
  site_res1[,1:2],
  list(site_res1$sample_env),
  mean
 )
 
# Some personal notes : now I understand the function
# This should give the same results
# SO it is an standard error of the weighted
# centroid with 95% confidence interva
#plot(lvm1$scores)
#ord<- vegan::ordiellipse(lvm1$scores, meta1$sample_env, display = "sites", 
#                   kind = "se", conf = 0.95, label = T)
#
#ord<-ordiellipse(sol, MyMeta$amt, display = "sites", 
#                   kind = "se", conf = 0.95, label = T)

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
 df_ell1 <- data.frame()
  for (g in levels(site_res1$sample_env)) {
    df_ell1 <- rbind(df_ell1,
                    cbind(as.data.frame(with(
                      site_res1[site_res1$sample_env == g, ],
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
        data = site_res1,
        aes(x = Factor1,
            y = Factor2,
            fill = sample_env,
            shape = sample_env),
        size = 3
     ) +
     geom_path(
        data = df_ell1,
        aes(x = Factor1,
            y = Factor2,
            color = sample_env),
        linewidth = 1,
        linetype = "dashed",
        show.legend = FALSE
     ) +
     scale_x_continuous(
        breaks = seq(-3, 3, 1),
        limits = c(-3.4, 3.4)
     ) +
     scale_y_continuous(
        breaks = seq(-3, 3, 1),
        limits = c(-3.4, 3.4)
     ) +
     scale_shape_manual(values = c(21, 24)) +
     scale_fill_manual(values = c("#E69F00",
                                  "#666666")) +
     scale_color_manual(values = c("#E69F00",
                                   "#666666")) +
     labs(fill = "Environment :",
          shape = "Environment :") +
     xlab("\nLatent variable 1") +
     ylab("Latent variable 2\n") +
     custom_theme

# ==================================================================
# ==================================================================





# ==================================================================
# 4. Plot the model for diet spiders
# ==================================================================

# Prepare plot data ------------------------------------------------

 alpha2 <- 0.7 # ou 2 comme dans l'exemple?
 site_res2 <- data.frame(lvm2$scores, meta2)
 site_res2$diet_treatment <- as.factor(site_res2$diet_treatment)
 sp_res2 <- data.frame(lvm2$loadings,
                      ASV = colnames(comm2))
 fam2 <- data.frame(ASV = rownames(taxa2),
                   family = taxa2[,5])
 sp_res2 <- merge(sp_res2, fam2, by = "ASV")



# Calculate ellipses for each environment --------------------------
 
 # Mean Factor axis for each environment
 lvm_mean2 <- aggregate(
  site_res2[,1:2],
  list(site_res2$diet_treatment),
  mean
 )
 

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
 df_ell2 <- data.frame()
  for (g in levels(site_res2$diet_treatment)) {
    df_ell2 <- rbind(df_ell2,
                    cbind(as.data.frame(with(
                      site_res2[site_res2$diet_treatment == g, ],
                      veganCovEllipse(cov.wt(
                        cbind(Factor1, Factor2),
                        wt = rep(1 / length(Factor1), length(Factor2))
                      )$cov,
                      center = c(mean(Factor1), mean(Factor2)))
                    )),
                    diet_treatment = g))
  }



# Compute the biplot and export ------------------------------------

 # Compute the plot
 plot2 <- ggplot() +
     geom_point(
        data = site_res2,
        aes(x = Factor1,
            y = Factor2,
            fill = diet_treatment,
            shape = diet_treatment),
        size = 3
     ) +
     geom_path(
        data = df_ell2,
        aes(x = Factor1, 
            y = Factor2,
            color = diet_treatment),
        linewidth = 1,
        linetype = "dashed",
        show.legend = FALSE
     ) +
     scale_x_continuous(
        breaks = seq(-3, 3, 1),
        limits = c(-3, 3)
     ) +
     scale_y_continuous(
        breaks = seq(-3, 3, 1),
        limits = c(-3, 3)
     ) +
     scale_shape_manual(values = c(23, 22, 25)) +
     scale_fill_viridis(discrete = TRUE,
                        option = "viridis") +
     scale_color_viridis(discrete = TRUE,
                         option = "viridis") +
     labs(fill = "Diet :",
          shape = "Diet :") +
     xlab("\nLatent variable 1") +
     ylab("Latent variable 2\n") +
     custom_theme

# ==================================================================
# ==================================================================





# ==================================================================
# 5. Export the plots
# ==================================================================
 
 # Folder path
 path <- file.path(
  getwd(),
  "outputs",
  "plots"
 )

 # Environment spiders
 ggsave(
  plot1,
  width = 12,
  height = 12,
  dpi = 300,
  units = c("cm"),
  file = file.path(
    path,
    "env-bac-lvm-bw.png"
  )
 )

 # Diet spiders
 ggsave(
  plot2,
  width = 12,
  height = 12,
  dpi = 300,
  units = c("cm"),
  file = file.path(
    path,
    "diet-bac-lvm-bw.png"
  )
 )

# ==================================================================
# ==================================================================