# ==================================================================

#                       Plot the copula models

# ==================================================================





# ==================================================================
# 1. Import libraries and data
# ==================================================================
 

# Import libraries -------------------------------------------------
 
 library(ggplot2)
 library(ggpubr)
 library(ecoCopula)



# Setup folder paths -----------------------------------------------
 
 # Setup paths to import data
 path1 <- file.path(getwd(), "data")
 path2 <- file.path(getwd(), "outputs", "fits-euk")



# Import data ------------------------------------------------------
 
 # Community matrices
 comm1 <- readRDS(
  file.path(
    path1,
    "data-clean-env",
    "comm-env-euk-bw.rds"
  )
 )
 
 comm2 <- readRDS(
  file.path(
    path1,
    "data-clean-env",
    "comm-env-euk-w.rds"
  )
 )


 # Metadata
 meta1 <- readRDS(
  file.path(
    path1,
    "data-clean-env",
    "metadata-env-euk-bw.rds"
  )
 )

 meta2 <- readRDS(
  file.path(
    path1,
    "data-clean-env",
    "metadata-env-euk-w.rds"
  )
 )


 # Taxa
 taxa1 <- readRDS(
  file.path(
    path1,
    "data-clean-env",
    "taxa-env-euk-bw.rds"
  )
 )
 
 taxa2 <- readRDS(
  file.path(
    path1,
    "data-clean-env",
    "taxa-env-euk-w.rds"
  )
 )


 # Load the models
 lvm1 <- readRDS(
  file.path(
    path2,
    "lvm-env-euk-bw.rds"
  )
 )

 lvm2 <- readRDS(
  file.path(
    path2,
    "lvm-env-euk-w.rds"
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



# Compute the biplot -----------------------------------------------

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
        limits = c(-3, 3)
     ) +
     scale_y_continuous(
        breaks = seq(-3, 3, 1),
        limits = c(-3, 3)
     ) +
     scale_shape_manual(values = c(21, 24)) +
     scale_fill_manual(values = c("#E69F00",
                                  "#666666")) +
     scale_color_manual(values = c("#E69F00",
                                   "#666666")) +
     labs(fill = "Environment :",
          shape = "Environment :") +
     xlab("Latent variable 1") +
     ylab("Latent variable 2") +
     custom_theme



# Add spider image to biplot ---------------------------------------
 
 # specifying the image path
 path3 <- file.path(getwd(), "data", "images")

 # function to open the file
 get_png <- function(filename) {
   grid::rasterGrob(png::readPNG(filename), interpolate = TRUE)
 }
 
 # read the file
 spid <- get_png(file.path(path3, "black-widow-spider.png"))

 # Add the file to the plot
 plot1 <- plot1 + 
 annotation_custom(
     spid, 
     xmin = 1.8, xmax = 3, 
     ymin = 1.8, ymax = 3)

# ==================================================================
# ==================================================================





# ==================================================================
# 4. Plot the model for env webs
# ==================================================================

# Prepare plot data ------------------------------------------------

 alpha2 <- 0.7 # ou 2 comme dans l'exemple?
 site_res2 <- data.frame(lvm2$scores, meta2)
 site_res2$sample_env <- as.factor(site_res2$sample_env)
 sp_res2 <- data.frame(lvm2$loadings,
                      ASV = colnames(comm2))
 fam2 <- data.frame(ASV = rownames(taxa2),
                   family = taxa2[,5])
 sp_res2 <- merge(sp_res2, fam2, by = "ASV")



# Calculate ellipses for each environment --------------------------
 
 # Mean Factor axis for each environment
 lvm_mean2 <- aggregate(
  site_res2[,1:2],
  list(site_res2$sample_env),
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
  for (g in levels(site_res2$sample_env)) {
    df_ell2 <- rbind(df_ell2,
                    cbind(as.data.frame(with(
                      site_res2[site_res2$sample_env == g, ],
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
 plot2 <- ggplot() +
     geom_point(
        data = site_res2,
        aes(x = Factor1,
            y = Factor2,
            fill = sample_env,
            shape = sample_env),
        size = 3
     ) +
     geom_path(
        data = df_ell2,
        aes(x = Factor1, 
            y = Factor2,
            color = sample_env),
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
     scale_shape_manual(values = c(21, 24)) +
     scale_fill_manual(values = c("#E69F00",
                                  "#666666")) +
     scale_color_manual(values = c("#E69F00",
                                   "#666666")) +
     labs(fill = "Environment :",
          shape = "Environment :") +
     xlab("Latent variable 1") +
     ylab("Latent variable 2") +
     custom_theme



# Add web image to biplot ------------------------------------------
 
 # specifying the image path
 path3 <- file.path(getwd(), "data", "images")

 # function to open the file
 get_png <- function(filename) {
   grid::rasterGrob(png::readPNG(filename), interpolate = TRUE)
 }
 
 # read the file
 web <- get_png(file.path(path3, "spider-web.png"))

 # Add the file to the plot
 plot2 <- plot2 + 
 annotation_custom(
     web, 
     xmin = 1.8, xmax = 3, 
     ymin = 1.8, ymax = 3)

# ==================================================================
# ==================================================================





# ==================================================================
# 5. Export the plots
# ==================================================================
 
# # Folder path
# path4 <- file.path(
#  getwd(),
#  "outputs",
#  "plots-euk"
# )
# 
# # Arrange one figure
# fig <- ggarrange(
#     plot1, NULL, plot2,
#     nrow = 1, ncol = 3,
#     labels = c("(A)", "", "(B)"),
#     common.legend = TRUE,
#     legend = "top",
#     hjust = -0.1, vjust = 1.5,
#     widths = c(1, 0.1, 1)
# )
#
# # Export figure
# ggexport(
#  fig,
#  filename = file.path(path4, "lvm-env-euk.png"),
#  width = 2500,
#  height = 1300,
#  res = 300
# )

# ==================================================================
# ==================================================================