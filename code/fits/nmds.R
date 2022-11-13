# ==================================================================

#                Compare web communities using NMDS

# ==================================================================





# ==================================================================
# 1. Import libraries and data
# ==================================================================


# Load libraries ---------------------------------------------------
 
 # Personal computer
 library(picante)
 library(ggplot2)



# Import data ------------------------------------------------------
 
 path <- file.path(getwd(),
                   "data",
                   "env-data-clean")

 # Community matrices
 comm_w <- readRDS(file.path(path, "env-bac-comm-w.rds"))
 
 # Metadata
 meta_w <- readRDS(file.path(path, "env-bac-metadata-w.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Compute NMDSs on the web data
# ==================================================================


 # Create Hellinger-transformed version of the community data
 comm_w_hel <- decostand(comm_w, method = "hellinger")
 
 # NMDS ordination based on Bray-Curtis distances
 nmds_w <- metaMDS(
  comm_w_hel,
  distance = "bray",
  trace = FALSE
 )

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Prepare a custom theme for the plots
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
# 4. Plot the NMDS results for webs
# ==================================================================
 

# Extract the nmds data -------------------------------------------

 # NMDS axes + environment
 dat <- data.frame(nmds_w$points[,1:2])
 dat <- cbind(dat, env = meta_w$sample_env)
 dat$env <- as.factor(dat$env)



# Prepare the plot ------------------------------------------------
 
 # Mean NMDS axis for each environment
 nmds_mean <- aggregate(
  dat[,1:2],
  list(dat$env),
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
df_ell <- data.frame()
 for (g in levels(dat$env)) {
   df_ell <- rbind(df_ell,
                   cbind(as.data.frame(with(
                     dat[dat$env == g, ],
                     veganCovEllipse(cov.wt(
                       cbind(MDS1, MDS2),
                       wt = rep(1 / length(MDS1), length(MDS1))
                     )$cov,
                     center = c(mean(MDS1), mean(MDS2)))
                   )),
                   env = g))
 }



# Compute the plot ------------------------------------------------

 plot1 <- ggplot() +
   geom_point(
    data = dat,
    aes(x = MDS1,
        y = MDS2,
        fill = env,
        shape = env),
    size = 3
   ) +
   geom_path(
     data = df_ell,
     aes(x = MDS1,
         y = MDS2,
         color = env),
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
   scale_shape_manual(
    values = c(21, 24)
   ) +
   scale_fill_manual(
    values = c("#E69F00",
               "#666666")
   ) +
   scale_color_manual(
    values = c("#E69F00",
               "#666666")
   ) +
   labs(
    fill = "Environment :",
    shape = "Environment :"
   ) +
   xlab("\nNMDS1") + ylab("NMDS2\n") +
   custom_theme

# ==================================================================
# ==================================================================





# ==================================================================
# 5. Save the plot
# ==================================================================
 
 # Folder
 path1 <- file.path(
  getwd(),
  "outputs",
  "plots"
 )
 
 # Export to folder
 ggsave(
  plot1,
  width = 12,
  height = 12,
  dpi = 300,
  units = c("cm"),
  file = file.path(path1, "env-bac-nmds-w.png")
 )



# ==================================================================
# ==================================================================