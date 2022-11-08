# ==================================================================

#                   Compare communities using NMDS

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
                   "env-folder",
                   "env-data",
                   "env-data-clean")

 # Community matrices
 comm_w <- readRDS(file.path(path, "env-bac-comm-w.rds"))
 comm_bw <- readRDS(file.path(path, "env-bac-comm-bw.rds"))
 
 # Metadata
 meta_w <- readRDS(file.path(path, "env-bac-metadata-w.rds"))
 meta_bw <- readRDS(file.path(path, "env-bac-metadata-bw.rds"))

# ==================================================================
# ==================================================================




# ==================================================================
# 2. Compute NMDSs on the web and spider data
# ==================================================================


 # Create Hellinger-transformed version of the community data
 comm_w_hel <- decostand(comm_w, method = "hellinger")
 comm_bw_hel <- decostand(comm_bw, method = "hellinger")
 
 # NMDS ordination based on Bray-Curtis distances
 nmds_w <- metaMDS(comm_w_hel, distance = "bray", trace = FALSE)
 nmds_bw <- metaMDS(comm_bw_hel, distance = "bray", trace = FALSE)
 
 plot(cca(comm_w_hel), type = "points")
 plot(cca(comm_bw_hel), type = "points")
 
 biplot(pcoa(vegdist(comm_w_hel,"bray")))
 biplot(pcoa(vegdist(comm_bw_hel,"bray")))
 
 obj1 <- pcoa(vegdist(comm_w_hel,"bray"))
 obj2 <- pcoa(vegdist(comm_bw_hel,"bray"))
# ==================================================================
# ==================================================================





# ==================================================================
# 3. Prepare a custom theme for the plots
# ==================================================================

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



# Calculate ellipses for each environment -------------------------
 
 # Mean NMDS axis for each environment
 nmds_mean <- aggregate(dat[,1:2], list(dat$env), mean)
 

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

# Draw the plot
 plot1 <- ggplot() +
   geom_point(data = dat,
              aes(
                x = MDS1,
                y = MDS2,
                color = env,
                shape = env
              ),
              size = 3,) +
   geom_path(
     data = df_ell,
     aes(x = MDS1, y = MDS2, color = env),
     size = 1,
     linetype = "dashed",
     show.legend = FALSE
   ) +
   scale_x_continuous(breaks = seq(-2, 2, 1),
                      limits = c(-2.8, 2.8)) +
   scale_y_continuous(breaks = seq(-2, 2, 1),
                      limits = c(-2.5, 2.5)) +
   scale_color_manual(values = c("#E69F00", "#666666")) +
   labs(color = "Environment :",
        shape = "Environment :") +
   xlab("\nNMDS1") + ylab("NMDS2\n") +
   custom_theme

# ==================================================================
# ==================================================================





# ==================================================================
# 5. Plot the NMDS results for spiders
# ==================================================================
 

# Extract the nmds data -------------------------------------------

 # NMDS axes + environment
 dat2 <- data.frame(nmds_bw$points[,1:2])
 dat2 <- cbind(dat2, env = meta_bw$sample_env)
 dat2$env <- as.factor(dat2$env)



# Calculate ellipses for each environment -------------------------
 
 # Mean NMDS axis for each environment
 nmds_mean2 <- aggregate(dat2[,1:2], list(dat2$env), mean)


 # Create a data frame of the ellipses results for ggplot
df_ell2 <- data.frame()
 for (g in levels(dat2$env)) {
   df_ell2 <- rbind(df_ell2,
                   cbind(as.data.frame(with(
                     dat2[dat2$env == g, ],
                     veganCovEllipse(cov.wt(
                       cbind(MDS1, MDS2),
                       wt = rep(1 / length(MDS1), length(MDS1))
                     )$cov,
                     center = c(mean(MDS1), mean(MDS2)))
                   )),
                   env = g))
 }



# Draw the plot ---------------------------------------------------


 plot2 <- ggplot() +
   geom_point(data = dat2,
              aes(
                x = MDS1,
                y = MDS2,
                color = env,
                shape = env
              ),
              size = 3,) +
   geom_path(
     data = df_ell2,
     aes(x = MDS1, y = MDS2, color = env),
     size = 1,
     linetype = "dashed",
     show.legend = FALSE
   ) +
   scale_x_continuous(
       breaks = seq(-1.5, 1.5, 0.5),
       limits = c(-1.5, 1.5)
   ) +
   scale_y_continuous(
       breaks = seq(-1.5, 1.5, 0.5),
       limits = c(-1.5, 1.5)
   ) +
   labs(color = "Environment :",
        shape = "Environment :") +
   scale_color_manual(values = c("#E69F00", "#666666")) +
   xlab("\nNMDS1") + ylab("NMDS2\n") +
   custom_theme

# ==================================================================
# ==================================================================





# ==================================================================
# 6. Save the plots
# ==================================================================
 
 # Foler
 path1 <- file.path(getwd(), "env-folder", "env-outputs")
 
 # Webs
 ggsave(plot1, file = file.path(path1, "nmds-plot-w.png"))

 # Spiders
 ggsave(plot2, file = file.path(path1, "nmds-plot-bw.png"))



# ==================================================================
# ==================================================================