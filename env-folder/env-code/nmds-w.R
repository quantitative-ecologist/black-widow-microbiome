# ==================================================================

#               Compare web communities using NMDS

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
 
 # Metadata
 meta_w <- readRDS(file.path(path, "env-bac-metadata-w.rds"))
 
# ==================================================================
# ==================================================================




# ==================================================================
# 2. Compute the NMDS on the web data
# ==================================================================


 # Create Hellinger-transformed version of the community data
 comm_w_hel <- decostand(comm_w, method = "hellinger")
 
 # NMDS ordination based on Bray-Curtis distances
 nmds_w <- metaMDS(comm_w_hel, distance = "bray", trace = FALSE)
 
# ==================================================================
# ==================================================================





# ==================================================================
# 3. Plot the NMDS results
# ==================================================================


 # Plot the NMDS
 ordiplot(nmds_w, cex = 0.6,
          type = "text",
          display = "sites")
 ordiellipse(nmds_w, display= "sites", meta_w$sample_env, label = TRUE)
 


# Extract the nmds data -------------------------------------------

 # NMDS axes + environment
 dat <- data.frame(nmds_w$points[,1:2])
 dat <- cbind(dat, env = meta_w$sample_env)
 dat$env <- as.factor(dat$env)



# Calculate ellipses for each environment -------------------------
 
 # Mean NMDS axis for each environment
 nmds_mean <- aggregate(dat[,1:2], list(dat$env), mean)
 
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
ggplot() +
    geom_point(
       data = dat,
       aes(x = MDS1, y = MDS2,
           color = env,
           shape = env),
       size = 3,
    ) +
    #geom_text(data = sp_res,
    #   aes(x = Factor1*alpha,
    #       y = Factor2*alpha,
    #       label = family),
    #   size = 2.5) +
    geom_path(
        data = df_ell,
        aes(x = MDS1, y = MDS2, color = env),
        size = 1,
        linetype = "dashed",
        show.legend = FALSE
    ) +
    scale_x_continuous(
        breaks = seq(-2, 2, 1),
        limits = c(-2.8, 2.8)
    ) +
    scale_y_continuous(
        breaks = seq(-2, 2, 1),
        limits = c(-2.5, 2.5)
    ) +
    scale_color_manual(
        values = c("#E69F00", "#666666")
    ) +
    xlab("\nNMDS1") + ylab("NMDS2\n") +                              
    theme_bw() + 
    theme(
        legend.position = "top",
        panel.grid = element_blank()
    )

# ==================================================================
# ==================================================================