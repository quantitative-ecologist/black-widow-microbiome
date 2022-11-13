# ==================================================================

#                   Plot the ASV richness plots

# ==================================================================




# ==================================================================
# 1. Import data, libraries, and model fits
# ==================================================================

# Libraries --------------------------------------------------------
 
 library(brms)
 library(data.table)
 library(ggplot2)
 library(cowplot)
 library(viridis)



# Data -------------------------------------------------------------

 # Folder path
 path <- file.path(getwd(), "outputs")

 fit1 <- readRDS(file.path(path, "fits", "env-bac-glm-bw.rds"))
 fit2 <- readRDS(file.path(path, "fits", "env-bac-glm-w.rds"))
 fit3 <- readRDS(file.path(path, "fits", "diet-bac-glm-bw.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Prepare the data to plot
# ==================================================================


# Extract the model draws ------------------------------------------
 
 # Draws into a data frame
 draws1 <- as_draws_df(fit1)
 draws2 <- as_draws_df(fit2)
 draws3 <- as_draws_df(fit3)
 
 # Keep columns of interest
 draws1 <- data.table(draws1[, c(1:3,5)])
 draws2 <- data.table(draws2[, c(1:4)])
 draws3 <- data.table(draws3[, c(1:4, 7:8)])


# Calculate the group means ----------------------------------------
 
 # Calculate differente between Intercept and b (Model 1)
 draws1[, ":=" 
    (mean_urban = b_Intercept + b_sample_envurban,
    sigma_urban = b_sigma_Intercept + b_sigma_sample_envurban)
 ]

 # Calculate differente between Intercept and b (Model 2)
 draws2[, ":=" 
    (mean_urban = b_Intercept + b_sample_envurban,
    sigma_urban = b_sigma_Intercept + b_sigma_sample_envurban)
 ]
 
 draws3[, ":=" 
    (mean_isopod = b_Intercept + b_diet_treatmentisopod,
     sigma_isopod = b_sigma_Intercept + b_sigma_diet_treatmentisopod,
     mean_notfed = b_Intercept + b_diet_treatmentnotfed,
     sigma_notfed = b_sigma_Intercept + b_sigma_diet_treatmentnotfed)
 ]


# Prepare the plot1 table ------------------------------------------
  
 # Reshape the table for easier computation
 draws1 <- draws1[, c(1,2,5,6)]
 setnames(draws1, "b_Intercept", "mean_desert")
 setnames(draws1, "b_sigma_Intercept", "sigma_desert")
 dat1 <- melt(
    draws1,
    variable.name = "coefficient",
    measure = patterns("_")
 )
 
 # Detransform sigma values
 dat1 <- dat1[
    coefficient %in% c("sigma_desert", "sigma_urban"),
    value := exp(value)
 ]

  # 95% CI functions
 lower_int <- function (x) {coda::HPDinterval(as.mcmc(x), 0.95)[1]}
 upper_int <- function (x) {coda::HPDinterval(as.mcmc(x), 0.95)[2]}

 # Calculate the posterior group means and CIs
 dat1 <- dat1[,
  .(mean = mean(value),
    up_ci = upper_int(value),
    lo_ci = lower_int(value)),
  by = "coefficient"
 ]
 
 # Add an env and type factor
 dat1[, env := as.factor(c(rep("desert",2), rep("urban",2)))]
 dat1[, type := as.factor(rep(c("mean", "sigma"), 2))]



# Prepare the plot2 table ------------------------------------------
  
 # Reshape the table for easier computation
 draws2 <- draws2[, c(1,2,5,6)]
 setnames(draws2, "b_Intercept", "mean_desert")
 setnames(draws2, "b_sigma_Intercept", "sigma_desert")
 dat2 <- melt(
    draws2,
    variable.name = "coefficient",
    measure = patterns("_")
 )
 
 # Detransform sigma values
 dat2 <- dat2[
    coefficient %in% c("sigma_desert", "sigma_urban"),
    value := exp(value)
 ]

 # Calculate the posterior group means and CIs
 dat2 <- dat2[,
  .(mean = mean(value),
    up_ci = upper_int(value),
    lo_ci = lower_int(value)),
  by = "coefficient"
 ]
 
 # Add an env and type factor
 dat2[, env := as.factor(c(rep("desert",2), rep("urban",2)))]
 dat2[, type := as.factor(rep(c("mean", "sigma"), 2))]



# Prepare the plot3 table ------------------------------------------
  
 # Reshape the table for easier computation
 draws3 <- draws3[, c(1,2,7:10)]
 setnames(draws3, "b_Intercept", "mean_cricket")
 setnames(draws3, "b_sigma_Intercept", "sigma_cricket")
 dat3 <- melt(
    draws3,
    variable.name = "coefficient",
    measure = patterns("_")
 )
 
 # Detransform sigma values
 dat3 <- dat3[
    coefficient %in% c("sigma_cricket", "sigma_isopod", "sigma_notfed"),
    value := exp(value)
 ]

 # Calculate the posterior group means and CIs
 dat3 <- dat3[,
  .(mean = mean(value),
    up_ci = upper_int(value),
    lo_ci = lower_int(value)),
  by = "coefficient"
 ]
 
 # Add an env and type factor
 dat3[, diet := as.factor(c(rep("cricket",2), rep("isopod",2), rep("not fed", 2)))]
 dat3[, type := as.factor(rep(c("mean", "sigma"), 3))]


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
     axis.title.y  = element_text(size = 13),
     axis.title.x  = element_text(size = 13),
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
# 4. Plot the differences between groups 
# ==================================================================


# Plot for env spiders ---------------------------------------------

 plot1 <- ggplot(dat1,
                 aes(x = type, y = mean,
                 fill = env,
                 shape = env)) +
      
      geom_pointrange(aes(ymin = lo_ci,
                          ymax = up_ci),
                      size = 1,
                      linewidth = 1,
                      position = position_dodge(width = 0.5)) +
      scale_shape_manual(values = c(21, 24)) +
      scale_fill_manual(values = c("#E69F00",
                                   "#666666")) +
      
      scale_y_continuous(breaks = seq(0, 6, 2),
                         limits = c(0, 6.8)) +
      scale_x_discrete(expand = c(1, 0)) +
      labs(fill = "Environment :",
           shape = "Environment :") +
      ylab("Log(posterior predicted richness)") +
      #xlab("\nParameter") +
      custom_theme +
      theme(axis.title.x = element_blank())



# Plot for env webs ------------------------------------------------

 plot2 <- ggplot(dat2,
                 aes(x = type, y = mean,
                 fill = env,
                 shape = env)) +
      
      geom_pointrange(aes(ymin = lo_ci,
                          ymax = up_ci),
                      size = 1,
                      linewidth = 1,
                      position = position_dodge(width = 0.5)) +
      scale_shape_manual(values = c(21, 24)) +
      scale_fill_manual(values = c("#E69F00",
                                   "#666666")) +
      
      scale_y_continuous(breaks = seq(0, 6, 2),
                         limits = c(0, 6.8)) +
      scale_x_discrete(expand = c(1, 0)) +
      labs(fill = "Environment :",
           shape = "Environment :") +
      ylab("Log(posterior predicted richness)") +
      #xlab("\nParameter") +
      custom_theme + 
      theme(axis.title.x = element_blank())



# Plot for diet spiders --------------------------------------------

 plot3 <- ggplot(dat3,
                 aes(x = type, y = mean,
                 fill = diet,
                 shape = diet)) +
      
      geom_pointrange(aes(ymin = lo_ci,
                          ymax = up_ci),
                      size = 1,
                      linewidth = 1,
                      position = position_dodge(width = 0.5)) +
      
      scale_shape_manual(values = c(23, 22, 25)) +
      scale_fill_viridis(discrete = TRUE,
                         option = "viridis") +
      #scale_y_continuous(breaks = seq(0, 6, 2),
      #                   limits = c(0, 6.8)) +
      scale_x_discrete(expand = c(1, 0)) +
      
      labs(fill = "Diet :",
           shape = "Diet :") +
      ylab("Log(posterior predicted richness)") +
      #xlab("\nParameter") +
      custom_theme + 
      theme(axis.title.x = element_blank())

# ==================================================================
# ==================================================================





# ==================================================================
# 5. Combine the plots as one figure
# ==================================================================


# Combine as one figure --------------------------------------------

 # Prepare the figure
 prow <- plot_grid(
    plot1 + theme(legend.position="none"),
    NULL,
    plot2 + theme(legend.position="none"),
    NULL,
    plot3 + theme(legend.position="none"),
    rel_widths = c(1, 0.08, 1, 0.08, 1),
    align = "hv",
    labels = c("(A)","", "(B)", "", "(C)"),
    hjust = -0.2,
    nrow = 1
 )
 
 # Create the two legends
 legend_a <- get_legend(plot1 + theme(legend.position="bottom"))
 legend_b <- get_legend(plot3 + theme(legend.position="bottom"))

 # Combine the plot with legends
 fig <- plot_grid(
   prow,
   plot_grid(legend_a, legend_b),
   ncol = 1, nrow = 2, rel_heights = c(1, .2)
 )

# Export in the outputs folder -------------------------------------
 
 ggexport(
    fig,
    filename = file.path(path, "plots", "ASV-richness-fig.png"),
    width = 3500,
    height = 1400,
    res = 300
 )
 
# ==================================================================
# ==================================================================