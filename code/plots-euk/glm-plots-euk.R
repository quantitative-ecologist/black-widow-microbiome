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
 library(ggpubr)



# Data -------------------------------------------------------------

 # Folder path
 path <- file.path(getwd(), "outputs")

 fit1 <- readRDS(file.path(path, "fits-euk", "glm-env-euk-bw.rds"))
 fit2 <- readRDS(file.path(path, "fits-euk", "glm-env-euk-w.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Prepare the data to plot
# ==================================================================


# Extract the model draws ------------------------------------------
 
 # Draws into a data frame
 draws1 <- as_draws_df(fit1)
 draws2 <- as_draws_df(fit2)
 
 # Keep columns of interest
 draws1 <- data.table(draws1[, c(1:3,5)])
 draws2 <- data.table(draws2[, c(1:3,5)])


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
      ylab("log(posterior predicted richness)") +
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
      ylab("log(posterior predicted richness)") +
      #xlab("\nParameter") +
      custom_theme + 
      theme(axis.title.x = element_blank())

# ==================================================================
# ==================================================================




# ==================================================================
# 5. Add web and spider icons to plots
# ==================================================================

 # Specifying the image path
 path1 <- file.path(getwd(), "data", "images")

 # function to open the file
 get_png <- function(filename) {
   grid::rasterGrob(png::readPNG(filename), interpolate = TRUE)
 }
 
 # read the file
 spider <- get_png(file.path(path1, "black-widow-spider.png"))
 web <- get_png(file.path(path1, "spider-web.png"))

 # Add the files to the plots
 plot1 <- plot1 + 
 annotation_custom(
     spider, 
     xmin = 2, xmax = 3, 
     ymin = 5.7, ymax = 6.7)

 plot2 <- plot2 + 
 annotation_custom(
     web, 
     xmin = 2, xmax = 3.2, 
     ymin = 5.5, ymax = 6.7)    

# ==================================================================
# ==================================================================





# ==================================================================
# 6. Combine the plots as one figure
# ==================================================================


# Combine as one figure --------------------------------------------

 # Prepare the figure
 fig <- ggarrange(
    plot1, NULL, plot2,
    nrow = 1, ncol = 3,
    labels = c("(A)", "", "(B)"),
    hjust = -0.1,
    vjust = 1.5,
    widths = c(1, 0.1, 1),
    common.legend = TRUE,
    legend.position = "top"
 )

 saveRDS(
    fig, 
    file = file.path(path, "plots-euk", "glm-env-euk-plot.rds")
 )


# Export in the outputs folder -------------------------------------
 
# ggexport(
#    fig[1],
#    filename = file.path(path, "plots-euk", "glm-env-euk.png"),
#    width = 2000,
#    height = 1200,
#    res = 300
# )

# ==================================================================
# ==================================================================