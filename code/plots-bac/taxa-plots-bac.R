# ==================================================================

#                Plot the taxa for each samples

# ==================================================================





# ==================================================================
# 1. Import the libraries and the data
# ==================================================================


# Libraries --------------------------------------------------------
 
 library(ggplot2)
 library(viridis)
 library(ggpubr)



# Community and taxa data ------------------------------------------
 
 # Folder paths
 path1 <- file.path(getwd(), "data", "data-clean-env")
 path2 <- file.path(getwd(), "data", "data-clean-diet")
 
 # Community matrices
 comm1 <- readRDS(file.path(path1, "comm-env-bac-bw.rds"))
 comm2 <- readRDS(file.path(path1, "comm-env-bac-w.rds"))
 comm3 <- readRDS(file.path(path2, "comm-diet-bac.rds"))
 
 # Taxa data
 tax1 <- readRDS(file.path(path1, "taxa-env-bac-bw.rds"))
 tax2 <- readRDS(file.path(path1, "taxa-env-bac-w.rds"))
 tax3 <- readRDS(file.path(path2, "taxa-diet-bac.rds"))
 
 # Metadata
 meta1 <- readRDS(file.path(path1, "metadata-env-bac-bw.rds"))
 meta2 <- readRDS(file.path(path1, "metadata-env-bac-w.rds"))
 meta3 <- readRDS(file.path(path2, "metadata-diet-bac.rds"))
 
# ==================================================================
# ==================================================================





# ==================================================================
# 2a. Plot the taxa for the env spiders --> Each sample
# ==================================================================


# Prepare data for the plot ----------------------------------------

 # community data aggregation at the phylum level :
 
 # take the sum of each phylum in each sample
 taxa_agg1 <- aggregate(t(comm1),
                        by = list(tax1[colnames(comm1),2]),
                        FUN = sum)
 
 # clean up resulting object
 rownames(taxa_agg1) <- taxa_agg1$Group.1
 data1 <- t(taxa_agg1[,-1])
 
 # convert abundances to relative abundances
 data1 <- data1/rowSums(data1)
 
 # remove rare families
 data1 <- data1[, colSums(data1) > 0.01]
 
 # now reshape phylum data to long format
 data1 <- reshape2::melt(data1)
 
 # rename columns
 colnames(data1)[1:2] <- c("sample_id", "bacteria_phylum")
 
 # Add env column
 meta1$sample_id <- as.factor(meta1$sample_id)
 data1 <- merge(
    data1,
    meta1[, c("sample_id", "sample_env")],
    by = "sample_id"
  )

 # Simplify sample names
 vec1 <- as.character(data1$sample_id)
 data1$sample_id <- as.character(strsplit(vec1,"-bac"))


# Compute the plot -------------------------------------------------

 # now we can plot phylum relative abundance per sample
 plot1a <- ggplot(data1,
                  aes(sample_id,
                      weight = value,
                      fill = bacteria_phylum)) +
            geom_bar(color = "black",
                     width = .7,
                     position = "fill") +
            ylab("\nRelative abundance (%)") +
            xlab("Sample ID") +
            #ggtitle("Samples") +
            labs(fill = "Phylum") +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_viridis_d(direction = -1L) +
            theme_classic() +
            theme(plot.title = element_text(face = "bold"),
                  axis.title.y = element_blank(),
                  legend.background = element_blank(),
                  legend.box.background = element_blank(),
                  panel.spacing = unit(1, "cm",
                                       data = NULL)) +
            coord_flip()

# ==================================================================
# ==================================================================





# ==================================================================
# 2b. Plot the taxa for the env spiders --> Each environment
# ==================================================================


# Prepare data for the plot ----------------------------------------
 
 # aggregate average family abundances per environment
 data_agg1 <- aggregate(
    data1$value,
    by = list(meta1[data1$sample_id,]$sample_env,
              data1$bacteria_phylum),
              FUN = mean
 )
 
 # rename columns
 colnames(data_agg1) <- c("env", "bacteria_phylum", "value")



# Compute the plot -------------------------------------------------

 # Now we can plot phylum relative abundance by environment
 plot1b <- ggplot(data_agg1,
                  aes(env,
                      weight = value,
                      fill = bacteria_phylum)) +
            geom_bar(color = "black",
                     width = .7,
                     position = "fill") +
            #ggtitle("Habitat") +
            ylab("\nRelative abundance (%)") +
            xlab("Environment") +
            labs(fill = "Bacteria phylum") +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_viridis_d(direction = -1L) +
            theme_classic() +
            theme(plot.title = element_text(face = "bold"),
                  axis.title.y = element_blank(),
                  legend.background = element_blank(),
                  legend.box.background = element_blank(),
                  panel.spacing = unit(1, "cm",
                                       data = NULL)) +
            coord_flip()

# ==================================================================
# ==================================================================





# ==================================================================
# 3a. Plot the taxa for the env webs --> Each sample
# ==================================================================


# Prepare data for the plot ----------------------------------------

 # community data aggregation at the phylum level :
 
 # take the sum of each phylum in each sample
 taxa_agg2 <- aggregate(t(comm2),
                        by = list(tax2[colnames(comm2),2]),
                        FUN = sum)
 
 # clean up resulting object
 rownames(taxa_agg2) <- taxa_agg2$Group.1
 data2 <- t(taxa_agg2[,-1])
 
 # convert abundances to relative abundances
 data2 <- data2/rowSums(data2)
 
 # remove rare families
 data2 <- data2[, colSums(data2) > 0.01]
 
 # now reshape phylum data to long format
 data2 <- reshape2::melt(data2)
 
 # rename columns
 colnames(data2)[1:2] <- c("sample_id", "bacteria_phylum")
 
 # Add env column
 meta2$sample_id <- as.factor(meta2$sample_id)
 data2 <- merge(
    data2,
    meta2[, c("sample_id", "sample_env")],
    by = "sample_id"
  )

 # Simplify sample id
 vec2 <- as.character(data2$sample_id)
 data2$sample_id <- as.character(strsplit(vec2, "-bac"))



# Compute the plot -------------------------------------------------

 # now we can plot phylum relative abundance per sample
 plot2a <- ggplot(data2,
                  aes(sample_id,
                      weight = value,
                      fill = bacteria_phylum)) +
            geom_bar(color = "black",
                     width = .7,
                     position = "fill") +
            #ggtitle("Samples") +
            ylab("\nRelative abundance (%)") +
            xlab("Sample ID") +
            labs(fill = "Phylum") +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_viridis_d(direction = -1L) +
            theme_classic() +
            theme(plot.title = element_text(face = "bold"),
                  axis.title.y = element_blank(),
                  legend.background = element_blank(),
                  legend.box.background = element_blank(),
                  panel.spacing = unit(1, "cm",
                                       data = NULL)) +
            coord_flip()

# ==================================================================
# ==================================================================





# ==================================================================
# 3b. Plot the taxa for the env webs --> Each environment
# ==================================================================


# Prepare data for the plot ----------------------------------------
 
 # aggregate average family abundances per environment
 data_agg2 <- aggregate(
    data2$value,
    by = list(meta2[data2$sample_id,]$sample_env,
              data2$bacteria_phylum),
              FUN = mean
 )
 
 # rename columns
 colnames(data_agg2) <- c("env", "bacteria_phylum", "value")



# Compute the plot -------------------------------------------------

 # Now we can plot phylum relative abundance by environment
 plot2b <- ggplot(data_agg2,
                  aes(env,
                      weight = value,
                      fill = bacteria_phylum)) +
            geom_bar(color = "black",
                     width = .7,
                     position = "fill") +
            #ggtitle("Habitat") +
            ylab("\nRelative abundance (%)") +
            xlab("Environment") +
            labs(fill = "Phylum") +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_viridis_d(direction = -1L) +
            theme_classic() +
            theme(plot.title = element_text(face = "bold"),
                  axis.title.y = element_blank(),
                  legend.background = element_blank(),
                  legend.box.background = element_blank(),
                  panel.spacing = unit(1, "cm",
                                       data = NULL)) +
            coord_flip()

# ==================================================================
# ==================================================================





# ==================================================================
# 4a. Plot the taxa for the diet spiders --> Each sample
# ==================================================================


# Prepare data for the plot ----------------------------------------

 # community data aggregation at the phylum level :
 
 # take the sum of each phylum in each sample
 taxa_agg3 <- aggregate(t(comm3),
                        by = list(tax3[colnames(comm3),2]),
                        FUN = sum)
 
 # clean up resulting object
 rownames(taxa_agg3) <- taxa_agg3$Group.1
 data3 <- t(taxa_agg3[,-1])
 
 # convert abundances to relative abundances
 data3 <- data3/rowSums(data3)
 
 # remove rare phylum
 data3 <- data3[, colSums(data3) > 0.01]
 
 # now reshape phylum data to long format
 data3 <- reshape2::melt(data3)
 
 # rename columns
 colnames(data3)[1:2] <- c("sample_id", "bacteria_phylum")
 
 # Add env column
 meta3$sample_id <- as.factor(meta3$sample_id)
 data3 <- merge(
    data3,
    meta3[, c("sample_id", "diet_treatment")],
    by = "sample_id"
  )


# Compute the plot -------------------------------------------------

 # now we can plot phylum relative abundance per sample
 plot3a <- ggplot(data3,
                  aes(sample_id,
                      weight = value,
                      fill = bacteria_phylum)) +
            geom_bar(color = "black",
                     width = .7,
                     position = "fill") +
            #ggtitle("Spiders in the diet experiment") +
            ylab("\nRelative abundance (%)") +
            xlab("Sample ID") +
            labs(fill = "Phylum") +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_viridis_d(direction = -1L) +
            theme_classic() +
            theme(plot.title = element_text(face = "bold"),
                  axis.title.y = element_blank(),
                  legend.key.size = unit(0.5, "cm"),
                  panel.spacing = unit(1, "cm",
                                       data = NULL)) +
            coord_flip()

# ==================================================================
# ==================================================================




# ==================================================================
# 4b. Plot the taxa for the diet spiders --> Each diet
# ==================================================================


# Prepare data for the plot ----------------------------------------
 
 # aggregate average family abundances per environment
 data_agg3 <- aggregate(
    data3$value,
    by = list(meta3[data3$sample_id,]$diet_treatment,
              data3$bacteria_phylum),
              FUN = mean
 )
 
 # rename columns
 colnames(data_agg3) <- c("env", "bacteria_phylum", "value")



# Compute the plot -------------------------------------------------

 # Now we can plot phylum relative abundance by environment
 plot3b <- ggplot(data_agg3,
                  aes(env,
                      weight = value,
                      fill = bacteria_phylum)) +
            geom_bar(color = "black",
                     width = .7,
                     position = "fill") +
            #ggtitle("Spiders in the diet experiment") +
            ylab("\nRelative abundance (%)") +
            xlab("Diet treatment") +
            labs(fill = "Phylum") +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_viridis_d(direction = -1L) +
            theme_classic() +
            theme(plot.title = element_text(face = "bold"),
                  axis.title.y = element_blank(),
                  legend.key.size = unit(0.5, "cm"),
                  panel.spacing = unit(1, "cm",
                                       data = NULL)) +
            coord_flip()

# ==================================================================
# ==================================================================





# ==================================================================
# 5. Export the plots as figures
# ==================================================================
 
 # Folder path
 path <- file.path(getwd(), "outputs", "plots-bac")



# Env figure -------------------------------------------------------
 
 # Env spiders
 fig1a <- ggarrange(
  plot1a, NULL, plot1b,
  nrow = 1, ncol = 3,
  common.legend = TRUE,
  legend = "right",
  widths = c(1, 0.05, 0.9)
 )

 # Env webs
 fig1b <- ggarrange(
  plot2a, NULL, plot2b,
  nrow = 1, ncol = 3,
  common.legend = TRUE,
  legend = "right",
  widths = c(1, 0.05, 0.9)
 )
 
 # Combine as 1 figure
 fig1c <- ggarrange(
   fig1a, fig1b,
   nrow = 2, ncol = 1
 )
 
 # Add labels to the figure
 fig1c <- annotate_figure(
   fig1c,
   top = text_grob(
     "Sample", hjust = 4, vjust = -0.5,
     face = "bold", size = "12"),
   left = text_grob(
     "Spiders", hjust = -2.4, vjust = 0,
     face = "bold", size = "12", rot = 90)
 )

 fig1c <- annotate_figure(
   fig1c,
   top = text_grob(
     "Habitat", hjust = -1.55, vjust = 1.25,
     face = "bold", size = "12"),
   left = text_grob(
     "Webs", hjust = 3.9, vjust = 1.55,
     face = "bold", size = "12", rot = 90)
 )

 # Export the figure
 ggexport(
   fig1c,
   filename = file.path(path, "taxa-env-bac.png"),
   width = 2500, height = 2000, res = 300
 )



# Diet figure ------------------------------------------------------

 # Diet as one figure
 fig2 <- ggarrange(
  plot3a, NULL, plot3b, NULL,
  nrow = 1, ncol = 4,
  common.legend = TRUE,
  legend = "right",
  widths = c(1, 0.05, 0.95, 0.01)
 )
 
 # Annotate the figure
 fig2 <- annotate_figure(
   fig2,
   top = text_grob(
     "Sample", hjust = 4.3, vjust = -0.5,
     face = "bold", size = "12")
 )

 fig2 <- annotate_figure(
   fig2,
   top = text_grob(
     "Diet", hjust = -3, vjust = 1.25,
     face = "bold", size = "12")
 )
 
 # Export the figure
 ggexport(
   fig2,
   filename = file.path(path, "taxa-diet-bac.png"),
   width = 2500, height = 1200, res = 300
 )

# ==================================================================
# ==================================================================