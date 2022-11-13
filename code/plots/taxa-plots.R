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
 path1 <- file.path(getwd(), "data", "env-data-clean")
 path2 <- file.path(getwd(), "data", "diet-data-clean")
 
 # Community matrices
 comm1 <- readRDS(file.path(path1, "env-bac-comm-bw.rds"))
 comm2 <- readRDS(file.path(path1, "env-bac-comm-w.rds"))
 comm3 <- readRDS(file.path(path2, "diet-bac-comm-bw.rds"))
 
 # Taxa data
 tax1 <- readRDS(file.path(path1, "env-bac-taxa-bw.rds"))
 tax2 <- readRDS(file.path(path1, "env-bac-taxa-w.rds"))
 tax3 <- readRDS(file.path(path2, "diet-bac-taxa-bw.rds"))
 
 # Metadata
 meta1 <- readRDS(file.path(path1, "env-bac-metadata-bw.rds"))
 meta2 <- readRDS(file.path(path1, "env-bac-metadata-w.rds"))
 meta3 <- readRDS(file.path(path2, "diet-bac-metadata-bw.rds"))
 
# ==================================================================
# ==================================================================





# ==================================================================
# 2a. Plot the taxa for the env spiders --> Each sample
# ==================================================================


# Prepare data for the plot ----------------------------------------

 # community data aggregation at the family level :
 
 # take the sum of each phylum in each sample
 taxa_agg1 <- aggregate(t(comm1),
                        by = list(tax1[colnames(comm1),5]),
                        FUN = sum)
 
 # clean up resulting object
 rownames(taxa_agg1) <- taxa_agg1$Group.1
 fam_data1 <- t(taxa_agg1[,-1])
 
 # convert abundances to relative abundances
 fam_data1 <- fam_data1/rowSums(fam_data1)
 
 # remove rare families
 fam_data1 <- fam_data1[, colSums(fam_data1) > 0.01]
 
 # now reshape phylum data to long format
 fam_data1 <- reshape2::melt(fam_data1)
 
 # rename columns
 colnames(fam_data1)[1:2] <- c("sample_id", "bacteria_family")
 
 # Add env column
 meta1$sample_id <- as.factor(meta1$sample_id)
 fam_data1 <- merge(
    fam_data1,
    meta1[, c("sample_id", "sample_env")],
    by = "sample_id"
  )



# Compute the plot -------------------------------------------------

 # now we can plot phylum relative abundance per sample
 plot1a <- ggplot(fam_data1,
                  aes(sample_id,
                      weight = value,
                      fill = bacteria_family)) +
            geom_bar(color = "black",
                     width = .7,
                     position = "fill") +
            ylab("\nRelative abundance (%)") +
            xlab("Sample ID") +
            labs(fill = "Bacteria family") +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_viridis_d(direction = -1L) +
            theme_classic() +
            theme(axis.title.y = element_blank(),
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
 fam_data_agg1 <- aggregate(
    fam_data1$value,
    by = list(meta1[fam_data1$sample_id,]$sample_env,
              fam_data1$bacteria_family),
              FUN = mean
 )
 
 # rename columns
 colnames(fam_data_agg1) <- c("env", "bacteria_family", "value")



# Compute the plot -------------------------------------------------

 # Now we can plot phylum relative abundance by environment
 plot1b <- ggplot(fam_data_agg1,
                  aes(env,
                      weight = value,
                      fill = bacteria_family)) +
            geom_bar(color = "black",
                     width = .7,
                     position = "fill") +
            ylab("\nRelative abundance (%)") +
            xlab("Environment") +
            labs(fill = "Bacteria family") +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_viridis_d(direction = -1L) +
            theme_classic() +
            theme(axis.title.y = element_blank(),
                  panel.spacing = unit(1, "cm",
                                       data = NULL)) +
            coord_flip()

# ==================================================================
# ==================================================================





# ==================================================================
# 3a. Plot the taxa for the env webs --> Each sample
# ==================================================================


# Prepare data for the plot ----------------------------------------

 # community data aggregation at the class level :
 
 # take the sum of each phylum in each sample
 taxa_agg2 <- aggregate(t(comm2),
                        by = list(tax2[colnames(comm2),3]),
                        FUN = sum)
 
 # clean up resulting object
 rownames(taxa_agg2) <- taxa_agg2$Group.1
 fam_data2 <- t(taxa_agg2[,-1])
 
 # convert abundances to relative abundances
 fam_data2 <- fam_data2/rowSums(fam_data2)
 
 # remove rare families
 fam_data2 <- fam_data2[, colSums(fam_data2) > 0.01]
 
 # now reshape phylum data to long format
 fam_data2 <- reshape2::melt(fam_data2)
 
 # rename columns
 colnames(fam_data2)[1:2] <- c("sample_id", "bacteria_class")
 
 # Add env column
 meta2$sample_id <- as.factor(meta2$sample_id)
 fam_data2 <- merge(
    fam_data2,
    meta2[, c("sample_id", "sample_env")],
    by = "sample_id"
  )



# Compute the plot -------------------------------------------------

 # now we can plot phylum relative abundance per sample
 plot2a <- ggplot(fam_data2,
                  aes(sample_id,
                      weight = value,
                      fill = bacteria_class)) +
            geom_bar(color = "black",
                     width = .7,
                     position = "fill") +
            ylab("\nRelative abundance (%)") +
            xlab("Sample ID") +
            labs(fill = "Bacteria class") +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_viridis_d(direction = -1L) +
            theme_classic() +
            theme(axis.title.y = element_blank(),
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
 fam_data_agg2 <- aggregate(
    fam_data2$value,
    by = list(meta2[fam_data2$sample_id,]$sample_env,
              fam_data2$bacteria_class),
              FUN = mean
 )
 
 # rename columns
 colnames(fam_data_agg2) <- c("env", "bacteria_class", "value")



# Compute the plot -------------------------------------------------

 # Now we can plot phylum relative abundance by environment
 plot2b <- ggplot(fam_data_agg2,
                  aes(env,
                      weight = value,
                      fill = bacteria_class)) +
            geom_bar(color = "black",
                     width = .7,
                     position = "fill") +
            ylab("\nRelative abundance (%)") +
            xlab("Environment") +
            labs(fill = "Bacteria class") +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_viridis_d(direction = -1L) +
            theme_classic() +
            theme(axis.title.y = element_blank(),
                  panel.spacing = unit(1, "cm",
                                       data = NULL)) +
            coord_flip()

# ==================================================================
# ==================================================================





# ==================================================================
# 4a. Plot the taxa for the diet spiders --> Each sample
# ==================================================================


# Prepare data for the plot ----------------------------------------

 # community data aggregation at the class level :
 
 # take the sum of each phylum in each sample
 taxa_agg3 <- aggregate(t(comm3),
                        by = list(tax3[colnames(comm3),3]),
                        FUN = sum)
 
 # clean up resulting object
 rownames(taxa_agg3) <- taxa_agg3$Group.1
 fam_data3 <- t(taxa_agg3[,-1])
 
 # convert abundances to relative abundances
 fam_data3 <- fam_data3/rowSums(fam_data3)
 
 # remove rare families
 fam_data3 <- fam_data3[, colSums(fam_data3) > 0.01]
 
 # now reshape phylum data to long format
 fam_data3 <- reshape2::melt(fam_data3)
 
 # rename columns
 colnames(fam_data3)[1:2] <- c("sample_id", "bacteria_class")
 
 # Add env column
 meta3$sample_id <- as.factor(meta3$sample_id)
 fam_data3 <- merge(
    fam_data3,
    meta3[, c("sample_id", "diet_treatment")],
    by = "sample_id"
  )



# Compute the plot -------------------------------------------------

 # now we can plot phylum relative abundance per sample
 plot3a <- ggplot(fam_data3,
                  aes(sample_id,
                      weight = value,
                      fill = bacteria_class)) +
            geom_bar(color = "black",
                     width = .7,
                     position = "fill") +
            ylab("\nRelative abundance (%)") +
            xlab("Sample ID") +
            labs(fill = "Bacteria class") +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_viridis_d(direction = -1L) +
            theme_classic() +
            theme(axis.title.y = element_blank(),
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
 fam_data_agg3 <- aggregate(
    fam_data3$value,
    by = list(meta3[fam_data3$sample_id,]$diet_treatment,
              fam_data3$bacteria_class),
              FUN = mean
 )
 
 # rename columns
 colnames(fam_data_agg3) <- c("env", "bacteria_class", "value")



# Compute the plot -------------------------------------------------

 # Now we can plot phylum relative abundance by environment
 plot3b <- ggplot(fam_data_agg3,
                  aes(env,
                      weight = value,
                      fill = bacteria_class)) +
            geom_bar(color = "black",
                     width = .7,
                     position = "fill") +
            ylab("\nRelative abundance (%)") +
            xlab("Diet treatment") +
            labs(fill = "Bacteria class") +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_viridis_d(direction = -1L) +
            theme_classic() +
            theme(axis.title.y = element_blank(),
                  panel.spacing = unit(1, "cm",
                                       data = NULL)) +
            coord_flip()

# ==================================================================
# ==================================================================





# ==================================================================
# 5. Export the plots as figures
# ==================================================================
 
 # Folder path
 path <- file.path(getwd(), "outputs", "plots")

 # Env spiders
 ggsave(plot1a, file = file.path(path, "env-bac-taxa-bw1.png"))
 ggsave(plot1b, file = file.path(path, "env-bac-taxa-bw2.png"))

 # Env webs
 ggsave(plot2a, file = file.path(path, "env-bac-taxa-w1.png"))
 ggsave(plot2b, file = file.path(path, "env-bac-taxa-w2.png"))
 
 # Diet spiders
 ggsave(plot3a, file = file.path(path, "diet-bac-taxa-bw1.png"))
 ggsave(plot3b, file = file.path(path, "diet-bac-taxa-bw2.png"))
# ==================================================================
# ==================================================================