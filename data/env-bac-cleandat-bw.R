# ==================================================================

#               Clean the spider data for analysis

# ==================================================================

# Code for sites in samples :
    # Urban samples :
    # UA = University of Arizona 
    # LO = Lowes

    # Desert samples :
    # DM = Dove mountain
    # CC = Chaos canyon

# W = web, VN = Black widow





# ==================================================================
# 1. Import packages and data
# ==================================================================


# Import libraries -------------------------------------------------

 library(picante)
 library(ggplot2)



# Import data ------------------------------------------------------

 # Folder path
 folder <- "./data"
 # Eventually, load the files from the OSF repo
 
 # Community data
 comm <- readRDS(
    file.path(folder, "env-data-raw",
              "env-bac-seqtabnochim.rds"))

  # Taxonomy data
 taxa_sp <- readRDS(
    file.path(folder, "env-data-raw",
              "env-bac-taxa-table2.rds"))
 # make sur order of ASVs match 
 taxo <- taxa_sp[colnames(comm),]
 rm(taxa_sp)

  # Raw metadata
 metadata <- read.csv(
    file.path(folder, "env-data-raw",
              "env-bac-metadata-raw.csv"),
              row.names = 1)
 # only spider samples
 metadata <- metadata[metadata$sample_type %in%
                      c("spider", "control"),]

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Assemble raw taxa and community data
# ==================================================================

# Rename ASVs from sequences to ASV number -------------------------
 
 # Create dataframe
 ASV_seq_info <- data.frame(ASV_name = paste0("ASV_",
                                              1:dim(comm)[2]),
                            ASV_sequence = colnames(comm))
 
 # Apply new names to comm and taxo tables
 colnames(comm) <- ASV_seq_info$ASV_name
 rownames(taxo) <- ASV_seq_info$ASV_name

 # Remove object
 rm(ASV_seq_info)



# Remove non target DNA --------------------------------------------
 
 # Inspect for ASVs other than bacteria in "domain"
 table(taxo[,"domain"]) # 2 eukaryota
 
 # Remove the eukaryotas
 taxo <- subset(taxo, taxo[,"domain"]!="Eukaryota")

 
 # Inspect if there are chloroplasts or mitochondria
 table(taxo[,"order"])["Chloroplast"]
 table(taxo[,"family"])["Mitochondria"]
 
 # Delete any mitochondria or chloroplast
 taxo <- subset(taxo,
                taxo[, "order"]!= "Chloroplast" &
                taxo[, "family"]!= "Mitochondria")


 # Apply changes to the community data
 comm <- comm[, rownames(taxo)]

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Data exploration
# ==================================================================



# Summary statistics -----------------------------------------------

 # Separate community to have only spider samples
 comm_bw <- comm[c(1:3, 7:10, 15:18, 23, 24:26), ]


 # Number of reads per sample
 rowSums(comm_bw)
 
 # Adjust the new number of reads to the metadata
 metadata$n_reads <- rowSums(comm_bw)

 # visualize log10 number of reads per sample
 hist(rowSums(comm_bw))
 hist(log10(rowSums(comm_bw)))
 # Most samples have around 35K to 100K reads

 # log10 of number of reads per ASV
 hist(log10(colSums(comm_bw)))



# Visualize the community ------------------------------------------

# PCA on Hellinger-transformed community data
comm_pca <- prcomp(decostand(comm_bw, "hellinger"))

# Prepare the ordination results for plotting
labs <- as.factor(metadata$sample_env)
score <- scores(comm_pca)[, 1:2]
score <- cbind(score, sample_env = metadata$sample_env)
score <- data.frame(score)
score[,1] <- as.numeric(score[,1])
score[,2]<- as.numeric(score[,2])
score[,3] <- as.factor(score[,3])
score$sample_id <- as.factor(rownames(score))
sp <- data.frame(comm_pca$rotation[, 1:2])
summ <- summary(comm_pca)$importance[2, 1:2]

# Plot the results
ggplot() +
  geom_point(data = sp,
             aes(x = PC1,
                 y = PC2),
             shape = 3,
             size = 1) +
  geom_text(data = score,
            aes(x = PC1,
                y = PC2,
                label = sample_id,
                color = sample_env),
            size = 3) +
  scale_x_continuous(breaks = seq(-1, 1, 0.5),
                     limits = c(-1, 1)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.5),
                     limits = c(-1, 1)) +
  scale_color_manual(values = c("black",
                               "#E69F00",
                               "#666666")) +
  labs(color = "Environment :") +
  xlab("\nPC1 (34.8%)") + ylab("PC2 (18.9%)\n") +
  theme_bw() + 
  theme(legend.position = "top",
        panel.grid = element_blank())

 # LO-VN-114-BAC is very similar to the negative control.

 # Number of sequences per sample mapped onto ordination axes
 ordisurf(comm_pca, rowSums(comm_bw),
          bubble = TRUE, cex = 2,
          main = "Library size (sequences/sample)")



# Check negative controls -----------------------------------------

 # Abundance of ASVs in negative control
 comm_bw["PCR-neg-CTRL-bac",][
     comm_bw["PCR-neg-CTRL-bac",]>0]
 # Some problems here. Some ASVs are above 100
 # See ASV_172, ASV_354, and ASV_399

 # Check the taxonomic identity of ASVs present in negative control
 taxo[names(comm_bw["PCR-neg-CTRL-bac",][
     comm_bw["PCR-neg-CTRL-bac",]>0]),]

# ==================================================================
# ==================================================================





# ==================================================================
# 4. Clean the community matrix
# ==================================================================


# Remove the neg control and similar sample ------------------------

 dim(comm_bw)

 # Remove negative control
 comm_sub <- comm_bw[-12,]

 # Remove sample too close to negative control
 comm_sub <- comm_sub[-8,]



# Remove very rare species -----------------------------------------

 # Remove ASVs that do not appear in any sample
 # These ASVs appear because of the webs
 comm_sub <- comm_sub[, colSums(comm_sub) > 0]
 # what is the dimension of the subset community data set?
 dim(comm_sub)


 # Remove ASVs that are excessively rare
 comm_sub <- comm_sub[, colSums(comm_sub) > 1]
 # what is the dimension of the subset community data set?
 dim(comm_sub)


 # subset metadata and taxonomy to match
 metadata_sub <- metadata[rownames(comm_sub),]
 metadata_sub$n_reads <- rowSums(comm_sub)
 taxo_sub <- taxo[colnames(comm_sub),]

 # descriptive stats for samples and ASVs
 # number of sequences per sample
 hist(rowSums(comm_sub))
 hist(log10(colSums(comm_sub)))



# PCA on the subset -----------------------------------------------

 # Inspect PCA for subset
 comm_sub_pca <- prcomp(decostand(comm_sub, "hellinger"))
 
 # Prepare the ordination results for plotting
 labs <- as.factor(metadata_sub$sample_env)
 score <- scores(comm_sub_pca)[, 1:2]
 score <- cbind(score, sample_env = metadata_sub$sample_env)
 score <- data.frame(score)
 score[,1] <- as.numeric(score[,1])
 score[,2]<- as.numeric(score[,2])
 score[,3] <- as.factor(score[,3])
 sp <- data.frame(comm_sub_pca$rotation[, 1:2])
 summ <- summary(comm_sub_pca)$importance[2, 1:2]
 
 # Plot the results
 ggplot() +
   geom_point(data = sp,
              aes(x = PC1,
                  y = PC2),
              shape = 3,
              size = 1) +
   geom_point(data = score,
              aes(x = PC1,
                  y = PC2,
                  color = sample_env,
                  shape = sample_env),
               size = 3) +
   stat_ellipse(data = score,
                geom = "polygon",
                aes(x = PC1,
                    y = PC2,
                    fill = sample_env,
                    color = sample_env),
                level = 0.95,
                alpha = 0.25,
                linetype = "dashed",
                linewidth = 1,
                show.legend = FALSE) +
   scale_x_continuous(breaks = seq(-1.5, 1.5, 0.5),
                      limits = c(-1.5, 1.6)) +
   scale_y_continuous(breaks = seq(-1.5, 1.5, 0.5),
                      limits = c(-1.5, 1.5)) +
   scale_fill_manual(values = c("#E69F00",
                                "#666666")) +
   scale_color_manual(name = "Environment :",
                      values = c("#E69F00",
                                 "#666666")) +
   labs(color = "Environment :",
        shape = "Environment :") +
   xlab("\nPC1 (42.0%)") + ylab("PC2 (28.6%)\n") +
   theme_bw() + 
   theme(legend.position = "top",
         panel.grid = element_blank())

# ==================================================================
# ==================================================================





# ==================================================================
# 6. Save outputs
# ==================================================================

 path <- file.path(folder,
                   "env-data-clean")
 
 # Save clean taxa table
 saveRDS(
  taxo_sub,
  file = file.path(path, "env-bac-taxa-bw.rds")
 )

 # Save clean community table
 saveRDS(
  comm_sub,
  file = file.path(path, "env-bac-comm-bw.rds")
 )

 # Save clean metadata
 saveRDS(
  metadata_sub,
  file = file.path(path, "env-bac-metadata-bw.rds")
 )

# ==================================================================
# ==================================================================