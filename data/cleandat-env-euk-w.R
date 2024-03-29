# ==================================================================

#          Clean the spider data for analysis of eukaryotes

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
    file.path(folder, "data-raw-env",
              "seqtabnochim-env-euk.rds"))

  # Taxonomy data
 taxa_sp <- readRDS(
    file.path(folder, "data-raw-env",
              "taxa-table-env-euk.rds"))
 # make sur order of ASVs match 
 taxo <- taxa_sp[colnames(comm),]
 rm(taxa_sp)

  # Raw metadata
 metadata <- read.csv(
    file.path(folder, "data-raw-env",
              "metadata-raw-env.csv"),
              row.names = 1)
 # only spider samples
 metadata <- metadata[metadata$sample_type %in%
                      c("web", "control"),]

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
 
 # Print column names
 colnames(taxo)

 # Inspect for ASVs other than eukaryota in "kingdom"
 table(taxo[,"kingdom"]) # only eukaryota
 
 # Inspect the phylum groups
 table(taxo[,"phylum"])
 # Remove unclassified eukaryotas
 taxo <- subset(taxo, taxo[, "phylum"] != "unclassified_Eukaryota")
 

 # Delete the Black Widows DNA
 table(taxo[,"genus"])
 taxo <- subset(taxo, taxo[, "genus"] != "Latrodectus")

 # Apply changes to the community data
 comm <- comm[, rownames(taxo)]

 # Test only with metazoans
 taxo_meta <- subset(taxo, taxo[, "phylum"] == "Metazoa")
 comm_meta <- comm[, rownames(taxo_meta)]
 comm_meta <- comm_meta[, colSums(comm_meta) > 0]
# ==================================================================
# ==================================================================





# ==================================================================
# 3. Data exploration
# ==================================================================



# Summary statistics -----------------------------------------------

 # Separate community to have only spider samples
 comm_w <- comm[c(4:6, 11:14, 19:22, 23, 27:29), ]
 comm_w_meta <- comm_meta[c(4:6, 11:14, 19:22, 23, 27:29), ]

 # Number of reads per sample
 rowSums(comm_w)
 

 # Adjust the new number of reads to the metadata
 metadata$n_reads_euk <- rowSums(comm_w)


 # visualize log10 number of reads per sample
 hist(rowSums(comm_w))
 hist(log10(rowSums(comm_w)))

 # log10 of number of reads per ASV
 hist(log10(colSums(comm_w)))



# Visualize the community ------------------------------------------

 # PCA on Hellinger-transformed community data
 comm_pca <- prcomp(decostand(comm_w, "hellinger"))
 #comm_pca <- prcomp(decostand(comm_w_meta, "hellinger"))

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
   xlab("\nPC1 (23.5%)") + ylab("PC2 (18.4%)\n") +
   theme_bw() + 
   theme(legend.position = "top",
         panel.grid = element_blank())


 # Number of sequences per sample mapped onto ordination axes
 ordisurf(comm_pca, rowSums(comm_w),
          bubble = TRUE, cex = 2,
          main = "Library size (sequences/sample)")



# Check negative controls -----------------------------------------

 # Abundance of ASVs in negative control
 comm_w["PCR-neg-CTRL-euc",][
     comm_w["PCR-neg-CTRL-euc",]>0]
 # Some problems here. Some ASVs are above 100
 # See ASV_44, ASV_47

 # Check the taxonomic identity of ASVs present in negative control
 taxo[names(comm_w["PCR-neg-CTRL-euc",][
     comm_w["PCR-neg-CTRL-euc",]>0]),]
 
 # It is mostly Sordariomycetes family

# ==================================================================
# ==================================================================





# ==================================================================
# 4. Clean the community matrix
# ==================================================================


# Remove the neg control and similar sample ------------------------
 
 # Check the dimension
 dim(comm_w)

 # Remove negative control
 comm_sub <- comm_w[-12,]



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
 metadata_sub <- metadata[-c(12),]
 metadata_sub$n_reads_euk <- rowSums(comm_sub)
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
 score$sample_id <- as.factor(rownames(score))
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
   geom_text(data = score,
             aes(x = PC1,
                 y = PC2,
                 label = sample_id,
                 color = sample_env),
             size = 3) +
   scale_x_continuous(breaks = seq(-1.5, 1.5, 0.5),
                      limits = c(-1.5, 1.5)) +
   scale_y_continuous(breaks = seq(-1.5, 1.5, 0.5),
                      limits = c(-1.5, 1.5)) +
   scale_fill_manual(values = c("#E69F00",
                                "#666666")) +
   scale_color_manual(name = "Environment :",
                      values = c("#E69F00",
                                 "#666666")) +
   labs(color = "Environment :",
        shape = "Environment :") +
   xlab("\nPC1 (27.3%)") + ylab("PC2 (19.8%)\n") +
   theme_bw() + 
   theme(legend.position = "top",
         panel.grid = element_blank())

# ==================================================================
# ==================================================================





# ==================================================================
# 6. Save outputs
# ==================================================================
 
 # Setup folder path
 path <- file.path(folder,
                   "data-clean-env")
 
 
 # Save clean taxa table
 saveRDS(
  taxo_sub,
  file = file.path(path, "taxa-env-euk-w.rds")
 )

 # Save clean community table
 saveRDS(
  comm_sub,
  file = file.path(path, "comm-env-euk-w.rds")
 )

 # Save clean metadata
 saveRDS(
  metadata_sub,
  file = file.path(path, "metadata-env-euk-w.rds")
 )

# ==================================================================
# ==================================================================