# ==================================================================

#          Clean the diet treated spider data for analysis

# ==================================================================





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
    file.path(folder, "data-raw-diet",
              "seqtabnochim-diet-bac.rds"))

  # Taxonomy data
 taxa_sp <- readRDS(
    file.path(folder, "data-raw-diet",
              "taxa-table2-diet-bac.rds"))
 # make sur order of ASVs match 
 taxo <- taxa_sp[colnames(comm),]
 rm(taxa_sp)

  # Raw metadata
 metadata <- read.csv(
    file.path(folder, "data-raw-diet",
              "metadata-raw-diet-bac.csv"),
              row.names = 1)

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
 table(taxo[,"domain"])

 
 # Inspect if there are chloroplasts or mitochondria
 table(taxo[,"order"])["Chloroplast"]
 table(taxo[,"family"])["Mitochondria"]


 # Delete unclassified ASVs at the phylum level
 table(taxo[, "phylum"])


 # Apply changes to the community data
 dim(comm)
 comm <- comm[, rownames(taxo)]
 dim(comm)
 # No changes were made

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Data exploration
# ==================================================================


# Summary statistics -----------------------------------------------

 # Number of reads per sample
 rowSums(comm)

 # visualize log10 number of reads per sample
 hist(rowSums(comm))
 hist(log10(rowSums(comm)))
 # Most samples have around 35K to 100K reads

 # log10 of number of reads per ASV
 hist(log10(colSums(comm)))
 # Distribution is lognormal



# Check negative controls -----------------------------------------

 # Abundance of ASVs in negative control
 comm["CTRL-PCR-neg-bac",][
     comm["CTRL-PCR-neg-bac",]>0]

 # Check the taxonomic identity of ASVs present in negative control
 taxo[names(comm["CTRL-PCR-neg-bac",][
     comm["CTRL-PCR-neg-bac",]>0]),]

 # There are no reads nor taxa in the negative control

 # Remove the negative control
 comm <- comm[-13,]



# Visualize the community ------------------------------------------
 
 # PCA on Hellinger-transformed community data
 comm_pca <- prcomp(decostand(comm, "hellinger"))
 
 # Prepare the ordination results for plotting
 labs <- as.factor(metadata[-13,]$diet_treatment)
 score <- scores(comm_pca)[, 1:2]
 score <- cbind(score, diet = metadata[-13,]$diet_treatment)
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
                 color = diet),
             size = 3) +
   scale_x_continuous(breaks = seq(-1, 1, 0.5),
                      limits = c(-1, 1)) +
   scale_y_continuous(breaks = seq(-1, 1, 0.5),
                      limits = c(-1, 1)) +
   scale_color_manual(values = c("black",
                                 "red",
                                 "blue")) +
   labs(color = "Diet :") +
   xlab("\nPC1 (25.0%)") + ylab("PC2 (17.2%)\n") +
   theme_bw() + 
   theme(legend.position = "top",
         panel.grid = element_blank())


 # Number of sequences per sample mapped onto ordination axes
 ordisurf(comm_pca, rowSums(comm),
          bubble = TRUE, cex = 2,
          main = "Library size (sequences/sample)")

# ==================================================================
# ==================================================================





# ==================================================================
# 4. Clean the community matrix
# ==================================================================


# Remove very rare species -----------------------------------------
 
 # Check dimensions
 dim(comm)

 # Remove ASVs that do not appear in any sample
 comm_sub <- comm[, colSums(comm) > 0]
 # what is the dimension of the subset community data set?
 dim(comm_sub)
 # No changes

 # Remove ASVs that are excessively rare
 comm_sub <- comm_sub[, colSums(comm_sub) > 1]
 # what is the dimension of the subset community data set?
 dim(comm_sub)
 # No changes


 # subset metadata and taxonomy to match
 metadata_sub <- metadata[rownames(comm_sub),]
 metadata_sub$n_reads <- rowSums(comm_sub)
 taxo_sub <- taxo[colnames(comm_sub),]
 # No changes in read count

 # descriptive stats for samples and ASVs
 # number of sequences per sample
 hist(rowSums(comm_sub))
 hist(log10(colSums(comm_sub)))
 
# ==================================================================
# ==================================================================





# ==================================================================
# 5. Save outputs
# ==================================================================
 
 # Setup folder path
 path <- file.path(
    folder,
    "data-clean-diet-bac"
 )
 

 # Save clean taxa table
 saveRDS(
  taxo_sub,
  file = file.path(path, "taxa-diet-bac.rds")
 )

 # Save clean community table
 saveRDS(
  comm_sub,
  file = file.path(path, "comm-diet-bac.rds")
 )

 # Save clean metadata
 saveRDS(
  metadata_sub,
  file = file.path(path, "metadata-diet-bac.rds")
 )

# ==================================================================
# ==================================================================