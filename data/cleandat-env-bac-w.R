# ==================================================================

#                   Clean the web data for analysis

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
              "seqtabnochim-env-bac.rds"))

  # Taxonomy data
 taxa_sp <- readRDS(
    file.path(folder, "data-raw-env",
              "taxa-table2-env-bac.rds"))
 # make sur order of ASVs match 
 taxo <- taxa_sp[colnames(comm),]
 rm(taxa_sp)

  # Raw metadata
 metadata <- read.csv(
    file.path(folder, "data-raw-env",
              "metadata-raw-env-bac.csv"),
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
 
 # Inspect for ASVs other than bacteria in "domain"
 table(taxo[,"domain"]) # 3 eukaryota
 
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

 # Separate community to have only web samples
 comm_w <- comm[c(4:6, 11:14, 19:22, 23, 27:29), ]

 # Number of reads per sample
 rowSums(comm_w)

 # Apply new read count to metadata
 metadata$n_reads_bac <- rowSums(comm_w)

 # visualize log10 number of reads per sample
 hist(rowSums(comm_w))
 hist(log10(rowSums(comm_w)))
 # Most samples have around 35K to 100K reads

 # log10 of number of reads per ASV
 hist(log10(colSums(comm_w)))
 # distribution is not lognormal



# Inspect rarefaction curves for web samples -----------------------

 # Prepare ploting options
 col <- c("black", "darkred", "forestgreen", "orange",
          "blue", "darkblue", "hotpink")
 lty <- c("solid", "dashed", "longdash", "dotdash")
 pars <- expand.grid(col = col, lty = lty,
                     stringsAsFactors = FALSE)
 samples <- data.frame(sample = as.factor(rownames(comm_w)))


 # Plot rarefaction curve over all samples (NO zoom)
 with(pars[1:15,],
    rarecurve(comm_w, step = 200, #sample = raremax,
              label = TRUE, col = pars$col,
              lty = pars$lty, cex = 0.7))
 with(samples,
    legend("topright", legend = levels(sample),
           col = pars$col, lty = pars$lty,
           cex = 0.7,
           bty = "n"))


 # Plot rarefaction curve over all samples (WITH zoom)
 with(pars[1:15,],
    rarecurve(comm_w, step = 200, #sample = raremax,
              label = FALSE, col = pars$col,
              lty = pars$lty, cex = 0.7,
              xlim = c(0, 30000)))
 with(samples,
    legend("topright", legend = levels(sample),
           col = pars$col, lty = pars$lty,
           cex = 0.7,
           bty = "n"))

 # Plateau seems to be reached around 15 000 reads



# Visualize the community before rarefying -------------------------

 # PCA on Hellinger-transformed community data
 comm_pca <- prcomp(decostand(comm_w, "hellinger"))

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
   xlab("\nPC1 (13.8%)") + ylab("PC2 (13.0%)\n") +
   theme_bw() + 
   theme(legend.position = "top",
         panel.grid = element_blank())

 # LO-VN-114-BAC is very similar to the negative control.

 # Number of sequences per sample mapped onto ordination axes
 ordisurf(comm_pca, rowSums(comm_w),
          bubble = TRUE, cex = 2,
          main = "Library size (sequences/sample)")



# Check negative controls -----------------------------------------

 # Abundance of ASVs in negative control
 comm_w["PCR-neg-CTRL-bac",][
     comm_w["PCR-neg-CTRL-bac",]>0]
 # Some problems here. Some ASVs are above 100
 # See ASV_169, ASV_352, and ASV_397 are high

 # Check the taxonomic identity of ASVs present in negative control
 taxo[names(comm_w["PCR-neg-CTRL-bac",][
     comm_w["PCR-neg-CTRL-bac",]>0]),]

# ==================================================================
# ==================================================================





# ==================================================================
# 4. Subset community
# ==================================================================


# Remove the neg control and similar sample ------------------------

 dim(comm_w)

 # Remove negative control
 comm_sub <- comm_w[-12,]

 # Remove sample too close to negative control
 comm_sub <- comm_sub[-8,]



# Remove very rare species -----------------------------------------
 

 # Remove ASVs that do not appear in any sample
 # Check dim before
 dim(comm_w)
 # Remove ASVs
 comm_sub <- comm_sub[, colSums(comm_sub) > 0]
 # what is the dimension of the subset community data set?
 dim(comm_sub)


 # Remove ASVs that are too rare
 comm_sub <- comm_sub[, colSums(comm_sub) > 1]
 # what is the dimension of the subset community data set?
 dim(comm_sub)


 # subset metadata and taxonomy to match
 metadata_sub <- metadata[rownames(comm_sub),]
 taxo_sub <- taxo[colnames(comm_sub),]

 # descriptive stats for samples and ASVs
 # number of sequences per sample
 hist(rowSums(comm_sub), breaks = 20)
 hist(log10(colSums(comm_sub)), breaks = 20)



# PCA on the subset -----------------------------------------------

 # Inspect PCA for subset
 comm_sub_pca <- prcomp(decostand(comm_sub,"hellinger"))
 
 # Prepare ordination results for the plots
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
              size = 1,
              shape = 3) +
   geom_point(data = score,
              aes(x = PC1,
                  y = PC2,
                  color = sample_env),
               shape = 17,
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
                size = 1,
                show.legend = FALSE) +
   scale_x_continuous(breaks = seq(-1, 1, 0.5),
                      limits = c(-1.2, 1.2)) +
   scale_y_continuous(breaks = seq(-1, 1, 0.5),
                      limits = c(-1.2, 1.2)) +
   scale_fill_manual(values = c("#E69F00",
                                 "#666666")) +
   scale_color_manual(name = "Environment :",
                      values = c("#E69F00",
                                 "#666666")) +
   xlab("\nPC1 (14.4%)") + ylab("PC2 (13.5%)\n") +
   theme_bw() + theme(legend.position = "top",
                      panel.grid = element_blank())

# ==================================================================
# ==================================================================





# ==================================================================
# 5. Data normalization
# ==================================================================


# Apply rarefaction ------------------------------------------------

 set.seed(0)
 # Randomly rarefy samples
 comm_rarfy <- rrarefy(comm_sub, sample = min(rowSums(comm_sub)))
 dim(comm_rarfy)

 # Remove any ASVs whose abundance is 0 after rarefaction
 comm_rarfy <- comm_rarfy[, colSums(comm_rarfy) > 0]
 dim(comm_rarfy)

 # Remove species that are too rare
 comm_rarfy <- comm_rarfy[, colSums(comm_rarfy) > 1]
 dim(comm_rarfy)

 # Match ASV taxonomy to rarefied community
 taxo_rarfy <- taxo_sub[colnames(comm_rarfy),]
 
 # Match rarefied community to metadata
 metadata_sub <- metadata[rownames(comm_rarfy),]
 metadata_sub$n_reads_bac_rarfy <- rowSums(comm_rarfy)



# Check the effect of rarefaction ---------------------------------

 richness_raw <- rowSums((comm_sub>0)*1)
 richness_rarfy <- rowSums((comm_rarfy>0)*1)
 
 plot(richness_rarfy ~ richness_raw,
      xlim = c(0, 1000), ylim = c(0,1000),
      pch = 16,
      xlab = "number of ASVs in raw data",
      ylab = "number of ASVs in rarefied data")
 abline(0:1000, 1:1000, lty = 2)
# ==================================================================
# ==================================================================





# ==================================================================
# 6. Save outputs
# ==================================================================

 path <- file.path(folder,
                   "data-clean-env")
 
 # Save clean taxa table
 saveRDS(
    taxo_rarfy,
    file = file.path(path, "taxa-env-bac-w.rds")
 )

 # Save clean community table
 saveRDS(
    comm_rarfy,
    file = file.path(path, "comm-env-bac-w.rds")
 )

 # Save clean metadata
 saveRDS(
    metadata_sub,
    file = file.path(path, "metadata-env-bac-w.rds")
 )

# ==================================================================
# ==================================================================