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
 folder <- "./env-folder/env-data"
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


 # Delete unclassified ASVs at the phylum level
 #table(taxo[, "phylum"])
 #taxo <- subset(taxo,
 #               taxo[, "phylum"]!= "unclassified_Bacteria")


 # Delete unclassified ASVs at the class level
 #table(taxo[, "class"])
 #vec <- as.character(taxo[,"class"])
 
 #taxo <- subset(taxo,
 #               taxo[, "class"] %in% unique(
 #                   grep("unclassified_",
 #                        vec,
 #                        invert = TRUE,
 #                        value = TRUE)))
 
 # Delete vec object
 #rm(vec)


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

# Delete samples with too few reads
 comm_bw <- comm_bw[rowSums(comm_bw)>300, ]
 metadata <- metadata[rownames(comm_bw),]
 metadata$n_reads <- rowSums(comm_bw)

 # visualize log10 number of reads per sample
 hist(rowSums(comm_bw))
 hist(log10(rowSums(comm_bw)))
 # Most samples have around 35K to 100K reads

 # log10 of number of reads per ASV
 hist(log10(colSums(comm_bw)))



# Inspect rarefaction curves for spider samples --------------------

 # Prepare ploting options
 col <- c("black", "darkred", "forestgreen", "orange",
          "blue", "darkblue", "hotpink")
 lty <- c("solid", "dashed", "longdash", "dotdash")
 pars <- expand.grid(col = col, lty = lty,
                     stringsAsFactors = FALSE)
 samples <- data.frame(sample = as.factor(rownames(comm_bw)))


 # Plot rarefaction curve over all samples (NO zoom)
 with(#pars[1:12,],
      pars[1:13,],
    rarecurve(comm_bw, step = 200, #sample = raremax,
              label = TRUE, col = pars$col,
              lty = pars$lty, cex = 0.7))
 with(samples,
    legend("topright", legend = levels(sample),
           col = pars$col, lty = pars$lty,
           bty = "n"))


 # Plot rarefaction curve over all samples (WITH zoom)
 with(#pars[1:12,],
      pars[1:13,],
    rarecurve(comm_bw, step = 200, #sample = raremax,
              label = TRUE, col = pars$col,
              lty = pars$lty, cex = 0.7,
              xlim = c(0, 5000)))
 with(samples,
    legend("topright", legend = levels(sample),
           col = pars$col, lty = pars$lty,
           bty = "n"))


# Plot rarefaction curve over all samples (WITH zoom)
 with(#pars[1:12,],
      pars[1:13,],
    rarecurve(comm_bw, step = 200, #sample = raremax,
              label = TRUE, col = pars$col,
              lty = pars$lty, cex = 0.7,
              xlim = c(0, 1000)))
 with(samples,
    legend("topright", legend = levels(sample),
           col = pars$col, lty = pars$lty,
           bty = "n"))

# Most samples seem to reach a plateau around 200
# Except for DM-VN-112-bac and LO-VN-116-bac (too low)



# Visualize the community before rarefying -------------------------

 # PCA on Hellinger-transformed community data
 comm_pca <- prcomp(decostand(comm_bw, "hellinger"))
 
 # Extract sample information for plotting
 labs <- metadata$sample_env
 
 # plot ordination results
 ordiplot(comm_pca, type = "points",
          display = "species",
          xlim = c(-1, 1),
          ylim = c(-1, 1))
 color <- c("#E69F00", "#666666")
 score <- scores(comm_pca)[, 1:2]
 
 text(score, rownames(score),
      col = color, cex = 0.8)
 ordiellipse(comm_pca,
             labs,
             label = TRUE, cex = 0.8, font = 4)

# We see that the sample LO-VN-114-BAC is very similar to the negative control.


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
# 4. Subset community
# ==================================================================


# Remove low sequence number samples -------------------------------

 dim(comm_bw)

 # take subset of communities with at least 1528 sequences
 #comm_sub <- comm_bw[rowSums(comm_bw)>=1528,]

 # Remove negative control
 comm_sub <- comm_bw[-12,]

 # Remove sample too close to negative control
 comm_sub <- comm_sub[-8,]


 # Remove ASVs that do not appear in any sample
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
 labs <- as.factor(metadata_sub$sample_env)
 
 comm_sub_pca <- prcomp(decostand(comm_sub, "hellinger"))
 
 # Prepare the ordination results for plotting
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
                linewidth = 1,
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
   xlab("\nPC1 (14.9%)") + ylab("PC2 (59.2%)\n") +
   theme_bw() + theme(legend.position = "top",
                      panel.grid = element_blank())

# ==================================================================
# ==================================================================





# ==================================================================
# 5. Data normalization
# ==================================================================


# Inspect rarefaction curves ---------------------------------------
 
 # Prepare ploting options
 col <- c("black", "darkred", "forestgreen", "orange",
          "blue", "darkblue", "hotpink")
 lty <- c("solid", "dashed", "longdash", "dotdash")
 pars <- expand.grid(col = col, lty = lty,
                     stringsAsFactors = FALSE)
 samples <- data.frame(sample = as.factor(rownames(comm_sub)))

 # Plot rarefaction curves
  with(pars[1:11,],
     rarecurve(comm_sub, step = 200,
               label = FALSE, col = pars$col,
               lty = pars$lty, cex = 0.7,
               xlim = c(0, max(rowSums(comm_sub)))))
  with(samples,
     legend("topright", legend = levels(sample),
            col = pars$col, lty = pars$lty,
            bty = "n"))

  with(pars[1:11,],
     rarecurve(comm_sub, step = 200,
               label = FALSE, col = pars$col,
               lty = pars$lty, cex = 0.7,
               xlim = c(0, 50000)))
  with(samples,
     legend("topright", legend = levels(sample),
            col = pars$col, lty = pars$lty,
            bty = "n"))

  with(pars[1:11,],
     rarecurve(comm_sub, step = 200,
               label = FALSE, col = pars$col,
               lty = pars$lty, cex = 0.7,
               xlim = c(0, 2000)))
  with(samples,
     legend("topright", legend = levels(sample),
            col = pars$col, lty = pars$lty,
            bty = "n"))


# Apply rarefaction 1000X ------------------------------------------

 # This function increases the speed of vegan::rrarefy
 vegan_rrarefy <- function(x, sample) {
   x <- as.matrix(x)
   if (!identical(all.equal(x, round(x)), TRUE))
     stop("function is meaningful only for integers (counts)")
   if (!is.integer(x))
     x <- round(x)
   if (ncol(x) == 1)
     x <- t(x)
   if (length(sample) > 1 && length(sample) != nrow(x))
     stop(gettextf("length of 'sample' and number of rows of 'x' do not match"))
   if (any(rowSums(x) < sample))
     warning("some row sums < 'sample' and are not rarefied")
   out <- apply(x, 1, \(y) .Call(vegan:::do_rrarefy, y, sample))
   rownames(out) <- colnames(x)
   t(out)
 }

 # This function runs rarefaction multiple times
 rarefy_vegan_multiSeeds <- function(mat, n, seed){
   mat <- mat[rowSums(mat) >= n,]
   mat <- mat[,colSums(mat) > 0]
   vapply(`names<-`(seed, seed), \(s) {
     set.seed(s)
     vegan_rrarefy(mat, n)
   }, matrix(numeric(), nrow(mat), ncol(mat)))
 }


 # Rarefy and resample 1000X
 set.seed(123)
 comm_rarfy <- rarefy_vegan_multiSeeds(mat = comm_sub,
                                       n = 1536,
                                       seed = 1:1000)
 
 comm_rarfy[1:10,1:211,1]

 # Calculate the mean abundance over all 1000 runs
 #comm_rarfy <- round(apply(comm_rarfy, c(1, 2), mean))
 #comm_rarfy <- comm_rarfy[, colSums(comm_rarfy) > 0]
#
#
 ## Match ASV taxonomy to rarefied community
 #taxo_rarfy <- taxo_sub[colnames(comm_rarfy),]
#
 ## Match rarefied community to metadata
 #metadata_sub <- metadata[rownames(comm_rarfy),]
 #metadata_sub$n_reads_rarfy <- rowSums(comm_rarfy)



# Check the effect of rarefaction ---------------------------------

 richness_raw <- rowSums((comm_sub>0)*1)
 richness_rarfy <- rowSums((comm_rarfy>0)*1)
 
 plot(richness_rarfy ~ richness_raw,
      xlim = c(0, 120), ylim = c(0,120),
      pch = 16,
      xlab = "number of ASVs in raw data",
      ylab = "number of ASVs in rarefied data")
 abline(0:120, 1:120, lty = 2)
 
 plot(richness_rarfy[-6] ~ richness_raw[-6],
      xlim = c(0, 30), ylim = c(0,30),
      pch = 16,
      xlab = "number of ASVs in raw data",
      ylab = "number of ASVs in rarefied data")
 abline(0:30, 1:30, lty = 2)
 
# ==================================================================
# ==================================================================





# ==================================================================
# 6. Save outputs
# ==================================================================

 path <- file.path(getwd(),
                   "env-folder",
                   "env-data",
                   "env-data-clean")
 
 # Save clean taxa table
 saveRDS(
  taxo_rarfy,
  file = file.path(path, "env-bac-taxa-bw.rds")
 )

 # Save clean community table
 saveRDS(
  comm_rarfy,
  file = file.path(path, "env-bac-comm-bw.rds")
 )

 # Save clean metadata
 saveRDS(
  metadata_sub,
  file = file.path(path, "env-bac-metadata-bw.rds")
 )

# ==================================================================
# ==================================================================