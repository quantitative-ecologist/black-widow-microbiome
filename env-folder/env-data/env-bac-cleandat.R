# ==================================================================

#        Synthetic data processing : part II - clean the data

# ==================================================================

# Code for sites in samples :
    # Urban samples :
    # UA = University of Arizona 
    # LO = Lowes

    # Desert samples :
    # DM = Dove mountain
    # CC = Chaos canyon

# W = web, VN = Black widow


# This script runs after prepdat to work


# ==================================================================
# 1. Import packages and data
# ==================================================================



# Import libraries -------------------------------------------------

 library(data.table)
 library(picante)



# Import data ------------------------------------------------------

 # Folder path
 folder <- "./env-folder/env-data"
 

  # Taxonomy data
 taxa_sp <- readRDS(
    file.path(folder, "env-data-raw",
              "env-bac-taxa-table2.rds"))
 # make sur order of ASVs match 
 taxo <- taxa_sp[colnames(comm),]
 rm(taxa_sp)

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
 table(taxo[, "phylum"])
 taxo <- subset(taxo,
                taxo[, "phylum"]!= "unclassified_Bacteria")


 # Delete unclassified ASVs at the class level
 table(taxo[, "class"])
 vec <- as.character(taxo[,"class"])
 
 taxo <- subset(taxo,
                taxo[, "class"] %in% unique(
                    grep("unclassified_",
                         vec,
                         invert = TRUE,
                         value = TRUE)))
 
 # Delete vec object
 rm(vec)


 # Apply changes to the community data
 comm <- comm[, rownames(taxo)]

# ==================================================================
# ==================================================================





# ==================================================================
# Summary statistics 
# ==================================================================


# Separate community
comm_w <- comm[c(4:6, 11:14, 19:23, 27:29),]
comm_vn <- comm[c(1:3, 7:10, 15:18, 23, 24:26),]


# Number of reads per sample
rowSums(comm_vn)
# Delete samples with too few reads
comm_vn <- comm_vn[rowSums(comm_vn)>700, ]

# visualize log10 number of reads per sample
hist(rowSums(comm_vn))
hist(log10(rowSums(comm_vn)))

# Most samples have around 35K to 100K reads

# log10 of number of reads per ASV
hist(log10(colSums(comm_vn)))




# Rarefaction on spider data ---------------------------------------

# Prepare ploting options
 col <- c("black", "darkred", "forestgreen", "orange",
          "blue", "darkblue", "hotpink")
 lty <- c("solid", "dashed", "longdash", "dotdash")
 pars <- expand.grid(col = col, lty = lty,
                     stringsAsFactors = FALSE)
 samples <- data.frame(sample = as.factor(rownames(comm_vn)))


 # Plot rarefaction curve over all samples (NO zoom)
 with(pars[1:12,],
    rarecurve(comm_vn, step = 200, #sample = raremax,
              label = TRUE, col = pars$col,
              lty = pars$lty, cex = 0.7))
 with(samples,
    legend("topright", legend = levels(sample),
           col = pars$col, lty = pars$lty,
           bty = "n"))

 # Plot rarefaction curve over all samples (WITH zoom)
 with(pars[1:12,],
    rarecurve(comm_vn, step = 200, #sample = raremax,
              label = TRUE, col = pars$col,
              lty = pars$lty, cex = 0.7,
              xlim = c(0, 50000)))
 with(samples,
    legend("topright", legend = levels(sample),
           col = pars$col, lty = pars$lty,
           bty = "n"))
# All the samples seem to reach a plateau around 5000
# Except for DM-VN-112-bac and LO-VN-116-bac (too low)



# Visualize community before rarefying ----------------------------

samples <- metadata[metadata$n_reads > 700,]
labs <- samples[samples$sample_type %in%
                c("spider", "control"),]$sample_env

# PCA on Hellinger-transformed community data
comm_pca <- prcomp(decostand(comm_vn, "hellinger"))

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
            samples[samples$sample_type %in%
                    c("spider", "control"),]$sample_env,
            label = TRUE, cex = 0.8, font = 4)

# Number of sequences per sample mapped onto ordination axes
ordisurf(comm_pca, rowSums(comm_vn),
         bubble = TRUE, cex = 2,
         main = "Library size (sequences/sample)")



# Check negative controls -----------------------------------------

# Abundance of ASVs in negative control
comm_vn["PCR-neg-CTRL-bac",][
    comm_vn["PCR-neg-CTRL-bac",]>0]
# Some problems here. Some ASVs are above 100
# See ASV_172, ASV_354, and ASV_399

# Check the taxonomic identity of ASVs present in negative control
taxo[names(comm_vn["PCR-neg-CTRL-bac",][
    comm_vn["PCR-neg-CTRL-bac",]>0]),]



# Remove low sequence number samples ------------------------------

dim(comm_vn)

# take subset of communities with at least 4000 sequences
comm.sub <- comm_vn[rowSums(comm_vn)>=1528,]
# also take subset of ASVs present in the remaining samples
comm.sub <- comm.sub[,apply(comm.sub,2,sum)>0]
# what is the dimension of the subset community data set?
dim(comm.sub)

# subset metadata and taxonomy to match
metadata.sub <- metadata[rownames(comm.sub),]
taxo.sub <- taxo[colnames(comm.sub),]

# descriptive stats for samples and ASVs
# number of sequences per sample
hist(rowSums(comm.sub))
hist(log10(colSums(comm.sub)))



# Inspect PCA for subset
samples <- metadata[metadata$n_reads >= 1528,]
labs <- samples[samples$sample_type %in%
                "spider",]$sample_env

comm.sub.pca <- prcomp(decostand(comm.sub,"hellinger"))

# plot ordination results
ordiplot(comm.sub.pca, type = "points",
         display = "species",
         xlim = c(-1, 1),
         ylim = c(-1, 1))
color <- c("#E69F00", "#666666")
score <- scores(comm.sub.pca)[, 1:2]

text(score, rownames(score),
     col = color, cex = 0.8)
ordiellipse(comm.sub.pca,
            samples[samples$sample_type %in%
                    "spider",]$sample_env,
            label = TRUE, cex = 0.8, font = 4)

# ==================================================================
# ==================================================================





# ==================================================================
# Data normalization
# ==================================================================

min(rowSums(comm.sub))

rarecurve(comm.sub, step=200, label=FALSE, xlim=c(0,8000))

# Apply rarefaction ------------------------------------------------

set.seed(0)
# Randomly rarefy samples
comm_rarfy <- rrarefy(comm.sub, sample = min(rowSums(comm.sub)))
# Remove any ASVs whose abundance is 0 after rarefaction
comm_rarfy <- comm_rarfy[,colSums(comm_rarfy)>1]
# Match ASV taxonomy to rarefied community
taxo_rarfy <- taxo.sub[colnames(comm_rarfy),]



# Check the effect of rarefaction ---------------------------------

richness_raw <- rowSums((comm.sub>0)*1)
richness_rarfy <- rowSums((comm_rarfy>0)*1)
plot(richness_rarfy~richness_raw,
     xlab = "number of ASVs in raw data",
     ylab = "number of ASVs in rarefied data")

# ==================================================================
# ==================================================================