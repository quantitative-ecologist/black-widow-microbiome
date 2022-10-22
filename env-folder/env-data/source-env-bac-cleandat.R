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





# ==================================================================
# 1. Import packages and data
# ==================================================================



# Import libraries -------------------------------------------------

 library(phyloseq)
 library(data.table)
 library(vegan)



# Import data ------------------------------------------------------

 # Folder path
 folder <- "./env-folder/env-data"
 

 # .rds objects
 tab <- readRDS(
    file.path(folder, "env-data-processed",
              "01_env-bac-phyloseq-data.rds"))

# ==================================================================
# ==================================================================






# ==================================================================
# 2. Prepare the table
# ==================================================================


# Table for spiders only -------------------------------------------
 
 # Remove samples with too few reads
 tab <- subset_samples(tab, !(sample_id %in% c("UA-VN-I-bac",
                                               "LO-VN-115-bac",
                                               "LO-VN-116-bac")))

 # Remove web samples
 tab <- subset_samples(tab, sample_type != "web")

 # Table as matrix
 sp <- data.frame(otu_table(tab))
 comm <- t(sp)

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Clean the data
# ==================================================================


# Remove taxa with 0 in any sample ---------------------------------
 
 # Delete columns with sum to 0
 comm <- comm[, colSums(comm) > 0]
 
 # Check the distribution of total reads per sample
 hist(rowSums(comm), breaks = 10)
 hist(log10(rowSums(comm)), breaks = 10)
 
 # Number of sequences per ASV
 hist(colSums(comm))
 hist(log10(colSums(comm)))
 


# Inspect rarefaction curves----------------------------------------
 
 # Prepare ploting options
 col <- c("black", "darkred", "forestgreen", "orange",
          "blue", "darkblue", "hotpink")
 lty <- c("solid", "dashed", "longdash", "dotdash")
 pars <- expand.grid(col = col, lty = lty,
                     stringsAsFactors = FALSE)
 samples <- data.frame(sample = as.factor(rownames(comm)))

 # Plot rarefaction curve over all samples
 with(pars[1:11,],
    rarecurve(comm, step = 200, #sample = raremax,
              label = TRUE, col = col, lty = lty, cex = 0.5))
 with(samples,
    legend("topright", legend = levels(sample),
           col = col, lty = lty, bty = "n"))
 
 # Check with different max values
 with(pars[1:11,],
     rarecurve(comm, step = 200, xlim = c(0, 50000),
               label = TRUE, col = col, lty = lty, cex = 0.5))
 with(samples,
     legend("topright", legend = levels(sample),
            col = col, lty = lty, bty = "n"))

# The cutoff should be around 18 059 (smallest sample)



# Cutoff samples and apply rarefaction------------------------------

# take subset of communities with at least 18059 samples
# THIS REMOVED ALMOST ALL MY URBAN SAMPLES
comm_sub <- comm[rowSums(comm) >= 18059,]

# also take subset of ASVs present in the remaining samples
# I DONT UNDERSTAND THIS STEP ******
comm_sub1 <- comm_sub[, apply(comm_sub, 2, sum) > 0]

# what is the dimension of the subset community data set?
dim(comm_sub)
dim(comm_sub1)

# PCA on Hellinger-transformed community data
comm.pca <- prcomp(decostand(comm, "hellinger"))
# plot ordination results
ordiplot(comm.pca, display = "sites", type = "text", cex = 0.5)


# ==================================================================
# ==================================================================




# Cette partie est un Test présentement.
# ==================================================================
# 10. Check for sequences in the negative controls
# ==================================================================
 
 # check l'étape du negatif controle
# remove negative control             
 track_tab <- track_tab[-(23), .(sample, nonchim)]


 # Subset negative control sequences
 neg <- subset(as.data.frame(t(seqtab.nochim)), `PCR-neg-CTRL-bac` != 0)
 

 # First, verify the distribution of sequences
 # Delete sequences from the negative control
 seqtab.nochim <- seqtab.nochim[-c(23),
                                -which(colnames(seqtab.nochim) %in%
                                       rownames(neg))
]

 saveRDS(seqtab.nochim,
         file = file.path(outputs, "env-bac-seqtab.nochim-clean.rds"))


 # Save the taxonomy file as a .csv
 write.csv(as.data.frame(taxid),
           file = file.path(outputs, "env-bac-ASVTax.csv"))
 

 # Save as a transposed ASV matrix to prevent loading bug (too many columns)
 write.csv(as.data.frame(t(seqtab.nochim)),
           file = file.path(outputs, "env-bac-ASVMatrix-t.csv"))

# ==================================================================
# ==================================================================