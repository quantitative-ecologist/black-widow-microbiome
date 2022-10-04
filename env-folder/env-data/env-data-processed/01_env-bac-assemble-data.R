# ==================================================================

#          Synthetic data processing with phyloseq : part I

# ==================================================================

# Code for sites in samples :
    # UA = University of Arizona,
    # DM = Dove mountain
    # CC = Chaos canyon,
    # LO = Lowes

# W = web, VN = Black widow





# ==================================================================
# 1. Import packages and data
# ==================================================================



# Import libraries -------------------------------------------------

 library(phyloseq)
 library(data.table)


# Import data ------------------------------------------------------

 # Folder path
 folder <- "./env-folder/env-data"
 

 # .rds objects
 seqtab.nochim <- readRDS(
    file.path(folder, "env-data-raw",
              "env-bac-seqtab-nochim-raw.rds"))
 
 taxid <- readRDS(
    file.path(folder, "env-data-raw",
              "env-bac-taxid.rds"))
 
 track_tab <- readRDS(
    file.path(folder, "env-data-raw",
              "env-bac-track-reads.rds"))
 
 setnames(track_tab, "nonchim", "nonchim_reads")

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Build the simplified phyloseq table
# ==================================================================


# Create data information for sample_data --------------------------

infos <- data.frame(
    sample_id = rownames(seqtab.nochim),
    sample_no = paste0("S", seq(1, 29)),
    sample_type = grepl("VN",
                        rownames(seqtab.nochim),
                        fixed = TRUE),
    sample_env = grepl(paste(c("DM", "CC"),
                             collapse = "|"),
                       rownames(seqtab.nochim)),
    nonchim_reads = track_tab$nonchim_reads)

infos$sample_type <- ifelse(infos$sample_type == TRUE,
                            "spider",
                            "web")
infos$sample_env <- ifelse(infos$sample_env == TRUE,
                           "desert",
                           "urban")
infos[23, 3] <- "control"
infos[23, 4] <- "control"

rownames(infos) <- rownames(seqtab.nochim)



# Build the simplified table ---------------------------------------
 
 ps <- phyloseq(otu_table(t(seqtab.nochim),
                          taxa_are_rows = TRUE),
                sample_data(infos),
                tax_table(as.matrix(taxid)))
 
 #ps1 <- phyloseq(otu_table(seqtab.nochim,
 #                         taxa_are_rows = FALSE),
 #                sample_data(infos),
 #                tax_table(as.matrix(taxid)))


 # Attribute short ASV names to the table
 taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Save the data
# ==================================================================



# Save the R object table as .rds ----------------------------------

# Save the phyloseq object to work with it
 saveRDS(ps, file.path(folder,
                       "env-data-processed",
                       "01_env-bac-phyloseq-data.rds"))



# Save the tables as .csv ------------------------------------------
 
 # Save the modified taxonomy table (Short Name)
 write.csv(as.data.frame(as(tax_table(ps), "matrix")),
           file = file.path(folder,
                            "env-data-processed", 
                            "01_env-bac-ASVTax-SN.csv"))
 
 # Save the modified transposed ASV matrix (Short Name)
 write.csv(as.data.frame(as(otu_table(ps), "matrix")),
           file = file.path(folder,
                            "env-data-processed", 
                            "01_env-bac-ASVMatrix-t-SN.csv"))

# ==================================================================
# ==================================================================