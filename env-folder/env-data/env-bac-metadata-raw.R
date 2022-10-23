# ==================================================================

#          Synthetic data processing : part I - raw metadata

# ==================================================================

# Code for sites in samples :
    # UA = University of Arizona,
    # DM = Dove mountain
    # CC = Chaos canyon,
    # LO = Lowes

# W = web, VN = Black widow


# Eventually download the data directly from OSF through API *****


# ==================================================================
# 1. Import data
# ==================================================================

 # Folder path
 folder <- "./env-folder/env-data"
 


# Load raw data ----------------------------------------------------

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

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Assemble raw metadata
# ==================================================================


# Create metadata --------------------------------------------------
 
 # Create dataframe
 metadata <- data.frame(
     sample_id = rownames(comm),
     sample_no = paste0("S", seq(1, 29)),
     sample_type = grepl("VN",
                         rownames(comm),
                         fixed = TRUE),
     sample_env = grepl(paste(c("DM", "CC"),
                              collapse = "|"),
                        rownames(comm)),
     sample_site = c(rep("CC", 6), rep("DM", 8),
                     rep("LO", 8), "control",
                     rep("UA", 6)),
     n_reads = rowSums(comm))
 
 # Add information
 metadata$sample_type <- ifelse(metadata$sample_type == TRUE,
                                "spider",
                                "web")
 metadata$sample_env <- ifelse(metadata$sample_env == TRUE,
                               "desert",
                               "urban")
 metadata[23, 3] <- "control"
 metadata[23, 4] <- "control"
 
 rownames(metadata) <- rownames(comm)

 
 # Save the data
 write.csv(metadata,
           file = file.path(folder,
                            "env-data-raw",
                            "env-bac-metadata-raw.csv"))

# ==================================================================
# ==================================================================