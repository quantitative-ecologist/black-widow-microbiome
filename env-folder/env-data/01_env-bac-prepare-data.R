# ==================================================================

#          Synthetic data processing : part I - Prepare data

# ==================================================================

# Code for sites in samples :
    # UA = University of Arizona,
    # DM = Dove mountain
    # CC = Chaos canyon,
    # LO = Lowes

# W = web, VN = Black widow


# Eventually download the data directly from OSF through API *****


# ==================================================================
# 1. Import packages and data
# ==================================================================


# Import libraries -------------------------------------------------

 library(phyloseq)
 library(data.table)



# Import data ------------------------------------------------------

 # Folder path
 folder <- "./env-folder/env-data"
 

 # load raw data :

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

 # Number of reads that made it through the pipeline
 reads_tab <- readRDS(
    file.path(folder, "env-data-raw",
              "env-bac-reads-tab.rds"))
 
 setnames(reads_tab, "nonchim", "nonchim_reads")

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Assemble metadata
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
     nonchim_reads = reads_tab$nonchim_reads)
 
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
 

 # Remove reads tab
 rm(reads_tab)

 # Save the data
 write.csv(metadata,
           file = file.path(folder,
                            "env-data-clean",
                            "env-bac-metadata.csv"))

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Assemble raw taxa and community data
# ==================================================================

# Rename ASVs from sequences to ASV number -------------------------
 
 # Create dataframe
 ASV_seq_info <- data.frame(ASV_name = paste0("ASV_",
                                              1:dim(comm)[2]),
                            ASV_sequence = colnames(comm))
 
 # Apply new names to comm and taxo tables
 colnames(comm) <- ASV_seq_info$ASV_name
 rownames(taxo) <- ASV_seq_info$ASV_name



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

 # Apply changes to the community data
 comm <- comm[, rownames(taxo)]

# ==================================================================
# ==================================================================