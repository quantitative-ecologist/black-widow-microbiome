# ==================================================================

#                Compare web communities using NMDS

# ==================================================================





# ==================================================================
# 1. Import libraries and data
# ==================================================================


# Load libraries ---------------------------------------------------
 
 # Personal computer
 library(picante)
 library(ggplot2)



# Import data ------------------------------------------------------
 
 path <- file.path(getwd(),
                   "data",
                   "data-clean-env")

 # Community matrices
 comm_w <- readRDS(file.path(path, "comm-env-bac-w.rds"))
 
 # Metadata
 meta_w <- readRDS(file.path(path, "metadata-env-bac-w.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Compute NMDSs on the web data
# ==================================================================


 # Create Hellinger-transformed version of the community data
 comm_w_hel <- decostand(comm_w, method = "hellinger")
 
 # NMDS ordination based on Bray-Curtis distances
 set.seed(123)
 nmds_w <- metaMDS(
  comm_w_hel,
  distance = "bray",
  trace = FALSE
 )
 
 # Save the output
 path1 <- file.path(getwd(), "outputs", "fits-bac")
 saveRDS(nmds_w, file = file.path(path1, "nmds-bac-w.rds"))

# ==================================================================
# ==================================================================