# ==================================================================

#         Compare communities using a latent variable model

# ==================================================================





# ==================================================================
# 1. Import libraries and data
# ==================================================================
 

# Import libraries -------------------------------------------------
 
 # For ncores
 options(mc.cores = parallel::detectCores())
 library(parallel)
 
 # For analysis
 library(boral)



# Import data ------------------------------------------------------
 
 # On Cedar
 path <- file.path(getwd(), "env-data-clean")

 # Community matrices
 #comm_bw <- readRDS(file.path(path, "env-bac-comm-bw.rds"))
 comm_w <- readRDS(file.path(path, "env-bac-comm-w.rds"))
 
 # Metadata
 #meta_bw <- readRDS(file.path(path, "env-bac-metadata-bw.rds"))
 meta_w <- readRDS(file.path(path, "env-bac-metadata-w.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# Fit the boral latent variable model
# ==================================================================

 # MCMC control options
 mcmc_control <- list(
     n.burnin = 500,
     n.iteration = 2000,
     n.thin = 2)
 
 # Fit the model
 lvm_w <- boral(
     y = comm_w,
     family = "negative.binomial",
     lv.control = list(num.lv = 2),
     row.eff = "fixed",
     mcmc.control = mcmc_control)

 # Save the output
 saveRDS(lvm_w, file = "lvm-boral-w.rds")

# ==================================================================
# ==================================================================