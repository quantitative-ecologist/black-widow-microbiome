# ==================================================================

#             Compare communities using copula models

# ==================================================================





# ==================================================================
# 1. Import libraries and data
# ==================================================================
 

# Import libraries -------------------------------------------------
 
 # For ncores (if there is parallel computing)
 #options(mc.cores = parallel::detectCores())
 #library(parallel)
 
 # For analysis
 library(ecoCopula)



# Import data ------------------------------------------------------

 # On my computer
 path1 <- file.path(getwd(), "data", "data-clean-env")

 # Community matrices
 comm_bw <- readRDS(file.path(path1, "comm-env-euk-bw.rds"))
 comm_w <- readRDS(file.path(path1, "comm-env-euk-w.rds"))
 
 # Metadata
 meta_bw <- readRDS(file.path(path1, "metadata-env-euk-bw.rds"))
 meta_w <- readRDS(file.path(path1, "metadata-env-euk-w.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Fit the models
# ==================================================================


# For env spiders --------------------------------------------------

 # fit marginal model
 fit1 <- stackedsdm(
    y = comm_bw,
    formula_X = ~ n_reads_euk,
    data = meta_bw,
    family = "negative.binomial"
 )
 
 # Inspect the residuals
 plot(
    residuals(fit1) ~ fitted(fit1),
    xlab = "Fitted values", 
    ylab = "Dunn-Smyth residuals"
 )
 abline(h = 0, col = "red")



# For env webs -----------------------------------------------------

 # fit marginal model
 fit2 <- stackedsdm(
    y = comm_w,
    formula_X = ~ n_reads_euk,
    # If parallelization is enabled
    #do_parallel = TRUE,
    #ncores = 48,
    data = meta_w,
    family = "negative.binomial"
 )
 
 # Inspect the residuals fit1
 plot(
    residuals(fit2) ~ fitted(fit2),
    xlab = "Fitted values", 
    ylab = "Dunn-Smyth residuals"
 )
 abline(h = 0, col = "red")

# ==================================================================
# ==================================================================




# ==================================================================
# 3. Fit the copula ordination
# ==================================================================


# For env spiders --------------------------------------------------
 
 # Fit
 lvm1 <- cord(fit1, seed = 123)



# For env webs -----------------------------------------------------
 
 # Fit
 lvm2 <- cord(fit2, seed = 123)

# ==================================================================
# ==================================================================




# ==================================================================
# 4. Save the outputs 
# ==================================================================

 # Setup the folder path
 path <- file.path(
   getwd(),
   "outputs",
   "fits-euk"
 )
 
 # Models for the env spiders
 saveRDS(fit1, file = file.path(path, "fit-env-euk-bw.rds"))
 saveRDS(lvm1, file = file.path(path, "lvm-env-euk-bw.rds"))

 
 # Models for the diet spiders
 saveRDS(fit2, file = file.path(path, "fit-env-euk-w.rds"))
 saveRDS(lvm2, file = file.path(path, "lvm-env-euk-w.rds"))

# ==================================================================
# ==================================================================