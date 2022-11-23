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
 path2 <- file.path(getwd(), "data", "data-clean-diet")

 # Community matrices
 comm_bw_env <- readRDS(file.path(path1, "comm-env-bac-bw.rds"))
 comm_bw_diet <- readRDS(file.path(path2, "comm-diet-bac.rds"))
 
 # Metadata
 meta_bw_env <- readRDS(file.path(path1, "metadata-env-bac-bw.rds"))
 meta_bw_diet <- readRDS(file.path(path2, "metadata-diet-bac.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Fit the models
# ==================================================================


# For env spiders --------------------------------------------------

 # fit marginal model
 fit1 <- stackedsdm(
    y = comm_bw_env,
    formula_X = ~ n_reads,
    # If parallelization is enabled
    #do_parallel = TRUE,
    #ncores = 48,
    data = meta_bw_env,
    family = "negative.binomial"
 )
 
 # Inspect the residuals
 plot(
    residuals(fit1) ~ fitted(fit1),
    xlab = "Fitted values", 
    ylab = "Dunn-Smyth residuals"
 )
 abline(h = 0, col = "red")



# For diet spiders -------------------------------------------------

 fit2 <- stackedsdm(
    y = comm_bw_diet,
    formula_X = ~ n_reads + spider_weight,
    # If parallelization is enabled
    #do_parallel = TRUE,
    #ncores = 48,
    data = meta_bw_diet,
    family = "negative.binomial"
 )

 # Inspect the residuals fit2
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



# For diet spiders -------------------------------------------------
 
 # Fit
 lvm2 <- cord(fit2, seed = 123)

# ==================================================================
# ==================================================================





# ==================================================================
# 5. Save the outputs 
# ==================================================================

 # Setup the folder path
 path <- file.path(
   getwd(),
   "outputs",
   "fits-bac"
 )
 
 # Models for the env spiders
 saveRDS(fit1, file = file.path(path, "fit-env-bac-bw.rds"))
 saveRDS(lvm1, file = file.path(path, "lvm-env-bac-bw.rds"))
 
 
 # Models for the diet spiders
 saveRDS(fit2, file = file.path(path, "fit-diet-bac-bw.rds"))
 saveRDS(lvm2, file = file.path(path, "lvm-diet-bac-bw.rds"))

# ==================================================================
# ==================================================================