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
 path1 <- file.path(getwd(), "data", "env-data-clean")
 path2 <- file.path(getwd(), "data", "diet-data-clean")

 # Community matrices
 comm_bw_env <- readRDS(file.path(path1, "env-bac-comm-bw.rds"))
 comm_bw_diet <- readRDS(file.path(path2, "diet-bac-comm-bw.rds"))
 
 # Metadata
 meta_bw_env <- readRDS(file.path(path1, "env-bac-metadata-bw.rds"))
 meta_bw_diet <- readRDS(file.path(path2, "diet-bac-metadata-bw.rds"))

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

 # fit marginal model
 fit2a <- stackedsdm(
    y = comm_bw_diet,
    formula_X = ~ n_reads,
    # If parallelization is enabled
    #do_parallel = TRUE,
    #ncores = 48,
    data = meta_bw_diet,
    family = "negative.binomial"
 )

 fit2b <- stackedsdm(
    y = comm_bw_diet,
    formula_X = ~ n_reads + spider_weight,
    # If parallelization is enabled
    #do_parallel = TRUE,
    #ncores = 48,
    data = meta_bw_diet,
    family = "negative.binomial"
 )
 
 # Inspect the residuals fit1
 plot(
    residuals(fit2a) ~ fitted(fit2a),
    xlab = "Fitted values", 
    ylab = "Dunn-Smyth residuals"
 )
 abline(h = 0, col = "red")

 # Inspect the residuals fit2
 plot(
    residuals(fit2b) ~ fitted(fit2b),
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
 lvm2a <- cord(fit2a, seed = 123)
 lvm2b <- cord(fit2b, seed = 123)

# ==================================================================
# ==================================================================




# ==================================================================
# 4. Save the outputs 
# ==================================================================

 # Setup the folder path
 path <- file.path(
   getwd(),
   "outputs",
   "fits"
 )
 
 # Models for the env spiders
 saveRDS(fit1, file = file.path(path, "env-bac-fit-bw.rds"))
 saveRDS(lvm1, file = file.path(path, "env-bac-lvm-bw.rds"))

 
 # Models for the diet spiders
 saveRDS(fit2a, file = file.path(path, "diet-bac-fit1-bw.rds"))
 saveRDS(lvm2a, file = file.path(path, "diet-bac-lvm1-bw.rds"))

 saveRDS(fit2b, file = file.path(path, "diet-bac-fit2-bw.rds"))
 saveRDS(lvm2b, file = file.path(path, "diet-bac-lvm2-bw.rds"))

# ==================================================================
# ==================================================================