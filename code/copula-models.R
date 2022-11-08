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
 meta_bw_diet <- readRDS(file.path(path2, "diet-bac-metadata-w.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Fit the models
# ==================================================================


# For spiders ------------------------------------------------------

 # fit marginal model
 fit1 <- stackedsdm(comm_bw_env, ~ n_reads,
                    # If parallelization is enabled
                    #do_parallel = TRUE,
                    #ncores = 48,
                    data = meta_bw_env,
                    family = "negative.binomial")
 
 # fit copula ordination 
 lvm1 <- cord(fit1, seed = 123)



# For webs ---------------------------------------------------------

 # fit marginal model
 fit2 <- stackedsdm(comm_bw_diet, ~ n_reads,
                    # If parallelization is enabled
                    #do_parallel = TRUE,
                    #ncores = 48,
                    data = meta_bw_diet,
                    family = "negative.binomial")
 
 # fit copula ordination 
 lvm2 <- cord(fit2, seed = 123)

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Save the outputs 
# ==================================================================

 # Setup the folder path
 path1 <- file.path(getwd(), "outputs")
 path2 <- 
 
 # Models for the env spiders
 saveRDS(fit1, file = file.path(path1, "env-bac-fit-bw.rds"))
 saveRDS(lvm1, file = file.path(path1, "env-bac-lvm-bw.rds"))

 
 # Models for the diet spiders
 saveRDS(fit2, file = file.path(path2, "diet-bac-fit-bw.rds"))
 saveRDS(lvm2, file = file.path(path2, "diet-bac-lvm-bw.rds"))

# ==================================================================
# ==================================================================