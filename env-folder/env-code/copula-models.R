# ==================================================================

#             Compare communities using copula model

# ==================================================================

# Notes CalculCan : 100 trop bas, 500 overhead, test 200



# ==================================================================
# 1. Import libraries and data
# ==================================================================
 

# Import libraries -------------------------------------------------
 
 # For ncores
 options(mc.cores = parallel::detectCores())
 library(parallel)
 
 # For analysis
 library(ecoCopula)



# Import data ------------------------------------------------------
 
 path <- file.path(getwd(), "env-data-clean")

 # Community matrices
 comm_bw <- readRDS(file.path(path, "env-bac-comm-bw.rds"))
 comm_w <- readRDS(file.path(path, "env-bac-comm-w.rds"))
 
 # Metadata
 meta_bw <- readRDS(file.path(path, "env-bac-metadata-bw.rds"))
 meta_w <- readRDS(file.path(path, "env-bac-metadata-w.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Fit the models
# ==================================================================


# For spiders ------------------------------------------------------

 # fit marginal model
 fit_bw <- stackedsdm(comm_bw, ~ 1,
                      data = meta_bw,
                      family = "negative.binomial",
                      do_parallel = TRUE,
                      ncores = 48)
 
 # fit copula ordination 
 bw_lv <- cord(fit_bw, seed = 123)

 # fit graphical model 
 bw_gr <- cgr(fit_bw, seed = 3)
 
 # Save the outputs
 saveRDS(fit_bw, file = "stackedsdm_bw.rds")
 saveRDS(bw_lv, file = "lvm_bw.rds")
 saveRDS(bw_gr, file = "grm_bw.rds")



# For webs ---------------------------------------------------------

 # fit marginal model
 fit_w <- stackedsdm(comm_w, ~ 1,
                     data = data.frame(meta_w$sample_env),
                     family = "negative.binomial",
                     do_parallel = TRUE,
                     ncores = 48)
 
 # fit copula ordination 
 w_lv <- cord(fit_w, n.samp = 200, seed = 123)

 # fit graphical model 
 #w_gr <- cgr(fit_w, n.samp = 200, seed = 3)

 # Save the outputs
 saveRDS(fit_w, file = "stackedsdm_w.rds")
 saveRDS(w_lv, file = "lvm_w.rds")
 #saveRDS(w_gr, file = "grm_w.rds")

# ==================================================================
# ==================================================================