# ==================================================================

#                 Compare ASV richness with glms

# ==================================================================




# ==================================================================
# 1. Libraries and data
# ==================================================================


# Libraries --------------------------------------------------------
 
 # Detect multicores
 options(mc.cores = parallel::detectCores())
 
 # Load libraries
 library(parallel)
 library(brms)
 library(data.table)



# Load the data ----------------------------------------------------

 # Path on Cedar
 path1 <- file.path(getwd(), "env-data-clean")
 path2 <- file.path(getwd(), "diet-data-clean")

 # Community matrices
 comm_bw_env <- readRDS(file.path(path1, "env-bac-comm-bw.rds"))
 comm_w_env <- readRDS(file.path(path1, "env-bac-comm-w.rds"))
 comm_bw_diet <- readRDS(file.path(path2, "diet-bac-comm-bw.rds"))
 
 # Metadata
 meta_bw_env <- readRDS(file.path(path1, "env-bac-metadata-bw.rds"))
 meta_w_env <- readRDS(file.path(path1, "env-bac-metadata-w.rds"))
 meta_bw_diet <- readRDS(file.path(path2, "diet-bac-metadata-bw.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Assemble the env spider data
# ==================================================================


# Prepare the env spider data --------------------------------------

 # Data frame from community matrix
 env_bw <- t(comm_bw_env)
 env_bw <- data.frame(env_bw)
 colnames(env_bw) <- rownames(comm_bw_env)

 # Add ASV name
 env_bw$ASV_id <- rownames(env_bw)

 # Reorder
 env_bw <- env_bw[,c(14,1:13)]

 # Transform
 env_bw <- data.table(env_bw)
 env_bw <- melt(env_bw,
                id.vars = "ASV_id",
                variable.name = "sample_id",
                value.name = "abundance")



# Compute the columns of interest ----------------------------------

 # Sum of abundances per sample
 env_bw[, sp_present := ifelse(abundance > 0, 1, 0)]
 env_bw[, 
        sp_richness := length(unique(ASV_id)),
        by = .(sample_id, sp_present)]
 
 # Extract the results in a synthetic table
 data1 <- unique(env_bw[sp_present == 1,
                        .(sample_id, sp_richness)])
 
 # Merge with the metadata to add covariates
 data1 <- merge(
   data1,
   meta_bw_env[,c(1,4,6)],
   by = "sample_id")



# Transform the data -----------------------------------------------

 # Log transform species richness
 data1$log_sp_richness <- log(data1$sp_richness)
 
 # Environment as factor
 data1$sample_env <- as.factor(data1$sample_env)

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Assemble the env web data
# ==================================================================


# Prepare the env web data -----------------------------------------

 # Data frame from community matrix
 env_w <- t(comm_w_env)
 env_w <- data.frame(env_w)
 colnames(env_w) <- rownames(comm_w_env)

 # Add ASV name
 env_w$ASV_id <- rownames(env_w)

 # Reorder
 env_w <- env_w[,c(14,1:13)]

 # Transform
 env_w <- data.table(env_w)
 env_w <- melt(env_w,
                id.vars = "ASV_id",
                variable.name = "sample_id",
                value.name = "abundance")



# Compute the columns of interest ----------------------------------

 # Sum of abundances per sample
 env_w[, sp_present := ifelse(abundance > 0, 1, 0)]
 env_w[, 
        sp_richness := length(unique(ASV_id)),
        by = .(sample_id, sp_present)]
 
 # Extract the results in a synthetic table
 data2 <- unique(env_w[sp_present == 1,
                        .(sample_id, sp_richness)])
 
 # Merge with the metadata to add covariates
 data2 <- merge(
   data2,
   meta_w_env[,c(1,4,6)],
   by = "sample_id")



# Transform the data -----------------------------------------------

 # Log transform species richness
 data2$log_sp_richness <- log(data2$sp_richness)
 
 # Environment as factor
 data2$sample_env <- as.factor(data2$sample_env)

# ==================================================================
# ==================================================================





# ==================================================================
# 4. Assemble the diet spider data
# ==================================================================


# Prepare the data -------------------------------------------------

 # Data frame from community matrix
 diet_bw <- t(comm_bw_diet)
 diet_bw <- data.frame(diet_bw)
 colnames(diet_bw) <- rownames(comm_bw_diet)

 # Add ASV name
 diet_bw$ASV_id <- rownames(diet_bw)

 # Reorder
 diet_bw <- diet_bw[,c(13,1:12)]

 # Transform
 diet_bw <- data.table(diet_bw)
 diet_bw <- melt(diet_bw,
                id.vars = "ASV_id",
                variable.name = "sample_id",
                value.name = "abundance")



# Compute the columns of interest ----------------------------------

 # Sum of abundances per sample
 diet_bw[, sp_present := ifelse(abundance > 0, 1, 0)]
 diet_bw[, 
        sp_richness := length(unique(ASV_id)),
        by = .(sample_id, sp_present)]
 
 # Extract the results in a synthetic table
 data3 <- unique(diet_bw[sp_present == 1,
                        .(sample_id, sp_richness)])
 
 # Merge with the metadata to add covariates
 data3 <- merge(
   data3,
   meta_bw_diet[,c(1,3, 4, 5)],
   by = "sample_id")



# Transform the data -----------------------------------------------

 # Log transform species richness
 data3$log_sp_richness <- log(data3$sp_richness)
 
 # Environment as factor
 data3$diet_treatment <- as.factor(data3$diet_treatment)
# ==================================================================
# ==================================================================





# ==================================================================
# 3. Fit the models
# ==================================================================


# Setup the model formulas -----------------------------------------
 
 # Gaussian model for env spiders
 form1 <- bf(
  log_sp_richness ~
    1 +
    sample_env +
    n_reads,
 sigma ~ 
    1 +
    sample_env +
    n_reads
 ) +
 gaussian()
 

 # Gaussian model for env webs
 form2 <- bf(
  log_sp_richness ~
    1 + 
    sample_env,
  sigma ~ 
    1 + 
    sample_env
 ) +
 gaussian()


 # Gaussian model for diet spiders
 form3 <- bf(
  log_sp_richness ~ 
    1 + 
    diet_treatment + 
    spider_weight  + 
    n_reads,
  sigma ~ 
    1 + 
    diet_treatment + 
    spider_weight + 
    n_reads
 ) + 
 gaussian()



# Setup the priors -------------------------------------------------
 
 # Gaussian priors for b and intercepts
 priors <- c(
   set_prior("normal(0, 2)", 
             class = "Intercept"),
   set_prior("normal(0, 2)", 
             class = "Intercept",
             dpar = "sigma"),
   set_prior("normal(0, 2)",
             class = "b"),
   set_prior("normal(0, 2)",
             class = "b",
             dpar = "sigma")
 )



# Run the models ---------------------------------------------------
 
 # Model for env spiders
 model1 <- brm(form1,
              warmup = 2000, 
              iter = 22000,
              thin = 80,
              chains = 4,
              seed = 123,
              init = 0,
              prior = priors,
              threads = threading(12),
              backend = "cmdstanr",
              control = list(
                adapt_delta = 0.99,
                max_treedepth = 12),
              data = data1)
 
 # Save the model output
 saveRDS(model1, file = "env-bac-glm-bw.rds")


 # Model for env webs
 model2 <- brm(form2,
              warmup = 2000, 
              iter = 22000,
              thin = 80,
              chains = 4,
              seed = 123,
              init = 0,
              prior = priors,
              threads = threading(12),
              backend = "cmdstanr",
              control = list(
                adapt_delta = 0.99,
                max_treedepth = 12),
              data = data2)
 
 # Save the model output
 saveRDS(model2, file = "env-bac-glm-w.rds")


 # Model for diet spiders
 model3 <- brm(form3,
              warmup = 5000, 
              iter = 65000,
              thin = 240,
              chains = 4,
              seed = 123,
              init = 0,
              prior = priors,
              threads = threading(12),
              backend = "cmdstanr",
              control = list(
                adapt_delta = 0.999,
                max_treedepth = 12),
              data = data3)
 
 # Save the model output
 saveRDS(model3, file = "diet-bac-glm-bw.rds")

# ==================================================================
# ==================================================================
