# ==================================================================

#          Compare ASV richness for env spiders and webs

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
 path <- file.path(getwd(), "data-clean-env")

 # Community matrices
 comm_bw <- readRDS(file.path(path, "comm-env-euk-bw.rds"))
 comm_w <- readRDS(file.path(path, "comm-env-euk-w.rds"))
 
 # Metadata
 meta_bw <- readRDS(file.path(path, "metadata-env-euk-bw.rds"))
 meta_w <- readRDS(file.path(path, "metadata-env-euk-w.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Assemble the env spider data
# ==================================================================


# Prepare the env spider data --------------------------------------

 # Data frame from community matrix
 env_bw <- t(comm_bw)
 env_bw <- data.frame(env_bw)
 colnames(env_bw) <- as.character(
    strsplit(rownames(comm_bw), "-euc")
 )

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
   meta_bw[,c(1,4,7)],
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
 env_w <- t(comm_w)
 env_w <- data.frame(env_w)
 colnames(env_w) <- as.character(
    strsplit(rownames(comm_w), "-euc")
 )

 # Add ASV name
 env_w$ASV_id <- rownames(env_w)

 # Reorder
 env_w <- env_w[,c(15,1:14)]

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
   meta_w[,c(1,4,7)],
   by = "sample_id")



# Transform the data -----------------------------------------------

 # Log transform species richness
 data2$log_sp_richness <- log(data2$sp_richness)
 
 # Environment as factor
 data2$sample_env <- as.factor(data2$sample_env)

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
    n_reads_euk,
 sigma ~ 
    1 +
    sample_env +
    n_reads_euk
 ) +
 gaussian()
 

 # Gaussian model for env webs
 form2 <- bf(
  log_sp_richness ~
    1 + 
    sample_env +
    n_reads_euk,
  sigma ~ 
    1 + 
    sample_env +
    n_reads_euk
 ) +
 gaussian()



# Setup the priors -------------------------------------------------
 
 # Gaussian priors for b and intercepts
 priors <- c(
   set_prior("normal(2, 1)", 
             class = "Intercept"),
   set_prior("normal(1, 2)", 
             class = "Intercept",
             dpar = "sigma"),
   set_prior("normal(0, 1)",
             class = "b"),
   set_prior("normal(0, 2)",
             class = "b",
             dpar = "sigma")
 )



# Run the models ---------------------------------------------------
 
 # Model for env spiders
 model1 <- brm(
    form1,
    warmup = 2000, 
    iter = 42000,
    thin = 160,
    chains = 4,
    seed = 123,
    init = 0,
    prior = priors,
    sample_prior = TRUE,
    threads = threading(12),
    backend = "cmdstanr",
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 12
    ),
    data = data1
 )
 
 # Save the model output
 saveRDS(model1, file = "glm-env-euk-bw.rds")


 # Model for env webs
 model2 <- brm(
    form2,
    warmup = 2000, 
    iter = 42000,
    thin = 160,
    chains = 4,
    seed = 123,
    init = 0,
    prior = priors,
    sample_prior = TRUE,
    threads = threading(12),
    backend = "cmdstanr",
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 12
    ),
    data = data2
 )
 
 # Save the model output
 saveRDS(model2, file = "glm-env-euk-w.rds")

# ==================================================================
# ==================================================================
