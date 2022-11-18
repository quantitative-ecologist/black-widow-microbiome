# ==================================================================

#              Compare ASV richness for diet spiders

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
 path <- file.path(getwd(), "diet-data-clean")

 # Community matrix
 comm_bw_diet <- readRDS(file.path(path, "diet-bac-comm-bw.rds"))
 
 # Metadata
 meta_bw_diet <- readRDS(file.path(path, "diet-bac-metadata-bw.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Assemble the diet spider data
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
 data <- unique(diet_bw[sp_present == 1,
                        .(sample_id, sp_richness)])
 
 # Merge with the metadata to add covariates
 data <- merge(
   data,
   meta_bw_diet[,c(1, 3, 4, 5)],
   by = "sample_id")



# Transform the data -----------------------------------------------

 # Log transform species richness
 data$log_sp_richness <- log(data$sp_richness)
 
 # Environment as factor
 data$diet_treatment <- as.factor(data$diet_treatment)

 # Transform covariates
 stdz <- function (x) {(x - mean(x)) / sd(x)}

 data[, Zspider_weight := lapply(.SD, stdz), .SDcols = "spider_weight"]
 data[, logn_reads := log(n_reads)]

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Fit the model
# ==================================================================


# Setup the model formula ------------------------------------------

 # Gaussian model for diet spiders
 form <- bf(
  log_sp_richness ~ 
    1 + 
    diet_treatment + 
    Zspider_weight  + 
    logn_reads
 ) + 
 gaussian()



# Setup the priors -------------------------------------------------

 # Gaussian priors for b and intercepts
 priors <- c(
   set_prior("normal(0, 2)", 
             class = "Intercept"),
   set_prior("normal(0, 2)",
             class = "b")
 )



# Run the model ----------------------------------------------------
 

 # Model for diet spiders
 model <- brm(form,
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
              data = data)
 
 # Save the model output
 saveRDS(model, file = "diet-bac-glm-bw.rds")

# ==================================================================
# ==================================================================