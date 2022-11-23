# ==================================================================

#            Plot shared ASVs between spiders and webs

# ==================================================================





# ==================================================================
# 1. Import the libraries and the data
# ==================================================================


# Libraries --------------------------------------------------------
 
 library(VennDiagram)
 library(ggpubr)
 library(data.table)




# Community and taxa data ------------------------------------------

 # Folder paths
 path1 <- file.path(getwd(), "data", "data-clean-env")
 
 # Community matrices
 comm1 <- readRDS(file.path(path1, "comm-env-bac-bw.rds"))
 comm2 <- readRDS(file.path(path1, "comm-env-bac-w.rds"))
 
 # Taxa data
 tax1 <- readRDS(file.path(path1, "taxa-env-bac-bw.rds"))
 tax2 <- readRDS(file.path(path1, "taxa-env-bac-w.rds"))
 
 # Metadata
 meta1 <- readRDS(file.path(path1, "metadata-env-bac-bw.rds"))
 meta2 <- readRDS(file.path(path1, "metadata-env-bac-w.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Arrange the data
# ==================================================================


# Step 1 - reshape the matrix --------------------------------------
 
 # Transpose the community matrices
 tab1 <- data.table(t(comm1), keep.rownames = T)
 setnames(tab1, "rn", "ASV")
 
 tab2 <- data.table(t(comm2), keep.rownames = T)
 setnames(tab2, "rn", "ASV")



# Step 2 - reshape the tables and count ASVs -----------------------
 
 # Long format
 tab1 <- melt(
     tab1,
     id.vars = "ASV",
     variable.name = "sample_id",
     value.name = "count")

 tab2 <- melt(
     tab2,
     id.vars = "ASV",
     variable.name = "sample_id",
     value.name = "count")
 

 # Merge env info
 tab1 <- merge(
     tab1,
     meta1[, c("sample_env", "sample_type", "sample_id")],
           by = "sample_id"
 )
 
 tab2 <- merge(
     tab2,
     meta2[, c("sample_env", "sample_type", "sample_id")],
           by = "sample_id"
 )
 
 
 # Bind the two tables together
 tab <- rbind(tab1, tab2)



# Step 3 - sum ASVs per group --------------------------------------

 tab[
   , 
   count_env_type := sum(count),
   by = .(ASV, sample_env, sample_type)
 ]

 tab[, presence := ifelse(count_env_type > 0, 1, 0)]



# Step 4 - Create a list of ASVs per group -------------------------
 
 # Extract only the ASVs that are present
 dat <- unique(
     tab[
         presence == 1, 
         .(ASV, sample_env, sample_type, presence)
     ]
 )


 # Create a list of vectors of ASV
 desert <- list(
     spiders = dat[
        sample_env == "desert" &
        sample_type == "spider",
        ASV
     ],
     webs = dat[
        sample_env == "desert" &
        sample_type == "web",
        ASV
     ]
 )

 urban <- list(
     spiders = dat[
        sample_env == "urban" &
        sample_type == "spider",
        ASV
     ],
     webs = dat[
        sample_env == "urban" &
        sample_type == "web",
        ASV
     ]
 )

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Produce the Venn diagrams
# ==================================================================


# Helper function to display Venn diagram --------------------------
 
 display_venn <- function(x, ...){
   library(VennDiagram)
   grid.newpage()
   venn_object <- venn.diagram(
    x,
    filename = NULL, 
    disable.logging = TRUE,
    ...
   )
   grid.draw(venn_object)
 }



# Produce the plots ------------------------------------------------

 # Plot for the desert environment
 dplot <- venn.diagram(
      desert,
      scaled = FALSE,
      category.names = c("Spiders", "Webs"),
      fill = c("red2", "gray80"),
      # Numbers
      cex = 1,
      cex.dist = c(0.03, 0.03),
      cex.pos = c(-1, 1),
      # Set names (titles)
      cat.cex = 1.2,
      cat.fontface = "bold",
      cat.default.pos = "text",
      cat.dist = c(0.08, 0.08),
      cat.pos = c(-1, 1),
      # Elements to show
      print.mode = c("raw", "percent"),
      # Disable outputs
      filename = NULL, 
      disable.logging = TRUE
 )


# Plot for the urban environment
 uplot <- venn.diagram(
      urban,
      scaled = FALSE,
      category.names = c("Spiders", "Webs"),
      fill = c("red2", "gray80"),
      # Numbers
      cex = 1,
      cex.dist = c(0.03, 0.03),
      cex.pos = c(-1, 1),
      # Set names (titles)
      cat.cex = 1.2,
      cat.fontface = "bold",
      cat.default.pos = "text",
      cat.dist = c(0.08, 0.08),
      cat.pos = c(-1, 1),
      # Elements to show
      print.mode = c("raw", "percent"),
      # Disable outputs
      filename = NULL, 
      disable.logging = TRUE
 )

# ==================================================================
# ==================================================================





# ==================================================================
# 4. Arrange as one figure and export
# ==================================================================

 # Folder path
 path <- file.path(getwd(), "outputs", "plots-bac")

 # Arrange
 fig <- ggarrange(
    dplot, uplot,
    nrow = 1, ncol = 2,
    labels = c("Desert samples", "Urban samples"),
    hjust = -0.7,
    vjust = 2.5
 )

 ggexport(
    fig,
    filename = file.path(path, "venn-env-bac.png"),
    width = 2000, height = 1000, res = 300
 )

# ==================================================================
# ==================================================================