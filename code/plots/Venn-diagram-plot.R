# ==================================================================

#            Plot shared ASVs between spiders and webs

# ==================================================================





# ==================================================================
# 1. Import the libraries and the data
# ==================================================================


# Libraries --------------------------------------------------------
 
 library(VennDiagram)
 library(data.table)



# Community and taxa data ------------------------------------------

 # Folder paths
 path1 <- file.path(getwd(), "data", "env-data-clean")
 
 # Community matrices
 comm1 <- readRDS(file.path(path1, "env-bac-comm-bw.rds"))
 comm2 <- readRDS(file.path(path1, "env-bac-comm-w.rds"))
 
 # Taxa data
 tax1 <- readRDS(file.path(path1, "env-bac-taxa-bw.rds"))
 tax2 <- readRDS(file.path(path1, "env-bac-taxa-w.rds"))
 
 # Metadata
 meta1 <- readRDS(file.path(path1, "env-bac-metadata-bw.rds"))
 meta2 <- readRDS(file.path(path1, "env-bac-metadata-w.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Arrange the data
# ==================================================================


# Step 1 - reshape the matrix --------------------------------------
 
 # Transpose the community matrix
 comm1
 tab <- data.table(t(comm1), keep.rownames = T)
 setnames(tab, "rn", "ASV")



# Step 2 - reshape the table and count ASVs ------------------------
 
 # Long format
 tab <- melt(
     tab,
     id.vars = "ASV",
     variable.name = "sample_id",
     value.name = "count")
 
 # Merge env info
 tab <- merge(
     tab,
     meta1[, c("sample_env", "sample_id")],
           by = "sample_id"
 )



# Step 3 - sum ASVs per environment --------------------------------

 tab[, count_env := sum(count), by = .(ASV, sample_env)]
 tab[, presence := ifelse(count_env > 0, 1, 0)]



# Step 4 - Create a list of ASVs per environment -------------------
 
 # Extract only the ASVs that are present
 tab1 <- unique(
     tab[
         presence == 1, 
         .(ASV, sample_env, presence)
     ]
 )

 # Create a list of vectors of ASV
 x <- list(
     desert = tab1[sample_env == "desert", ASV],
     urban = tab1[sample_env == "urban", ASV]
 )

# ==================================================================
# ==================================================================





# ==================================================================
# 3. Produce the Venn diagram
# ==================================================================


# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(
    x,
    category.names = c("Desert", "Urban"),
    fill = c("#E69F00", "#666666"),
    # Numbers
    cex = 1.5,
    cex.dist = c(0.03, 0.03),
    cex.pos = c(-1, 1),
    # Set names
    cat.cex = 1.8,
    cat.fontface = "bold",
    cat.default.pos = "text",
    cat.dist = c(0.05, 0.05),
    cat.pos = c(-1, 1),
    # Element to show
    print.mode = c("raw", "percent")
)

# ==================================================================
# ==================================================================