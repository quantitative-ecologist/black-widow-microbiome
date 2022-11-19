# ==================================================================

#     Plot shared ASVs between spiders with varying diets

# ==================================================================





# ==================================================================
# 1. Import the libraries and the data
# ==================================================================


# Libraries --------------------------------------------------------
 
 library(VennDiagram)
 library(viridis)
 library(ggpubr)
 library(data.table)




# Community and taxa data ------------------------------------------

 # Folder paths
 path <- file.path(getwd(), "data", "data-clean-diet")
 
 # Community matrices
 comm <- readRDS(file.path(path, "comm-diet-bac-bw.rds"))
 
 # Metadata
 meta <- readRDS(file.path(path, "metadata-diet-bac-bw.rds"))

# ==================================================================
# ==================================================================





# ==================================================================
# 2. Arrange the data
# ==================================================================


# Step 1 - reshape the matrix --------------------------------------
 
 # Transpose the community matrices
 tab <- data.table(t(comm), keep.rownames = T)
 setnames(tab, "rn", "ASV")



# Step 2 - reshape the tables and count ASVs -----------------------
 
 # Long format
 tab <- melt(
     tab,
     id.vars = "ASV",
     variable.name = "sample_id",
     value.name = "count")
 

 # Merge env info
 tab <- merge(
     tab,
     meta[, c("diet_treatment", "sample_id")],
           by = "sample_id"
 )



# Step 3 - sum ASVs per group --------------------------------------

 tab[
   , 
   count_diet := sum(count),
   by = .(ASV, diet_treatment)
 ]

 tab[, presence := ifelse(count_diet > 0, 1, 0)]



# Step 4 - Create a list of ASVs per group -------------------------
 
 # Extract only the ASVs that are present
 dat <- unique(
     tab[
         presence == 1, 
         .(ASV, diet_treatment, presence)
     ]
 )


 # Create a list of vectors of ASV
 spiders <- list(
     not_fed = dat[
        diet_treatment == "not fed",
        ASV
     ],
     isopod = dat[
        diet_treatment == "isopod",
        ASV
     ],
     cricket = dat[
        diet_treatment == "cricket",
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
 splot <- venn.diagram(
      spiders,
      scaled = FALSE,
      category.names = c("Not fed", "Isopod", "Cricket"),
      fill = c("#FDE725FF", "#1F968BFF", "#440154FF"),
      # Numbers
      cex = 1,
      #cex.dist = c(0.03, 0.03),
      #cex.pos = c(-1, 1),
      # Set names (titles)
      cat.cex = 1.2,
      cat.fontface = "bold",
      #cat.default.pos = "text",
      #cat.dist = c(0.08, 0.08),
      #cat.pos = c(-1, 1),
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
    splot,
    nrow = 1, ncol = 1
 )

 ggexport(
    fig,
    filename = file.path(path, "venn-diet-bac.png"),
    width = 1200, height = 1200, res = 300
 )

# ==================================================================
# ==================================================================