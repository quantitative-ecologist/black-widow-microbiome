
library(ggpubr)

path1 <- file.path(getwd(), "code", "plots-bac")
path2 <- file.path(getwd(), "code", "plots-euk")

#fig1 <- readRDS(file.path(path1, "glm-env-bac-plot.rds"))
#fig2 <- readRDS(file.path(path2, "glm-env-euk-plot.rds"))

# These files have to be run in this order
source(file.path(path1, "glm-plots-bac.R"))
source(file.path(path2, "glm-plots-euk.R"))


# Combine into one figure
fig1 <- ggarrange(
    plot1, NULL, plot2,
    plot4, NULL, plot5,
    nrow = 2, ncol = 3,
    widths = c(1, 0.1, 1, 1, 0.1, 1),
    common.legend = TRUE,
    legend = "top",
    labels = c("(A)", "", "(B)", "(C)", "", "(D)"),
    hjust = c(-0.1, -0.1),
    vjust = c(1.5, 1.5)
)

 # Export the figure
 ggexport(
    fig1,
    filename = file.path(path, "figures", "figure1.png"),
    width = 2000,
    height = 2400,
    res = 300
 )
