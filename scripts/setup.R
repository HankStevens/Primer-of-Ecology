# Shared setup for every chapter.
#
# Unlike bookdown -- which rendered the whole book in one persistent R
# session, so a single library() call in index.Rmd stayed loaded for every
# later chapter -- Quarto renders each .qmd file in its own fresh R
# session. So every chapter needs to load what it uses. Rather than repeat
# the same library() list in all 15 files, each chapter's first chunk just
# does source("scripts/setup.R").
#
# Add a package here once and every chapter picks it up on next render.

x <- c("bbmle", "bipartite", "data.table", "dagitty", "diagram", "DiagrammeR", "DiagrammeRsvg",
       "igraph", "kableExtra", "knitr", "lattice",
       "lavaan", "magrittr", "MARSS", "patchwork", "phaseR",
       "pracma",
       "primer", "reshape2", "rARPACK", "rsvg", "scatterplot3d", "semPlot",
       "tidyverse", "tinytable", "tufte", "untb", "vegan", "xtable")

invisible(lapply(x, library, character.only = TRUE))

theme_set(theme_minimal() +
theme(panel.grid.major = element_line(colour = "darkgrey"),
      panel.grid.minor = element_line(colour = "darkgrey")) )
theme(panel.grid.major = element_line(colour = "darkgrey"))
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
