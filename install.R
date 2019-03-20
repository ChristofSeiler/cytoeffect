devtools::install_github("ChristofSeiler/cytoeffect")
devtools::install_github("RGLab/ggcyto", ref="trunk")
pkgs_needed = c("devtools","tidyverse","magrittr","SummarizedExperiment",
                "ggthemes","cowplot","RColorBrewer","broom","hexbin",
                "intergraph","igraph","ggnetwork","ggcorrplot","MASS",
                "parallel","dplyr","knitr","dagitty","ggdag")
source("http://bioconductor.org/biocLite.R")
biocLite(pkgs_needed)
