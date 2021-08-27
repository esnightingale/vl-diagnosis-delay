################################################################################
# Setup packages and environment for data setup
################################################################################
################################################################################

# Dev version of ggmap
# devtools::install_github("dkahle/ggmap", force = TRUE)
# devtools::install_github('joenomiddlename/PStestR', dependencies=T, build_vignettes=T)
# remotes::install_github("gfalbery/ggregplot")

# Packages required
packages <- c("tidyverse","lubridate","gridExtra","sf","spdep","rgdal","ggspatial",
              "ggmap","patchwork","scales","here", "gstat", "variosig", "biscale",
              "INLA","inlabru", "inlatools", "ggregplot","mapr", "mgcv")

# Check and install if necessary
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# ggmap::register_google(key = "AIzaSyDBz7TtWW69ZGDHmThLgSzEuVEXrEt8pEQ", write = TRUE)

# Set default theme for plotting
theme_set(theme_minimal())

# Source helper functions
list.files(here::here("code","utils"), full.names = TRUE) %>% purrr::walk(source)

