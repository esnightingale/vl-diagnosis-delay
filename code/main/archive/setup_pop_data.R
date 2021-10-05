devtools::install_github('wpgp/wopr')
library(wopr)

# Retrieve the WOPR data catalogue
catalogue <- getCatalogue()

# Select files from the catalogue by subsetting the data frame
selection <- subset(catalogue,
                    country == 'NGA' &
                      category == 'Population' & 
                      version == 'v1.2')

# Download selected files
downloadData(selection)