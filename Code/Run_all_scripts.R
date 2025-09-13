# Clear all objects from the R environment
rm(list = ls())

# Set the path to the folder containing your R scripts
code_folder <- "Code"


script_files <- list.files(path = code_folder,
                           pattern = "^[0-9a-zA-Z]+\\..*\\.R$",
                           full.names = TRUE)

# Use gtools::mixedsort() for robust natural sorting
# This handles files like 1.xxx.R, 2.yyy.R, and 10.zzz.R correctly
if (!requireNamespace("gtools", quietly = TRUE)) {
  install.packages("gtools")
}
library(gtools)

script_files_ordered <- gtools::mixedsort(script_files)

# Loop through each ordered file and source it
for (file in script_files_ordered) {
  message("Sourcing file: ", file)
  source(file)
}