
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(signal)
  library(tidyverse)
  library(ggsci)
  library(readr)
  

# --- 1. Define your .txt file path ---
file_path <- "Data/20190416 layout ATPase mutants and Mn.txt" 

# --- 2. Read the relevant data from the .txt file ---
raw_data <- read_delim(
  file_path,
  delim = "\t", 
  skip = 1,     # Skip the very first row (Excel's row 1 with the '#' comment)
  col_names = FALSE,
  trim_ws = TRUE, # Trim leading/trailing whitespace
  show_col_types = FALSE, # Suppress messages about column types
  col_types = readr::cols(.default = "c") # IMPORTANT: Read ALL columns as character to preserve your factors
)


# --- 3. Identify the rows containing the descriptors using the '>' marker ---
# Look for rows where the first column (R's index 1, Excel's A) contains ">".
descriptor_rows_indices <- which(raw_data[[1]] == ">")

# Initialize an empty list to store the intermediate melted data frames for each block.
all_melted_blocks_long_format <- list()
block_counter <- 0

# --- 4. Loop through each identified descriptor block and process it (melt it) ---
for (i in seq_along(descriptor_rows_indices)) {
  block_counter <- block_counter + 1
  
  current_descriptor_r_row <- descriptor_rows_indices[i]
  
  # Extract the descriptor text from column B (index 2 in raw_data)
  # Clean the descriptor name for use as a column name (replace spaces/special chars)
  descriptor_value_raw <- as.character(raw_data[current_descriptor_r_row, 2])
  descriptor_value_clean <- str_replace_all(descriptor_value_raw, "[^[:alnum:]]", "_") %>% # Replace non-alphanumeric with underscore
    str_replace_all("__+", "_") %>% # Replace multiple underscores with one
    str_trim() %>% # Trim leading/trailing underscores/spaces
    str_to_lower() # Convert to lowercase for consistent column names
  
  # Ensure the descriptor doesn't start with a number if it becomes purely numeric after cleaning
  if (grepl("^[0-9]", descriptor_value_clean)) {
    descriptor_value_clean <- paste0("desc_", descriptor_value_clean)
  }
  
  # Determine the start and end rows (in R's dataframe indexing) for the well plate data.
  data_start_r_row <- current_descriptor_r_row + 1
  data_end_r_row <- current_descriptor_r_row + 16
  
  if (data_end_r_row > nrow(raw_data)) {
    warning(paste("Data block for descriptor '", descriptor_value_raw, "' ends prematurely. Adjusting end row."))
    data_end_r_row <- nrow(raw_data)
  }
  
  well_data_block <- raw_data[data_start_r_row:data_end_r_row, 2:25]
  
  colnames(well_data_block) <- as.character(1:24) # Assign column numbers
  
  row_letters <- LETTERS[1:nrow(well_data_block)]
  well_data_block$row_letter <- row_letters
  
  # --- 5. Melt (pivot longer) the extracted data block ---
  melted_block <- well_data_block %>%
    pivot_longer(
      cols = -row_letter,
      names_to = "column_number",
      values_to = descriptor_value_clean # Assign the clean descriptor name directly as the value column name
    ) %>%
    mutate(
      well = paste0(row_letter, column_number) # Create the well identifier
    ) %>%
    select(well, all_of(descriptor_value_clean)) # Select only well and the new descriptor column
  
  all_melted_blocks_long_format[[block_counter]] <- melted_block
}

# --- 6. Combine all melted blocks ---
# Start with the first block to set up the 'well' column
f <- all_melted_blocks_long_format[[1]]

# Join subsequent blocks by 'well' to add them as new columns
if (length(all_melted_blocks_long_format) > 1) {
  for (j in 2:length(all_melted_blocks_long_format)) {
    f <- full_join(f, all_melted_blocks_long_format[[j]], by = "well")
  }
}

# --- 7. View the final data frame structure and some of its data ---
print(head(f))
print(tail(f))
print(summary(f))
print(str(f)) # This should now show character columns for your descriptor data





#read data
file_path_data <- "Data/20190416 data ATPase mutants and Mn.txt" 

data<-read.table(file_path_data, header = TRUE, sep="\t", dec = ".")






