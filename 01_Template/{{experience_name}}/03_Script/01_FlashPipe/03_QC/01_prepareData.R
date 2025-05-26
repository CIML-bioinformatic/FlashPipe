## ####################################################
## This script aims to load and prepare the data
##  for the downstream analysis
## ####################################################

## @knitr prepare_data

# Retrieve the name of the plate
# Put them in a list, and will be used in all our code.
PLATES_LIST <- list.dirs( PATH_ZUMIS_OUTPUT, full.names = FALSE, recursive = FALSE)
PLATES_LIST <- PLATES_LIST[PLATES_LIST != ""]
PLATES_LIST <- gsub("/$", "", PLATES_LIST)
# Removes the zUMIs_output directory, which has no place in this list.
if (any(PLATES_LIST == "zUMIs_output")){
  PLATES_LIST = PLATES_LIST[PLATES_LIST != "zUMIs_output"]
}

# Load the configuration file
config_file_yaml <- yaml.load_file(PATH_CONFIG_FILE)
origin_plate_name_list = gsub( " ", "", strsplit( config_file_yaml$plate_names[[1]], fixed= TRUE, split = ",")[[1]])

# Compare the observed list of plates and the user provided list of plates
filtered_origin_plate_name_list <- origin_plate_name_list[origin_plate_name_list %in% PLATES_LIST]

# Order the plate list by the user provided order
PLATES_LIST = factor( PLATES_LIST, levels = filtered_origin_plate_name_list)

cat("<BR><b>Filtered plates to analyse=</b>", paste(filtered_origin_plate_name_list, collapse="; "))
rm(config_file_yaml)
rm(origin_plate_name_list)
rm(filtered_origin_plate_name_list)

########################################### Verification of every input and potential error/warning ###########################################

### Verification for 02_formatCountTable

# Initialize an error list
plate_errors_cell_barcode_well <- c()

# Load the barcode-to-well mapping file
CELL_BARCODE_WELL_DF <- read.csv(file = PATH_WELL_BARCODE_FILE)
# Check if some column name are in the barcode cell with well
if (!"WellID" %in% colnames(CELL_BARCODE_WELL_DF)){
  plate_errors_cell_barcode_well <- c(plate_errors_cell_barcode_well, paste("\nERROR : There is no column named 'WellID' in the barcode_well file. \nDirectory :", PATH_WELL_BARCODE_FILE))
}
if (!"BarcodeSequence" %in% colnames(CELL_BARCODE_WELL_DF)){
  plate_errors_cell_barcode_well <- c(plate_errors_cell_barcode_well, paste("\nERROR : There is no column named 'BarcodeSequence' in the barcode_well file. \nDirectory :", PATH_WELL_BARCODE_FILE))
}
if (length(CELL_BARCODE_WELL_DF[[WELL_ID]]) > 96 | length(CELL_BARCODE_WELL_DF$BarcodeSequence) > 96){
  plate_errors_cell_barcode_well <- c(plate_errors_cell_barcode_well, paste("\nERROR : The number of WellID in the cell_barcode_well file is too high, number of wells detected :", length(CELL_BARCODE_WELL_DF$WellID)))
}
if (length(CELL_BARCODE_WELL_DF[[WELL_ID]]) != length(CELL_BARCODE_WELL_DF$BarcodeSequence)){
  plate_errors_cell_barcode_well <- c(plate_errors_cell_barcode_well, paste("\nERROR : The number of WellID and BarcodeSequence differs. -- Number of wells :", length(CELL_BARCODE_WELL_DF$WellID), " -- Number of barcode :", length(CELL_BARCODE_WELL_DF$BarcodeSequence)))
}
if (!is(CELL_BARCODE_WELL_DF[[WELL_ID]], "character")){
  plate_errors_cell_barcode_well <- c(plate_errors_cell_barcode_well, "\nERROR : The column WellID is not of type character in the cell_barcode file")
}
if (!is(CELL_BARCODE_WELL_DF$BarcodeSequence, "character")){
  plate_errors_cell_barcode_well <- c(plate_errors_cell_barcode_well, "\nERROR : BarcodeSequence column is not of type character in the cell_barcode file")
}

# Display errors if any were found 
if (length(plate_errors_cell_barcode_well) > 0) {
  cat("\n••••• ERROR DETECTED IN CELL BARCODE WELL FILE (01_reference)•••••\n")
  cat(plate_errors_cell_barcode_well, "\n")
}

# Initialize an error list
plate_errors_zUMIs_output <- c()

for (plate_name in PLATES_LIST) {
  # Define the directory path for the current plate
  plate_directory <- file.path(PATH_ZUMIS_OUTPUT, plate_name)
  
  # Paths to files
  rds_file <- file.path(plate_directory, "zUMIs_output", "expression", paste0(plate_name, ".dgecounts.rds"))
  gene_names_file <- file.path(plate_directory, "zUMIs_output", "expression", paste0(plate_name, ".gene_names.txt"))
  
  if (!file.exists(rds_file)){
    plate_errors_zUMIs_output <- c(plate_errors_zUMIs_output, paste("\nERROR : RDS file or gene_names.txt file for plate", plate_name, "does not exist at specified location :", rds_file))
  }
  if (!file.exists(gene_names_file)){
    plate_errors_zUMIs_output <- c(plate_errors_zUMIs_output, paste("\nERROR : Gene_names.txt file for plate", plate_name, "does not exist at specified location : ", gene_names_file))
  }
  
  file_rds <- readRDS(rds_file)
  if (!is(file_rds, "list")){
    plate_errors_zUMIs_output <- c(plate_errors_zUMIs_output, paste("The RDS file for", plate_name,"is not of the correct type (other than a list) please check zUMIs output results"))
  }
  
  # Load the gene names file (tab-delimited)
  gene_names_df <- read.table(gene_names_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Check the type of data inside
  if (!is(gene_names_df, "data.frame")){
    plate_errors_zUMIs_output <- c(plate_errors_zUMIs_output, paste("\nERROR : Gene_names.txt is not of type data.frame for", plate_name ,"please verify the data"))
  }
  if (!is(gene_names_df$gene_id, "character")){
    plate_errors_zUMIs_output <- c(plate_errors_zUMIs_output, paste("\nERROR : gene_id column is not of type character in the Gene_names.txt for", plate_name))
  }
  if (!is(gene_names_df$gene_name, "character")){
    plate_errors_zUMIs_output <- c(plate_errors_zUMIs_output, paste("\nERROR : gene_name column is not of type character in the Gene_names.txt for", plate_name))
  }
}

# Display errors if any were found 
if (length(plate_errors_zUMIs_output) > 0) {
  cat("\n••••• ERROR DETECTED IN zUMIs OUTPUT•••••\n")
  cat(plate_errors_zUMIs_output, "\n")
}


rm(plate_directory)
rm(rds_file)
rm(gene_names_file)
rm(file_rds)
rm(gene_names_df)

### Verification for 03_computeERCCUMIPercent

### Verification for 04_computeERCCAccuracy

ERCC_EXPECTED_DF <- read.csv(PATH_ERCC_CONCENTRATION)

# Initialize an error list
plate_errors_ERCC <- c()

# Check the type of data inside
if (!is.data.frame(ERCC_EXPECTED_DF)){
  plate_errors_ERCC <- c(plate_errors_ERCC, "\nERROR : The ERCC file (in 01_Reference/00_Experiment), is not a dataframe. Check that this is the correct file.")
}
if (!"ERCC_ID" %in% colnames(ERCC_EXPECTED_DF)){
  plate_errors_ERCC <- c(plate_errors_ERCC, paste("\nERROR : There is no column named 'ERCC_ID' in the barcode_well file. \nDirectory :", PATH_ERCC_CONCENTRATION))
}
if (!"Concentration_attomoles_ul" %in% colnames(ERCC_EXPECTED_DF)){
  plate_errors_ERCC <- c(plate_errors_ERCC, paste("\nERROR : There is no column named 'Concentration_attomoles_ul' in the barcode_well file. \nDirectory :", PATH_ERCC_CONCENTRATION))
}
if (!is(ERCC_EXPECTED_DF, "data.frame")){
  plate_errors_ERCC <- c(plate_errors_ERCC, "\nERROR : BarcodeSequence column is not of type data.frame in the ERRC csv file")
}
if (!is(ERCC_EXPECTED_DF$ERCC_ID, "character")){
  plate_errors_ERCC <- c(plate_errors_ERCC, "\nERROR : BarcodeSequence column is not of type character in the ERRC csv file")
}
if (!is(ERCC_EXPECTED_DF$Concentration_attomoles_ul, "numeric")){
  plate_errors_ERCC <- c(plate_errors_ERCC, "\nERROR : BarcodeSequence column is not of type numeric in the ERRC csv file")
}

# Display errors if any were found 
if (length(plate_errors_ERCC) > 0) {
  cat("\n••••• ERROR DETECTED IN ERCC FILE•••••\n")
  cat(plate_errors_ERCC, "\n")
}

### Verification for 05_computeRNAUMI

### Verification for 06_computeFeatureNumber

### Verification for 07_computeIndexSort
if (INDEXSORT == TRUE){
  # Initialize an error list
  plate_errors_indexsort <- c()
  
  columns_names = ""
  index_sort_colnames = c()
  for (plate_name in levels(PLATES_LIST)) {
    path_index_sort_file = paste(PATH_EXPERIMENT_RAWDATA, "/01_IndexSort/", plate_name, "_indexsort.csv", sep = "")
    index_sort_df = read.csv(path_index_sort_file)
    
    for (categorial_term_set in CATEGORIAL_TERM_SET){
      if (!categorial_term_set %in% colnames(index_sort_df)){
        plate_errors_indexsort <- c(plate_errors_indexsort, paste("\nERROR : There is no column named ", categorial_term_set," (the name here is the one provided in the config user) in the indexsort file for plate : ", plate_name))
      }
    }
    
    if (length(columns_names)>1){
      for (name_column in colnames(index_sort_df)){
        if (!name_column %in% columns_names){
          plate_errors_indexsort <- c(plate_errors_indexsort, paste("\nERROR : The column name between",plate_name_verification, "and", plate_name,"are different in indexsort file. Correct the name column between them"))
        }
      } 
    }
    columns_names = colnames(index_sort_df)
    plate_name_verification = plate_name
    
    if (!WELL_ID %in% colnames(index_sort_df)){
      plate_errors_indexsort <- c(plate_errors_indexsort, paste("\nERROR : There is no column named 'WellID' in the indexsort file for plate : ", plate_name))
    }
    
    if (!is(index_sort_df[[WELL_ID]], "character")){
      plate_errors_indexsort <- c(plate_errors_indexsort, paste("\nERROR : The column WellID is not of type character in the indexsort file for plate : ", plate_name))
    }
    
    # Check the length of index_sort_df$Well.ID
    well_count = length(index_sort_df[[WELL_ID]])
    if (well_count < 96) {
      cat(paste("WARNING : Some wells are not supplied in plate:", plate_name, "\nTotal number of wells present:", well_count))
    } else if (well_count > 96) {
      plate_errors_indexsort <- c(plate_errors_indexsort, paste("\nERROR : Too many wells/lines in plate index sort file:", plate_name, "\nTotal number of wells/lines:", well_count))
    }
    
    # Check if the plates have the same column names
    if( length( index_sort_colnames) >0){
      if( sum( sort( index_sort_colnames) == sort( names( index_sort_df))) != ncol( index_sort_df) ){
        plate_errors_indexsort <- c(plate_errors_indexsort, paste( "\nERROR : The column names of the plates are not equals :\n Previous plates used:", paste( sort( index_sort_colnames), collapse = ";"), "\nwhereas plate", plate_name, "uses :", paste( sort( names( index_sort_df)), collapse = ";")))
      }
    }else{
      index_sort_colnames = names( index_sort_df)
    }
  }
  
  # Display errors if any were found 
  if (length(plate_errors_indexsort) > 0) {
    cat("\n••••• ERROR DETECTED IN INDEXSORT FILE•••••\n")
    cat(plate_errors_indexsort, "\n")
  }
  
  rm(path_index_sort_file)
  rm(index_sort_df)
  rm(categorial_term_set)
  rm(columns_names)
  rm(plate_name_verification)
  rm(index_sort_colnames)
}

### Verification for 08_displayMetaData

# Initialize an error list
plate_errors_metadata <- c()

## Check before run the process.
# Check names between config file and sheets and library (both present in GSF).
sheets_gsf = excel_sheets(PATH_GSF)
plate_sheet_verify = read_excel(PATH_GSF, sheet = "General Info")
plate_sheet_verify = as.data.frame(plate_sheet_verify)
for (plate_name in levels(PLATES_LIST)){
  if (any(plate_sheet_verify == plate_name) != TRUE){
    cat(paste("WARNING : The plate", plate_name, "given in the config file doesn't match the one given in the 'General Info' tab of the GSF."))
  }
  if (any(plate_name == sheets_gsf) != TRUE){
    plate_errors_metadata <- c(plate_errors_metadata, paste("\nERROR : The plate", plate_name, "given in the config file does not match the GSF sheets."))
  }
}

# Check if the column name are the same in all the different plate file
columns_names = ""
for (plate_name in levels(PLATES_LIST)){
  plate_sheet = read_excel(PATH_GSF, sheet = plate_name)
  if (length(columns_names)>1){
    for (name_column in colnames(plate_sheet)){
      if (!name_column %in% columns_names){
        plate_errors_metadata <- c(plate_errors_metadata, paste("\nERROR : The column name between",plate_name_verification, "and", plate_name,"in the GSF file (sheet) are different. Correct the name column between them"))
      }
    } 
  }
  columns_names = colnames(plate_sheet)
  plate_name_verification = plate_name
  if (!WELL_ID %in% colnames(plate_sheet)){
    plate_errors_metadata <- c(plate_errors_metadata, paste0("\nERROR : WELL_ID (", WELL_ID, ") column is not found in the column name of the GSF file.\nCheck the variable in analysisParams.R or you're column in the GSF file."))
  }
}

# Display errors if any were found 
if (length(plate_errors_metadata) > 0) {
  cat("\n••••• ERROR DETECTED IN METADATA FILE (gsf)•••••\n")
  cat(plate_errors_metadata, "\n")
}

rm(sheets_gsf)
rm(plate_sheet_verify)
rm(plate_sheet)
rm(plate_name_verification)

# Regular expression to capture plate names beginning with “P” followed by digits and a “_”.
# plate_pattern <- "^P\\d{1,2}_.+" 

# Extract the names of plates in sheets_gsf that correspond to this format
# plates_gsf <- sheets_gsf[grepl(plate_pattern, sheets_gsf)]
# 
# # Verify that every plate in plates_gsf is present in plate_sheet_verify
# for (plate_name in plates_gsf) {
#   if (!any(plate_name == plate_sheet_verify)) {
#     warning("The plate", plate_name ," is not present in the libriaire 'General Info' of the GSF.")
#   }
# }

### Verification for 09_computeTrust4
if (BCR == TRUE & TCR == TRUE){
  # Initialize an error list
  plate_errors_trust4 <- c()
  
  columns_names = ""
  for (plate_name in levels(PLATES_LIST)) {
    plate_directory <- file.path(PATH_TRUST4_OUTPUT, plate_name)
    tsv_air_file <- list.files(plate_directory, pattern = paste0("^", plate_name, "_barcode_airr.tsv$"), full.names = TRUE)
    if (!file.exists(tsv_air_file)) {
      plate_errors_trust4 <- c(plate_errors_trust4, paste("\nERROR : There is no barcode_airr.tsv file (TRUST4 output) for", plate_name))
    }
    file_air_tsv <- read.csv(tsv_air_file, sep = '\t')
    if (length(file_air_tsv) <= 1){
      plate_errors_trust4 <- c(plate_errors_trust4, paste("\nERROR : Le fichier barcode_airr.tsv file (TRUST4 output) for", plate_name, "is empty"))
    }
    if (!is(file_air_tsv, "data.frame")){
      plate_errors_trust4 <- c(plate_errors_trust4, paste("\nERROR : The barcode_airr.tsv file (TRUST4 output) is not a dataframe for", plate_name))
    }
    if (length(columns_names)>1){
      for (name_column in colnames(file_air_tsv)){
        if (!name_column %in% columns_names){
          plate_errors_trust4 <- c(plate_errors_trust4, paste("\nERROR : The column name between",plate_name_verification, "and", plate_name,"in the barcode_airr.tsv file (TRUST4 output) are different. Correct the name column between them."))
        }
      } 
    }
    columns_names = colnames(file_air_tsv)
    plate_name_verification = plate_name
    
    if (!C_CALL %in% colnames(file_air_tsv)){
      plate_errors_trust4 <- c(plate_errors_trust4, paste0("\nERROR : C_CALL (", C_CALL, ") column is not found in the column name of barcode_airr.tsv file (TRUST4 output).\nCheck the variable in analysisParams.R or you're column in the tsv file."))
    }
    if (!V_CALL %in% colnames(file_air_tsv)){
      plate_errors_trust4 <- c(plate_errors_trust4, paste0("\nERROR : V_CALL (", V_CALL, ") column is not found in the column name of barcode_airr.tsv file (TRUST4 output).\nCheck the variable in analysisParams.R or you're column in the tsv file."))
    }
    if (!J_CALL %in% colnames(file_air_tsv)){
      plate_errors_trust4 <- c(plate_errors_trust4, paste0("\nERROR : J_CALL (", J_CALL, ") column is not found in the column name of barcode_airr.tsv file (TRUST4 output).\nCheck the variable in analysisParams.R or you're column in the tsv file."))
    }
    if (!PRODUCTIVE %in% colnames(file_air_tsv)){
      plate_errors_trust4 <- c(plate_errors_trust4, paste0("\nERROR : PRODUCTIVE (", PRODUCTIVE, ") column is not found in the column name of barcode_airr.tsv file (TRUST4 output).\nCheck the variable in analysisParams.R or you're column in the tsv file."))
    }
    if (!CELL_ID %in% colnames(file_air_tsv)){
      plate_errors_trust4 <- c(plate_errors_trust4, paste0("\nERROR : CELL_ID (", CELL_ID, ") column is not found in the column name of barcode_airr.tsv file (TRUST4 output).\nCheck the variable in analysisParams.R or you're column in the tsv file."))
    }
  }
  
  # Display errors if any were found 
  if (length(plate_errors_trust4) > 0) {
    cat("\n••••• ERROR DETECTED IN OUTPUT TRUST4•••••\n")
    cat(plate_errors_trust4, "\n")
  }
  rm(plate_directory)
  rm(tsv_air_file)
  rm(plate_name)
  rm(columns_names)
  rm(name_column)
  rm(file_air_tsv)
  rm(plate_name_verification)
  
  if (length(plate_errors_zUMIs_output) > 0 | length(plate_errors_ERCC) > 0 | length(plate_errors_cell_barcode_well) > 0 | length(plate_errors_indexsort) > 0 | length(plate_errors_metadata) > 0 | length(plate_errors_trust4) > 0){
    stop("Some error has been detecting. Stop analysis")
  }
  rm(plate_errors_zUMIs_output)
  rm(plate_errors_ERCC)
  rm(plate_errors_cell_barcode_well)
  rm(plate_errors_indexsort)
  rm(plate_errors_metadata)
  rm(plate_errors_trust4) 
}