## ####################################################
## The purpose of this script is to prepare raw data 
## for use in quality control and future analyses.
## ####################################################

## @knitr prepare_data

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# General description of the code :
# 1. Retrieve plate data from yaml file and verify
# 2. Vérification all data produced by tools (script by script)
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

###########################################################
# ## 1. Retrieve plate data from yaml file and verify #####
###########################################################

# Retrieves the names of directories in the zUMIs directory, which contains only the names of plates (automatically built by Copy).
## Retrieves the list of directories (names of "plates") in the zUMIs output folder.
PLATES_LIST <- list.dirs( PATH_ZUMIS_OUTPUT, full.names = FALSE, recursive = FALSE)

## Deletes empty elements (in case list.dirs returns an empty string)
PLATES_LIST <- PLATES_LIST[PLATES_LIST != ""]

## Deletes any terminal slash ("/") at the end of directory names
PLATES_LIST <- gsub("/$", "", PLATES_LIST)

# Removes the “zUMIs_output” entry from the list, which does not represent a “plate” to be analyzed
if (any(PLATES_LIST == "zUMIs_output")){
  PLATES_LIST = PLATES_LIST[PLATES_LIST != "zUMIs_output"]
}

# Load the configuration file
config_file_yaml <- yaml.load_file(PATH_CONFIG_FILE)

# Checking information: Retrieve the plate name information given in the FlashPipe configuration file, and the plate names obtained from the directories in zUMIs.
## Retrieves the list of “plate” names supplied by the user from the YAML file
## Removes spaces and splits string using comma as separator
origin_plate_name_list = gsub( " ", "", strsplit( config_file_yaml$plate_names[[1]], fixed= TRUE, split = ",")[[1]])
## Filters config plate names to keep only those actually present in the zUMIs folder (and therefore the plates necessarily present in the analysis)
filtered_origin_plate_name_list <- origin_plate_name_list[origin_plate_name_list %in% PLATES_LIST]
# Warning if there is a difference between them
if (length(filtered_origin_plate_name_list) != length(origin_plate_name_list) | length(filtered_origin_plate_name_list) != length(PLATES_LIST) ){
  cat("WARNING : A difference in the number of plates between the configuration file (FlashPipe) and the directories created for the analysis was detected.")
}

# Reorders the PLATES_LIST list according to the order given by the user in the configuration file
PLATES_LIST = factor( PLATES_LIST, levels = filtered_origin_plate_name_list)

# Displays the list of plates selected for analysis
cat("<BR><b>Filtered plates to analyse=</b>", paste(filtered_origin_plate_name_list, collapse="; "))

# Environment cleanup: removes temporary objects
rm(config_file_yaml)
rm(origin_plate_name_list)
rm(filtered_origin_plate_name_list)

########################################################################
# ## 2. Vérification all data produced by tools (script by script) #####
########################################################################

##################################################
# ## Verification for code in 02_formatCountTable
##################################################

# We store the various warning in variables in order to display them all together instead of stopping the code at each warning.

# Initialize an error list
plate_errors_cell_barcode_well <- c()

# Loads the correspondence file between barcodes and wells
CELL_BARCODE_WELL_DF <- read.csv(file = PATH_WELL_BARCODE_FILE)

# Checks that the ‘WellID’ column is present in the file
if (!"WellID" %in% colnames(CELL_BARCODE_WELL_DF)){
  plate_errors_cell_barcode_well <- c(plate_errors_cell_barcode_well, paste("ERROR : There is no column named 'WellID' in the barcode_well file. \nDirectory :", PATH_WELL_BARCODE_FILE))
}
# Checks that the ‘BarcodeSequence’ column is present in the file
if (!"BarcodeSequence" %in% colnames(CELL_BARCODE_WELL_DF)){
  plate_errors_cell_barcode_well <- c(plate_errors_cell_barcode_well, paste("ERROR : There is no column named 'BarcodeSequence' in the barcode_well file. \nDirectory :", PATH_WELL_BARCODE_FILE))
}
# Check that the number of lines in WellID and BarcodeSequence does not exceed 96 (corresponding to the 96 wells of a standard plate).
if (length(CELL_BARCODE_WELL_DF[[COLUMN_HEADER_WELL_ID]]) > 96 | length(CELL_BARCODE_WELL_DF$BarcodeSequence) > 96){
  plate_errors_cell_barcode_well <- c(plate_errors_cell_barcode_well, paste("ERROR : The number of WellID in the cell_barcode_well file is too high, number of wells detected :", length(CELL_BARCODE_WELL_DF$WellID)))
}
# Checks that the number of wells (WellID) corresponds exactly to the number of barcode sequences
if (length(CELL_BARCODE_WELL_DF[[COLUMN_HEADER_WELL_ID]]) != length(CELL_BARCODE_WELL_DF$BarcodeSequence)){
  plate_errors_cell_barcode_well <- c(plate_errors_cell_barcode_well, paste("ERROR : The number of WellID and BarcodeSequence differs.\n-- Number of wells :", length(CELL_BARCODE_WELL_DF$WellID), "\n-- Number of barcode :", length(CELL_BARCODE_WELL_DF$BarcodeSequence)))
}
# Checks that the WellID column is of character type
if (!is(CELL_BARCODE_WELL_DF[[COLUMN_HEADER_WELL_ID]], "character")){
  plate_errors_cell_barcode_well <- c(plate_errors_cell_barcode_well, "ERROR : The column WellID is not of type character in the cell_barcode file")
}
# Checks that the BarcodeSequence column is of character type
if (!is(CELL_BARCODE_WELL_DF$BarcodeSequence, "character")){
  plate_errors_cell_barcode_well <- c(plate_errors_cell_barcode_well, "ERROR : BarcodeSequence column is not of type character in the cell_barcode file")
}

# Display errors if any were found 
if (length(plate_errors_cell_barcode_well) > 0) {
  # Creates a header
  header <- "\n============== ERROR(S) DETECTED IN CELL BARCODE FILE (01_reference) FROM SCRIPT '01_prepareData' ==============\n"
  # Numbering of all errors found in this section
  error_list <- paste0(seq_along(plate_errors_cell_barcode_well), ". ", plate_errors_cell_barcode_well, collapse = "\n")
  # Final message
  full_message <- paste0(header, error_list, "\n================================================================================================================\n")
  # Displays and stops
  stop(full_message)
}

# Initialize an error list
plate_errors_zUMIs_output <- c()

# Check files and parameters for zUMIs execution
for (plate_name in PLATES_LIST) {
  # Define the directory path for the current plate
  plate_directory <- file.path(PATH_ZUMIS_OUTPUT, plate_name)
  
  # Builds path to .rds file containing expression counts
  rds_file <- file.path(plate_directory, "zUMIs_output", "expression", paste0(plate_name, ".dgecounts.rds"))
  # Builds the path to the text file containing gene names
  gene_names_file <- file.path(plate_directory, "zUMIs_output", "expression", paste0(plate_name, ".gene_names.txt"))
  
  # Checks that the RDS file exists, otherwise adds an error
  if (!file.exists(rds_file)){
    plate_errors_zUMIs_output <- c(plate_errors_zUMIs_output, paste("ERROR : RDS file for plate", plate_name, "does not exist at specified location :", rds_file))
  }
  # Checks that the gene name file exists, otherwise adds an error
  if (!file.exists(gene_names_file)){
    plate_errors_zUMIs_output <- c(plate_errors_zUMIs_output, paste("ERROR : Gene_names.txt file for plate", plate_name, "does not exist at specified location :", gene_names_file))
  }
  # Load RDS file (expression counts)
  file_rds <- readRDS(rds_file)
  # Checks that the RDS file is a list (expected type for zUMIs results)
  if (!is(file_rds, "list")){
    plate_errors_zUMIs_output <- c(plate_errors_zUMIs_output, paste("The RDS file for", plate_name,"is not of the correct type (other than a list) please check zUMIs output results."))
  }
  # Load the gene names file (tab-delimited)
  gene_names_df <- read.table(gene_names_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Checks that the object loaded is a data.frame
  if (!is(gene_names_df, "data.frame")){
    plate_errors_zUMIs_output <- c(plate_errors_zUMIs_output, paste("ERROR : Gene_names.txt is not of type data.frame for", plate_name ,"please verify the data."))
  }
  # Checks that the "gene_id" column is of character type
  if (!is(gene_names_df$gene_id, "character")){
    plate_errors_zUMIs_output <- c(plate_errors_zUMIs_output, paste0("ERROR : gene_id column is not of type character in the Gene_names.txt for ", plate_name, "."))
  }
  # Checks that the "gene_name" column is of character type
  if (!is(gene_names_df$gene_name, "character")){
    plate_errors_zUMIs_output <- c(plate_errors_zUMIs_output, paste0("ERROR : gene_name column is not of type character in the Gene_names.txt for ", plate_name, "."))
  }
}

# Display errors if any were found 
if (length(plate_errors_zUMIs_output) > 0) {
  # Creates a header
  header <- "\n============== ERROR(S) DETECTED IN zUMIs OUTPUT FROM SCRIPT '01_prepareData' ==============\n"
  # Numbering of all errors found in this section
  error_list <- paste0(seq_along(plate_errors_zUMIs_output), ". ", plate_errors_zUMIs_output, collapse = "\n")
  # Final message
  full_message <- paste0(header, error_list, "\n============================================================================================\n")
  # Displays and stops
  stop(full_message)
}

# Environment cleanup: removes temporary objects
rm(plate_directory)
rm(rds_file)
rm(gene_names_file)
rm(file_rds)
rm(gene_names_df)

###############################################
# ## Verification for 03_computeERCCUMIPercent
###############################################

# Nothing to be done

#############################################
# ## Verification for 04_computeERCCAccuracy
#############################################

# We store the various warning in variables in order to display them all together instead of stopping the code at each warning.
ERCC_EXPECTED_DF <- read.csv(PATH_ERCC_CONCENTRATION)

# Initialize an error list
plate_errors_ERCC <- c()

# Checks that the ERCC_EXPECTED_DF object is a data.frame
if (!is.data.frame(ERCC_EXPECTED_DF)){
  plate_errors_ERCC <- c(plate_errors_ERCC, "ERROR : The ERCC file (in 01_Reference/00_Experiment), is not a dataframe. Check that this is the correct file.")
}
# Checks that the ‘ERCC_ID’ column exists in the data.frame
if (!"ERCC_ID" %in% colnames(ERCC_EXPECTED_DF)){
  plate_errors_ERCC <- c(plate_errors_ERCC, paste("ERROR : There is no column named 'ERCC_ID' in the barcode_well file. \nDirectory :", PATH_ERCC_CONCENTRATION))
}
# Checks that the ‘Concentration_attomoles_ul’ column exists in the data.frame
if (!"Concentration_attomoles_ul" %in% colnames(ERCC_EXPECTED_DF)){
  plate_errors_ERCC <- c(plate_errors_ERCC, paste("ERROR : There is no column named 'Concentration_attomoles_ul' in the barcode_well file. \nDirectory :", PATH_ERCC_CONCENTRATION))
}
# Checks that the ‘ERCC_ID’ column 
if (!is(ERCC_EXPECTED_DF$ERCC_ID, "character")){
  plate_errors_ERCC <- c(plate_errors_ERCC, "ERROR : BarcodeSequence column is not of type character in the ERRC csv file.")
}
# Checks that the ‘Concentration_attomoles_ul’ column is numeric
if (!is(ERCC_EXPECTED_DF$Concentration_attomoles_ul, "numeric")){
  plate_errors_ERCC <- c(plate_errors_ERCC, "ERROR : BarcodeSequence column is not of type numeric in the ERRC csv file.")
}

# Display errors if any were found 
if (length(plate_errors_ERCC) > 0) {
  # Creates a header
  header <- "\n============== ERROR(S) DETECTED IN ERCC FILE FROM SCRIPT '01_prepareData' ==============\n"
  # Numbering of all errors found in this section
  error_list <- paste0(seq_along(plate_errors_ERCC), ". ", plate_errors_ERCC, collapse = "\n")
  # Final message
  full_message <- paste0(header, error_list, "\n=========================================================================================\n")
  # Displays and stops
  stop(full_message)
}

#######################################
# ## Verification for 05_computeRNAUMI
#######################################

# Nothing to be done

##############################################
# ## Verification for 06_computeFeatureNumber
##############################################

# Nothing to be done

##########################################
# ## Verification for 07_computeIndexSort
##########################################

# We store the various warning in variables in order to display them all together instead of stopping the code at each warning.
plate_errors_indexsort <- c()
if (PARAMS_INDEXSORT == TRUE){
  # Initialize an error list
  SKIPPED_INDEXSORT_PLATES <- c()
  columns_names <- ""
  index_sort_colnames_number <- c()
  global_column_indexsort <- NULL
  acc = 0
  for (plate_name in levels(PLATES_LIST)) {
    path_index_sort_file = paste(PATH_EXPERIMENT_RAWDATA, "/01_IndexSort/", plate_name, "_indexsort.csv", sep = "")
    
    # Verify that the IndexSort file is not empty.
    if (file.info(path_index_sort_file)$size > 1) {
      # Detects the type of spacer in the csv file (based only on “,” or “;”)
      separator <- detect_sep_csv(path_index_sort_file)
      # Read the file with the correct separator present in the csv file.
      index_sort_df <- read.csv(path_index_sort_file, sep = separator)
      
      # We've come to store the names of the columns present in the plates, in case they differ (the aim is to add 0s to the others, to continue the analysis).
      current_columns <- colnames(index_sort_df)
      # Check if the plate has missing columns or more than the global list
      missing_in_plate <- setdiff(global_column_indexsort, current_columns)
      new_in_plate     <- setdiff(current_columns, global_column_indexsort)
      
      if (plate_name == PLATES_LIST[1]){
        # Add new columns to the global list
        for (col in new_in_plate) {
          global_column_indexsort <- c(global_column_indexsort, col)
        }
      } else {
        if (length(new_in_plate) > 0) {
          acc =+ 1
          for (col in new_in_plate) {
            global_column_indexsort <- c(global_column_indexsort, col)
          }
        }
        if (length(missing_in_plate) > 0){
          acc =+ 1
        }
      }
      if (acc > 0){
        cat("\n \nDifferences were detected between the columns of the IndexSort plates. The analysis will continue by adding the missing column(s) between the plates with 0's.\n \n ")
      }
      
      for (column_name_in_not_fluorescent in NOT_FLUORESCENT){
        if (!column_name_in_not_fluorescent %in% colnames(index_sort_df)){
          plate_errors_indexsort <- c(plate_errors_indexsort, paste("ERROR : There is no column named", column_name_in_not_fluorescent,"(the name here is the one provided in the config user) in the indexsort file for plate :", plate_name))
        }
      }
      
      if (length(columns_names)>1){
        for (name_column in colnames(index_sort_df)){
          if (!name_column %in% columns_names){
            plate_errors_indexsort <- c(plate_errors_indexsort, paste0("ERROR : The column name between ",plate_name_verification, " (Column : ", setdiff(columns_names, colnames(index_sort_df)), ") and ", plate_name, " (Column : ", setdiff(colnames(index_sort_df), columns_names), ") are different in indexsort file. Please correct the column name between the file."))
          }
        } 
      }
      
      
      if (!COLUMN_HEADER_WELL_ID %in% colnames(index_sort_df)){
        plate_errors_indexsort <- c(plate_errors_indexsort, paste("ERROR : There is no column named 'WellID' in the indexsort file for plate :", plate_name))
      }
      
      if (!is(index_sort_df[[COLUMN_HEADER_WELL_ID]], "character")){
        plate_errors_indexsort <- c(plate_errors_indexsort, paste("ERROR : The column WellID is not of type character in the indexsort file for plate :", plate_name))
      }
      
      # Check the length of index_sort_df$Well.ID
      well_count = length(index_sort_df[[COLUMN_HEADER_WELL_ID]])
      if (well_count < 96) {
        cat(paste("WARNING : Some wells are not supplied in plate:", plate_name, "--- Total number of wells present:", well_count, "\n"))
      } else if (well_count > 96) {
        plate_errors_indexsort <- c(plate_errors_indexsort, paste0("ERROR : Too many wells/lines in plate index sort file : ", plate_name, ".\nTotal number of wells/lines : ", well_count))
      }
      
      # Check if the plates have the same column names
      if( length( index_sort_colnames_number) >0){
        if( index_sort_colnames_number != ncol( index_sort_df) ){
          plate_errors_indexsort <- c(plate_errors_indexsort, paste( "ERROR : The column names of the plates are not equals. For the plate", plate_name_verification,"there is", index_sort_colnames_number, "column and for the plate", plate_name, "there is", ncol( index_sort_df), "column."))
        }
      }else{
        index_sort_colnames_number = ncol( index_sort_df)
      }
      columns_names = colnames(index_sort_df)
      plate_name_verification = plate_name
    } else {
      # If the IndexSort file is empty, we store it in a list for further analysis (to avoid bugs).
      SKIPPED_INDEXSORT_PLATES <- c(SKIPPED_INDEXSORT_PLATES, plate_name)
      warning("The FACS data (IndexSort) of plate", plate_name, "are empty. However, the analysis continues without this data.")
    }
  }
  
  # Affiche les erreurs s'il y en a
  if (length(plate_errors_indexsort) > 0) {
    # Creates a header
    header <- "\n============== ERROR(S) DETECTED IN INDEXSORT FILE FROM SCRIPT '01_prepareData' ==============\n"
    # Numbering of all errors found in this section
    error_list <- paste0(seq_along(plate_errors_indexsort), ". ", plate_errors_indexsort, collapse = "\n")
    # Final message
    full_message <- paste0(header, error_list, "\n==============================================================================================\n")
    # Displays and stops
    stop(full_message)
  }
  
  # Environment cleanup: removes temporary objects
  rm(path_index_sort_file)
  rm(index_sort_df)
  rm(column_name_in_not_fluorescent)
  rm(columns_names)
  rm(plate_name_verification)
  rm(index_sort_colnames)
  rm(separator)
  rm(acc)
  rm(new_in_plate)
  rm(missing_in_plate)
  rm(current_columns)
}

#########################################
# ## Verification for 08_displayMetaData
#########################################

# We store the various warning in variables in order to display them all together instead of stopping the code at each warning.
plate_errors_metadata <- c()
if (PARAMS_METADATA){
  # Initialize an error list
  SKIPPED_PLATES_METADATA <- c()
  
  ## Check before to run the process.
  # Check names between config file and sheets and library (both present in GSF).
  sheets_gsf = excel_sheets(PATH_GSF)
  plate_sheet_verify = read_excel(PATH_GSF, sheet = "General Info")
  plate_sheet_verify = as.data.frame(plate_sheet_verify)
  for (plate_name in levels(PLATES_LIST)){
    if (!any(plate_sheet_verify == plate_name)){
      cat(paste("WARNING : The plate", plate_name, "given in the config file doesn't match the one given in the 'General Info' tab of the GSF."))
    }
    if (!any(plate_name == sheets_gsf)){
      plate_errors_metadata <- c(plate_errors_metadata, paste("ERROR : The plate", plate_name, "given in the config file does not match the GSF sheets."))
    }
  }
  
  # Check if the column name are the same in all the different plate file
  columns_names_previous_plate <- ""
  for (plate_name in levels(PLATES_LIST)){
    plate_sheet = read_excel(PATH_GSF, sheet = plate_name)
    
    if (length(plate_sheet) == 0){
      # If the IndexSort file is empty, we store it in a list for further analysis (to avoid bugs).
      SKIPPED_PLATES_METADATA <- c(SKIPPED_PLATES_METADATA, plate_name)
      warning("The sheet Metadata of plate", plate_name, "are empty. However, the analysis continues without this data.")
    } else {
      if (length(columns_names_previous_plate)>1){
        for (name_column in colnames(plate_sheet)){
          if (!name_column %in% columns_names_previous_plate){
            plate_errors_metadata <- c(plate_errors_metadata, paste("ERROR in plate ", plate_name ," : The column name between", plate_name_verification, " (Column : ", setdiff(columns_names_previous_plate, colnames(index_sort_df)), ") and", plate_name, " : The column name between",plate_name_verification, " (Column : ", setdiff(colnames(index_sort_df), columns_names_previous_plate), ") in the GSF file (sheet) are different. Correct the name column between them"))
          }
        } 
      }
      columns_names_previous_plate = colnames(plate_sheet)
      plate_name_verification = plate_name
      if (!COLUMN_HEADER_WELL_ID %in% colnames(plate_sheet)){
        plate_errors_metadata <- c(plate_errors_metadata, paste0("ERROR in plate ", plate_name ," : WellID (", COLUMN_HEADER_WELL_ID, ") column is not found in the column name of the GSF file.\nCheck the variable in analysisParams.R or you're column in the GSF file."))
      }
    } 
  }
  
  # Display errors if any were found 
  if (length(plate_errors_metadata) > 0) {
    # Creates a header
    header <- "\n============== ERROR(S) DETECTED IN METADATA FILE (gsf) FROM SCRIPT '01_prepareData' ==============\n"
    # Numbering of all errors found in this section
    error_list <- paste0(seq_along(plate_errors_metadata), ". ", plate_errors_metadata, collapse = "\n")
    # Final message
    full_message <- paste0(header, error_list, "\n===================================================================================================\n")
    # Displays and stops
    stop(full_message)
  }
  
  # Environment cleanup: removes temporary objects
  rm(sheets_gsf)
  rm(plate_sheet_verify)
  rm(plate_name_verification)
}

#######################################
# ## Verification for 09_computeTrust4
#######################################

# We store the various warning in variables in order to display them all together instead of stopping the code at each warning.
if (PARAMS_BCR | PARAMS_TCR){
  # Initialize an error list
  plate_errors_trust4 <- c()
  
  columns_names = ""
  for (plate_name in levels(PLATES_LIST)) {
    plate_directory <- file.path(PATH_TRUST4_OUTPUT, plate_name)
    tsv_air_file <- list.files(plate_directory, pattern = paste0("^", plate_name, "_barcode_airr.tsv$"), full.names = TRUE)
    if (!file.exists(tsv_air_file)) {
      plate_errors_trust4 <- c(plate_errors_trust4, paste("ERROR : There is no barcode_airr.tsv file (TRUST4 output) for", plate_name))
    }
    file_air_tsv <- read.csv(tsv_air_file, sep = '\t')
    if (length(file_air_tsv) <= 1){
      plate_errors_trust4 <- c(plate_errors_trust4, paste("ERROR : Le fichier barcode_airr.tsv file (TRUST4 output) for", plate_name, "is empty"))
    }
    if (!is(file_air_tsv, "data.frame")){
      plate_errors_trust4 <- c(plate_errors_trust4, paste("ERROR : The barcode_airr.tsv file (TRUST4 output) is not a dataframe for", plate_name))
    }
    if (length(columns_names)>1){
      for (name_column in colnames(file_air_tsv)){
        if (!name_column %in% columns_names){
          plate_errors_trust4 <- c(plate_errors_trust4, paste("ERROR : The column name between",plate_name_verification, "and", plate_name,"in the barcode_airr.tsv file (TRUST4 output) are different. Correct the name column between them."))
        }
      } 
    }
    columns_names = colnames(file_air_tsv)
    plate_name_verification = plate_name
    
    if (!COLUMN_HEADER_C_CALL %in% colnames(file_air_tsv)){
      plate_errors_trust4 <- c(plate_errors_trust4, paste0("ERROR : C_CALL (named :", COLUMN_HEADER_C_CALL, ") column is not found in the column name of barcode_airr.tsv file (TRUST4 output).\nCheck the variable in analysisParams.R or you're column in the tsv file."))
    }
    if (!COLUMN_HEADER_V_CALL %in% colnames(file_air_tsv)){
      plate_errors_trust4 <- c(plate_errors_trust4, paste0("ERROR : V_CALL (named :", COLUMN_HEADER_V_CALL, ") column is not found in the column name of barcode_airr.tsv file (TRUST4 output).\nCheck the variable in analysisParams.R or you're column in the tsv file."))
    }
    if (!COLUMN_HEADER_J_CALL %in% colnames(file_air_tsv)){
      plate_errors_trust4 <- c(plate_errors_trust4, paste0("ERROR : J_CALL (named :", COLUMN_HEADER_J_CALL, ") column is not found in the column name of barcode_airr.tsv file (TRUST4 output).\nCheck the variable in analysisParams.R or you're column in the tsv file."))
    }
    if (!COLUMN_HEADER_PRODUCTIVE %in% colnames(file_air_tsv)){
      plate_errors_trust4 <- c(plate_errors_trust4, paste0("ERROR : PRODUCTIVE (named :", COLUMN_HEADER_PRODUCTIVE, ") column is not found in the column name of barcode_airr.tsv file (TRUST4 output).\nCheck the variable in analysisParams.R or you're column in the tsv file."))
    }
    if (!COLUMN_HEADER_CELL_ID %in% colnames(file_air_tsv)){
      plate_errors_trust4 <- c(plate_errors_trust4, paste0("ERROR : COLUMN_HEADER_CELL_ID (named :", COLUMN_HEADER_CELL_ID, ") column is not found in the column name of barcode_airr.tsv file (TRUST4 output).\nCheck the variable in analysisParams.R or you're column in the tsv file."))
    }
  }
  
  # Display errors if any were found 
  if (length(plate_errors_trust4) > 0) {
    # Creates a header
    header <- "\n============== ERROR(S) DETECTED IN OUTPUT OF TRUST4 FROM SCRIPT '01_prepareData' ==============\n"
    # Numbering of all errors found in this section
    error_list <- paste0(seq_along(plate_errors_trust4), ". ", plate_errors_trust4, collapse = "\n")
    # Final message
    full_message <- paste0(header, error_list, "\n================================================================================================\n")
    # Displays and stops
    stop(full_message)
  }
  
  # Environment cleanup: removes temporary objects
  rm(plate_directory)
  rm(tsv_air_file)
  rm(plate_name)
  rm(columns_names)
  rm(name_column)
  rm(file_air_tsv)
  rm(plate_name_verification)
  rm(plate_errors_zUMIs_output)
  rm(plate_errors_ERCC)
  rm(plate_errors_cell_barcode_well)
  rm(plate_errors_indexsort)
  rm(plate_errors_metadata)
  rm(plate_errors_trust4) 
}