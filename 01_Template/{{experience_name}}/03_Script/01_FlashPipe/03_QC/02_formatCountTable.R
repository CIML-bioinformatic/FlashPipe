## ####################################################
## This script aims to re organise the output of zUMIs, 
## and Trust4 for the comprehension and the analysis.
## He will create 3 csv file for each plate.
## And prepare list of data for the quality control.
## ####################################################

## @knitr format_count_table

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# General description of the code :
# 1. Preprocess Trust4 output data for analysis
# 2. Preprocessing of zUMIs output data for analysis with preprocessed file output
# 3. Recover zUMIs data for Seurat object creation
# 4. Retrieve metadata information (gsf) for seurat object metadata
# 5. Retrieve sort index information (FACS data) for seurat object metadata
# 6. Create seurat object with sorted data (zUMIs and metadata)
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

########################################################
# ## 1. Preprocess Trust4 output data for analysis #####
########################################################

# Read the Trust 4 output file and sort the TCR and BCR data into different lists for use in quality control later.
if (PARAMS_BCR | PARAMS_TCR){
  # Initialize lists that will be used in subsequent quality checks, simplifying analysis by storing data that has already been checked.
  BCR_df_list = list()
  TCR_df_list = list()
  for (plate_name in PLATES_LIST) {
    # Reading output tsv files from Airrflow (if params selected)
    if (PARAMS_TOOLS_BCR_TCR_ANALYSIS == "airrflow"){
      PATH_AIRRFLOW_OUTPUT_TEMP <- sub("/+$", "", PATH_AIRRFLOW_OUTPUT)
      plate_directory <- file.path(PATH_AIRRFLOW_OUTPUT_TEMP, "trust4")
    }
    # Reading output tsv files from trust4
    if (PARAMS_TOOLS_BCR_TCR_ANALYSIS == "trust4"){
      plate_directory <- file.path(PATH_TRUST4_OUTPUT, plate_name)
    }else{
      plate_directory <- PATH_TRUST4_OUTPUT
    }
    tsv_air_file <- list.files(plate_directory, pattern = paste0("^", plate_name, "_barcode_airr.tsv$"), full.names = TRUE)
    file_air_tsv <- read.csv(tsv_air_file, sep = '\t')
    
    if (PARAMS_BCR){
      # Separate information from data containing “IG” which corresponds to the BCR in the “v_call” column
      bcr_df <- file_air_tsv[grepl("IG", file_air_tsv[[COLUMN_HEADER_V_CALL]]), ]
      # Storing sorted Trust4 data from BCRs in the list.
      BCR_df_list[[ plate_name]] = bcr_df
    }
    
    if (PARAMS_TCR){
      # Separate information from data containing “TR” which corresponds to the TCR in the “v_call” column
      tcr_df <- file_air_tsv[grepl("TR", file_air_tsv[[COLUMN_HEADER_V_CALL]]), ] 
      # Stores sorted Trust4 data from TCRs in the list.
      TCR_df_list[[ plate_name]] = tcr_df
    }
  }
}

###########################################################################################
# ## 2. Preprocessing of zUMIs output data for analysis with preprocessed file output #####
###########################################################################################

# Initialize empty lists to store the dataframes for ERCC and RNA count data for each plate 
ERCC_count_df_list = list()
RNA_count_df_list = list()
# Create an empty dataframe to store the RNA UMI results for the Seurat object.
df_ARN_object_seurat <- NULL

# Loop that recovers data from zUMIs to separate them from ERCC and RNA data, and simplify the data for future analysis.
# But also to retrieve zUMis data from UMI RNAs for the Seurat object.
for (plate_name in PLATES_LIST) {
  # Define the directory path for the current plate
  plate_directory <- file.path(PATH_ZUMIS_OUTPUT, plate_name)
  
  # Check if the RDS file for the current plate exists
  rds_file <- file.path(plate_directory, "zUMIs_output", "expression", paste0(plate_name, ".dgecounts.rds"))
  # Load the RDS file
  file_rds <- readRDS(rds_file)
  
  # Check if the gene_names.txt file for the current plate exists
  gene_names_file <- file.path(plate_directory, "zUMIs_output", "expression", paste0(plate_name, ".gene_names.txt"))
  # Load the gene names file (tab-delimited)
  gene_names_df <- read.table(gene_names_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Access counting data from the RDS file
  data_all <- file_rds$umicount$exon$all
  data_all_dense <- as.matrix(data_all)
  
  # Create a dataframe for ERCC rows only
  df_ERCC <- data_all_dense[grepl("ERCC", rownames(data_all_dense)), ]
  
  # Create a dataframe for gene data excluding ERCC rows
  gene_id_data <- data_all_dense[!grepl("ERCC", rownames(data_all_dense)), ]
  
  # Define the output directory for the current plate
  plate_output_directory <- file.path(PATH_ANALYSIS_OUTPUT, plate_name, '02_formatCountTable')
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(plate_output_directory)) {
    dir.create(plate_output_directory, recursive = TRUE)
  }
  
  # Save the wellbarcode and geneid before the modification of the gene name.
  write.csv(gene_id_data, file.path(plate_output_directory, paste0(plate_name, "_RNA_geneId_wellBarcode_count.csv")))
  
  # 3rd Dataset: gene_name_data.csv with gene_ids replaced by gene_names
  # Create a mapping from gene_id to gene_name
  gene_id_to_name <- setNames(gene_names_df$gene_name, gene_names_df$gene_id)
  
  # Store and replace Ensembl names with gene names.
  new_row_names = gene_id_to_name[rownames(gene_id_data)]
  names( new_row_names) = rownames(gene_id_data)
  # Check data and presence of NA. If present, we keep the original name.
  if (any( is.na( new_row_names))){
    na_index_set = which( is.na( new_row_names))
    new_row_names[ na_index_set] = names( new_row_names)[ which( is.na( new_row_names))] 
  }
  
  # Replace the gene_ids in the rownames with the corresponding gene_names
  rownames(gene_id_data) <- new_row_names
  
  # Create a dataframe with the new gene names as rownames
  gene_name_data <- data.frame(gene_id_data)
  # Copy gene_id_data values and modify them to suit Seurat object
  gene_name_data_object_seurat <- as.data.frame(gene_id_data)
  
  # Replace column names with WellID using the barcode-to-well mapping from the cell_barcode_well_df
  # Create a mapping from BarcodeSequence to WellID
  barcode_to_well <- setNames(CELL_BARCODE_WELL_DF[[COLUMN_HEADER_WELL_ID]], CELL_BARCODE_WELL_DF$BarcodeSequence)
  
  # Update the column names in gene_name_data and df_ERCC to use WellID
  colnames(gene_name_data) <- barcode_to_well[colnames(gene_name_data)]
  colnames(df_ERCC) <- barcode_to_well[colnames(df_ERCC)]
  
  # Store the ERCC and RNA dataframes in their respective lists for the current plate
  ERCC_count_df_list[[ plate_name]] = df_ERCC
  RNA_count_df_list[[ plate_name]] = gene_name_data
  
  # Save the ERCC, and RNA gene name data as CSV files in the plate output directory
  ERCC_count_table_output_file = file.path(plate_output_directory, paste0(plate_name, "_ERCC_erccId_wellName_count.csv"))
  write.csv(df_ERCC, ERCC_count_table_output_file)
  write.csv(gene_name_data, file.path(plate_output_directory, paste0(plate_name, "_RNA_geneSymbol_wellName_count.csv")))
  
  ###########################################################
  # ## 3. Recover zUMIs data for Seurat object creation ##### 
  ###########################################################
  
  # Rename columns with the name of the plate and well together (Example: Plate_1_B7) for the object seurat
  colnames(gene_name_data_object_seurat) <- paste0(plate_name, "_", barcode_to_well[colnames(gene_name_data_object_seurat)])
  
  # Merge with global dataframe
  if (is.null(df_ARN_object_seurat)) {
    df_ARN_object_seurat <- gene_name_data_object_seurat
  } 
  else {
    # Merge by aligning rows by gene (row.names), keeping all columns
    df_ARN_object_seurat <- merge(df_ARN_object_seurat, gene_name_data_object_seurat,
                                  by = "row.names", all = TRUE)
    
    # Restore rownames after merging
    rownames(df_ARN_object_seurat) <- df_ARN_object_seurat$Row.names
    df_ARN_object_seurat$Row.names <- NULL
  }
}

# Checks the presence of ERCC data in the various plates
ERCC_lengths <- sapply(ERCC_count_df_list[PLATES_LIST], length)
ERCC <- any(ERCC_lengths > 0)
if (!ERCC) {
  warning("No ERCC was found in the RDS files of the individual plates. The analysis continues without ERCC data.")
}

# Replace the NAs present in the Seurat object dataframe with 0.
# This remains coherent because these are genes that have not been found, and therefore not expressed (even if the genes are not the base plates).
df_ARN_object_seurat[is.na(df_ARN_object_seurat)] <- 0

############################################################################
# ## 4. Retrieve metadata information (gsf) for seurat object metadata #####
############################################################################
if (PARAMS_METADATA){
  df_metadata_object_seurat <- NULL
  for (plate_name in PLATES_LIST) {
    # Skip plates marked as empty
    if (plate_name %in% SKIPPED_PLATES_METADATA) {
      next
    }
    plate_sheet <- read_excel(PATH_GSF, sheet = plate_name)
    # Delete rows where all columns are empty
    well_col_index <- which(names(plate_sheet) == COLUMN_HEADER_WELL_ID)
    other_cols <- plate_sheet[, -well_col_index, drop = FALSE]
    keep_rows <- apply(other_cols, 1, function(row) any(!is.na(row) & row != ""))
    plate_sheet <- plate_sheet[keep_rows, ]
    
    # Change tible to data.frame
    plate_sheet = as.data.frame( plate_sheet)
    # Create a new column containing the new column names( Plate_Name_WellID), then set the rows as columns for the seurat object.
    plate_sheet$Plate_Well_ID <- paste0(plate_name, "_", plate_sheet[[COLUMN_HEADER_WELL_ID]])
    row.names(plate_sheet) = plate_sheet$Plate_Well_ID
    # Delete the WellID column and Plate_Well_ID, which no longer provides information for the seurat object.
    plate_sheet[[COLUMN_HEADER_WELL_ID]] = NULL
    plate_sheet$Plate_Well_ID = NULL
    # We transpose the data used to create the seurat object.
    df_metadata_object_seurat_transform = data.table::transpose(plate_sheet)
    names( df_metadata_object_seurat_transform) = row.names( plate_sheet)
    row.names( df_metadata_object_seurat_transform) = names( plate_sheet)
    
    # Merge with global dataframe
    if (is.null(df_metadata_object_seurat)) {
      df_metadata_object_seurat <- df_metadata_object_seurat_transform
    } else {
      # Merge by aligning rows by gene (row.names), keeping all columns
      df_metadata_object_seurat <- merge(df_metadata_object_seurat, df_metadata_object_seurat_transform,
                                         by = "row.names", all = TRUE)
    }
  } 
  row.names(df_metadata_object_seurat) = df_metadata_object_seurat$Row.names
  df_metadata_object_seurat$Row.names = NULL
}

####################################################################################
# ## 5. Retrieve sort index information (FACS data) for seurat object metadata #####
####################################################################################
if (PARAMS_INDEXSORT){
  df_indexsort_object_seurat <- NULL
  for (plate_name in PLATES_LIST) {
    # Skip plates marked as empty
    if (plate_name %in% SKIPPED_INDEXSORT_PLATES) {
      next
    }
    path_index_sort_file <- file.path(PATH_EXPERIMENT_RAWDATA, "01_IndexSort", paste0(plate_name, "_indexsort.csv"))
    # Automatic separator detection
    separator <- detect_sep_csv(path_index_sort_file)
    # Read the file with the correct separator
    index_sort_df <- read.csv(path_index_sort_file, sep = separator, stringsAsFactors = FALSE)
    
    # Filter empty rows (except well column)
    well_col_index <- which(names(index_sort_df) == COLUMN_HEADER_WELL_ID)
    other_cols <- index_sort_df[ , -well_col_index, drop = FALSE]
    keep_rows <- apply(other_cols, 1, function(row) any(!is.na(row) & row != ""))
    index_sort_df <- index_sort_df[keep_rows, ]
    
    # Change tible to data.frame
    index_sort_df = as.data.frame( index_sort_df)
    # Create a new column containing the new column names( Plate_Name_WellID), then set the rows as columns for the seurat object.
    index_sort_df$Plate_Well_ID <- paste0(plate_name, "_", index_sort_df[[COLUMN_HEADER_WELL_ID]])
    row.names(index_sort_df) = index_sort_df$Plate_Well_ID
    # Delete the WellID column and Plate_Well_ID, which no longer provides information for the seurat object.
    index_sort_df[[COLUMN_HEADER_WELL_ID]] = NULL
    index_sort_df$Plate_Well_ID = NULL
    # We transpose the data used to create the seurat object.
    df_indexsort_object_seurat_transform = data.table::transpose(index_sort_df)
    names( df_indexsort_object_seurat_transform) = row.names( index_sort_df)
    row.names( df_indexsort_object_seurat_transform) = names( index_sort_df)
    
    # Merge with global dataframe
    if (is.null(df_indexsort_object_seurat)) {
      df_indexsort_object_seurat <- df_indexsort_object_seurat_transform
    } else {
      # Merge by aligning rows by gene (row.names), keeping all columns
      df_indexsort_object_seurat <- merge(df_indexsort_object_seurat, df_indexsort_object_seurat_transform,
                                          by = "row.names", all = TRUE)
    }
  } 
  row.names(df_indexsort_object_seurat) = df_indexsort_object_seurat$Row.names
  df_indexsort_object_seurat$Row.names = NULL
}

########################################################################
# ## 6. Create seurat object with sorted data (zUMIs and metadata) #####
########################################################################

## Select and create the metadata type for the seurat object, with or without index sort, and with or without metadata.
# Merge the 2 dataframes (containing fluorescence (from FACS data) and metadata) into 1, if the option have metatada and indexsort file set to TRUE.
if (PARAMS_METADATA & PARAMS_INDEXSORT){
  
  # Make sure columns are in the same order
  common_columns <- intersect(colnames(df_metadata_object_seurat), colnames(df_indexsort_object_seurat))
  
  # Reorder columns and rbind
  DATAFRAME_FOR_METADATA_SEURAT <- rbind(df_metadata_object_seurat[, common_columns],
                                         df_indexsort_object_seurat[, common_columns])
  
  # Environment cleanup: removes temporary objects
  rm(df_metadata_object_seurat)
  rm(df_metadata_object_seurat_transform)
  rm(path_index_sort_file)
  rm(separator)
  rm(index_sort_df)
  rm(df_indexsort_object_seurat_transform)
  rm(df_indexsort_object_seurat)
  rm(common_columns)
}
# Create the seurat object with only Metadata, depends on the option selected.
if (PARAMS_METADATA &! PARAMS_INDEXSORT){
  DATAFRAME_FOR_METADATA_SEURAT = df_metadata_object_seurat
  
  # Environment cleanup: removes temporary objects
  rm(df_metadata_object_seurat)
}
# Create the seurat object with only IndexSort, depends on the option selected.
if (PARAMS_INDEXSORT &! PARAMS_METADATA){
  DATAFRAME_FOR_METADATA_SEURAT = df_indexsort_object_seurat
  
  # Environment cleanup: removes temporary objects
  if( exists( "path_index_sort_file")) rm("path_index_sort_file")
  if( exists( "separator")) rm("separator")
  if( exists( "index_sort_df")) rm("index_sort_df")
  if( exists( "df_indexsort_object_seurat_transform")) rm("df_indexsort_object_seurat_transform")
  if( exists( "df_indexsort_object_seurat")) rm("df_indexsort_object_seurat")
}

# Seurat object create with no normalized or scaled data, and the cell and feature minimums are set to 0, to leave the choice open afterwards.
FINAL_SEURAT_OBJECT = CreateSeuratObject( df_ARN_object_seurat, project = EXPERIMENT_NAME, assay = "RNA",
                                          min.cells = 0, min.features = 0, names.field = 1,
                                          names.delim = "_", meta.data = DATAFRAME_FOR_METADATA_SEURAT)

PATH_DIRECTORY_ANALYSIS <- file.path(PATH_FLASHPIPE_OUTPUT, '04_Analysis')
# Registers the seurat object.
saveRDS( FINAL_SEURAT_OBJECT, file = file.path( PATH_DIRECTORY_ANALYSIS, "complete_seurat_object.RDS"))

# Environment cleanup: removes temporary objects
if( exists( "FINAL_SEURAT_OBJECT")) rm("FINAL_SEURAT_OBJECT")
if( exists( "plate_sheet")) rm("plate_sheet")
if( exists( "df_ARN_object_seurat")) rm("df_ARN_object_seurat")
if( exists( "gene_name_data_object_seurat")) rm("gene_name_data_object_seurat")
if( exists( "well_col_index")) rm("well_col_index")
if( exists( "other_cols")) rm("other_cols")
if( exists( "keep_rows")) rm("keep_rows")
