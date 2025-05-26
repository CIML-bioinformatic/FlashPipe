## ####################################################
## This script aims to re organise the output of zUMIs, 
## for the comprehension and the analysis.
## He will create 3 csv file for each plate.
## ####################################################

## @knitr format_count_table

# Initialize empty lists to store the dataframes for ERCC and RNA count data for each plate
ERCC_count_df_list = list()
RNA_count_df_list = list()

# Loop over the list of plates to process each one
for (plate_name in PLATES_LIST) {
  # Define the directory path for the current plate
  plate_directory <- file.path(PATH_ZUMIS_OUTPUT, plate_name)
  
  # Check if the RDS file for the current plate exists
  rds_file <- file.path(plate_directory, "zUMIs_output", "expression", paste0(plate_name, ".dgecounts.rds"))
  
  # Check if the gene_names.txt file for the current plate exists
  gene_names_file <- file.path(plate_directory, "zUMIs_output", "expression", paste0(plate_name, ".gene_names.txt"))
  
  # Load the RDS file
  file_rds <- readRDS(rds_file)
  
  # Load the gene names file (tab-delimited)
  gene_names_df <- read.table(gene_names_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Access the 'all' data under the umicount$exon object from the RDS file
  data_all <- file_rds$umicount$exon$all
  
  # Convert sparse matrix to dense matrix
  data_all_dense <- as.matrix(data_all)
  
  # Create a dataframe for ERCC rows only
  df_ERCC <- data_all_dense[grepl("ERCC", rownames(data_all_dense)), ]
  
  # Create a dataframe for gene data excluding ERCC rows
  gene_id_data <- data_all_dense[!grepl("ERCC", rownames(data_all_dense)), ]
  
  # 3rd Dataset: gene_name_data.csv with gene_ids replaced by gene_names
  # Create a mapping from gene_id to gene_name
  gene_id_to_name <- setNames(gene_names_df$gene_name, gene_names_df$gene_id)
  
  # Replace the gene_ids in the rownames with the corresponding gene_names
  rownames(data_all_dense) <- gene_id_to_name[rownames(data_all_dense)]
  
  # Create a dataframe with the new gene names as rownames
  gene_name_data <- data.frame(data_all_dense)
  
  # Replace column names with WellID using the barcode-to-well mapping from the cell_barcode_well_df
  # Create a mapping from BarcodeSequence to WellID
  barcode_to_well <- setNames(CELL_BARCODE_WELL_DF[[WELL_ID]], CELL_BARCODE_WELL_DF$BarcodeSequence)
  
  # Update the column names in gene_name_data and df_ERCC to use WellID
  colnames(gene_name_data) <- barcode_to_well[colnames(gene_name_data)]
  colnames(df_ERCC) <- barcode_to_well[colnames(df_ERCC)]
  
  # Define the output directory for the current plate
  plate_output_directory <- file.path(PATH_ANALYSIS_OUTPUT, plate_name)
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(plate_output_directory)) {
    dir.create(plate_output_directory)
  }
  
  # Store the ERCC and RNA dataframes in their respective lists for the current plate
  ERCC_count_df_list[[ plate_name]] = df_ERCC
  RNA_count_df_list[[ plate_name]] = gene_name_data
  
  # Save the ERCC, RNA gene ID, and RNA gene name data as CSV files in the plate output directory
  ERCC_count_table_output_file = file.path(plate_output_directory, paste0(plate_name, "_ERCC_erccId_wellName_count.csv"))
  write.csv(df_ERCC, ERCC_count_table_output_file)
  write.csv(gene_id_data, file.path(plate_output_directory, paste0(plate_name, "_RNA_geneId_wellBarcode_count.csv")))
  write.csv(gene_name_data, file.path(plate_output_directory, paste0(plate_name, "_RNA_geneSymbol_wellName_count.csv")))
}

