## ####################################################
## This script aims to load and send back plot and data
## from the trust4 output, for the BCR analysis.
## ####################################################

## @knitr computeTrust4_bcr

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# General description of the code :
# 1. Pre-processing of BCR data provided by Trust4 for Isotypes
# 2. General graphics of all Isotype plates
# 3. Plate-by-plate graphs for Isotypes
# 4. Pre-processing of BCR data supplied by Trust4 for string composition
# 5. General graph of all plates for string composition
# 6. Plate-by-plate graph for chain composition
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# ######################################################################
# ## 1. Pre-processing of BCR data provided by Trust4 for Isotypes #####
# ######################################################################

if (PARAMS_BCR){
  cat("\n  \n")
  cat("## BCR {.tabset .tab-fade} \n\n")
  ### Code for isotypes class.
  # Initialize empty lists to store the dataframes
  all_isotype_class_bcr <- data.frame()
  for (plate_name in levels(PLATES_LIST)) {
    BCR_df = BCR_df_list[[ plate_name]]
    
    # Extract the data from 'c_call' and 'cell_id' columns from dataframe
    df_class_by_cell_id_bcr <- BCR_df[, c(COLUMN_HEADER_C_CALL, COLUMN_HEADER_CELL_ID)]
    
    # Filter lines that do not contain 'IGKC' or 'IGLC' in the 'c_call' column
    df_class_by_cell_id_bcr <- df_class_by_cell_id_bcr[!grepl("IGKC|IGLC", df_class_by_cell_id_bcr[,COLUMN_HEADER_C_CALL]), ]
    
    # There may be multiple cell_IDs, but the aim is to study Cell_IDs only once.
    # detect them and work only on them. (Do after in the for loop).
    unique_cell_ids <- unique(df_class_by_cell_id_bcr[[COLUMN_HEADER_CELL_ID]])
    
    # List to store results for each plate
    result_vect <- c()
    cell_id_vect <- c()
    
    for (cell_id in unique_cell_ids) {
      # Extract rows matching cell_id
      subset_df <- df_class_by_cell_id_bcr[df_class_by_cell_id_bcr[[COLUMN_HEADER_CELL_ID]] == cell_id, ]
      
      # Filter empty c_call values (if c_call is an empty string, we don't add it)
      non_empty_calls <- subset_df[[COLUMN_HEADER_C_CALL]][subset_df[[COLUMN_HEADER_C_CALL]] != ""]
      
      # Extract the part before the star (*)
      non_empty_calls <- sub("\\*.*", "", non_empty_calls)
      
      # If several non-empty values, we add them
      if (length(non_empty_calls) > 0) {
        # If only one non-empty value, we add it directly
        if (length(non_empty_calls) == 1) {
          result_vect <- c(result_vect, non_empty_calls)
          cell_id_vect <- c(cell_id_vect, rep(cell_id, length(non_empty_calls)))
        } else {
          # If several non-empty values, paste them (with underscore)
          result_vect <- c(result_vect, paste(unique(non_empty_calls), collapse = "_"))
          cell_id_vect <- c(cell_id_vect, rep(cell_id, length(non_empty_calls)))
        }
      }
    }
    # Create a dataframe to associate each class with its cell_id
    result_df <- data.frame(cell_id = cell_id_vect,
                            c_call_class = result_vect)
    
    if (length(result_df != 0)){
      # Add the plate to the dataframe, to add it to the main dataframe (all_isotype_class_bcr)
      result_df[ , COLUMN_HEADER_PLATE_NAME] = plate_name
      
      # If there is a difference between the barcodes in TRUST4 and the barcode file supplied, they are accumulated.
      erroneous_barcodes = setdiff( result_df[[COLUMN_HEADER_CELL_ID]], CELL_BARCODE_WELL_DF$BarcodeSequence)
      
      # If erroneous_barcodes is greater than 0, then warning (because there are barcodes in trust4 not found in the supplied bacrode file).
      if (length(erroneous_barcodes) > 0){
        warning(paste("WARNING : There is a discripency between the list of Well Barcode in the TRUST4 results and the list of Well Barcode from the provider.\nSome barcodes present in the TRUST4 results are not found in the provided cell barcode file in Plate :", plate_name, "\n Barcodes not found:", paste( erroneous_barcodes, collapse=";")))
      }
      
      # Creation of a WellID column containing the well id assigned to the barcode (for plate plots)
      result_df[COLUMN_HEADER_WELL_ID] = CELL_BARCODE_WELL_DF$WellID[match(result_df[[COLUMN_HEADER_CELL_ID]], CELL_BARCODE_WELL_DF$BarcodeSequence)]
      # Creation of a dataframe containing all results for each plate
      all_isotype_class_bcr = rbind(all_isotype_class_bcr, result_df)
    }
  }
  
  # ##################################################
  # ## 2. General graphics of all Isotype plates #####
  # ##################################################
  
  # Order the plate names by the order provided by the user
  all_isotype_class_bcr[ , COLUMN_HEADER_PLATE_NAME] = factor(all_isotype_class_bcr[ , COLUMN_HEADER_PLATE_NAME], levels = levels(PLATES_LIST))
  
  # Define a factor for each class, then apply the colors between classes (which are fixed) to the factor for better understanding.
  categorical_v = unique(all_isotype_class_bcr$c_call_class)
  categorical_v = sort(categorical_v)
  all_isotype_class_bcr$c_call_class = factor(all_isotype_class_bcr$c_call_class, levels = categorical_v)
  color_categorical = GENERAL_COLOR_PANEL[ 1:length( levels( all_isotype_class_bcr$c_call_class))]
  names( color_categorical) = levels( all_isotype_class_bcr$c_call_class)
  
  # Display stacked barchat for all classes across all plates
  cat("### Isotype class {.tabset .tab-fade} \n\n")
  plot_stacked_barchat <- ggplot(all_isotype_class_bcr,  aes( x = .data[[COLUMN_HEADER_PLATE_NAME]])) + 
    geom_bar(aes(fill=.data[[COLUMN_HEADER_CLASS_CALL]])) + theme_light() +
    scale_fill_manual(values = color_categorical)
  print(plot_stacked_barchat)
  
  # ##############################################
  # ## 3. Plate-by-plate graphs for Isotypes #####
  # ##############################################
  
  cat("\n  \n")
  # Plate plot of each class of each plate
  for (plate_name in levels(PLATES_LIST)) {
    cat("#### ", plate_name, " \n\n")
    class_plate = all_isotype_class_bcr[ which( all_isotype_class_bcr[ , COLUMN_HEADER_PLATE_NAME] == plate_name),]
    
    plate_plot_class <- simple_plate_plot(
      data = class_plate,
      position = COLUMN_HEADER_WELL_ID,
      value = COLUMN_HEADER_CLASS_CALL,
      plate_size = 96,
      plate_type = "round",
      title = paste("Isotype class of", plate_name)
    )
    print(plate_plot_class)
    cat("\n  \n")
  }
  
  # ################################################################################
  # ## 4. Pre-processing of BCR data supplied by Trust4 for string composition #####
  # ################################################################################
  
  ### Count for the productivity Heavy and light chain.
  ## This part works in the same way as above, but not on the same data.
  all_productivity_by_cell_id_bcr = data.frame()
  for (plate_name in levels(PLATES_LIST)) {
    BCR_df = BCR_df_list[[ plate_name]]
    
    # Only keeps COLUMN_HEADER_V_CALL, COLUMN_HEADER_J_CALL, COLUMN_HEADER_PRODUCTIVE and COLUMN_HEADER_CELL_ID in the dataset.
    df_productivity_by_cell_id <- BCR_df[, c(COLUMN_HEADER_V_CALL, COLUMN_HEADER_J_CALL, COLUMN_HEADER_PRODUCTIVE, COLUMN_HEADER_CELL_ID)]
    unique_cell_ids <- unique(df_productivity_by_cell_id$cell_id)
    
    # List to store values
    count_chain_composition_vect <- c()
    cell_id_vect <- c()
    
    for (cell_id in unique_cell_ids) {
      subset_df <- df_productivity_by_cell_id[df_productivity_by_cell_id$cell_id == cell_id, ]
      chain_components <- c()
      
      # For the entire contents of the df subset
      for (number_row in 1:nrow(subset_df)) {
        # Stores line data
        v_call_value <- subset_df[[COLUMN_HEADER_V_CALL]][number_row]
        j_call_value <- subset_df[[COLUMN_HEADER_J_CALL]][number_row]
        productive_value <- subset_df[[COLUMN_HEADER_PRODUCTIVE]][number_row]
        
        # Check if it's productive and if it contains IGH or anything else.
        if (grepl("IGH", v_call_value) | grepl("IGH", j_call_value)) {
          if (productive_value) {
            chain_components <- c(chain_components, "PH")
          } else {
            chain_components <- c(chain_components, "NPH")
          }
        }
        
        # Check the productivity of is set to True, and if it's a chain that contain K (for Light Chain Kappa).
        if (grepl("IGK", v_call_value) | grepl("IGK", j_call_value)) {
          if (productive_value) {
            chain_components <- c(chain_components, "PK")
          } else {
            chain_components <- c(chain_components, "NPK")
          }
        }
        
        if (grepl("IGL", v_call_value) | grepl("IGL", j_call_value)) {
          if (productive_value) {
            chain_components <- c(chain_components, "PL")
          } else {
            chain_components <- c(chain_components, "NPL")
          }
        }
      }
      
      # Normalize string
      chaine_caractere <- normalize_chain(chain_components)
      count_chain_composition_vect <- c(count_chain_composition_vect, chaine_caractere)
      cell_id_vect <- c(cell_id_vect, cell_id)
    }
    result_df <- data.frame(cell_id = cell_id_vect, chain_composition = count_chain_composition_vect)
    result_df[ , COLUMN_HEADER_PLATE_NAME]  <- plate_name
    result_df[COLUMN_HEADER_WELL_ID] <- CELL_BARCODE_WELL_DF$WellID[match(result_df$cell_id, CELL_BARCODE_WELL_DF$BarcodeSequence)]
    all_productivity_by_cell_id_bcr = rbind(all_productivity_by_cell_id_bcr, result_df)
  }
  
  # ##############################################################
  # ## 5. General graph of all plates for string composition #####
  # ##############################################################
  
  # Order plate names according to the order provided by the user.
  all_productivity_by_cell_id_bcr[, COLUMN_HEADER_PLATE_NAME] = factor(all_productivity_by_cell_id_bcr[, COLUMN_HEADER_PLATE_NAME], levels = levels(PLATES_LIST))
  
  # Define a factor for each string composition class, then apply fixed colors to these classes.
  categorical_v = unique(all_productivity_by_cell_id_bcr$chain_composition)
  categorical_v = sort(categorical_v)
  all_productivity_by_cell_id_bcr$chain_composition = factor(all_productivity_by_cell_id_bcr$chain_composition, levels = categorical_v)
  
  # Define colors for each level of the chain composition
  color_categorical = GENERAL_COLOR_PANEL[1:length(levels(all_productivity_by_cell_id_bcr$chain_composition))]
  names(color_categorical) = levels(all_productivity_by_cell_id_bcr$chain_composition)
  cat("### Chain composition {.tabset .tab-fade} \n\n")
  cat("<b>Lexique : </b>\n\n")
  cat("P = Productive | NP = Non-productive | H = Heavy chain | L = Light chain lambda | K = Light chain kappa")
  cat(" \n\n")
  
  # Display a stacked barchat graph for all string composition classes on all plates
  plot_stacked_barchat <- ggplot(all_productivity_by_cell_id_bcr, aes(x = .data[[COLUMN_HEADER_PLATE_NAME]])) +
    geom_bar(aes(fill = .data[[COLUMN_HEADER_CHAIN_COMPOSITION]])) +
    theme_light() +
    scale_fill_manual(values = color_categorical)
  print(plot_stacked_barchat)
  cat("\n  \n")
  
  # ######################################################
  # ## 6. Plate-by-plate graph for chain composition #####
  # ######################################################
  
  # For each plate, display a graph of chain composition
  for (plate_name in levels(PLATES_LIST)) {
    cat("#### ", plate_name, " \n\n")
    
    # Filter data for current plate
    class_plate = all_productivity_by_cell_id_bcr[which(all_productivity_by_cell_id_bcr[, COLUMN_HEADER_PLATE_NAME] == plate_name), ]
    plate_plot_class <- simple_plate_plot(
      data = class_plate,
      position = COLUMN_HEADER_WELL_ID,
      value = COLUMN_HEADER_CHAIN_COMPOSITION,
      plate_size = 96,
      plate_type = "round",
      title = paste("Chain composition of", plate_name))
    print(plate_plot_class)
    cat("\n  \n")
  }
}
