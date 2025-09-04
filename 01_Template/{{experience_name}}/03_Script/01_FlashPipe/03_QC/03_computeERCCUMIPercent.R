## ####################################################
## This script aims to create a plot for the
## UMI ERRC in percent, for each well and each plate.
## ####################################################

## @knitr compute_ERCC_UMI_percent

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# General description of the code :
# 1. Verification of processed zUMI data for ERCC and RNA UMI
# 2. Construction of dataframe containing ERCC data set for each plate
# 3. Construction of general graph containing all ERCC data from all plates
# 4. Construction of plate-by-plate graphs of ERCC data
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# ####################################################################
# ## 1. Verification of processed zUMI data for ERCC and RNA UMI #####
# ####################################################################

# ••••••••••••••••••••••••••••• Verification (must be done here as the data is retrieved from the 02_formatCountTable) •••••••••••••••••••••••••••••
# Initialize an error list
if (ERCC) {
  plate_errors_ERCC_UMI_percent <- c()
  
  # For each plate, retrieve the information from ERCC UMI and RNA UMI and compute the percentage of ERCC UMI in each well
  for (plate_name in levels(PLATES_LIST)) {
    # Load ERCC observation data
    ERCC_observed_df <- data.frame( ERCC_count_df_list[[plate_name]])
    # Read the RNA_Count file (required for the UMI plot)
    RNA_count_df = RNA_count_df_list[[ plate_name]]
    
    # Verify list of wells in ERCC and RNA
    ercc_wells = sort( names( ERCC_observed_df))
    rna_wells = sort( names( RNA_count_df))
    if( sum( ercc_wells == rna_wells) != length( ercc_wells)){
      plate_errors_ERCC_UMI_percent <- c(plate_errors_ERCC_UMI_percent, paste( "The plate", plate_name, "has different wells between ERCC and RNA :\nERCCwells=", paste( ercc_wells, collapse =";"), "\nRNAWells=", paste( rna_wells, collapse =";")))
    }
    if (length(ERCC_observed_df) > 96){
      plate_errors_ERCC_UMI_percent <- c(plate_errors_ERCC_UMI_percent, paste("The plate", plate_name, "has too much well in the ERCC file -- Number of well detected :", length(ERCC_observed_df)))
    }
    if (length(RNA_count_df) > 96){
      plate_errors_ERCC_UMI_percent <- c(plate_errors_ERCC_UMI_percent, paste("The plate", plate_name, "has too much well in the RNA file -- Number of well detected :", length(ERCC_observed_df)))
    }
  }
  # Display errors if any were found 
  if (length(plate_errors_ERCC_UMI_percent) > 0) {
    cat("\n••••• ERROR DETECTED IN SOME OUTPUT FROM SCRIPT '03_computeERCCAccuracy' •••••\n")
    cat(plate_errors_ERCC_UMI_percent, "\n")
    stop("Some error has been detecting. Stop analysis")
  }
  
  # Environment cleanup: removes temporary objects
  rm(plate_errors_ERCC_UMI_percent)
  rm(ERCC_observed_df)
  rm(plate_name)
  rm(RNA_count_df)
  rm(ercc_wells)
  rm(rna_wells)
  # ••••••••••••••••••••••••••••• End of verification •••••••••••••••••••••••••••••
  
  # #############################################################################
  # ## 2. Construction of dataframe containing ERCC data set for each plate #####
  # #############################################################################
  
  # Create an empty df for the UMI ERCC percentage
  all_ERCC_ratio_by_well_df <- data.frame(stringsAsFactors = FALSE)
  
  # Initialize the list of all the dataframe for each plates
  ERCC_ratio_by_well_df_list = list()
  
  # For each plate, retrieve the information from ERCC UMI and RNA UMI and compute the percentage of ERCC UMI in each well
  for (plate_name in levels(PLATES_LIST)) {
    # Load ERCC observation data
    ERCC_observed_df <- data.frame( ERCC_count_df_list[[plate_name]])
    # Read the RNA_Count file (required for the UMI plot)
    RNA_count_df = RNA_count_df_list[[ plate_name]]
    
    # Verify list of wells in ERCC and RNA
    ercc_wells = sort( names( ERCC_observed_df))
    rna_wells = sort( names( RNA_count_df))
    
    # Create a dataframe with the sum of UMI for ERCC and RNA  
    summary_df = data.frame( RNA.sum = colSums(RNA_count_df)[ ercc_wells],
                             ERCC.sum = colSums(ERCC_observed_df)[ ercc_wells])
    
    # Create an empty data.frame to store the final results with ratio of ERCC UMI
    ERCC_ratio_by_well_df <- data.frame( Plate_Name = plate_name, WellID = ercc_wells)
    
    # Compute the percentage of ERCC UMI in each well
    ERCC_ratio_by_well_df$UMI_ERCC_Percent <- summary_df[ercc_wells, "ERCC.sum"] / (summary_df[ ercc_wells, "RNA.sum"] + summary_df[ercc_wells, "ERCC.sum"])
    
    # Here, we set up a df that collect every values for the UMI_ERCC_Percent and the well for each plate to make a dot plot with all the value on one graph.
    all_ERCC_ratio_by_well_df <- rbind(all_ERCC_ratio_by_well_df, ERCC_ratio_by_well_df)
    
    # We store in a list the values obtained for a plate (by assigning it to the plate name), for later reuse.
    ERCC_ratio_by_well_df_list[[ plate_name]] = ERCC_ratio_by_well_df
    
  }
  
  # Order the plate names by the order provided by the user
  all_ERCC_ratio_by_well_df[ , COLUMN_HEADER_PLATE_NAME] = factor(all_ERCC_ratio_by_well_df[ , COLUMN_HEADER_PLATE_NAME], levels = levels(PLATES_LIST))
  
  # ##################################################################################
  # ## 3. Construction of general graph containing all ERCC data from all plates #####
  # ##################################################################################
  
  # Print the summary of UMI ERCC for all plates
  plot_UMI_ERCC_Percent <- ggplot(all_ERCC_ratio_by_well_df, aes(x = .data[[COLUMN_HEADER_PLATE_NAME]], y = UMI_ERCC_Percent)) + 
    geom_violin(aes(fill = .data[[COLUMN_HEADER_PLATE_NAME]])) + 
    geom_jitter( width = 0.2) +
    stat_summary(fun = "mean", geom = "crossbar", colour = "red", linewidth = 0.5, width = 0.08, aes(linetype = "Mean", group = 1)) + 
    stat_summary(fun = "median", geom = "crossbar", colour = "orange", linewidth = 0.5, width = 0.08, aes(linetype = "Median", group = 1)) +
    scale_y_continuous() +  
    scale_linetype_manual(values = c("Mean" = "solid", "Median" = "solid")) + 
    labs(linetype = "Statistics") + 
    ggtitle( "Violin plot of UMI ERCC percentage for all plates") +
    xlab("Plate name") + ylab("ERCC UMI percentage") + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(plot_UMI_ERCC_Percent)
  
  # ##############################################################
  # ## 4. Construction of plate-by-plate graphs of ERCC data #####
  # ##############################################################
  
  # Compute the minimum and maximum values of ERCC to set a common range on figures
  umi_ercc_min = min( all_ERCC_ratio_by_well_df$UMI_ERCC_Percent)
  umi_ercc_max = max( all_ERCC_ratio_by_well_df$UMI_ERCC_Percent)
  
  # Print the ERCC UMI per plate in plate plots
  for (plate_name in levels(PLATES_LIST)) {
    
    cat("\n \n")
    cat("\n#### ", plate_name, "\n\n")
    
    # Graph for each well for one plate.
    plates_plot_UMI_ERCC_percent <- plate_plot(
      data = ERCC_ratio_by_well_df_list[[ plate_name]] ,
      position = WellID,
      value = UMI_ERCC_Percent,
      plate_size = 96,
      plate_type = "round",
      limits = c(umi_ercc_min, umi_ercc_max),
      title = paste( "Plate", plate_name, ":  ERCC UMI percentage")
    )
    plates_plot_UMI_ERCC_percent = plates_plot_UMI_ERCC_percent + labs(fill = "ERCC UMI percent")
    print( plates_plot_UMI_ERCC_percent)
    
    # Save the plot to file
    plate_output_directory <- file.path(PATH_ANALYSIS_OUTPUT, plate_name, '03_computeERCCUMIPercent')
    
    # Create the output directory if it doesn't exist
    if (!dir.exists(plate_output_directory)) {
      dir.create(plate_output_directory)
    }
    
    plot_filename_UMI_ERCC_Percent <- file.path(plate_output_directory, paste0("PlatePlot_UMI_ERCC_Percent_", plate_name, ".png"))
    ggsave(plot_filename_UMI_ERCC_Percent, plot = plates_plot_UMI_ERCC_percent)
    
  }
}
