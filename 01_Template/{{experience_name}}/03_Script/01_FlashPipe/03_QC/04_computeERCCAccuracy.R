## ####################################################
## This script aims to create a plot for the 
## ERRC correlation (Accuracy), for each well and each plate.
## ####################################################

## @knitr compute_ERCC_accuracy

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# General description of the code :
# 1. Calculation of the number of molecules for ERCC data
# 2. Calculation of correlation between expected and observed ERCC data
# 3. Construction of general ERCC correlation graph
# 4. Construction of ERCC correlation scatter plots
# 5. Construction of individual plate graphs for ERCC correlation
# 6. Retrieve graphs of minimum, maximum and mean ERCC values to display graphs in HTML report.
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# ################################################################
# ## 1. Calculation of the number of molecules for ERCC data #####
# ################################################################

if (ERCC){
  # Create an empty dataframe to store the ERCC correlation results
  all_ERCC_correlation <- data.frame(stringsAsFactors = FALSE)
  
  # Calculate the number of molecules in the expected data
  ERCC_EXPECTED_DF$NbrOfMolecules <- ERCC_EXPECTED_DF$Concentration_attomoles_ul * UL_BY_WELL * ATTOMOLE_TO_MOLECULE
  
  # Set the row names for both dataframes using ERCC_ID for the expected data
  row.names(ERCC_EXPECTED_DF) <- ERCC_EXPECTED_DF$ERCC_ID
  ERCC_EXPECTED_DF$ERCC_ID <- NULL
  
  # Initialize the list of all the dataframe for each plates
  all_ercc_correlation_by_plate_list = list()
  
  # ##############################################################################
  # ## 2. Calculation of correlation between expected and observed ERCC data #####
  # ##############################################################################
  
  # For each plate, and each well retrieve the observed ERCC UMI and compare them to
  # the expected ERCC UMI with a Pearson correlation test.
  # Show the results by plate and for all plates
  for (plate_name in levels(PLATES_LIST)) {
    # Load the observed ERCC data
    ERCC_observed_df <- data.frame( ERCC_count_df_list[[plate_name]])
    
    # Check the common ERCC IDs between the two dataframes
    common_ERCC <- intersect(row.names(ERCC_EXPECTED_DF), row.names(ERCC_observed_df))
    
    # Select only the ERCCs that are present in both dataframes
    ERCC_sort <- ERCC_EXPECTED_DF[common_ERCC, ]
    
    # Vectors to store correlation values and counts
    ERCC_correlation = c()
    ERCC_length_number = c()
    
    # Loop through each well to calculate correlation
    for (well_id in colnames(ERCC_observed_df)) {
      
      # Order the ERCCs in the same order as those present
      ERCC_well_observed = ERCC_observed_df[ row.names(ERCC_sort), well_id]
      
      # Calculate the correlation between observed and expected values
      corr_value <- cor(ERCC_well_observed, ERCC_sort['NbrOfMolecules'], method = "pearson")
      ERCC_correlation = c(ERCC_correlation, corr_value)
      ERCC_length_number = c(ERCC_length_number, length(ERCC_well_observed))
    }
    
    # Store the correlation results and lengths in a dataframe
    names(ERCC_correlation) = colnames(ERCC_observed_df)  
    names(ERCC_length_number) = colnames(ERCC_observed_df)
    ERCC_correlation <- as.data.frame(ERCC_correlation)
    ERCC_length_number <- as.data.frame(ERCC_length_number)
    ERCC_correlation <- merge(x=ERCC_correlation, y= ERCC_length_number, by = "row.names", all = TRUE)
    
    # Rename the Row.names column to WellID
    names(ERCC_correlation)[names(ERCC_correlation) == "Row.names"] <- COLUMN_HEADER_WELL_ID
    
    # ##########################################################
    # ## 3. Construction of general ERCC correlation graph #####
    # ##########################################################
    
    # Stores ERCC data sorted by plate name for subsequent display in a plate_plot
    all_ercc_correlation_by_plate_list[[ plate_name]] = ERCC_correlation
    
    # For the general plot with all points
    plate_correlation <- data.frame(Plate_Name = rep(plate_name, nrow(ERCC_correlation)),
                                    ERCC_correlation = ERCC_correlation)
    
    # Add to our global dataframe, to later use it for a plot outside the loop
    all_ERCC_correlation <- rbind(all_ERCC_correlation, plate_correlation)
  }
  
  # Rename the columns of the ERCC_Correlation dataframe
  colnames(all_ERCC_correlation) <- gsub("^ERCC_correlation.", "", colnames(all_ERCC_correlation))
  
  # Order the plate names by the order provided by the user
  all_ERCC_correlation[, COLUMN_HEADER_PLATE_NAME] = factor(all_ERCC_correlation[, COLUMN_HEADER_PLATE_NAME], levels = levels(PLATES_LIST))
  
  # Dot plot of all plates for each ERCC.
  plot_all_correlation <- ggplot(all_ERCC_correlation, aes(x = .data[[COLUMN_HEADER_PLATE_NAME]], y = ERCC_correlation)) + 
    geom_violin(aes(fill = .data[[COLUMN_HEADER_PLATE_NAME]])) + 
    geom_jitter(width = 0.2) +
    stat_summary(fun = "mean", geom = "crossbar", colour = "red", linewidth = 0.5, width = 0.08, aes(linetype = "Mean", group = 1)) + 
    stat_summary(fun = "median", geom = "crossbar", colour = "orange", linewidth = 0.5, width = 0.08, aes(linetype = "Median", group = 1)) +
    scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
    scale_linetype_manual(values = c("Mean" = "solid", "Median" = "solid")) + 
    labs(linetype = "Statistics") + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Plate name") + ylab("ERCC pearson correlation") +
    ggtitle("Violin Plot of ERCC accuracy for all plates")
  print(plot_all_correlation)
  
  # ##########################################################
  # ## 4. Construction of ERCC correlation scatter plots #####
  # ##########################################################
  
  # Print the ERCC UMI per plate in plate plots
  for (plate_name in levels(PLATES_LIST)) {
    # Load the observed ERCC data
    ERCC_observed_df <- data.frame( ERCC_count_df_list[[plate_name]])
    
    # Check the common ERCC IDs between the two dataframes
    common_ERCC <- intersect(row.names(ERCC_EXPECTED_DF), row.names(ERCC_observed_df))
    
    # Select only the ERCCs that are present in both dataframes
    ERCC_sort <- ERCC_EXPECTED_DF[common_ERCC, ]
    
    # Create the output path
    plate_output_directory <- file.path(PATH_ANALYSIS_OUTPUT, plate_name, '04_computeERCCAccuracy')
    path_scatter_plot <- file.path(plate_output_directory, "Scatter_Plot_by_Well/Scatter_Plot_by_Well_log_ERCC")
    path_scatter_log_plot <- file.path(plate_output_directory, "Scatter_Plot_by_Well/Scatter_Plot_by_Well_ERCC")
    
    # Create the output directory if it doesn't exist
    if (!dir.exists(path_scatter_plot)) {
      dir.create(path_scatter_plot, recursive = TRUE)
    }
    
    # Create the output directory if it doesn't exist
    if (!dir.exists(path_scatter_log_plot)) {
      dir.create(path_scatter_log_plot, recursive = TRUE)
    }
    
    # Loop through each well to calculate correlation and generate plots
    scatter_plot_list = list()
    for (well_id in colnames(ERCC_observed_df)) {
      
      # Order the ERCCs in the same order as those present
      ERCC_well_observed = ERCC_observed_df[ row.names(ERCC_sort), well_id]
      
      # Generate scatter plot for each well
      plot_filename_log <- file.path(path_scatter_log_plot, paste0("ScatterPlot_log_", well_id, "_ERCC.png"))
      plot_filename <- file.path(path_scatter_plot, paste0("ScatterPlot_", well_id, "_ERCC.png"))
      
      # Scatter plot without log
      plot <- ggplot(data.frame(Observed = ERCC_well_observed, Expected = ERCC_sort$NbrOfMolecules), aes(x = Expected , y = Observed)) +
        geom_point() +
        labs(title = paste("Scatter plot for well", well_id),
             x = "Expected Values",
             y = "Observed Values") + geom_smooth(method=lm, se=FALSE, col = 'red')
      
      # Get all the plot for each well_id
      scatter_plot_list[[ well_id]] = plot
      
      # Scatter plot with log
      plot_log <- ggplot(data.frame(Observed = ERCC_well_observed, Expected = ERCC_sort$NbrOfMolecules), aes(x = log10( Expected +1) , y = log10( Observed + 1))) +
        geom_point() +
        labs(title = paste("Scatter plot for well", well_id),
             x = "Expected Values",
             y = "Observed Values")
      
      # Save the plots
      ggsave(plot_filename_log, plot = plot_log)
      ggsave(plot_filename, plot = plot)
    }
    
    # ########################################################################
    # ## 5. Construction of individual plate graphs for ERCC correlation #####
    # ########################################################################
    
    cat("\n \n")
    cat("\n#### ", plate_name, "{.tabset .tab-fade} \n\n")
    
    # Retrieve ERCC data for each plate and display it plate by plate
    ERCC_correlation = all_ercc_correlation_by_plate_list[[plate_name]]
    
    # Save the ERCC correlation results to a CSV file
    write.csv(ERCC_correlation, file.path(plate_output_directory, "ERCC_Correlation_Pearson.csv"))
    
    # Generate a plate-wide plot
    plates_plot_ERCC <- plate_plot(
      data = ERCC_correlation,
      position = WellID,
      value = ERCC_correlation,
      plate_size = 96,
      plate_type = "round",
      limits = c(0,1),
      title = paste("Plate", plate_name, " : ERCC pearson correlation")
    )
    plates_plot_ERCC = plates_plot_ERCC + labs(fill = "ERCC pearson correlation")
    print( plates_plot_ERCC)
    
    # Save the plot to file
    plate_output_directory <- file.path(PATH_ANALYSIS_OUTPUT, plate_name)
    plot_filename_ERRC <- file.path(plate_output_directory, paste0("ScatterPlot_", plate_name, "ERCC_Correlation.png"))
    ggsave(plot_filename_ERRC, plot = plates_plot_ERCC)
    
    # ######################################################################################################
    # ## 6. Retrieve graphs of minimum, maximum and mean ERCC values to display graphs in HTML report. #####
    # ######################################################################################################
    
    # Recovery of maximum ERRC correlation
    ERCC_max <- which.max(ERCC_correlation$ERCC_correlation)
    ERCC_max_df <- ERCC_correlation[ERCC_max,]
    ERCC_max_df[ , COLUMN_HEADER_STATISTIC] <- "Maximum"
    
    # Minimum ERRC correlation recovery
    ERCC_min <- which.min(ERCC_correlation$ERCC_correlation)
    ERCC_min_df <- ERCC_correlation[ERCC_min,]
    ERCC_min_df[ , COLUMN_HEADER_STATISTIC] <- "Minimum"
    
    # Recovery of the average ERRC correlations and selection of the closest line
    # If NAs are present in the data, then there may be errors in the code. So we still want results. We ignore them (because other values are present).
    ERCC_mean <- mean(ERCC_correlation$ERCC_correlation, na.rm = TRUE)
    ERCC_mean_closest <- Closest(ERCC_correlation$ERCC_correlation, a = ERCC_mean, na.rm = TRUE)
    ERCC_mean_closest_df <- ERCC_correlation[ERCC_correlation$ERCC_correlation == ERCC_mean_closest & !is.na(ERCC_correlation$ERCC_correlation), ]
    ERCC_mean_closest_df[ , COLUMN_HEADER_STATISTIC] <- "Mean"
    
    # Creation of a dataframe containing these 3 values (min, mean, max)
    statistic_df <- data.frame(stringsAsFactors = FALSE)
    statistic_df <- rbind(ERCC_min_df, ERCC_mean_closest_df, ERCC_max_df)
    
    # Display scatter plots of minimum, mean_closest and maximum ERCC correlation.
    for (name_stat in (statistic_df[ , COLUMN_HEADER_STATISTIC])){
      cat("\n \n")
      cat("\n#####", name_stat, "\n\n")
      cat( "\nERCC Correlation for well ID", statistic_df$WellID[statistic_df[ , COLUMN_HEADER_STATISTIC] == name_stat], ":", statistic_df$ERCC_correlation[statistic_df[ , COLUMN_HEADER_STATISTIC] == name_stat], "\n \n")
      print( scatter_plot_list[[  statistic_df$WellID[statistic_df[ , COLUMN_HEADER_STATISTIC] == name_stat]]])
      
    }}
  
}