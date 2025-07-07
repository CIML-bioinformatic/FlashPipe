## ####################################################
## This script aims to plot information on the data
## provide in Index Sort file.
## ####################################################

## @knitr compute_index_sort

if (INDEXSORT){
  cat("\n  \n")
  cat("## FACS {.tabset .tab-fade} \n\n")
  
  # Creation of the dataframe that will receive all the data from all the plates
  all_index_sort_df <- data.frame()
  
  # Loop that recovers all the data present in the files (for each plate) for subsequent graphing.
  for (plate_name in levels(PLATES_LIST)) {
    path_index_sort_file = paste(PATH_EXPERIMENT_RAWDATA, "/01_IndexSort/", plate_name, "_indexsort.csv", sep = "")
    # Detects the type of spacer in the csv file (based only on "," or ";")
    separator <- detect_sep_csv(path_index_sort_file)
    # Read the file with the correct separator present in the csv file.
    index_sort_df <- read.csv(path_index_sort_file, sep = separator)
    
    ## Step to keep only the wells that contain values in the columns. The others are deleted to avoid overloading the graph.
    # Position of WellID column
    well_col_index <- which(names(index_sort_df) == WELL_ID)
    other_cols <- index_sort_df[ , -well_col_index, drop = FALSE]
    # Delete rows where all other columns are NA or empty ("")
    keep_rows <- apply(other_cols, 1, function(row) {
      any(!is.na(row) & row != "")
    })
    
    # Lines to remove, with name for the warning and prevention
    if (nrow(index_sort_df[!keep_rows, ]) > 0){
      well_remove = index_sort_df[!keep_rows, ]
      cat("Empty wells were detected in the IndexSort file for the plate :", plate_name)
      cat("\n")
      cat("As a result, the wells :", well_remove[[WELL_ID]], " were removed for IndexSort analysis.")
      cat("\n\n")}
    # Filter lines to keep
    index_sort_df <- index_sort_df[keep_rows, ]
    
    # Apply a calcul to the non lineary column
    for (column_name in colnames(index_sort_df)) {
      if( class(index_sort_df[ , column_name]) == "character"){
        index_sort_df[[column_name]] = trimws( index_sort_df[[column_name]] )
      }
      if (!(startsWith(tolower(column_name), "lin_"))) {
        if (is.numeric(index_sort_df[[column_name]]) && !column_name %in% CATEGORIAL_TERM_SET) {
          index_sort_df[[column_name]] = asinh(index_sort_df[[column_name]] / 50)}}
      }
    
    # Adds current plate data to dataframe (accumulates information)
    index_sort_df[ , PLATE_NAME] = plate_name
    all_index_sort_df = rbind(all_index_sort_df, index_sort_df)
  }
  
  # Order the plate names by the order provided by the user
  all_index_sort_df[ , PLATE_NAME] = factor(all_index_sort_df[ , PLATE_NAME], levels = levels(PLATES_LIST))
  cat("<b> The non-linear data were transformed by the application of: asinh(x/50) </b>")
  cat("\n  \n")
  
  # Loops over categorical terms given as user input in config file.
  for (categorical_value in CATEGORIAL_TERM_SET){
    # Retrieves the set of (unique) names present in the category (if it's not numeric) and their attribute a factor
    if (!is.numeric(all_index_sort_df[[categorical_value]])){
      categorical_v = unique(all_index_sort_df[[categorical_value]])
      categorical_v = sort(categorical_v)
      all_index_sort_df[,categorical_value] = factor(all_index_sort_df[,categorical_value], levels = categorical_v)
    }
    # Each unique name in the column (factor) is assigned a color (to preserve colors between graphs).
    color_categorical = GENERAL_COLOR_PANEL[ 1:length( levels( all_index_sort_df[[categorical_value]]))]
    names( color_categorical) = levels( all_index_sort_df[[categorical_value]])
    
    # If the column contains numerical data, a violin plot is displayed (containing all the plates).
    if (is.numeric(all_index_sort_df[[categorical_value]])){
      cat("### ", categorical_value ," {.tabset .tab-fade} \n\n")
      # Print violin plot for column with numerical value
      violin_plot_numerical_value <- ggplot(all_index_sort_df, aes(x = .data[[PLATE_NAME]], y = .data[[categorical_value]])) + 
        geom_violin( aes( fill = .data[[PLATE_NAME]])) + 
        geom_jitter( width = 0.2) +
        stat_summary(fun = "mean", geom = "crossbar", colour = "red", linewidth = 0.5, width = 0.08, aes(linetype = "Mean", group = 1)) + 
        stat_summary(fun = "median", geom = "crossbar", colour = "orange", linewidth = 0.5, width = 0.08, aes(linetype = "Median", group = 1)) +
        scale_y_continuous() +  
        scale_linetype_manual(values = c("Mean" = "solid", "Median" = "solid")) + 
        labs(linetype = "Statistics") + 
        ggtitle( paste("Dotplot of ", categorical_value ," for all plates.")) + 
        xlab("Plate name") + 
        ylab(categorical_value) + 
        theme_light() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(violin_plot_numerical_value)
    } else {
      # If data not numeric, a stacked barchat is displayed.
      cat("### ", categorical_value," {.tabset .tab-fade} \n\n")
      plot_stacked_barchat <- ggplot(all_index_sort_df,  aes( x = .data[[PLATE_NAME]])) + 
        geom_bar(aes(fill=.data[[categorical_value]])) + theme_light() +
        scale_fill_manual(values = color_categorical) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(plot_stacked_barchat)
    }
    
    # Plate plot of all plates on target column
    cat("\n  \n")
    for (plate_name in levels(PLATES_LIST)) {
      cat("#### ", plate_name, " \n\n")
      index_sort_df = all_index_sort_df[ which( all_index_sort_df[ , PLATE_NAME] == plate_name),]
      
      plate_plot_plates <- simple_plate_plot(
        data = index_sort_df,
        position = WELL_ID,
        value = categorical_value,
        plate_size = 96,
        plate_type = "round",
        # limits = c(min_val, max_val),
        title = paste("Plate", plate_name, " : ", categorical_value))
      print(plate_plot_plates)
      cat("\n  \n")
    }
  }
  
  
  # The code does the same as before for all other columns.
  # This part is not included in the loop to allow these tabs to be displayed in a separate tab called “Index Sort”.
  cat("### Index Sort {.tabset .tab-fade} \n\n")
  for (colname in colnames(all_index_sort_df)) {
    if (!colname %in% CATEGORIAL_TERM_SET && colname != WELL_ID && colname != PLATE_NAME) {
      cat("#### ", colname, " {.tabset .tab-fade} \n\n")
      violin_plot_colname_index_sort <- ggplot(all_index_sort_df, aes(x = .data[[PLATE_NAME]], y = .data[[colname]])) + 
        geom_violin(aes( fill = .data[[PLATE_NAME]])) + 
        geom_jitter(width = 0.2) +
        stat_summary(fun = "mean", geom = "crossbar", colour = "red", linewidth = 0.5, width = 0.08, aes(linetype = "Mean", group = 1)) + 
        stat_summary(fun = "median", geom = "crossbar", colour = "orange", linewidth = 0.5, width = 0.08, aes(linetype = "Median", group = 1)) +
        scale_linetype_manual(values = c("Mean" = "solid", "Median" = "solid")) + 
        labs(linetype = "Statistics") + 
        xlab("Plate name") +
        ggtitle(paste("Dotplot of : ", colname ," for all plates")) +
        theme_light() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(violin_plot_colname_index_sort)
      cat("\n\n")
      
      limits_min = min( all_index_sort_df[[colname]])
      limits_max = max( all_index_sort_df[[colname]])
      
      for (plate_name in levels(PLATES_LIST)) {
        cat("##### ", plate_name, " {.tabset .tab-fade} \n\n")
        index_sort_df = all_index_sort_df[ which( all_index_sort_df[ , PLATE_NAME] == plate_name),]
        
        plate_plot_col_name <- simple_plate_plot(
          data = index_sort_df,
          position = WELL_ID,
          value = colname,
          plate_size = 96,
          plate_type = "round",
          limits = c(limits_min, limits_max),
          title = paste("Plate", plate_name, " : ", colname))
        print(plate_plot_col_name)
        cat("\n\n")
      }
      cat("\n\n")
    }
  }
}
