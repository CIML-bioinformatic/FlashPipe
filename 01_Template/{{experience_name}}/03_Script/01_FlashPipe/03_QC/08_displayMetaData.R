## ####################################################
## This script aims to plot information of the GSF file provide.
## ####################################################

## @knitr display_metadata

### Browse the information in the GSF.
### The information retrieved will be positioned in the Metadata tab.
### If the data found is numeric, then we'll have violin plot followed by plate plot (of each plate).
### For values other than numeric (character), we'll have a stacked barchat followed by plate plot (for each plate).

if (METADATA){
  cat("\n  \n")
  cat("## MetaData {.tabset .tab-fade} \n\n")
  
  # Retrieves data from each plate and combines it into a single dataset.
  all_plate_sheet_gsf = data.frame()
  for (plate_name in levels(PLATES_LIST)){
    plate_sheet = read_excel(PATH_GSF, sheet = plate_name)
    
    ## Step to keep only the wells that contain values in the columns. The others are deleted to avoid overloading the graph.
    # Position of WellID column
    well_col_index <- which(names(plate_sheet) == WELL_ID)
    other_cols <- plate_sheet[ , -well_col_index, drop = FALSE]
    # Delete rows where all other columns are NA or empty ("")
    keep_rows <- apply(other_cols, 1, function(row) {
      any(!is.na(row) & row != "")
    })
    
    # Lines to remove, with name for the warning and prevention
    if (nrow(plate_sheet[!keep_rows, ]) > 0){
      well_remove = plate_sheet[!keep_rows, ]
      cat("Empty wells were detected in the MetaData file for the plate :", plate_name)
      cat("\n")
      cat("As a result, the wells :", well_remove[[WELL_ID]], " were removed for Metadata analysis.")
      cat("\n\n")
    }
    
    # Filter lines to keep
    plate_sheet <- plate_sheet[keep_rows, ]
    plate_sheet[ , PLATE_NAME] = plate_name
    all_plate_sheet_gsf = rbind(all_plate_sheet_gsf, plate_sheet)
  }
  
  # The Plate_name column is assigned the plate_list factors to maintain colors between different graphs.
  all_plate_sheet_gsf[[PLATE_NAME]] = factor(all_plate_sheet_gsf[[PLATE_NAME]], levels = levels(PLATES_LIST))
  
  # Check the data types present in the general df to display the plots according to their types.
  for (column_name in colnames(all_plate_sheet_gsf)){
    # If not numeric, then factors are assigned to the data to give them a fixed color for the graphs.
    if (!is.numeric(all_plate_sheet_gsf[[column_name]])){
      categorical_v = unique(all_plate_sheet_gsf[[column_name]])
      categorical_v = sort(categorical_v)
      all_plate_sheet_gsf[[column_name]] = factor(all_plate_sheet_gsf[[column_name]], levels = categorical_v)
    }
    # Each unique name in the column (factor) is assigned a color (to preserve colors between graphs).
    color_categorical = GENERAL_COLOR_PANEL[ 1:length( levels( all_plate_sheet_gsf[[column_name]]))]
    names( color_categorical) = levels( all_plate_sheet_gsf[[column_name]])
    
    # Plot the appropriate graph for each column (if numeric = Violin Plot; if other = Stacked Barchat)
    if (is.numeric(all_plate_sheet_gsf[[column_name]]) && column_name != WELL_ID){
      cat("### ", column_name ," {.tabset .tab-fade} \n\n")
      violin_plot_numerical_value <- ggplot(all_plate_sheet_gsf, aes(x = .data[[PLATE_NAME]], y = .data[[column_name]])) + 
        geom_violin( aes( fill = .data[[PLATE_NAME]])) + 
        # Specificity: If there are different values in the column, a jitter is applied (otherwise none).
        {if (length(unique(all_plate_sheet_gsf[[column_name]])) != 1){
          geom_jitter( width = 0.2)
        } else {}} +
        stat_summary(fun = "mean", geom = "crossbar", colour = "red", linewidth = 0.5, width = 0.08, aes(linetype = "Mean", group = 1)) + 
        stat_summary(fun = "median", geom = "crossbar", colour = "orange", linewidth = 0.5, width = 0.08, aes(linetype = "Median", group = 1)) +
        scale_y_continuous() +  
        scale_linetype_manual(values = c("Mean" = "solid", "Median" = "solid")) + 
        labs(linetype = "Statistics") + 
        ggtitle( paste("Dotplot of ", column_name ," for all plates.")) + 
        xlab("Plate name") + 
        ylab(column_name) + 
        theme_light() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(violin_plot_numerical_value)
    } else if (is.factor(all_plate_sheet_gsf[[column_name]]) && column_name != WELL_ID && column_name != PLATE_NAME){
      cat("### ", column_name," {.tabset .tab-fade} \n\n")
      plot_stacked_barchat <- ggplot(all_plate_sheet_gsf,  aes( x = .data[[PLATE_NAME]])) + 
        geom_bar(aes(fill=.data[[column_name]])) + 
        theme_light() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_manual(values = color_categorical)
      print(plot_stacked_barchat)
    }
    
    # Plate plot for each column (after general plots).
    if (column_name != WELL_ID && column_name != PLATE_NAME){
      cat("\n  \n")
      for (plate_name in levels(PLATES_LIST)) {
        cat("#### ", plate_name, " \n\n")
        plate_sheet <- all_plate_sheet_gsf[ which( all_plate_sheet_gsf[[PLATE_NAME]] == plate_name),]
        # Plate_sheet is not a df, but a tible. To plot it with the function, you need a dataframe.
        plate_sheet_df <- as.data.frame(plate_sheet)
        
        plate_plot_plates <- simple_plate_plot(
          data = plate_sheet_df,
          position = WELL_ID,
          value = column_name,
          plate_size = 96,
          plate_type = "round",
          title = paste("Plate", plate_name, " : ", column_name))
        print(plate_plot_plates)
        cat("\n  \n")
      }
    }
  }
}