# ##################################################
# Global declarations and libraries for the analysis
# ##################################################

######## R Libraries
######## Add here all the 'library( ...)' required for the analysis
library( pander)
library( digest)
library( ggplot2)
library( ggplate)
library( DescTools)
library( yaml)
library( UpSetR)
library( tidyr)
library( viridis)
library( readxl)

######## DEFINITION OF FUNCTIONS #########################
######## Add here all the definition of functions

# Function to sort and concatenate string components
normalize_chain <- function(components) {
  sorted_components <- sort(components)
  return(paste(sorted_components, collapse = "_"))
}

# Function plate plot modify to be able to use it inside a for loop
simple_plate_plot = function (data, position, value, label, plate_size = 96, plate_type = "square", 
                              colour, limits, title, title_size, show_legend = TRUE, legend_n_row, 
                              label_size, silent = TRUE, scale) 
{
  if (missing(scale)) {
    scale <- min((graphics::par("fin")[1]/5.572917), (graphics::par("fin")[2]/3.177083))
  }
  if (!silent) {
    message(paste0("width: ", round(graphics::par("fin")[1], 
                                    digits = 3), " height: ", round(graphics::par("fin")[2], 
                                                                    digits = 3)))
    message(paste0("scaling factor: ", round(scale, digits = 3)))
  }
  
  if (is.numeric(data[, value])) {
    if (missing(limits)) {
      min_val <- min(data[, value])
      max_val <- max(data[, value])
      n_distinct_values <- length(unique(data[, value]))
      if (n_distinct_values == 1 & is.numeric(min_val)) {
        max_val <- min_val + abs(min_val)
      }
    }
    else {
      min_val <- ifelse(is.na(limits[1]), min(data[, value]), limits[1])
      max_val <- ifelse(is.na(limits[2]), max(data[, value]), limits[2])
    }
    if (missing(colour)) {
      viridis_colours <- "placeholder"
      utils::data("viridis_colours", envir = environment())
      fill_colours <- viridis_colours
      colfunc <- scales::gradient_n_pal(viridis_colours, 
                                        values = c(min_val, max_val))
      data_colours <- colfunc(data[, value])
    }
    else {
      fill_colours <- colour
      colfunc <- scales::gradient_n_pal(colour, values = c(min_val, max_val))
      data_colours <- colfunc(data[, value])
    }
  }
  if (!is.numeric(data.frame( data)[, value])) {
    if (missing(colour)) {
      fill_colours = GENERAL_COLOR_PANEL[ 1:length( levels( data[[value]]))]
      names( fill_colours) = levels( data[[value]])
      data_colours = fill_colours[ as.integer( data[[value]])]
    }
    else {
      if (length(colour) < length(unique(data[, value]))) {
        stop("There are more categories in the \"value\" column than provided colours. Please add more colours to the \"colour\" argument!")
      }
      fill_colours <- colour
      data_colours <- purrr::map_chr(as.data.frame(data[, value]), .f = ~{
        fill_colours[which(.x == unique(data[, value]))]
      })
    }
  }
  if (!missing(label)) {
    hcl <- farver::decode_colour(data_colours, "rgb", "hcl")
    label_col <- ifelse(hcl[, "l"] > 50, "#000000", "#FFFFFF")
  }
  else {
    label_col <- ""
  }
  max_label_length <- max(nchar(as.character(unique(data[, value]))))
  MORELETTERS <- c(LETTERS, "AA", "AB", "AC", "AD", "AE", "AF")
  # data_prep <- dplyr::mutate(dplyr::ungroup(data), row = stringr::str_extract({{position}}, pattern = "[:upper:]+"), 
  #                            col = as.numeric(stringr::str_extract({{position}}, pattern = "\\d+")), 
  #                            row_num = as.numeric(match(.data$row, MORELETTERS)), 
  #                            colours = data_colours, 
  #                            label_colours = label_col)
  data_prep <- dplyr::mutate(dplyr::ungroup(data), 
                             row = stringr::str_extract(.data[[position]], pattern = "[:upper:]+"), 
                             col = as.numeric(stringr::str_extract(.data[[position]], pattern = "\\d+")), 
                             row_num = as.numeric(match(.data$row, MORELETTERS)), 
                             colours = data_colours, 
                             label_colours = label_col)
  
  if (!is.numeric(data[, value])) {
    data_prep <- dplyr::mutate(data_prep,!!value := forcats::fct_inorder(data_prep[[value]]))
  }
  if (show_legend) {
    label_is_numeric <- is.numeric(data_prep[, value])
  }
  else {
    label_is_numeric <- TRUE
  }
  if (plate_size != 96) {
    stop("The number plate given is different than 96.")
  }
  if (plate_size == 96) {
    n_cols <- 12
    n_rows <- 8
    size <- 9.5 * scale
    min_x_axis <- 0.75
    max_x_axis <- n_cols + 0.25
    min_y_axis <- 0.75
    max_y_axis <- n_rows + 0.25
    text_size <- 8 * scale
    legend_text_size <- text_size
    legend_title_size <- text_size
    title_size_preset <- 12 * scale
    legend_size <- size
    stroke_width <- 0.6 * scale
    if (show_legend) {
      size <- (9.5 - max((max_label_length - 9) * 0.15, 
                         0)) * scale
      legend_size <- 5 * scale
    }
  }
  data_prep <- dplyr::mutate(data_prep, row_num = abs(.data$row_num - 
                                                        (n_rows + 1)))
  if (missing(label_size)) {
    label_size_scaled <- size/3
  }
  else {
    label_size_scaled <- label_size * scale
  }
  if (missing(title_size)) {
    title_size <- title_size_preset
  }
  else {
    title_size <- title_size * scale
  }
  plate_layout <- expand.grid(columns = seq(1, n_cols), rows = seq(1, 
                                                                   n_rows))
  if (plate_type == "round") {
    shape <- 21
  }
  if (plate_type == "square") {
    shape <- 22
  }
  if (!missing(title)) {
    plot_title <- title
  }
  else {
    plot_title <- paste0(plate_size, "-well Plate Layout")
  }
  
  plot <- ggplot2::ggplot(data_prep, ggplot2::aes(x = col, y = .data$row_num)) + 
    ggplot2::geom_point(data = plate_layout, ggplot2::aes(x = .data$columns, y = .data$rows), color = "grey90", size = size, shape = shape) + 
    ggplot2::geom_point(ggplot2::aes(fill = data[, value]), size = size, shape = shape, stroke = stroke_width) + 
    ggplot2::coord_fixed(ratio = ((n_cols + 1)/n_cols)/((n_rows + 1)/n_rows), 
                         xlim = c(min_x_axis, max_x_axis), 
                         ylim = c(min_y_axis, max_y_axis)) +
    ggplot2::scale_y_continuous(breaks = seq(1, n_rows), labels = rev(MORELETTERS[1:n_rows])) +
    ggplot2::scale_x_continuous(breaks = seq(1, n_cols), position = "top")+ 
    ggplot2::labs(title = plot_title, x = "", y = "") + ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = text_size, face = "bold"), axis.ticks.x = ggplot2::element_blank(), 
                   axis.ticks.y = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = legend_text_size), 
                   legend.title = ggplot2::element_text(size = legend_title_size), 
                   plot.title = ggplot2::element_text(size = title_size), 
                   axis.title.x = ggplot2::element_blank(), panel.border = ggplot2::element_rect(linewidth = stroke_width))
  
  if( !is.numeric(data_prep[, value])){
    plot = plot + ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = legend_size)))
    plot = plot + ggplot2::scale_fill_manual(values = fill_colours)
  }
  else{
    plot = plot + ggplot2::scale_fill_gradientn(colors = fill_colours, limits = c(min_val, max_val), guide = ggplot2::guide_colorbar(ticks.linewidth = max(0.5 * scale, 0.2)))
  }
  
  plot = plot + labs(fill=value) 
  plot
}

######## Generate a datatable summarizing values
# For environments (parameters), all values/variables are shown

showSimpleDT = function( dataToShow, rownames = TRUE, colnames = "Value")
{
  valuesDF = NULL;
  
  if(is.environment( dataToShow))
  {
    # Extract values from environment as character strings
    envValues = sapply(lapply(dataToShow, function(x) {ifelse(is.null(x), "NULL", x)}), paste, collapse = ", ");
    # Sort them by name and convert to data.frame
    valuesDF = data.frame("Value" = envValues[order(names(envValues))]);
  } else
  {
    valuesDF = dataToShow;
  }
  
  # Create a datatable with collected information
  # Create datatable
  datatable( as.data.frame(valuesDF),
             class = "compact",
             rownames = rownames,
             colnames = colnames,
             options = list(dom = "<'row'rt>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                            autoWidth = FALSE,
                            columnDefs = list( # Center all columns
                              list( targets = "_all",
                                    className = 'dt-center')),
                            orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                            ordering = FALSE,
                            paging = FALSE, # Disable pagination (show all)
                            processing = TRUE,
                            scrollCollapse = TRUE,
                            scroller = TRUE,  # Only load visible data
                            scrollX = TRUE,
                            scrollY = "525px",
                            stateSave = TRUE));
}


