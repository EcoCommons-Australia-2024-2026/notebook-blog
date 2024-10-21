# Load necessary libraries
library(shiny)
library(ggplot2)
library(terra)   # For handling rasters
library(sf)      # For handling shapefiles
library(dplyr)

# Define the folder paths where the rasters and CSV files are stored
folder_path_rasters <- "qld_3species_Marxan/QLD_feature/"
folder_path_csvs <- "qld_3species_Marxan/model_evaluation/"
output_dir <- "qld_3species_Marxan/QLD_feature/Marxan_feature_input/"

QLD_Unit <- "qld_3species_Marxan/QLD_plannningunits/cost-surface-template.shp"  #This cost-surface-template was prepared by the Marxan Mapp with a resolution of 189 Km2, which is the highest resolution Marxan Mapp can give at this scale.

QLD_Unit  <- st_read(QLD_Unit)

st_crs(QLD_Unit)
QLD_Unit  <- st_simplify(QLD_Unit , dTolerance = 0.01) 


# Get a list of all .tif files and CSV files in the folder
raster_files <- list.files(path = folder_path_rasters, pattern = "\\.tif$", full.names = TRUE)
csv_files <- list.files(path = folder_path_csvs, pattern = "\\.csv$", full.names = TRUE)

# Extract the species names from the file names (removing the folder path and .tif/.csv extension)
species_names <- tools::file_path_sans_ext(basename(raster_files))

# Read all raster files in one go using lapply
raster_list <- lapply(raster_files, rast)  # Use rast() from terra for reading rasters

# Transform the raster CRS to match the vector CRS and apply masking in one step
raster_list <- lapply(raster_list, function(r) {
  r_transformed <- project(r, crs(vect(QLD_Unit)))
  mask(r_transformed, vect(QLD_Unit))
})

# Prepare a named list of rasters
species_rasters <- setNames(raster_list, species_names)
species_csvs <- setNames(csv_files, species_names)

# Define UI for the application
ui <- fluidPage(
  titlePanel("Interactive TSS-based threshold setting for probability of presence for multiple Species"),
  
  # Use a loop to create a row for each species
  lapply(species_names, function(species) {
    fluidRow(
      column(3, 
             h4(paste("Species:", species)),
             sliderInput(paste0("tss_value_", species), 
                         "Select TSS Value:", 
                         min = 0, max = 1, value = 0.5, step = 0.01),
             actionButton(paste0("run_analysis_", species), "Run Species Analysis"),
             br(),
             textOutput(paste0("tpr_tnr_", species))
      ),
      
      column(4, 
             plotOutput(paste0("plot_", species), width = "600px")
      ),
      
      column(5, 
             plotOutput(paste0("species_plot_", species))
      )
    )
  })
)

# Define server logic
server <- function(input, output, session) {
  
  selected_raster <- function(species) {
    species_rasters[[species]]
  }
  
  species_eval_data <- function(species) {
    csv_path <- species_csvs[[species]]
    
    if (!file.exists(csv_path)) {
      showNotification(paste("CSV file for", species, "not found!"), type = "error")
      return(NULL)
    }
    
    eval_data <- read.csv(csv_path)
    
    if (!all(c("tpr", "tnr", "tpv") %in% names(eval_data))) {
      showNotification(paste("Required columns missing in CSV for", species), type = "error")
      return(NULL)
    }
    
    if (nrow(eval_data) == 0) {
      showNotification(paste("No data found in CSV for", species), type = "error")
      return(NULL)
    }
    
    eval_data$tss <- round(eval_data$tpr + eval_data$tnr - 1, 3)
    return(eval_data)
  }
  
  lapply(species_names, function(species) {
    eval_data <- species_eval_data(species)
    
    if (!is.null(eval_data)) {
      min_tss <- min(eval_data$tss, na.rm = TRUE)
      max_tss <- max(eval_data$tss, na.rm = TRUE)
      
      updateSliderInput(session, paste0("tss_value_", species), 
                        min = min_tss, 
                        max = max_tss, 
                        value = max_tss,
                        step = 0.01)
    }
    
    observeEvent(input[[paste0("tss_value_", species)]], {
      if (!is.null(eval_data)) {
        row <- which.min(abs(eval_data$tss - input[[paste0("tss_value_", species)]]))
        
        tpr <- eval_data$tpr[row]
        tnr <- eval_data$tnr[row]
        
        output[[paste0("tpr_tnr_", species)]] <- renderText({
          paste0("TPR (Sensitivity): ", round(tpr, 3), 
                 ", TNR (Specificity): ", round(tnr, 3))
        })
      }
    })
    
    output[[paste0("plot_", species)]] <- renderPlot({
      if (is.null(eval_data)) return(NULL)
      
      ggplot(eval_data, aes(x = tpv)) +
        geom_line(aes(y = tpr, colour = "TPR"), linewidth = 1) +
        geom_line(aes(y = tnr, colour = "TNR"), linewidth = 1) +
        geom_line(aes(y = tss, colour = "TSS"), linewidth = 1) +
        geom_vline(xintercept = eval_data$tpv[which.min(abs(eval_data$tss - input[[paste0("tss_value_", species)]]))],
                   linetype = "dotted", color = "red", linewidth = 1) +
        labs(title = paste("Sensitivity, Specificity, and TSS for", species),
             x = "Threshold Probability Value",
             y = "Value") +
        scale_colour_manual(values = c("TPR" = "blue", "TNR" = "green", "TSS" = "red")) +
        theme_minimal()
    })
    
    observeEvent(input[[paste0("run_analysis_", species)]], {
      species_shp <- process_species(selected_raster(species), QLD_Unit, species, output_dir, input[[paste0("tss_value_", species)]])
      
      output[[paste0("species_plot_", species)]] <- renderPlot({
        ggplot() +
          geom_sf(data = QLD_Unit, fill = NA, color = "grey") +
          geom_sf(data = species_shp, aes(fill = feature), color = NA) +
          scale_fill_viridis_c(option = "plasma") +
          labs(title = paste("Species Distribution for", species),
               x = "Longitude", y = "Latitude") +
          theme_minimal()
      })
    })
  })
}

process_species <- function(raster_data, planning_unit, species_name, output_dir, tss_threshold) {
  raster_data_transformed <- project(raster_data, crs(vect(planning_unit)))
  extracted_values <- extract(raster_data_transformed, vect(planning_unit), fun = mean, na.rm = TRUE)
  names(planning_unit)[names(planning_unit) == "cost"] <- "feature"
  planning_unit$feature <- extracted_values[, 2]
  
  QLD_species <- subset(planning_unit, feature >= tss_threshold)
  shapefile_base <- file.path(output_dir, species_name)
  st_write(QLD_species, paste0(shapefile_base, ".shp"), delete_layer = TRUE)
  
  return(QLD_species)
}

# Run the application 
shinyApp(ui = ui, server = server)
