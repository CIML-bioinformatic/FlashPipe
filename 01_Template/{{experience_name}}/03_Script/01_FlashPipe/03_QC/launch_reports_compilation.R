# ####################################################################
# This script :
#   * define the folder the code will work in
#   * load the variables values from the parameter files
#   * launch the compilation of the report (if required)
# ####################################################################

options(future.globals.maxSize= 1048576000)

library( knitr)
library( rmarkdown)
library( funr)
library( DT)

## ......................................................................
## Define the working folder according to the different execution mode
## ......................................................................

EXECUTION_MODE_SNAKEMAKE = "Snakemake mode"
EXECUTION_MODE_RSCRIPT = "Rscript mode"
EXECUTION_MODE_RSTUDIO = "Rstudio mode"
EXECUTION_MODE_RSESSION = "R session mode"

WORKING_DIR = NULL

if( exists( "snakemake")){
  WORKING_DIR = snakemake@scriptdir
  EXECUTION_MODE = EXECUTION_MODE_SNAKEMAKE
}else{
  WORKING_DIR = tryCatch( dirname( sys.script()), error = function( e){ NA})
  if( !is.na( WORKING_DIR) && nchar( WORKING_DIR) > 0){
    EXECUTION_MODE = EXECUTION_MODE_RSCRIPT
  }else{
    WORKING_DIR = dirname( rstudioapi::getSourceEditorContext()$path)
    cat("\nWORKING_DIR=", WORKING_DIR)
    if( !is.na( WORKING_DIR) && nchar( WORKING_DIR) > 0){
      EXECUTION_MODE = EXECUTION_MODE_RSTUDIO
    }else{
      WORKING_DIR = getwd()
      EXECUTION_MODE = EXECUTION_MODE_RSESSION
      cat("WARNING: YOU ARE RUNNING IN RSESSION MODE. ENSURE YOUR WORKING DIR IS CORRECTLY SET.\nWorking dir=", WORKING_DIR)
    }
  }
}

cat("\nINFO :: EXECUTION MODE =", EXECUTION_MODE)
cat("\nINFO :: EXECUTION IN DIRECTORY =", WORKING_DIR)

## ......................................................
## Load the variable values from the parameters files
## ......................................................
# Define an environment that will contain parameters
paramsEnv = new.env();

# Assign the WORKING_DIR to the paramsEnv
assign( "WORKING_DIR" , WORKING_DIR , env = paramsEnv )

# Load file defining the project parameters
projectParamsFilePath = file.path( WORKING_DIR, "../projectParams.R");
if(file.exists(projectParamsFilePath)) {
  source( projectParamsFilePath, local = paramsEnv);
} else {
  warning("The file defining the Project parameters is missing.");
}

# Load file defining the expriment parameters
experimentParamsFilePath = file.path( WORKING_DIR, "../experimentParams.R");
if(file.exists(experimentParamsFilePath)) {
  source( experimentParamsFilePath, local = paramsEnv);
} else {
  warning("The file defining the Experiment parameters is missing.");
}

# Load file defining sample parameters
sampleParamsFilePath = file.path( WORKING_DIR, "../sampleParams.R");
if(file.exists(sampleParamsFilePath)) {
  source( sampleParamsFilePath, local = paramsEnv);
} else {
  warning("The file defining the Sample parameters is missing.");
}

# Load file defining analysis parameters
analysisParamsFilePath = file.path( WORKING_DIR, "analysisParams.R");
if(file.exists(analysisParamsFilePath)) {
  source( analysisParamsFilePath, local = paramsEnv);
} else {
  warning("The file defining the Analysis parameters is missing.");
}

# Keep trace of the initial analysis output path
assign( "ORIGINAL_PATH_ANALYSIS_OUTPUT" , paramsEnv$PATH_ANALYSIS_OUTPUT , env = paramsEnv )

## ......................................................
## Clean the local environment and 
## load the variable values in the local environment
## ......................................................

# Clean the global Environment (but not the paramsEnv)
list_variables = ls()
list_variables = list_variables[ list_variables != "paramsEnv"]
rm( list = list_variables)

# Assign loaded values in paramsEnv to current environment
invisible( lapply( ls( paramsEnv, all.names = TRUE), function(x, envir)
{ 
  assign( x = x, 
          value = get( x, pos = paramsEnv), 
          pos = envir)
}, 
environment()));

## ......................................................................................
## Execute the report compilation only in Rscript mode and Snakemake mode
## In Rstudio mode and R session mode, the report compilation must be launched manually
## ......................................................................................
# if( EXECUTION_MODE == EXECUTION_MODE_SNAKEMAKE || EXECUTION_MODE == EXECUTION_MODE_RSCRIPT){
rmarkdown::render( input = file.path( WORKING_DIR, "Report.Rmd"),
                   output_dir = PATH_ANALYSIS_OUTPUT,
                   output_file  = paste0( SCIENTIFIC_PROJECT_NAME, "_", EXPERIMENT_NAME, "_", ANALYSIS_STEP_NAME, ".html"),
                   quiet = FALSE)
# }

