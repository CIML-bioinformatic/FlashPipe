import os
import subprocess
import sys
from optparse import OptionParser
import pandas as pd
import shutil

# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# General description of the code :
# 1. Retrieves the paths required for the project (existing and the one who will be created).
# 2. Loads and reads the configuration file created by the user (Config_FlashPipe).
# 3. Retrieves plate names and FastQ files to iterate over each plate and each read1 and read2.
# 4. Selects reference files according to species defined in config.
# 5. Creates a dictionary with the values from the config file for the copier.yml file.
# 6. Modifies the templave config file (copier.yml) with the new values.
# 7. Runs Copier on the template to create the necessary directories and files.
# 8. Copy Index Sorting and GSF files.
# 9. Calls the function to generate the file containing the analysis information for Airrflow.
# 10. Check if every file is correctly copied/created by Copier or not.
# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 1. Retrieves the paths for the project
# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

print("Retrieving options")

# Argument parser for the project and the experience name.
parser = OptionParser()
parser.add_option("-p", "--project_name", dest="project_name", help="Name of your project", metavar="PROJECT_NAME")
parser.add_option("-e", "--experience_name", dest="experience_name", help="Name of your experience", metavar="EXPERIENCE_NAME")
parser.add_option("-w", "--working_dir", dest="working_dir", help="Your working directory", metavar="WORKING_DIR")
parser.add_option("-t", "--template_path", dest="template_path", help="Name of your project", metavar="TEMPLATE_PATH")
    
# Parsing arguments
(options, args) = parser.parse_args()
# Check that both arguments have been supplied
if not options.project_name or not options.experience_name or not options.template_path :
    print("Error: You must provide a project name, an experience name, a working_dir and a template path.")
    parser.print_help()


# •••••••••••••••••
PATH_EXPERIENCE = options.working_dir
PATH_TEMPLATE = options.template_path
PATH_CONTAINER_TEMPLATE = os.path.join(PATH_TEMPLATE, "{{experience_name}}/02_Container/")
PATH_PROJECT_FLASH_PIPE = os.path.dirname(PATH_EXPERIENCE)
PATH_TMP_FLASHPIPE = "/tmp/FlashPipe/"
PATH_TMP_FLASHPIPE_TEMPLATE = os.path.join(PATH_TMP_FLASHPIPE, "01_Template")
PATH_RAWDATA = os.path.join(PATH_EXPERIENCE, '00_RawData/')
PATH_REFERENCE = os.path.join(PATH_EXPERIENCE, '01_Reference/')
PATH_CONTAINER = os.path.join(PATH_EXPERIENCE, '02_Container/')
PATH_WORKFLOW = os.path.join(PATH_EXPERIENCE, '04_Workflow/')
PATH_SNAKEMAKE = os.path.join(PATH_WORKFLOW, '01_snakemake/')
PATH_OUTPUT = os.path.join(PATH_EXPERIENCE, '05_Output/')
PATH_RNA = os.path.join(PATH_RAWDATA, '00_RNA/')
PATH_INDEX_SORTING = os.path.join(PATH_RAWDATA, '01_IndexSort/')
PATH_EXPERIMENT_REFERENCE = os.path.join(PATH_REFERENCE, '00_Experiment/')
PATH_ZUMIS_REFERENCE = os.path.join(PATH_REFERENCE, '01_zUMIs/')
PATH_AIRRFLOW_REFERENCE = os.path.join(PATH_REFERENCE, '02_airrflow/')
PATH_OUTPUT_FLASHPIPE = os.path.join(PATH_OUTPUT, '01_FlashPipe/')
PATH_OUTPUT_ZUMIS = os.path.join(PATH_OUTPUT_FLASHPIPE, '01_zUMIs/')
PATH_OUTPUT_AIRRFLOW = os.path.join(PATH_OUTPUT_FLASHPIPE, '02_airrflow/')
PATH_OUTPUT_TRUST4 = os.path.join(PATH_OUTPUT_FLASHPIPE, '02_trust4/')
PATH_OUTPUT_QC = os.path.join(PATH_OUTPUT_FLASHPIPE, '03_QC/')
PATH_OUTPUT_ANALYSIS = os.path.join(PATH_OUTPUT_FLASHPIPE, '04_Analysis/')

# •••••••••••••••••

# •••••••••••••••••
### Create variable for the code 
EXPERIENCE_NAME_FILE = "experience_name"
PLATE_NAMES_FILE = "plate_names"
GSF_FILE = "gsf_file"
FASTQ_DIRECORIES_FILE = "fastq_directories"
INDEX_SORT_FILE = "index_sort"
NOT_FLUORESCENT = "not_fluorescent"
SPECIES_FILE = "species"
STAR_INDEX_FILE = "star_index"
GTF_FILE = "gtf_file"
OUTDIR_TEMP_FILE = "outdir_temp"
TRUST4_IMGT_BCR_TCR = "trust4_imgt_BCR_TCR"
TRUST4_IMGT_VDJ = "trust4_imgt_VDJ"
METHOD = "method_analysis"
INDEXSORT = "index_sort_analysis"
BCR = "bcr_repertoire_analysis"
TCR = "tcr_repertoire_analysis"
METADATA = "metadata_analysis"
TOOLS_BCR_TCR_ANALYSIS = "tools_bcr_tcr_analysis"
# •••••••••••••••••

# Import the custom library and functions
sys.path.append(os.path.join(PATH_EXPERIENCE, "03_Script/01_FlashPipe/00_organizeStructure/"))
from create_folder_structure_function import open_file_yml, verify_empty_values_config_file, verify_name_experience_path_and_config_file, verify_separator_in_config_file, verify_file_exist, replace_value_config_template_copier, verify_parameters, prepare_fastq_symlinks_and_paths, verify_method, copy_index_sort, copy_gsf, airrflow_parameter, generate_airrflow_samplesheet

# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 2. Loading the configuration file (config FlashPipe)
# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

print("Reading config file")

# Define the path to the config file
file_config_flashpipe_path = os.path.join(PATH_REFERENCE, 'config_FlashPipe.yml')

# Checks the existence of the file (i.e. checks that the project and experience names are correct) (for more detail check the file that contain function)
verify_file_exist(file_config_flashpipe_path)

# Reading the config file
file_config_flashpipe = open_file_yml(file_config_flashpipe_path)

# Check all values in the config file, to make sure any section is not empty
verify_empty_values_config_file(file_config_flashpipe)

# Check that the name of the experiment provided in the path matches that given in the config file (more details in the function).
verify_name_experience_path_and_config_file(PATH_EXPERIENCE, file_config_flashpipe.get(EXPERIENCE_NAME_FILE))

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 3. Retrieves plate names and FastQ files to iterate over each plate and each read1 and read2.
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# Checking the separators in the various options in the config file.
str_split = [":", ";", "/", ".", "?"]
verify_separator_in_config_file(file_config_flashpipe.get(PLATE_NAMES_FILE), str_split)
verify_separator_in_config_file(file_config_flashpipe.get(NOT_FLUORESCENT), str_split)

# Recovers plate names from configuration and cleans them up (deletes spaces).
plates_list = file_config_flashpipe.get(PLATE_NAMES_FILE).split(',')
plates_list = [plate.strip() for plate in plates_list]

#PATH_CONFIG_AIRRFLOW = '/mnt/DOSI/MASTER_TEMP/CB2M/Test/Test_Airrflow_v1/01_Reference/config_FlashPipe.yml'
#CONFIG_FILE_AIRRFLOW = open_file_yml(PATH_CONFIG_AIRRFLOW)
#CONFIG_FILE_AIRRFLOW.get(BCR)

# Check if the parameter are the type of True or False / Single-cell or Mini Bulk for the analysis
print("•••••Parameters selection•••••")
verify_method(file_config_flashpipe.get(METHOD), file_config_flashpipe.get(INDEXSORT))
indexsort = verify_parameters(file_config_flashpipe.get(INDEXSORT), INDEXSORT)
bcr = verify_parameters(file_config_flashpipe.get(BCR), BCR)
tcr = verify_parameters(file_config_flashpipe.get(TCR), TCR)
metadata = verify_parameters(file_config_flashpipe.get(METADATA), METADATA)
print("••••••••••••••••••••••••••••••")

# Retrieve category data to be set aside for certain analysis data.
not_fluorescent_list = file_config_flashpipe.get(NOT_FLUORESCENT).split(',')
not_fluorescent_list = [categorial_term.strip() for categorial_term in not_fluorescent_list]

# Recover reads 1 and 2 from fastQ (to put them in the copier.yml file, in dictionary form) and create the symbolic link.
fastq_files_read1, fastq_files_read2 = prepare_fastq_symlinks_and_paths(file_config_flashpipe.get(FASTQ_DIRECORIES_FILE), 
                                                                        PATH_RNA, 
                                                                        plates_list)

# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 4. Selects reference files according to species defined in config.
# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

species = file_config_flashpipe.get(SPECIES_FILE)
star_index = file_config_flashpipe.get(STAR_INDEX_FILE).get(species)
gtf_file = file_config_flashpipe.get(GTF_FILE).get(species)
trust4_imgt_BCR_TCR = file_config_flashpipe.get(TRUST4_IMGT_BCR_TCR).get(species)
trust4_imgt_VDJ = file_config_flashpipe.get(TRUST4_IMGT_VDJ).get(species)

### Define the exact position of barcode and ERCC file for analysis (more precisely : zUMIs)
ERCC_PATH = os.path.join(PATH_EXPERIMENT_REFERENCE, 'ERCC_concentration.csv')
BARCODE_WELL_PATH = os.path.join(PATH_EXPERIMENT_REFERENCE, 'cell_barcode_well.csv')

# Define a barcode path (without well ID) for zUMIs (this file is create down).
BARCODE_PATH = os.path.join(PATH_EXPERIMENT_REFERENCE, "cell_barcode.txt")
GSF_FILE_PATH = os.path.join(PATH_EXPERIMENT_REFERENCE, os.path.basename(file_config_flashpipe.get(GSF_FILE)))

# Select the parameter for the analysis for Airrflow
CLONAL_PARAMETER_AIRRFLOW = airrflow_parameter(file_config_flashpipe)

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 5. Creates a dictionary with the values from the config file for the copier.yml file
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# Assemble the parameter values ina single dictionary used for the template config file
values_config_flashpipe = {
    'experience_name': file_config_flashpipe.get(EXPERIENCE_NAME_FILE),
    'barcode_file' : BARCODE_PATH,
    'fastq_files_read1': ','.join(fastq_files_read1),
    'fastq_files_read2': ','.join(fastq_files_read2), 
    'gtf_file': gtf_file,
    'trust4_imgt_VDJ' : trust4_imgt_VDJ,
    'trust4_imgt_BCR_TCR' : trust4_imgt_BCR_TCR,
    'plate_names': ','.join(plates_list),
    'star_index': star_index,
    'gsf_file': GSF_FILE_PATH,
    'index_sort': file_config_flashpipe.get(INDEX_SORT_FILE),
    'species': species,
    'ercc_concentration_file': ERCC_PATH,
    'outdir_temp' : file_config_flashpipe.get(OUTDIR_TEMP_FILE),
    'project_name' : options.project_name,
    'template_path' : PATH_TEMPLATE,
    'path_project' : PATH_PROJECT_FLASH_PIPE,
    'not_fluorescent' : ','.join(not_fluorescent_list),
    'method_analysis' : file_config_flashpipe.get(METHOD),
    'index_sort_analysis' : indexsort,
    'bcr_repertoire_analysis' : bcr,
    'tcr_repertoire_analysis' : tcr,
    'metadata_analysis' : metadata,
    'airrflow_or_trust4' : file_config_flashpipe.get(TOOLS_BCR_TCR_ANALYSIS),
    'path_output_airrflow' : PATH_OUTPUT_AIRRFLOW,
    'clonal_parameter_airrflow' : CLONAL_PARAMETER_AIRRFLOW
}

# Copier le tempalte dans le tmp pour permettre de modifier le fichier copier.yml (et d'éviter les conflits sur un même fichier)
shutil.copytree(PATH_TEMPLATE, PATH_TMP_FLASHPIPE_TEMPLATE, ignore=shutil.ignore_patterns("02_Container"))

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 6. Modifies the templave config file (copier.yml) with the new values.
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# Define the path to the template config file and replace the values with the custom ones
FILE_COPIER_TERMPLATE_YML = os.path.join(PATH_TMP_FLASHPIPE_TEMPLATE, 'copier.yml')
replace_value_config_template_copier(FILE_COPIER_TERMPLATE_YML, values_config_flashpipe)

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 7. Runs Copier on the template to create the necessary directories and files.
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# Create the Copier command and launch it
print("Create directories and files structure, might take few minutes...")
cmd = f'copier copy -f {PATH_TMP_FLASHPIPE_TEMPLATE} {PATH_PROJECT_FLASH_PIPE}'
subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

# Copy the Container that was forbidden to import (in the template) into the project
shutil.copytree(PATH_CONTAINER_TEMPLATE, PATH_CONTAINER, dirs_exist_ok=True)

# Delete the temporary template located in the tmp.
shutil.rmtree(PATH_TMP_FLASHPIPE)

# Create a file barcode (without well ID) for zUMIs.
barcode_well = pd.read_csv(BARCODE_WELL_PATH)
barcode_well["BarcodeSequence"].to_csv(BARCODE_PATH, sep=' ', index=False, header=False)

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 8. Copy Index Sorting files and gsf files.
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# Copy Index Sorting files (if not FALSE, which means that the section has been filled (and not empty))
if file_config_flashpipe.get(INDEX_SORT_FILE) != "FALSE" :
    copy_index_sort(file_config_flashpipe.get(INDEX_SORT_FILE), PATH_INDEX_SORTING, plates_list=plates_list)
    
# Copy GSF file
copy_gsf(file_config_flashpipe.get(GSF_FILE), PATH_EXPERIMENT_REFERENCE)

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 9. Calls the function to generate the file containing the analysis information for Airrflow.
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# Path to the .tsv file already copied and msi in place in the project
CSV_AIRRFLOW_PATH = os.path.join(PATH_AIRRFLOW_REFERENCE, "assembled_samplesheet.tsv")

generate_airrflow_samplesheet(path_airrflow_reference=CSV_AIRRFLOW_PATH, plates_list=plates_list,
                              fastq_files_read1=fastq_files_read1, fastq_files_read2=fastq_files_read2,
                              species=species, file_config_flashpipe=file_config_flashpipe)

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 10. Check if every file is correctly copied/created by Copier.
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
verify_file_exist(ERCC_PATH)
verify_file_exist(BARCODE_WELL_PATH)
verify_file_exist(BARCODE_PATH)
verify_file_exist(os.path.join(PATH_EXPERIENCE, '03_Script/01_FlashPipe/03_QC/analysisParams.R'))
verify_file_exist(os.path.join(PATH_EXPERIENCE, '03_Script/01_FlashPipe/projectParams.R'))
verify_file_exist(os.path.join(PATH_EXPERIMENT_REFERENCE, os.path.basename(file_config_flashpipe.get(GSF_FILE))))
verify_file_exist(os.path.join(PATH_WORKFLOW, '01_snakemake/snakefile.yaml'))
verify_file_exist(CSV_AIRRFLOW_PATH)

print("Success : The project structure is now in place.")



