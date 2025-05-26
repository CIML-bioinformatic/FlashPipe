import os
import subprocess
import sys
from optparse import OptionParser
import pandas as pd

# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# General description of the code :
# 1. Retrieves the paths required for the project (existing and tthe one who will be create).
# 2. Loads and reads the configuration file used by the user (Config_FlashPipe).
# 3. Retrieves plate names and FastQ files to iterate over each plate and each read1 and read2.
# 4. Selects files according to species defined in config.
# 5. Creates a dictionary with the values from the config file for the copier.yml file.
# 6. Modifies the copier.yml file with the new values.
# 7. Runs the template to create the necessary directories and files.
# 8. Creates FASTQ symbolic links.
# 9. Copy Index Sorting and GSF.
# 10. Check if every file is correctly import by copier or not.
# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# •••••••••••••••••
### 1. Retrieves the paths for the project

# Argument parser for the project and the experience name.
parser = OptionParser()
parser.add_option("-p", "--project_name", dest="project_name", help="Name of you're project", metavar="PROJECT_NAME")
parser.add_option("-e", "--experience_name", dest="experience_name", help="Name of you're experience", metavar="EXPERIENCE_NAME")
parser.add_option("-w", "--working_dir", dest="working_dir", help="You're working directory", metavar="WORKING_DIR")
parser.add_option("-t", "--template_path", dest="template_path", help="Name of you're project", metavar="TEMPLATE_PATH")
    
# Parsing arguments
(options, args) = parser.parse_args()
# Check that both arguments have been supplied
if not options.project_name or not options.experience_name or not options.template_path :
    print("Error: You must provide a project name, an experience name, an working_dir and an template path.")
    parser.print_help()


PATH_EXPERIENCE = options.working_dir
#PATH_EXPERIENCE = "/mnt/DOSI/MASTER_TEMP/CB2M/Test/Analysis_MCMV-infected-cells/250108_VH00228_244_AAFVYNVM5_Bulk/"
PATH_COPIER_TEMPLATE = options.template_path
#PATH_COPIER_TEMPLATE = "/mnt/DOSI/MASTER_TEMP/CB2M/Project/01_FlashPipe/01_Template"
PATH_PROJECT_FLASH_PIPE = os.path.dirname(PATH_EXPERIENCE)

PATH_RAWDATA = os.path.join(PATH_EXPERIENCE, '00_RawData/')
PATH_REFERENCE = os.path.join(PATH_EXPERIENCE, '01_Reference/')
PATH_WORKFLOW = os.path.join(PATH_EXPERIENCE, '04_Workflow/')
PATH_SNAKEMAKE = os.path.join(PATH_WORKFLOW, '01_snakemake/')
PATH_OUTPUT = os.path.join(PATH_EXPERIENCE, '05_Output/')
PATH_RNA = os.path.join(PATH_RAWDATA, '00_RNA/')
PATH_INDEX_SORTING = os.path.join(PATH_RAWDATA, '01_IndexSort/')
PATH_EXPERIMENT_REFERENCE = os.path.join(PATH_REFERENCE, '00_Experiment/')
PATH_ZUMIS_REFERENCE = os.path.join(PATH_REFERENCE, '01_zUMIs/')
PATH_TRUST4_REFERENCE = os.path.join(PATH_REFERENCE, '02_trust4/')
PATH_OUTPUT_FLASHPIPE = os.path.join(PATH_OUTPUT, '01_FlashPipe/')
PATH_OUTPUT_ZUMIS = os.path.join(PATH_OUTPUT_FLASHPIPE, '01_zUMIs/')
PATH_OUTPUT_TRUST4 = os.path.join(PATH_OUTPUT_FLASHPIPE, '02_trust4/')
PATH_OUTPUT_QC = os.path.join(PATH_OUTPUT_FLASHPIPE, '03_QC/')
PATH_OUTPUT_ANALYSIS = os.path.join(PATH_OUTPUT_FLASHPIPE, '04_Analysis/')
# •••••••••••••••••

# •••••••••••••••••
### Get variable for the code 
EXPERIENCE_NAME_FILE = "experience_name"
PLATE_NAMES_FILE = "plate_names"
GSF_FILE = "gsf_file"
FASTQ_DIRECORIES_FILE = "fastq_directories"
INDEX_SORT_FILE = "index_sort"
CATEGORIAL_TERM_SET_FILE = "categoriale_term_set"
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
# •••••••••••••••••

# Necessary to import the library
sys.path.append(os.path.join(PATH_EXPERIENCE, "03_Script/01_FlashPipe/00_organizeStructure/"))
from create_folder_structure_function import open_file_yml, verify_empty_values_config_file, verify_name_experience_path_and_config_file, verify_separator_in_config_file, verify_file_exist, replace_value_config_template_copier, verify_parameters, verify_method, fastq_read_files, create_symlinks, copy_index_sort, copy_gsf

print("Reading config file")
### 2. Loading the configuration file (config FlashPipe)
file_config_flashpipe_path = os.path.join(PATH_REFERENCE, 'config_FlashPipe.yml')

# Checks the existence of the file (i.e. checks that the project and experience names are correct) (for more detail check the file that contain function)
verify_file_exist(file_config_flashpipe_path)

# Reading the config file
file_config_flashpipe = open_file_yml(file_config_flashpipe_path)

# Check all values in the config file, to make sure any section is not empty
verify_empty_values_config_file(file_config_flashpipe)

# Check that the name of the experiment provided in the path matches that given in the config file (more details in the function).
verify_name_experience_path_and_config_file(PATH_EXPERIENCE, file_config_flashpipe.get(EXPERIENCE_NAME_FILE))

### 3. Plate configuration, categoriale term set and FastQ files

# Checking the separators in the various options in the config file.
str_split = [":", ";", "/", ".", "?"]
verify_separator_in_config_file(file_config_flashpipe.get(PLATE_NAMES_FILE), str_split)
verify_separator_in_config_file(file_config_flashpipe.get(CATEGORIAL_TERM_SET_FILE), str_split)

# Recovers plate names from configuration and cleans them up (deletes spaces).
plates_list = file_config_flashpipe.get(PLATE_NAMES_FILE).split(',')
plates_list = [plate.strip() for plate in plates_list]

# Check if the parameter are the type of True or False / Single-cell or Mini Bulk for the analysis
print("•••••Paramater Selection•••••")
verify_method(file_config_flashpipe.get(METHOD), file_config_flashpipe.get(INDEXSORT))
indexsort = verify_parameters(file_config_flashpipe.get(INDEXSORT), INDEXSORT)
bcr = verify_parameters(file_config_flashpipe.get(BCR), BCR)
tcr = verify_parameters(file_config_flashpipe.get(TCR), TCR)
metadata = verify_parameters(file_config_flashpipe.get(METADATA), METADATA)
print("•••••••••••••••••••••••••••••")

# Retrieve category data to be set aside for certain analysis data.
categoriale_term_set_list = file_config_flashpipe.get(CATEGORIAL_TERM_SET_FILE).split(',')
categoriale_term_set_list = [categorial_term.strip() for categorial_term in categoriale_term_set_list]

# Recover reads 1 and 2 from fastQ (to put them in the copier.yml file, in dictionary form).
fastq_files_read1, fastq_files_read2 = fastq_read_files(file_config_flashpipe, plates_list, PATH_RNA)

### 4. Selecting files by species
species = file_config_flashpipe.get(SPECIES_FILE)
star_index = file_config_flashpipe.get(STAR_INDEX_FILE).get(species)
gtf_file = file_config_flashpipe.get(GTF_FILE).get(species)
trust4_imgt_BCR_TCR = file_config_flashpipe.get(TRUST4_IMGT_BCR_TCR).get(species)
trust4_imgt_VDJ = file_config_flashpipe.get(TRUST4_IMGT_VDJ).get(species)

### Define the exact position of barcode and ercc file for analysis (more precisely : zUMIs)
ERCC_PATH = os.path.join(PATH_EXPERIMENT_REFERENCE, 'ERCC_concentration.csv')
BARCODE_WELL_PATH = os.path.join(PATH_EXPERIMENT_REFERENCE, 'cell_barcode_well.csv')
# Define a barcode path (without well ID) for zUMIs (this file is create down).
BARCODE_PATH = os.path.join(PATH_EXPERIMENT_REFERENCE, "cell_barcode.txt")
GSF_FILE_PATH = os.path.join(PATH_EXPERIMENT_REFERENCE, os.path.basename(file_config_flashpipe.get(GSF_FILE)))

### 5. Creation of a dictionary retrieving/containing values from the Config_Flashpipe file
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
    'template_path' : PATH_COPIER_TEMPLATE,
    'path_project' : PATH_PROJECT_FLASH_PIPE,
    'categoriale_term_set' : ','.join(categoriale_term_set_list),
    'method_analysis' : file_config_flashpipe.get(METHOD),
    'index_sort_analysis' : indexsort,
    'bcr_repertoire_analysis' : bcr,
    'tcr_repertoire_analysis' : tcr,
    'metadata_analysis' : metadata
}

### 6. Modify the copier.yml file by adding the values present in the “values_config_flashpipe” dictionary.
FILE_COPIER_TERMPLATE_YML = os.path.join(PATH_COPIER_TEMPLATE, 'copier.yml')
replace_value_config_template_copier(FILE_COPIER_TERMPLATE_YML, values_config_flashpipe)

### 7. Execute the copier.yml file, and consequently the template, to create the necessary structure and files.
Path_template = PATH_COPIER_TEMPLATE
Path_Output = PATH_PROJECT_FLASH_PIPE
cmd = f'copier copy -f {Path_template} {Path_Output}'
print("Launch copier, might take few minutes")
subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

#create a file barcode (without well ID) for zUMIs.
barcode_well = pd.read_csv(BARCODE_WELL_PATH)
barcode_well["BarcodeSequence"].to_csv(BARCODE_PATH, sep=' ', index=False, header=False)

### 8. Creating symbolic links for FASTQ files
create_symlinks(file_config_flashpipe.get(FASTQ_DIRECORIES_FILE), PATH_RNA, plates_list)

### 9. Copy Index Sorting files and gsf.
# Call function to copy Index Sorting files (if not FALSE, which means that the section has been filled (and not empty))
if file_config_flashpipe.get(INDEX_SORT_FILE) != "FALSE" :
    copy_index_sort(file_config_flashpipe.get(INDEX_SORT_FILE), PATH_INDEX_SORTING, plates_list=plates_list)
# Call function to copy GSF
copy_gsf(file_config_flashpipe.get(GSF_FILE), PATH_EXPERIMENT_REFERENCE)

### 10. Check if every file is correctly import by copier or not.
verify_file_exist(ERCC_PATH)
verify_file_exist(BARCODE_WELL_PATH)
verify_file_exist(BARCODE_PATH)
verify_file_exist(os.path.join(PATH_EXPERIENCE, '03_Script/01_FlashPipe/03_QC/analysisParams.R'))
verify_file_exist(os.path.join(PATH_EXPERIENCE, '03_Script/01_FlashPipe/projectParams.R'))
verify_file_exist(os.path.join(PATH_EXPERIMENT_REFERENCE, os.path.basename(file_config_flashpipe.get(GSF_FILE))))
verify_file_exist(os.path.join(PATH_WORKFLOW, '01_snakemake/snakefile.yaml'))



print("Success : The project structure is now in place.")

