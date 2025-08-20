import yaml
import os
import re
import shutil
import sys

# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# General code description :
# 1. open_file_yml : opens a YAML file.
# 2. Function replace_value_config_template_copier : Replaces values present in copier.yml to define project directories and files.
# 3. Identify the path of each read (1 or 2) for each plate in the “template” config file (copier.yml) and create the symlink.
# 4. Function Copy_index_sort : Copies the Index Sorting (of the plates supplied in the FlashPipe Config file) and 
# 5. Function copy_gsf : Copies the GSF file into the new directories.
# 6. Function verify_method : Checks whether the parameter for the method is correct (single-cell or minibulk).
# 7. Function verify_parameter : Check if the provided paramter has the right type of value.
# 8. Function verify_name_experience_path_and_config_file : Checks if the experiment name given in the config file matches the one retrieved by the directory path
# 9. Function verify_file_exist : Checks if the provided file exists.
# 10. Function verify_separator_in_config_file : Check that the separators in the file are "," and not something else.
# 11. Function verify_empty_values_config_file : Checks that all mandarity sections in the user config file (yaml) are not empty.
# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 1. open_file_yml function: opens a YAML file and returns its contents.
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
def open_file_yml(file_config_user_path):
    '''
    Opens and loads a YAML file.
    
    file_config_user_path: Path to the YAML config file
    '''
    with open(file_config_user_path, 'r') as file: 
        return yaml.safe_load(file)

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 2. Function replace_value_config_template_copier: Replaces default values in copier.yml with those in the dictionary.
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
def replace_value_config_template_copier(config_path, dict_new_values):
    '''
    Replaces default values in copier.yml file with those provided in dict_new_values.
    
    config_path: Path to file copier.yml.
    dict_new_values : Dictionary with new values to be inserted in copier.yml.
    '''
    data = open_file_yml(config_path)

    # Retrieves keys and values from the dictionary (containing FlashPipe config values)
    for key, value in dict_new_values.items():
        # If the key (example: “barcode”) is found in the copier.yml file (template config)
        if key in data:
            # The variable to be applied to the key (e.g. “barcode”) will be replaced by the value (values) present in our dictionary.
            data[key]['default'] = value

    # Save changes to file
    with open(config_path, 'w') as file:
        yaml.dump(data, file, default_flow_style=False, allow_unicode=True)

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 3. Identify the path of each read (1 or 2) for each plate in the “template” config file (copier.yml) and create the symlink.
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
def prepare_fastq_symlinks_and_paths(source_dir, output_dir, plates_list):
    """
    Identifies fastq files (R1 and R2) by plate, creates symbolic links, and returns the paths to these links for each read (for the template file).

    dir_entry: Source directory containing the files.
    output_dir: Output directory where symbolic links will be created.
    plates_list: List of plates to be taken into account when creating links.
    """
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Initialize array of Read1 and Read2 fastq files
    fastq_files_read1 = []
    fastq_files_read2 = []

    # Regex to identify FastQ Read1 and Read2 files.
    r1_pattern = re.compile(r'^(?P<plate_name>.*)_.*_R1_.*\.fastq.gz$')
    r2_pattern = re.compile(r'^(?P<plate_name>.*)_.*_R2_.*\.fastq.gz$')

    # Dictionary to track files found by plate
    found_reads = {plate: {'R1': None, 'R2': None} for plate in plates_list}

    # Browse files in the source directory
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            if not file.endswith('.fastq.gz'):
                continue

            file_path = os.path.join(root, file)
            link_path = os.path.join(output_dir, file)

            r1_match = r1_pattern.match(file)
            r2_match = r2_pattern.match(file)

            matched_plate_name = None
            read_type = None

            if r1_match:
                matched_plate_name = r1_match.group('plate_name')
                read_type = 'R1'
            elif r2_match:
                matched_plate_name = r2_match.group('plate_name')
                read_type = 'R2'

            if matched_plate_name in plates_list and read_type:
                found_reads[matched_plate_name][read_type] = link_path

                # Creating symbolic link
                try:
                    os.symlink(file_path, link_path)
                except FileExistsError:
                    print(f"The {link_path} symbolic link already exists.") # Link_name n'existe pas dans cette fonction...
                except Exception as e:
                    print(f"Error creating {link_path} symbolic link: {e}") # Link_name n'existe pas dans cette fonction...

    # Check that each plate has an R1 and an R2
    for plate, reads in found_reads.items():
        if not reads['R1'] or not reads['R2']:
            print(f"ERROR : R1 or R2 files missing for '{plate}' plate")
            sys.exit()

        fastq_files_read1.append(reads['R1'])
        fastq_files_read2.append(reads['R2'])

    return fastq_files_read1, fastq_files_read2

# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 4. copy_index_sort function: copies the Index Sorting (from the plates supplied in the FlashPipe Config file)
# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
def copy_index_sort(source_dir, destination_dir, plates_list=None):
    '''
    Copies files to a destination directory.
    If plates_list is supplied, only the specified plates will be copied.

    source_dir: Source directory containing the files to be copied.
    destination_dir: Destination directory for copied files.
    plates_list: List of plates to be processed for index_sort files.
    '''
    # Create destination directory if it doesn't exist
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    # Definition of regular expression for .csv files with a plate name followed by something
    csv_pattern = re.compile(r'^(?P<plate_name>.*)_.*\.csv$')

    # If plates_list is provided, only the .csv files corresponding to the plates are copied
    if plates_list:
        for plate in plates_list:
            found = False  # Flag to check if at least one file was found for this plate
            for file_name in os.listdir(source_dir):
                # Check if the file matches the pattern and the plate name
                match = csv_pattern.match(file_name)
                if match and match.group('plate_name') == plate:
                    found = True
                    # Retrieves plate IndexSort file (path)
                    csv_file = os.path.join(source_dir, file_name)
                    # Create the output path, to make sure it doesn't already exist.
                    dest_file = os.path.join(destination_dir, file_name)
                    if os.path.exists(dest_file):
                        print(f"Skipped: File already exists in destination : {dest_file}")
                        continue
                    elif os.path.exists(csv_file):
                        try:
                            shutil.copy2(csv_file, destination_dir)
                        except Exception as e:
                            print(f"Error when copying file {csv_file}: {e}")
                    else:
                        print(f"The file {csv_file} does not exist.")
            if not found:
                print(f"ERROR: No indexsort file found for plate : {plate}")
                sys.exit()

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 5. copy_gsf : Copies the GSF file into the new directories.
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

def copy_gsf(source_file, destination_dir):
    '''
    Copies a single file (typically a GSF file) to a destination directory.
    If a file with the same name already exists in the target directory, it will be renamed
    by adding the suffix "_old" before the new file is copied.

    Arguments :
    - source_file : path to the source file to be copied.
    - destination_dir : destination directory.
    '''
    if not os.path.exists(source_file) :
        print("WARNING: The GSF file supplied does not exist.")
        sys.exit()
    elif os.path.isfile(source_file):
        # Recovers only the file name from the source path.
        # Rebuilds the complete path in the destination directory.
        file_name = os.path.basename(source_file)
        file_output_template = os.path.join(destination_dir, file_name)
        # Checks that the source file is not already located in the destination directory.
        if file_output_template != source_file:
            # If a file with the same name already exists in the destination,
            # rename it by adding the suffix “_old” before overwriting the file.
            if os.path.exists(file_output_template):
                base, ext = os.path.splitext(file_output_template)
                old_file = base + "_old" + ext
                os.rename(file_output_template, old_file)

            shutil.copy(source_file, file_output_template)
        # If the source file is already present in the same location with the same name,
        # no action is required (avoids copying a file onto itself).
        elif os.path.exists(file_output_template):
            print( "WARNING: Using already existing GSF file")
    else:
        print("ERROR:", source_file, "is not a file.")
        sys.exit()

# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 6. Function verify_method : Checks whether the parameter for the method is correct (single-cell or minibulk).
# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

def verify_method(parameter, indexsort):
    '''
    Checks whether the parameter for the method is correct (single-cell or minibulk).

    parameter: The value contained in the config file (in this case it will be either single-cell or minibulk).
    indexsort: Here the function is specific to indexsort (as it requires its base value to check for bugs).
    '''
    # If it's not single-cell or minibulk then there's an error
    if parameter != "single-cell" and parameter != "minibulk" :
        print("ERROR : The method (single-cell or minibulk) is erroneous (check the config file for the correct term).")
        sys.exit()
    else : 
        print("You selected method: ", parameter)    
        
    # If method is minibulk, check the status of the indexsort parameter to ensure coherence
    if parameter == "minibulk" : 
        # If indexsort is set to true (then error, in mini bulk it can't have an index sort)
        if indexsort == True : 
            indexsort = False
            print("You set the paramater IndexSort to 'yes', but you selected 'minibulk' in method. The value is now for IndexSort nas been changed to 'no'")
            return indexsort

# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 7. Check if the provided paramter has the right type of value
# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
def verify_parameters(parameter, name_params):
    '''
    Checks whether the type of parameter entered is correct (yes or no).

    parameter: The value contained in the config file (in this case it will be either yes or no).
    name_params: Name of the variable to be checked, to be used for output (print). 
    '''
    # Checks whether True or False (in config file 'yes' or 'no')
    if parameter != True and parameter != False :
        print("ERROR : The parameter", name_params, "is erroneous. It should be 'yes' or 'no' (check the config file to fix the term).")
        sys.exit()
    else :
        # If the value is True, then this means it's selected
        if parameter :
            print("You selected ", name_params, " (Set to ", parameter, ")", sep = '')
            # We return it with 'TRUE' and not 'True' because in R, TRUE is recognized, which facilitates analysis afterwards.
            return "TRUE"
        # Value set to false, then this means it's not selected
        else : 
            print("You do not select ", name_params, " (Set to ", parameter, ")", sep = '')
            # We return it with 'FALSE' and not 'False', because in R it's TRUE that's recognized, which makes it easier to analyze afterwards.
            return "FALSE"

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 8. Checks if the experiment name given in the config file matches the one retrieved by the directory path
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
def verify_name_experience_path_and_config_file(path_experience, name_experience_config):
    '''
    Checks if the experiment name given in the config file matches the one retrieved by the directory path

    path_experience: Path to experiment folder (Example: mnt/DOSI/MASTER_TEMP/project_name/experience_name)
    name_experience_config: Name of experiment given in config (Config_flashpipe) file by user
    '''
    if os.path.basename(path_experience) != name_experience_config :
        print("ERROR : Path to experiment (",os.path.basename(path_experience) ,") is different from that given in config file (", name_experience_config, ").", sep = "")
        sys.exit()

# ••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 9. Checks if the provided file exists.
# ••••••••••••••••••••••••••••••••••••••••••••••••••••
def verify_file_exist(file_path):
    '''
    Checks if the provided file exists.

    file_path: Path to file (check for existence)
    '''
    if not os.path.isfile(file_path) :
        print("ERROR: The file does not exist (", file_path ,"). (Can be cause by the project or experience name).", sep = "")
        sys.exit()

# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 10. Check that the separators in the file are "," and not something else.
# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
def verify_separator_in_config_file(file_config_flashpipe_option, str_split):
    '''
    Check that the separators in the file are "," and not something else.

    file_config_flashpipe_option: Category to be selected in the config file to be checked (example: plate_names)
    str_split: List containing possible separators set by the user (example: ; / | etc...)
    '''
    if any(char in file_config_flashpipe_option for char in str_split):
        print("ERROR: The separator you put in the config file (", file_config_flashpipe_option, ") are not good. \nPut ',' to separate the different names.", sep="")
        sys.exit()

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# ## 11. Checks that all mandarity sections in the user config file (yaml) are not empty.
# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
def verify_empty_values_config_file(file_config_flashpipe):
    '''
    Checks that all mandarity sections in the user config file (yaml) are not empty.
    
    file_config_flashpipe: The yaml file to be checked
    '''
    acc = 0
    skip_keys = set()

    # If index_sort_analysis is set to False, the checks for certain keys are ignored.
    if file_config_flashpipe.get("index_sort_analysis") is False:
        skip_keys.update(["index_sort", "not_fluorescent"])
        
        # Force value to False (just to have a string, because if you don't specify copy doesn't work)
        file_config_flashpipe["index_sort"] = "FALSE"
        file_config_flashpipe["not_fluorescent"] = "FALSE"

    for key in file_config_flashpipe:
        if key in skip_keys:
            continue
        if file_config_flashpipe[key] is None:
            acc += 1
            print("ERROR: Section (", key ,") is empty in the config file. \nPlease provide at least 1 elements.", sep = "")

    if acc > 0:
        sys.exit()

