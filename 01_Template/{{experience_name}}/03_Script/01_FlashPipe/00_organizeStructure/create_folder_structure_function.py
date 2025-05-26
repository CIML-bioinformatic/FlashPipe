import yaml
import os
import re
import shutil
import sys

# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
# General code description :
# 1. open_file_yml function: opens a YAML file.
# 2. Function replace_value_config_template_copier : Replaces values present in copier.yml to define project directories and files.
# 3. Function fastq_read_files: Orders and creates FASTQ file paths (symbolic links) according to plate name and read type (1 or 2). (The symbolic links will be created later, this function allows you to place the paths in the necessary config files).
# 4. Function create_symlinks : Creates symbolic links for FastQ files associated with specified plates. 
# 5. Copy_files function: copies the Index Sorting (of the plates supplied in the FlashPipe Config file) and the GSF File into the new directories.
# ••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

### 1. open_file_yml function: opens a YAML file and returns its contents.
def open_file_yml(file_config_user_path):
    '''
    Opens and loads a YAML file.
    '''
    with open(file_config_user_path, 'r') as file: 
        return yaml.safe_load(file)

### 2. Function replace_value_config_template_copier: Replaces default values in copier.yml with those in the dictionary.
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

### 3. Fastq_read_files function: orders and creates FASTQ file paths (symbolic links) according to plate names and read type (1 or 2). 
def fastq_read_files(file_config_flashpipe, plates_list, PATH_SYM_OUTPUT) :
    '''
    Set the path of each read (1 or 2) for each plate in the “template” config file (copier.yml).
    
    file_config_flashpipe : Open yml file of the config (user).
    plates_list: List of plates to be taken into account when creating links.
    PATH_SYM_OUTPUT: Path to directory containing future symbolic link positions.
    '''
    
    fastq_files_read1, fastq_files_read2 = [], []
    
    # Regex to identify FastQ Read1 and Read2 files.
    r1_pattern = re.compile(r'^(?P<plate_name>.*)_.*_R1_.*\.fastq.gz$')
    r2_pattern = re.compile(r'^(?P<plate_name>.*)_.*_R2_.*\.fastq.gz$')
    
    ## Browses plates and files to associate the correct Read1 and Read2 with the plates.
    # Browses the plate names in the configuration file.
    for plate in plates_list:
        # Browse files in fastq directory (provided by user)
        for file in os.listdir(file_config_flashpipe.get('fastq_directories')):
            # Checks files for regex.
            if r1_pattern.match(file):
                matched_plate_name = r1_pattern.match(file).group('plate_name')
                if matched_plate_name == plate:
                    # Adds the path of future fastq symbolic links.
                    fastq_files_read1.append(os.path.join(PATH_SYM_OUTPUT, file))

            if r2_pattern.match(file):
                matched_plate_name = r2_pattern.match(file).group('plate_name')
                if matched_plate_name == plate:
                    fastq_files_read2.append(os.path.join(PATH_SYM_OUTPUT, file))
                    
    return fastq_files_read1, fastq_files_read2

### 4. create_symlinks function: Creates symbolic links for files associated with the specified plates.
def create_symlinks(dir_entry, output_dir, plates_list):
    '''
    Creates symbolic links for FastQ files associated with plates specified in plates_list.
    
    dir_entry: Source directory containing the files.
    output_dir: Output directory where symbolic links will be created.
    plates_list: List of plates to be taken into account when creating links.
    '''
    # Create output directory if none exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Regex to identify FastQ Read1 and Read2 files.
    r1_pattern = re.compile(r'^(?P<plate_name>.*)_.*_R1_.*\.fastq.gz$')
    r2_pattern = re.compile(r'^(?P<plate_name>.*)_.*_R2_.*\.fastq.gz$')
    
    # Dictionnaire pour suivre les fichiers trouvés par plaque
    found_files_per_plate = {plate: False for plate in plates_list}
    
    # Browse files in the source directory
    for root, dirs, files in os.walk(dir_entry):
        for file in files:
            if file.endswith('.fastq.gz'):  # Filtering on FastQ files
                r1_match = r1_pattern.match(file)
                r2_match = r2_pattern.match(file)

                matched_plate_name = None
                if r1_match:
                    matched_plate_name = r1_match.group('plate_name')
                elif r2_match:
                    matched_plate_name = r2_match.group('plate_name')

                if matched_plate_name and matched_plate_name in plates_list:
                    found_files_per_plate[matched_plate_name] = True
                    source_file = os.path.join(root, file)
                    link_name = os.path.join(output_dir, file)

                    # Creating symbolic link
                    try:
                        os.symlink(source_file, link_name)
                    except FileExistsError:
                        print(f"The {link_name} symbolic link already exists.")
                    except Exception as e:
                        print(f"Error creating {link_name} symbolic link: {e}")

    # Checks which plates have no associated files
    for plate, found in found_files_per_plate.items():
        if not found:
            print(f"ERROR: No fastq files found for the plate : {plate}")
            sys.exit()

### 5. copy_index_sort function: copies the Index Sorting (from the plates supplied in the FlashPipe Config file) and the GSF File into the new directories.
def copy_index_sort(source_dir, destination_dir, plates_list=None):
    '''
    Copies files to a destination directory.
    If plates_list is supplied, only the specified plates will be copied.

    source_dir: Source directory containing the files to be copied.
    destination_dir: Destination directory for copied files.
    plates_list: List of plates to be processed for index_sort files. Optional parameter. If provided, then specifically for index_sort.
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
                    csv_file = os.path.join(source_dir, file_name)
                    if os.path.exists(csv_file):
                        try:
                            shutil.copy2(csv_file, destination_dir)
                        except Exception as e:
                            print(f"Error when copying file {csv_file}: {e}")
                    else:
                        print(f"The file {csv_file} does not exist.")
            if not found:
                print(f"ERROR: No indexsort file found for plate : {plate}")
                sys.exit()
                        

def copy_gsf(source_file, destination_dir):
    '''
    Single file will be copied (in our case, it will be used to copy the GSF).
    (Can give a directory, or the link)

    source_file: Source directory containing the files to be copied.
    destination_dir: Destination directory for copied files.
    '''
    if os.path.isfile(source_file) :
        file_name = os.path.basename(source_file)
        shutil.copy(source_file, os.path.join( destination_dir, file_name))
    else :
        print("ERROR:", source_file, "is not a file.")
        sys.exit()

def verify_method(parameter, indexsort):
    '''
    Checks whether the type of parameter entered is correct (single-cell or minibulk).

    parameter: The value contained in the config file (in this case it will be either single-cell or minibulk).
    indexsort: Here the function is specific to indexsort (as it requires its base value to check for bugs).
    '''
    # Si le paramètre est minibulk
    if parameter == "minibulk" : 
        # If indexsort is set to true (then error, in mini bulk it can't have an index sort)
        if indexsort == True : 
            indexsort = "FALSE"
            print("You set the paramater IndexSort to 'yes', but you selected 'minibulk' in method. The value is now, fixed to 'no'")
            return indexsort
    # If it's not single-cell or minibulk then there's an error
    if parameter != "single-cell" and parameter != "minibulk" :
        print("ERROR : The method (single-cell or minibulk) is erroneous (check the config file for the correct term).")
        sys.exit()
    else : 
        print("You selected", parameter)

def verify_parameters(parameter, name_params):
    '''
    Checks whether the type of parameter entered is correct (yes or no).

    parameter: The value contained in the config file (in this case it will be either yes or no).
    name_params: Name of the variable to be checked, to be used for output (print). 
    '''
    # Checks whether True or False (in config file 'yes' or 'no')
    if parameter != True and parameter != False :
        print("ERROR : The paramater", name_params, "is erroneous (check the config file for the correct term).")
        sys.exit()
    else :
        # If the value is True, then this means it's selected
        if parameter :
            print("You selected ", name_params, " (Set to ", parameter, ")", sep = '')
            # We return it with 'TRUE' and not 'True' because in R, TRUE is recognized, which facilitates analysis afterwards.
            return "TRUE"
        # Special case: if the parameter is already set to FALSE, it's because the sort index was set by the previous function (because of the mini bulk).
        elif parameter == "FALSE":
            print("IndexSort (FACS) is fixed to FALSE")
        # Value set to false, then this means it's not selected
        else : 
            print("You not selected ", name_params, " (Set to ", parameter, ")", sep = '')
            # We return it with 'FALSE' and not 'False', because in R it's TRUE that's recognized, which makes it easier to analyze afterwards.
            return "FALSE"

def verify_name_experience_path_and_config_file(path_experience, name_experience_config):
    '''
    Checks if the experiment name given in the config file matches the one retrieved by the directory path

    path_experience: Path to experiment folder (Example: mnt/DOSI/MASTER_TEMP/project_name/experience_name)
    name_experience_config: Name of experiment given in config (Config_flashpipe) file by user
    '''
    if os.path.basename(path_experience) != name_experience_config :
        print("ERROR : Path to experiment (",os.path.basename(path_experience) ,") is different from that given in config file (", name_experience_config, ").", sep = "")
        sys.exit()

def verify_file_exist(file_path):
    '''
    Checks if the config file exists (this command can also be used to check the existence of other files).

    file_path: Path to file (check for existence)
    '''
    if not os.path.isfile(file_path) :
        print("ERROR: The file does not exist (", file_path ,"). (Can be cause by the project or experience name).", sep = "")
        sys.exit()

def verify_separator_in_config_file(file_config_flashpipe_option, str_split):
    '''
    Check that the separators in the file are , and not something else.

    file_config_flashpipe_option: Category to be selected in the config file to be checked (example: plate_names)
    str_split: List containing possible separators set by the user (example: ; / | etc...)
    '''
    if any(char in file_config_flashpipe_option for char in str_split):
        print("ERROR: The separator you put in the config file (", file_config_flashpipe_option, ") are not good. \nPut ',' to separate the different names.", sep="")
        sys.exit()


def verify_empty_values_config_file(file_config_flashpipe):
    '''
    Checks that all sections in the supplied (yaml) file are not empty.
    
    file_config_flashpipe: The yaml file to be checked
    '''
    acc = 0
    skip_keys = set()

    # If index_sort_analysis is set to False, the checks for certain keys are ignored.
    if file_config_flashpipe.get("index_sort_analysis") is False:
        skip_keys.update(["index_sort", "categoriale_term_set"])
        
        # Force value to False (just to have a string, because if you don't specify copy doesn't work)
        file_config_flashpipe["index_sort"] = "FALSE"
        file_config_flashpipe["categoriale_term_set"] = "FALSE"

    for key in file_config_flashpipe:
        if key in skip_keys:
            continue
        if file_config_flashpipe[key] is None:
            acc += 1
            print("ERROR: Section (", key ,") is empty in the config file. \nPlease provide at least 1 elements.", sep = "")

    if acc > 0:
        sys.exit()

