#!/usr/bin/env python
import os


######################################################
# Argument parsing functions
######################################################
def parse_string(s, default=None):
    s = s.strip()
    if s is not None and s != "":
        return s
    else:
        return default


def parse_list(s, delimiter=' ', default=None):
    s = s.strip()
    return s.split(delimiter)


def parse_boolean(s, default=None):
    s = s.strip()
    if s.isnumeric():
        return int(s) == 1
    else:
        return default


def parse_float(s, default=None):
    s = s.strip()
    if s.isnumeric():
        return float(s)
    else:
        return default


def parse_int(s, default=None):
    s = s.strip()
    if s.isnumeric():
        return int(s)
    else:
        return default


def make_sane_string(s):
    return ''.join([char for char in s.strip() if char.isalpha() or char.isdigit() or char == ' ']).rstrip()


######################################################
# Starting config creation
######################################################
arguments = dict()

print("Creating a config file for selex-ngs-prep workflow.")
print(
    "Please answer the following questions. If you don't know an answer, you can press enter without an answer to opt for the default value.")

# Experiment Name
while not arguments.get("experiment"):
    arguments["experiment"] = parse_string(input("What's your SELEX experiment's unique identifying name?\n"))
    if not arguments["experiment"]:
        print("There is no default value for the experiment's name.")

# Input Directory
arguments["input_dir"] = parse_string(input("In which directory are the raw files FASTQ-files stored? (default: "
                                            "./input_dir)\n"), "./input_dir")

# Output Directory
arguments["output_dir"] = parse_string(
    input("In which directory to write output files to? (default: './output_{}/)\n".format(make_sane_string(arguments["experiment"]))),
    './output_{}'.format(make_sane_string(arguments["experiment"]))
)

# Check input directory for FASTQ-Files
if os.path.exists(arguments["input_dir"]):
    fastq_files = [f for f in os.listdir(arguments["input_dir"]) if
                   os.path.isfile(os.path.join(arguments["input_dir"], f))
                   and (f.lower().endswith(".fastq") or f.lower().endswith(".fq"))]
    fastq_files.sort()

    print("Files found in the input directory '{}'".format(arguments["input_dir"]))
    print('\n'.join(fastq_files))
else:
    print("Input directory '{}' does not exist.".format(arguments["input_dir"]))

# Ask for round names
while not arguments.get("round_order"):
    arguments["round_order"] = parse_list(
        input(
            "Please provide the list of SELEX rounds in sequential order, separated by a space. (e.g. R0 R2 R4 R6)\n"),
        delimiter=" ")
    if not arguments.get("round_order"):
        print("Please specify the round order of your FASTQ-files!")

# Ask if the user wants to set trim_delimiter and fwd/rev pattern automatically.
auto_detect_pattern = None
while auto_detect_pattern is None:
    auto_detect_pattern = parse_boolean(
        input("Should I try setting file lookup patterns automatically (recommended)? 1 for yes | 0 for no.\t"))

if auto_detect_pattern:
    # Finding trimming delimiter
    delimiters = {}
    for round_name in arguments["round_order"]:
        for file in fastq_files:
            if file.startswith(round_name):
                delimiter = file[len(round_name)]
                delimiters[delimiter] = delimiters.get(delimiter, 0) + 1
                break

    max_delimiter = None
    max_delimiter_count = 0
    for delimiter, count in delimiters.items():
        if count > max_delimiter_count:
            max_delimiter = delimiter
            max_delimiter_count = count

    # Find files of first round to determine forward and reverse read
    round_name = arguments["round_order"][0]
    paired_end_names = []
    for file in fastq_files:
        if file.startswith(round_name):
            paired_end_names.append(file)
    paired_end_names.sort()

    assert len(paired_end_names) == 2
    assert len(paired_end_names[0]) == len(paired_end_names[1])

    file1 = paired_end_names[0]
    file2 = paired_end_names[1]
    start = -1
    end = -1
    for i in range(0, len(file1)):  # Look for the first part of the strings that are not equal
        if file1[i] != file2[i]:
            start = i
        if start > -1 and file1[i] == file2[i]:
            end = i
            break
    pattern = file1[end:len(file1)]

    file1_is_fwd = parse_boolean(
        input(
            "Considering these files:\n{file1} and {file2}\nIs {file1} the forward read? 1 for yes | 0 for no.\t".format(
                file1=file1,
                file2=file2
            )))
    if file1_is_fwd:
        fwd = file1[start:end]
        rev = file2[start:end]
    else:
        fwd = file2[start:end]
        rev = file1[start:end]

    fastq_pattern = "*{{{fwd},{rev}}}{pattern}".format(fwd=fwd, rev=rev, pattern=pattern)
    arguments["fastq_pattern"] = fastq_pattern
    arguments["trim_delimiter"] = max_delimiter

    auto_detect_success = parse_boolean(
        input("Detected search pattern was '{}'. Is this pattern okay? 1 for yes | 0 for no.\t".format(fastq_pattern)))
    auto_detect_success = auto_detect_success and parse_boolean(input(
        "Detected trimming delimiter was '{}'. Is this trim delimiter okay? 1 for yes | 0 for no.\t".format(
            max_delimiter)))

if not auto_detect_pattern or not auto_detect_success:
    # Ask for FASTQ-Pattern
    while not arguments.get("fastq_pattern"):
        arguments["fastq_pattern"] = parse_string(input(
            "Please provide a FASTQ-search pattern, including a wild-card. For example: \nIf you had R0_S1_100_1.fastq as forward read and R0_S1_100_2.fastq as reverse read, your search pattern should look like this: \'*{1,2}.fastq\'\n"))
        if not arguments["fastq_pattern"]:
            print("Please specify a FASTQ-search pattern!")

    # Ask for trim delimiter
    while not arguments.get("trim_delimiter"):
        arguments["trim_delimiter"] = parse_string(input(
            "Please provide a file name trimming delimiter. If you had R0_S0.fastq, the delimiter would be _.\n"))
        if not arguments["trim_delimiter"]:
            print("Please specify a file trimmming delimiter!")

while not arguments.get("primer5"):
    arguments["primer5"] = parse_string(input("Please provide the forward primer (5'-primer)."))
    if not arguments["primer5"]:
        print("Please specify a forward primer!")

while not arguments.get("primer3"):
    arguments["primer3"] = parse_string(input("Please provide the reverse primer (3'-primer)."))
    if not arguments["primer3"]:
        print("Please specify a reverse primer!")

while not arguments.get("random_region"):
    arguments["random_region"] = parse_int(input("How long is the random region of your aptamers?"))
    if not arguments["random_region"]:
        print("Please specify an exact random region length!")
arguments["random_region_min"] = parse_int(
    input("How long should the shortest accepted random region be? (default: random region exact length)"),
    arguments["random_region"])
arguments["random_region_max"] = parse_int(
    input("How long should the longest accepted random region be? (default: random region exact length)"),
    arguments["random_region"])

default_advanced = None
while not default_advanced:
    default_advanced = parse_boolean(input(
        "Will ask for advanced settings now (max. error rate, quality score threshold, merging options).\nGo with default options for them? 1 for yes | 0 for no"))

if not default_advanced:
    arguments["trim.max_error_rate"] = parse_float(input(
        "Primer trimming: Please specify the max allowed error rate expected in flanking primers. (default: 0.2)"), 0.2)
    arguments["filter.min_phred_quality"] = parse_int(
        input("Quality Filter: Please specify the minimum average aptamer quality required. (default: 30)"), 30)
    arguments["merge.base_correction"] = parse_boolean(input(
        "Paired-end read merging: Should erroneous base calls be corrected during merging if possible? 1 for yes | 0 for no (default: yes)"),
        1)
    arguments["merge.max_mismatches"] = parse_int(
        input("Paired-end-read merging: How many mismatches should be allowed for successful merging? (default: 1)"), 1)

# Finished reading config parameters.
print("")
config_file_name = "./{}.config".format(make_sane_string(arguments["experiment"]))
with open(config_file_name, "w") as config_file:
    config_file.write("params {\n")
    for key, value in arguments.items():
        config_file.write("\t")
        if isinstance(value, list):
            config_file.write('{key} = ["'.format(key=key))
            config_file.write('","'.join(value))
            config_file.write('"]\n')
        elif isinstance(value, int):
            config_file.write('{key} = {value}\n'.format(key=key, value=value))
        elif isinstance(value, str):
            config_file.write('{key} = "{value}"\n'.format(key=key, value=value))
        elif isinstance(value, bool):
            if value:
                config_file.write('{key} = true\n')
            else:
                config_file.write('{key} = false\n')
    config_file.write("}\n")

print("Done.")
print("Your config file was saved to '{}'.".format(config_file_name))
