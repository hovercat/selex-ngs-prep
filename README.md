# selex-ngs-prep
*Selex-ngs-prep* is a NextFlow workflow made for data preparation and quality assessment of next generation sequencing (NGS) files resulting from SELEX experiments.

The pipeline works with demultiplexed files in FASTQ format.
Sequences in the FASTQ files are expected to consist of a 5'-primer, random region, and 3'-primer.

The workflow localizes and trims off flanking regions, filters the sequences by quality and merges paired-end reads.
It also provides extensive information about the SELEX experiment's success and the quality of the NGS run.

FASTA files resulting from *selex-ngs-prep* are required by other SELEX analysis pipelines we developed, it can be used as a standalone analysis workflow to gain first insights in the experiment's data though.


## Prerequisites and Installation
### Conda Environment
A conda-like python environment manager is required.
The workflow was tested using conda 4.9.0, though any conda-like environment manager should work.

Info on the installation of conda can be found here: [conda.io](https://docs.conda.io/)

New packages can be installed to the base environment, which is activated by default.
We suggest to create a new environment to avoid conflicting packages.
#### Optional: Creating a new environment
Create a new environment in which the required packages will be installed.
Be sure to always activate this environment before using the workflow, as seen below.
```bash
# Creating a new environment
conda create -n selex-pipelines
# Activateing the new environment
conda activate selex-pipelines
```

### NextFlow Workflow Manager
Install the latest version of the bioinfomatic workflow manager NextFlow on your device in the environment of your choice.
NextFlow can be obtained from the channel Bioconda.
```bash
# Installing NextFlow
conda install -c bioconda nextflow
```

## Output
<!--- The workflow will put created new files into the folder *output/*.
Preprocessed FASTA files are in put into *output/preprocessed*.
Discarded sequences are put into *output/discarded*.
Plots and charts concerning the success of the SELEX experiment (selex round enrichment, nucleotide distribution along aptamers, nucleotide distribution over SELEX rounds) are put into *output/selex_analysis*.
Plots and charts concerning the quality of the next generation sequencing is put into *output/ngs_quality*. -->

The workflow generates a vast number of different files.

### Preprocessed FASTA/FASTQ Files

The most important output files are the preprocessed (trimmed, filtered, merged) FASTA and FASTQ files, containing the aptamer random regions.
These files should be used for further analysis.

### Aptamer CSV Lists

Two files called *selex.aptamers.csv* and *selex.aptamers.rpm.csv* are created.
They contain one row for every unique random region encountered, along with read counts for every SELEX round.

Values in the rpm-file (rpm is reads per million) are based on the percentage a sequences takes of a SELEX round, multiplicated by 10^6. 

### SELEX Enrichment

Plots and CSV files are created, in which random regions are binned logarithmically by base 2 and base 10.

An example plot: ![SELEX Enrichment](readme_images/selex_success.log2.csv.png)

### Nucleotide Distribution

### NGS Quality


## Usage
### Directory Structure
The workflow is embedded in a directory structure.
By default new directories are created containing resulting files and a cache.

    .
    ├── bin/                # R and Python scripts for analysis
    ├── input_data/         # Directory used by default for input of raw FASTQ files
    ├── output_YOUREXPERIMENT/             
    │   ├── analysis.ngs_quality/       # HTML-Files and plots on sequencing quality
    │   ├── analysis.nt_distribution/   # HTML-Files and plots on nucleotide distribution on aptamers
    │   ├── analysis.preprocessing/     # Plots and charts on loss during preprocessing
    │   ├── analysis.selex_success/     # Plots and charts on aptamers enrichment
    │   ├── fastx*                      # Directories with preprocessed FASTQ/FASTA files
    │   ├── selex.aptamers.csv          # CSV chart containing all aptamers in the experiment
    │   └── selex.aptamers.rpm.csv      # CSV chart containing all aptamers in the experiment, counted in reads per million
    ├── work/               # cache dir containing intermediate files and temporary conda environments
    |   └── ...
    ├── selex-ngs-prep.nf   # The data preparation workflow
    ├── nextflow.config     # Fallback config file (must not be changed)
    ├── create_config.py    # Config creation wizard
    ├── LICENSE
    └── README.md
    
### Execution
Before execution a config file has to be created in which details about the SELEX experiment have to be included.
Below is the command which can be used to execute the workflow.
```bash
nextflow run selex-ngs-prep.nf -c YOUR_SELEX.config
```

Output and Input directory can be changed in the config file.

If you are working on a cluster you probably don't want NextFlow to place intermediate files in the working directory.
The working directory can be changed as seen below:
```bash
nextflow run selex-ngs-prep.nf -c YOUR_SELEX.config -w /scratch/xyz/work
```

## Configuration
selex-ngs-prep **needs** you to define a config file containing information about your experiment.

For a sucessful selex-ngs-prep only a handful of settings have to be set, as seen below.
Create a new text file (or modify one from *config/*) and save it as *YOUR_EXPERIMENT_NAME.config*.
Make sure to choose a proper name for you config file.
```groovy
params {
    fastq_pattern = "*R{1,2}_001.fastq"     // wildcard pattern to find FASTQ files
    round_order = ["R0", "R1", "R2", "R3", "R4"]
    
    // Length of random region
    random_region = 40
    
    // Forward and reverse primers (flanking regions)
    primer5 = "TAGGGAAGAGAAGGACATATGAT"
    primer3 = "TTGACTAGTACATGACCACTTGA"
}
```
The config parameter 'fastq_pattern' has to be a wildcard pattern.
In the case of "`*R{1,2}_001.fastq`" the workflow will look for forward of the pattern `'*R1_001.fastq'` and reverse reads of the pattern `'*R2_001.fastq'`.

### Default Config
The default config-file with all options is shown here:
```groovy
// ngs-prep.config
params {
    input_dir = "./input_data"          // optional
    output_dir = "./output"             // optional
    fastq_pattern = null
    selex_rounds = ["R0", "R1", "R2", "R3", "R4"]    // optional, if null then alphabetical order is used
    
    random_region = null
    random_region_max_deviation = 3     // optional
    
    primer5 = null
    primer3 = null      
    
    trim.max_error_rate = 0.2           // maximum error rate in the primers to be recognized
    filter.min_phred_quality = 30       // default average sequence read quality to include sequence
    merge.base_correction = true        // enables base correctio if paired-end reads mismatch and there are major quality differences
    merge.max_mismatches = 1            // allowed number of mismatches for merging
}


```


## License
    
