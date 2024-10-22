#!/usr/bin/env nextflow
"""
========================================================
Groovy Helper Functions
========================================================    
"""
def reverse_complement(String s) {
    complement(s.reverse());
}

def complement(String s) {
    def acgt_map = [
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a"
    ];

    char[] sc = new char[s.length()];
    for (int i = 0; i < s.length(); i++) {
        sc[i] = acgt_map[s[i]];
    }
    new String(sc);
}

def remove_all_extensions(String s) {
    s.substring(0, s.indexOf("."));
}

def trim_fn(String s) {    
    if (params.trim_filenames && s.indexOf(params.trim_delimiter) > 0) {
        return s.substring(0, s.indexOf(params.trim_delimiter));
    } else {
        return s;
    }
}

def filter_selex_rounds(String round_name) {
    if (params.round_order == null || params.round_order == "" || params.round_order.size() == 0) return true;
    // // building regex to check files if they match the rounds specified in YOUR_SELEX
    // round_regex = "(" + params.round_order.join('|') + ")" + params.trim_delimiter + ".*"; // TODO replace trim_delimiter with autodetect
    // return s.matches(round_regex);
    
    return params.round_order.contains(round_name);
}

def get_round_id(String round_name) {
    if (params.round_order == null || params.round_order == "" || params.round_order.size() == 0) return 0;
    else return params.round_order.indexOf(round_name);
}

"""
========================================================
Make output directory if it doesn't exist
========================================================    
"""
dir_output = file(params.output_dir)
if (!dir_output.exists()) {
    if (!dir_output.mkdir()) {
        println("Couldn't create output directory.")
    }
}


"""
========================================================
Read FASTQ-files from specified input directory and order by round_order (if present)
========================================================    
"""
fastq_files = Channel.fromFilePairs(params.input_dir + "/" + params.fastq_pattern, checkIfExists:true, type: "file")
fastq_files
    .map { it -> [get_round_id(trim_fn(it[0])), trim_fn(it[0]), it[1][0], it[1][1]	] }
    .filter { filter_selex_rounds(it[1]) }
    .view()	
    .into { fastq_files_ngs_quality; fastq_files_preprocess }
    

(fastq_files_fwd, fastq_files_rev) = fastq_files_ngs_quality.separate(2) { it -> [[it[0], it[2]], [it[0], it[3]]] }
fastq_files_fwd_sorted = fastq_files_fwd
    .toSortedList( { a, b -> a[0] <=> b[0] } )
    .collect { it -> return it.collect { it[1] } } // double collect here, single did not work
    
fastq_files_rev_sorted = fastq_files_rev
    .toSortedList( { a, b -> a[0] <=> b[0] } )
    .collect { it -> return it.collect { it[1] } } // double collect here, single did not work
    

"""
========================================================
NGS quality analysis
========================================================
"""
process ngs_quality_raw {
    echo true
    //conda 'conda-forge:r-argparse conda-forge:r-here conda-forge:r-biocmanager conda-forge:r-rmarkdown conda-forge:r-knitr conda-forge:r-ggplot2'

    publishDir "${params.output_dir}/analysis.ngs_quality",
        mode: "copy"
    
    input:
        file(fastq_fwd) from fastq_files_fwd_sorted
        file(fastq_rev) from fastq_files_rev_sorted
    output:
        file("*.html") into ngs_quality_raw_html
        file("raw_plots/") into ngs_quality_raw_png
    script:
    """
        ngs_quality_report.R \
            --title '${params.experiment}: Raw Files' \
            --out_html ./ngs_quality_raw.html \
            --out_dir ./raw_plots \
            --fwd_reads ${fastq_fwd} \
            --rev_reads ${fastq_rev} 
    """
}


"""
========================================================
Data preparation
========================================================
"""
process trim {
    //conda 'bioconda::cutadapt'
    
    input:
        tuple val(round_id), val(round_name), file(fastq_fwd), file(fastq_rev) from fastq_files_preprocess
    output:
        tuple val(round_id), val(round_name), file("${fastq_fwd.simpleName}.trim.fastq"), file("${fastq_rev.simpleName}.trim.fastq") into trim_keep        
        tuple val(round_id), val(round_name), file("${fastq_fwd.baseName}.fastq"), file("${fastq_fwd.baseName}.trim.fastq") into preprocessing_analysis_trim
        
    script:
    """
        # hacked way to get cutadapt to output the correct log name
        mv ${fastq_fwd} fwd.${fastq_fwd}
        mv ${fastq_rev} ${fastq_fwd}
        
        
        mkdir discarded_fastq
        cutadapt \
                --action=trim \
                -e ${params.trim.max_error_rate} \
                -g ^${params.primer5}...${params.primer3} \
                -G ^${reverse_complement(params.primer3)}...${reverse_complement(params.primer5)} \
                --report=full \
                \
                --untrimmed-output discarded_fastq/${fastq_fwd.simpleName}.no_primers.trim.fastq \
                --untrimmed-paired-output discarded_fastq/${fastq_rev.simpleName}.no_primers.trim.fastq \
                \
                --output ${fastq_fwd.simpleName}.trim.fastq \
                --paired-output ${fastq_rev.simpleName}.trim.fastq \
                fwd.${fastq_fwd} ${fastq_fwd} > ${fastq_fwd.simpleName}.cutadapt.log
                
    #            mv ${fastq_fwd.simpleName}.trim.fastq ${fastq_fwd.simpleName}.fastq
    #            mv ${fastq_rev.simpleName}.trim.fastq ${fastq_rev.simpleName}.fastq
    
        mv ${fastq_fwd} ${fastq_rev}
        mv fwd.${fastq_fwd} ${fastq_fwd}
    """
}


process filter {
    //conda 'bioconda::fastp'
    
    input:
        tuple val(round_id), val(round_name), file(fastq_fwd), file(fastq_rev) from trim_keep
    output:
        tuple val(round_id), val(round_name), file("${fastq_fwd.baseName}.filter.fastq"), file("${fastq_rev.baseName}.filter.fastq") into filter_keep
        tuple val(round_id), file("${fastq_fwd.baseName}.filter.fastq") into preprocessing_analysis_filter
    script:
    """
        fastp -i ${fastq_fwd} -I ${fastq_rev} \
            -o ${fastq_fwd.baseName}.filter.fastq -O ${fastq_rev.baseName}.filter.fastq \
            --disable_adapter_trimming \
            --average_qual $params.filter.min_phred_quality \
            --json=${fastq_fwd.baseName}.filter.fastp.json
    """
}

process merge {
    //conda 'bioconda::fastp'
    publishDir "${params.output_dir}/fastq.prepped.all_lengths",
        pattern: '*merge.fastq',
        saveAs: { filename -> "${trim_fn(remove_all_extensions(filename))}.fastq"},
        mode: "copy"
    publishDir "${params.output_dir}/fasta.prepped.all_lengths",
        pattern: '*merge.fasta',
        saveAs: { filename -> "${trim_fn(remove_all_extensions(filename))}.fasta"},
        mode: "copy"
    
    input:
        tuple val(round_id), val(round_name), file(fastq_fwd), file(fastq_rev) from filter_keep
    output:
        tuple val(round_id), val(round_name), file("${fastq_fwd.baseName}.merge.fastq") into fastp_keep
        tuple val(round_id), val(round_name), file("${fastq_fwd.baseName}.merge.fasta") into fasta_keep
        tuple val(round_id), file("${fastq_fwd.baseName}.merge.fastq") into preprocessing_analysis_merge
    script:
    """
#        mkdir discarded_fastq
#        touch discarded_fastq/${fastq_fwd.baseName}.failed_filter_merge.fastq
#        touch discarded_fastq/${fastq_rev.baseName}.failed_filter_merge.fastq
        
        fastp -i ${fastq_fwd} -I ${fastq_rev} \
            --merge \
            --merged_out=${fastq_fwd.baseName}.merge.fastq \
            --disable_adapter_trimming \
            --json=${fastq_fwd.simpleName}.fastp.json
        
        # remove quality information to convert fastq to fasta
        sed -n 'p;n;p;n;n' ${fastq_fwd.baseName}.merge.fastq | sed 's/@M/>M/g' > ${fastq_fwd.baseName}.merge.fasta
    """
}


process restrict_random_region_length {
    //conda 'conda-forge:pandas'
    publishDir "${params.output_dir}/fastq.prepped",
        pattern: '*fastq',
        //saveAs: { filename -> "${trim_fn(remove_all_extensions(filename))}.fastq"},
        mode: "copy"
    publishDir "${params.output_dir}/fasta.prepped",
        pattern: '*fasta',
        //saveAs: { filename -> "${trim_fn(remove_all_extensions(filename))}.fasta"},
        mode: "copy"
        
    input:
        tuple val(round_id), val(round_name), file(fastq_preprocessed) from fastp_keep
    output:
        tuple val(round_id), val(round_name), file("${round_name}.fastq") into fastq_preprocessed
        tuple val(round_id), val(round_name), file("${round_name}.fasta") into fasta_preprocessed
        tuple val(round_id), val(round_name), file("${round_name}.fasta") into fasta_preprocessed_nt_distribution
        tuple val(round_id), val(round_name), file("${round_name}.fastq") into fastq_preprocessed_ngs_analysis
        
        tuple val(round_id), file("${round_name}.fastq") into preprocessing_analysis_restricted_length
    script:
    """
        restrict_random_region_length.py -i ${fastq_preprocessed} -o ${round_name}.fastq \
            --min ${params.random_region_min} --mino ${round_name}.too_short.fastq \
            --max ${params.random_region_max} --maxo ${round_name}.too_long.fastq
        
        # remove quality information to convert fastq to fasta
        sed -n 'p;n;p;n;n' ${round_name}.fastq | sed 's/@M/>M/g' > ${round_name}.fasta
    """
}

"""
========================================================
NGS Quality Analysis of Preprocessed Files
========================================================
"""

fastq_preprocessed_ngs_analysis_sorted = fastq_preprocessed_ngs_analysis
    .toSortedList( { a, b -> a[0] <=> b[0] } )
    .collect { it -> return it.collect { it[2] } } // double collect here, single did not work

process ngs_quality_analysis_preprocessed {
    echo true
    //conda 'conda-forge:r-argparse conda-forge:r-here conda-forge:r-biocmanager conda-forge:r-rmarkdown conda-forge:r-knitr conda-forge:r-ggplot2'


    publishDir "${params.output_dir}/analysis.ngs_quality",
        mode: "copy"
    
    input:
        file(fastq_file) from fastq_preprocessed_ngs_analysis_sorted
    output:
        file("*.html") into ngs_quality_html_preprocessed
        file("prepped_plots") into ngs_quality_png_preprocessed
    script:
    """
        ngs_quality_report.R  \
            --title '${params.experiment}: Preprocessed Files' \
            --out_html ./ngs_quality_preprocessed.html \
            --out_dir ./prepped_plots \
            --fwd_reads ${fastq_file}
    """
}

"""
========================================================
Preprocessing Analysis
========================================================
"""
preprocessing_analysis_trim
    .join(preprocessing_analysis_filter)
    .join(preprocessing_analysis_merge)
    .join(preprocessing_analysis_restricted_length)
    .set{ preprocessing_analysis }

process preprocessings_analysis {
    input:
        tuple val(round_id), val(round_name), file("raw.fastq"), file("trim.fastq"), file("filter.fastq"), file("merge.fastq"), file("restricted_length.fastq") from preprocessing_analysis
    output:
        tuple val(round_id), file("${round_id}.csv") into preprocessing_analysis_csv_lines
    script:
    """        
        raw_reads=\$((\$(cat raw.fastq | wc -l)/4))
        trim_reads=\$(( \$(cat trim.fastq | wc -l)/4))
        filtered_reads=\$((\$(cat filter.fastq | wc -l)/4))
        merged_reads=\$((\$(cat merge.fastq | wc -l)/4))
        restricted_length_reads=\$((\$(cat restricted_length.fastq | wc -l)/4))


        touch ${round_id}.csv
        echo -n $round_id >> ${round_id}.csv
        echo -n '\t' >> ${round_id}.csv
        echo -n $round_name >> ${round_id}.csv
        echo -n '\t' >> ${round_id}.csv
        echo -n \$raw_reads >> ${round_id}.csv
        echo -n '\t' >> ${round_id}.csv
        echo -n \$trim_reads >> ${round_id}.csv
        echo -n '\t' >> ${round_id}.csv
        echo -n \$filtered_reads >> ${round_id}.csv
        echo -n '\t' >> ${round_id}.csv
        echo -n \$merged_reads >> ${round_id}.csv
        echo -n '\t' >> ${round_id}.csv
        echo -n \$restricted_length_reads >> ${round_id}.csv
        echo >> ${round_id}.csv # new line
    """
}

preprocessing_analysis_csv_sorted = preprocessing_analysis_csv_lines.toSortedList( { a -> a[0] } ).transpose().last().collect()
process preprocessings_analysis_combine_csv {
    //conda 'conda-forge:r-argparse conda-forge:r-dplyr conda-forge:r-tidyr conda-forge:r-ggplot2 conda-forge:r-viridis'
    publishDir "${params.output_dir}/analysis.preprocessing/",
        pattern: "preprocessing*",
        mode: "copy"
    
    input:
        file(f) from preprocessing_analysis_csv_sorted
    output:
        file("preprocessing*.csv") into prep_analysis_csv
        file("preprocessing*.png") into prep_analysis_plots
    script:
    """
        touch preprocessing_analysis.csv
        echo 'round_id\tround_name\traw\ttrim\tfilter\tmerge\trestrict_size' >> preprocessing_analysis.csv
        cat $f >> preprocessing_analysis.csv
        preprocessing_analysis.R -i preprocessing_analysis.csv -o preprocessing.png -p preprocessing.perc.png
    """
}



