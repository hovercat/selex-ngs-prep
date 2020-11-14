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
    if (params.trim_filenames) {
        return s.substring(0, s.indexOf(params.trim_delimiter));
    } else {
        return s;
    }
}


fastq_files = Channel.fromFilePairs(params.input_dir + "/" + params.fastq_pattern, checkIfExists:true, type: "file")
fastq_files
    .into { fastq_files_ngs_quality; fastq_files_preprocess }

"""
========================================================
NGS Quality Analysis
========================================================
"""
process ngs_quality {
    conda 'bioconda::fastp'
    publishDir 'output/ngs_quality',
        pattern: "*.html",
        mode: "copy"
    
    input:
        tuple val(identifier), file(files) from fastq_files_ngs_quality
    output:
        file("*.json") into multiqc_afterqc
        file("*.html") into ngs_quality_html
    script:
    fastq_fwd = files[0]
    fastq_rev = files[1]
    """
        #after.py -1 $fastq_fwd -2 $fastq_rev --qc_only
        #fastp -1 $fastq_fwd -2 $fastq_rev --qc_only
        #mv QC/* .
        
        ngs_quality.R -1 $fastq_fwd -2 $fastq_rev
    """
}


"""
========================================================
Data preparation
========================================================
"""
process trim {
    conda 'bioconda::cutadapt'
    publishDir 'output/',
        pattern: 'discarded_fastq/*.fastq',
        mode: "copy"
    
    input:
        tuple val(identifier), file(files) from fastq_files_preprocess
    output:
        tuple val(identifier), file("${fastq_fwd.simpleName}.fastq"), file("${fastq_rev.simpleName}.fastq") into trim_keep
        file("discarded_fastq/*") into trim_discard
        file("${fastq_fwd.simpleName}.cutadapt.log") into multiqc_cutadapt
        
    script:
    fastq_fwd = files[0]
    fastq_rev = files[1]
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
            --minimum-length ${params.random_region - params.random_region_max_deviation} \
            --maximum-length ${params.random_region + params.random_region_max_deviation} \
            \
            --untrimmed-output discarded_fastq/${fastq_fwd.simpleName}.missing_primers.trim.fastq \
            --untrimmed-paired-output discarded_fastq/${fastq_rev.simpleName}.missing_primers.trim.fastq \
            --too-short-output discarded_fastq/${fastq_fwd.simpleName}.too_short.trim.fastq \
            --too-short-paired-output discarded_fastq/${fastq_rev.simpleName}.too_short.trim.fastq \
            --too-long-output discarded_fastq/${fastq_fwd.simpleName}.too_long.trim.fastq \
            --too-long-paired-output discarded_fastq/${fastq_rev.simpleName}.too_long.trim.fastq \
            \
            --output ${fastq_fwd.simpleName}.trim.fastq \
            --paired-output ${fastq_rev.simpleName}.trim.fastq \
            fwd.${fastq_fwd} ${fastq_fwd} > ${fastq_fwd.simpleName}.cutadapt.log
            
            mv ${fastq_fwd.simpleName}.trim.fastq ${fastq_fwd.simpleName}.fastq
            mv ${fastq_rev.simpleName}.trim.fastq ${fastq_rev.simpleName}.fastq
    """
}


process filter_merge {
    conda 'bioconda::fastp'
    publishDir 'output/',
        pattern: 'discarded_fastq/*.fastq',
        mode: "copy"
    
    input:
        tuple val(identifier), file(fastq_fwd), file(fastq_rev) from trim_keep
    output:
        tuple val(identifier), file("${fastq_fwd.baseName}.filter.merge.fastq") into fastp_keep
        file("discarded_fastq/*") into fastp_discarded
        file("${fastq_fwd.simpleName}.fastp.json") into multiqc_fastp_filter
    script:
    """
        mkdir discarded_fastq
        touch discarded_fastq/${fastq_fwd.baseName}.failed_filter_merge.fastq
        touch discarded_fastq/${fastq_rev.baseName}.failed_filter_merge.fastq
        
        fastp -i ${fastq_fwd} -I ${fastq_rev} \
            --merge \
            --merged_out=${fastq_fwd.baseName}.filter.merge.fastq \
            --disable_adapter_trimming \
            --average_qual $params.filter.min_phred_quality \
            --json=${fastq_fwd.simpleName}.fastp.json
            
            #-o ${fastq_fwd.baseName}.filter.merge.fastq \
            #-O ${fastq_rev.baseName}.filter.merge.fastq \
#             --unpaired1=discarded_fastq/${fastq_fwd.baseName}.low_quality.filter.fastq \
 #           --unpaired2=discarded_fastq/${fastq_rev.baseName}.low_quality.filter.fastq \
    """
}


/*process merge {
    conda 'bioconda::fastp'
    publishDir 'output/',
        pattern: 'discarded_fastq/*.fastq',
        mode: "copy"
        
    input:
        tuple val(identifier), file(fastq_fwd), file(fastq_rev) from filter_keep
    output:
        file("${fastq_fwd.baseName}.merge.fastq") into merge_keep
        file("discarded_fastq/*") into merge_discard
        file("${fastq_fwd.simpleName}.merge.fastp.json") into multiqc_fastp_merge
    script:
    """
        mkdir discarded_fastq
        fastp -i ${fastq_fwd} -I ${fastq_rev} \
            -o discarded_fastq/${fastq_fwd.baseName}.failed.merge.fastq \
            -O discarded_fastq/${fastq_rev.baseName}.failed.merge.fastq \
            --merge \
            --merged_out=${fastq_fwd.baseName}.merge.fastq \
            --disable_adapter_trimming \
            --json=${fastq_fwd.simpleName}.merge.fastp.json
    """
}*/


process publish_preprocessed_files {
    publishDir 'output/preprocessed_fastq',
        pattern: '*.fastq',
        saveAs: { filename -> "${trim_fn(remove_all_extensions(filename))}.fastq"},
        mode: "copy"
    
    publishDir 'output/preprocessed_fasta',
        pattern: '*.fasta',
        saveAs: { filename -> "${trim_fn(remove_all_extensions(filename))}.fasta"},
        mode: "copy"
        
    input:
        tuple val(identifier), file(fastq_preprocessed) from fastp_keep
    output:
        file(fastq_preprocessed) into fasta_preprocessed
        file("${fastq_preprocessed.baseName}.fasta") into fastq_preprocessed
    script:
    """
        # remove quality information to convert fastq to fasta
        sed -n 'p;n;p;n;n' $fastq_preprocessed | sed 's/@M/>M/g' > ${fastq_preprocessed.baseName}.fasta
    """
}



"""
========================================================
MultiQC of Preprocessing
========================================================
"""
//multiqc_afterqc
 //   .join(multiqc_cutadapt)
  //  .join(multiqc_fastp_filter)
   // .join(multiqc_fastp_merge)
    //.set { multiqc_input }


process multiqc {
    conda 'conda_envs/multiqc.yaml'
    
    input:
        //file(afterqc) from multiqc_afterqc.collect()
        file(cutadapt) from multiqc_cutadapt.collect()
        file(fastp_filter) from multiqc_fastp_filter.collect()
        //file(fastp_merge) from multiqc_fastp_merge.collect()
        //tuple val(identifier), file(afterqc), file("cutadapt.log"), file("fastp_filter.json"), file("fastp_merge.json") from multiqc_input
    output:
        file("*") into multiqc_output
    script:
    """
        multiqc .
    """
}

"""
========================================================
Analysing SELEX Success
========================================================
"""
