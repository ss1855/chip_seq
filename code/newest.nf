#!/usr/bin/env nextflow

params.fastq_reads       = "${baseDir}/fastq_files/*.fastq.gz"
params.output_dir        = "${baseDir}/fastqc_results"
params.trim_dir          = "${baseDir}/trimmed"
//change the index without decoy
params.genome_index ="${baseDir}/genome_index/GRCh38_noalt_decoy_as"
log.info '''
==================================
        CHIP-SEQ Pipeline

==================================
'''

// Step 1: Load FASTQ files and sample ID
reads_ch = Channel
    .fromPath(params.fastq_reads, checkIfExists: true)
    .map { file -> tuple(file.getBaseName().replace('.fastq', ''), file) }
//reads_ch.view()

// Step 2: FastQC on raw reads

process quality_check {
    tag { id }
    cpus 4
    publishDir params.output_dir, mode: 'copy'

    input:
    tuple val(id), path (reads)

    output:
    tuple val(id), path("*.zip")

    script:
    """
    fastqc $reads
    """
}

// Step 3: Check for adapter contamination (warn/fail)
process check_adapter_flag {
    tag { id }
    cpus 1

    input:
    tuple val(id), path (zipfile)

    output:
    stdout

    script:
    """
    if unzip -p $zipfile */fastqc_data.txt | grep -E '>>Adapter Content|>>Per base sequence quality' | grep -E 'pass' > /dev/null; then
        echo -e "$id"
    fi
    """
}


// Step 4: Trim flagged reads
process trimming {
    tag { id }
    cpus 4
    publishDir params.trim_dir, mode: 'move'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_trimmed.fq.gz")

    script:
    """
    trim_galore --cores 4 $reads --output_dir . --basename ${id}_trimmed
    """
}

// Step 5: FastQC on trimmed reads{optional not done in this}
process trimmed_quality_check {
    tag { trimmed.getBaseName() }
    cpus 4
    publishDir params.output_dir, mode: 'copy'

    input:
    path trimmed

    output:
    path "*_fastqc.*"

    script:
    """
    fastqc $trimmed
    """
}

process alignment {
    tag {id}
    cpus 8
    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}.sam")

    script:
    """
    echo 'Aligning $reads with Bowtie2'
    bowtie2 -x ${params.genome_index} -U $reads -S ${id}.sam --threads 8
    """
}


workflow {

    // Step 1: FastQC
    log.info "Performing quality check..."
    qc_out = quality_check(reads_ch)


    // Step 2: Find reads needing trimming
    log.info "Gathering passed_ids..."
    passed_ids = check_adapter_flag(qc_out)
    .map { it.trim() }
    .filter { it }

    passed_ids_set = passed_ids.collect()

    // ✅ Get reads that failed QC (not in passed list)
    reads_failed = reads_ch
        .filter { id, _ -> !passed_ids_set.contains(id) }

    //reads_failed.view()
    reads_passed = reads_ch
        .filter { id, _ -> passed_ids_set.contains(id) }
    //reads_passed.view()

    trim_reads = trimming(reads_failed)
    //trimmed_quality = trimmed_quality_check(trim_reads)
    //trim_reads.view()
    log.info "Merging trimmed + passed reads..."
    cleaned_reads = trim_reads.mix(reads_passed)
    cleaned_reads.view()
    log.info "Permorming alignment using bowtie2..."
    aligned = alignment(cleaned_reads)
    aligned.view()


}
