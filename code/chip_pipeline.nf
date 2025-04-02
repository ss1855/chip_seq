params.fastq_reads = "$baseDir/fastq_files/*.fastq.gz"
params.output_dir = "$baseDir/fastqc_results"
params.trim_dir = "$baseDir/trimmed"

//println "dir: $params.fastq_reads"


log.info '''
==================================
        CHIP-SEQ Pipeline
==================================
'''

reads_ch = Channel.fromPath(params.fastq_reads, checkIfExists: true)
reads_ch.view()

process quality_check {
    tag "{$reads}"
    cpus 4
    publishDir params.output_dir, mode: 'copy'

    input:
    path reads

    output:
    path "*_fastqc.*"

    script:
    """
    fastqc $reads
    """
}
 process trimming{
 	tag "{$reads}"
  cpus 4
  publishDir params.trim_dir, mode:'move'
  input:
  path reads
  output:
  path "*_trimmed.*"
  script:
  """
  trim_galore --cores 4  $reads

  """
 	}

workflow {
    quality_check(reads_ch)
    trimming(reads_ch)
}
