ChIP-Seq Nextflow Pipeline
This Nextflow pipeline performs a streamlined ChIP-Seq analysis workflow, including quality control, adapter contamination detection, conditional trimming, alignment, and duplicate marking.

Overview
This pipeline supports:

Quality control of raw reads using FastQC

Conditional trimming with Trim Galore, only if adapter contamination is detected

Alignment with Bowtie2

Samtools

Deduplication using Picard 

Requirements
Nextflow

FastQC

Trim Galore

Bowtie2

Samtools

Picard Tools

