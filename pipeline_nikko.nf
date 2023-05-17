//Pipline to analyse NGS data (single read) from fetching the data of the NCBI to the assembly
// Created by Nikolas Reppert   
// Creation date 09.05.2023

// Comment to Bilal: The user repository of the sra-tools is set to: SRA_data
nextflow.enable.dsl = 2

params.accession = null
params.outdir = "SRA_data"
params.kmerlen = 71
params.cov_cutoff = 0

process PREFETCH {
    //storeDir is needed if data are downloaded and not created by a process (for that you need publishDir)
    storeDir "${params.outdir}"

    container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"

    input:
    val accession
    // the output needs to be into an already existing directory (here: SRA_data/sra)
    // emit is like an adress to the output data for the workflow, is more imported if there are multiple output data
    output:
    path "sra/${accession}.sra" , emit: dotsra

    //in script use the $variable the ${variable} does not seem to work? only in " " you need ${} like f-String in python
    script:
    """
    prefetch $accession
    """
}

process CONVERTTOFASTQ {
    //here it needs to be storeDir as well, I dont know why exactly
    storeDir "${params.outdir}/sra/${accession}/fastq"
    
    container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"

//input variable is defined by workflow
    input:
    val accession
    path srafile
// *.fastq ==> everything that comes out of fastq-dump saved as fastq file
    output:
    path "*.fastq" , emit: fastqfileadress

// split file is actually only needed if the sra data are of a paired read (rev and for sequence by synthesis),
// which is than put into two different fastq files
    script:
    """
    fastq-dump --split-files $srafile
    """
}

process QUALITYCHECK{
    storeDir "${params.outdir}/sra/${params.accession}/fastqc"
// use next singularity image for fastqc
    container "https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0"

    input:
    path fastqfile
// works this way, I do not 100% understand why
    output:
    path "*.html", emit: fastQC_html
    path "*.zip", emit: fastQC_zip

    script:
    // fastqc produces an zip and an html file with name accession_fastqc
    """
    fastqc $fastqfile 
    """
}

process TRIMM{
    // I guess storDir makes the file be stored only in this directory and not in work as well
    storeDir "${params.outdir}/sra/${params.accession}/fastp"
   // use next singularity image for fastqc
    container "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0"

    input:
    path fastqfile


    output:
    path "*.fastq", emit: fastpfile
    path "*.html", emit: report

    script:
    """
    fastp -i $fastqfile -o after_fastp.fastq -h after_fastp.html
    """
}

process QUALITYCHECK2{
    storeDir "${params.outdir}/sra/${params.accession}/fastqc_after_fastp"
// use next singularity image for fastqc
    container "https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0"

    input:
    path fastq_trimm
// works this way, I do not 100% understand why
    output:
    path "*.html", emit: fastQC_html
    path "*.zip", emit: fastQC_zip

    script:
    // fastqc produces an zip and an html file with name accession_fastqc
    """
    fastqc $fastq_trimm
    """
}

process MULTIQC{
    storeDir "${params.outdir}/sra/${params.accession}/multiqc"

    container "https://depot.galaxyproject.org/singularity/multiqc%3A1.13a--pyhdfd78af_1"

    input:
    path report_list
   

    output:
    path "*"

    script:
    """
    multiqc .
    """

}

process velveth{
    storeDir "${params.outdir}/sra/${params.accession}/velvet"

    container "https://depot.galaxyproject.org/singularity/velvet%3A1.2.10--h7132678_5"

    input:
    path fastq_trimm

    output:
    path "hash"

    script:
    """
    velveth hash ${params.kmerlen} -fastq ${fastq_trimm}
    """
}

process velvetg{
    storeDir "${params.outdir}/sra/${params.accession}/velvet"

    container "https://depot.galaxyproject.org/singularity/velvet%3A1.2.10--h7132678_5"

    input:
    path hash_data

    output:
    path "hash"

    script:
    """
    velvetg hash -cov_cutoff ${params.cov_cutoff}
    """
}

workflow{
    // the emit "adresses" go to the process they were used in
    srafile=PREFETCH(params.accession).dotsra
    fastqfile=CONVERTTOFASTQ(params.accession, srafile).fastqfileadress
    //there would be a function like flatten8) and collect() necessary, if data of paired read were used
    fastqcfile_flat=QUALITYCHECK(fastqfile)
    fastq_trimm= TRIMM(fastqfile).fastpfile
    fastqc_after_trim_flat=QUALITYCHECK2(fastq_trimm)
    // the multiqc function needs a list of channels, so here I creat an empty channel
    // and add the two fastqc channels by "concat"
    reportQC= Channel.empty()
    reportQC_2=reportQC.concat(fastqcfile_flat,fastqc_after_trim_flat)
    report_list = reportQC_2.collect()
    MULTIQC(report_list)
    hash_data= velveth(fastq_trimm)
    velvetg(hash_data)
}