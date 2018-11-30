#!/usr/bin/env nextflow

/*
===============================================================
SLUBI/rnaseq
===============================================================
 RNA-Seq Analysis Pipeline. Started November 2018.
 #### Homepage / Documentation
https://...
 #### Authors
Oskar Karlsson-Lindsjö <Oskar.E.Karlsson@slu.se>
Juliette Hayer <juliette.hayer@slu.se>
---------------------------------------------------------------
*/


// Configurable variables
params.name = false
params.project = false
params.readsPath = false
params.output = 'results'
params.ref = false
params.gff = false
params.samples = false
params.trim_qual = 20


// Validate inputs
if ( params.ref ){
    fasta = file(params.ref)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.ref}"
}
else {
    exit 1, "No reference transcriptome specified!"
}
if ( params.gff ){
    gff = file(params.gff)
    if( !gff.exists() ) exit 1, "GFF file not found: ${params.gff}"
}
else {
    exit 1, "No genome annotation file (GFF3) specified!"
}
if ( params.samples ){
    samples = file(params.samples)
    if( !samples.exists() ) exit 1, "Samples file not found: ${params.samples}"
}
else {
    exit 1, "No samples/conditions tab file specified!"
}


// The reference genome file 
transcriptome_file = file(params.ref)
// the annotation file
gff_file = file(params.gff)
// samples and conditions file
samples_file = file(params.samples)


/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.readsPath + '*_{R1,R2}_001.fastq.gz', size: 2, flat: true)
    .ifEmpty { error "Cannot find any reads matching: ${params.readsPath}" }
    .into { read_pairs_multiqc_raw ; read_pairs_fastp_raw }

/*
* Trimming quality and adapters with fastp
*/
process fastp {
    publishDir "$params.output/qc", mode: 'copy'
    input:
        set val(id), file(read1), file(read2) from read_pairs_fastp_raw

    output:
        set val(id), file("${id}_R1_trimmed.fastq"), file("${id}_R2_trimmed.fastq") into trimmed

    script:
        """
        fastp -i $read1 -o ${id}_R1_trimmed.fastq -I $read2 -O ${id}_R2_trimmed.fastq \
            --detect_adapter_for_pe --qualified_quality_phred $params.trim_qual \
            --cut_by_quality5 --cut_by_quality3 --cut_mean_quality $params.trim_qual \
            --html ${id}_fastp_report.html
        """
}
trimmed.into {trimmed_reads_pe ; trimmed_pe_trinity ; trimmed_pe_salmon}

process fastqc {
    input:
        set val(id), file(read1), file(read2) from read_pairs_multiqc_raw
        set val(id), file(trimmed_read1), file(trimmed_read2) from trimmed_reads_pe

    output:
        file "*_fastqc.{zip,html}" into fastqc_results

    script:
        """
        fastqc -t 4 $read1 $read2 $trimmed_read1 $trimmed_read2
        """
}

process multiqc {
    publishDir "$params.output/qc", mode: 'copy'
    input:
        file 'fastqc/*' from fastqc_results.collect()

    output:
        file 'multiqc_report.html'

    script:
        """
        multiqc .
        """
}

process salmon {
    publishDir params.output, mode: 'copy'
    input:
        file fasta from transcriptome_file
        set val(id), file(trimmed_read1), file(trimmed_read2) from trimmed_pe_salmon

    output:
        file("salmon/${id}") into salmon_quant

    script:
    """
    salmon index -t $fasta -i punivalens_trans_index
    salmon quant -i punivalens_trans_index --libType A \
          -1 ${trimmed_read1} -2 ${trimmed_read2} -o salmon/${id}
    """
}

process deseq2 {
    publishDir params.output, mode: 'copy'
    input:
        file gff from gff_file
        file samples_cond from samples_file
        set val(id), file(quant) from salmon_quant

    output:
        file("deseq2/") into deseq2

    script:
    """
    #!/usr/bin/env Rscript

    library(tximport)
    library(GenomicFeatures)
    library(readr)

    # importing annotation from the gff3 file
    txdb <- makeTxDbFromGFF($gff)
    k <- keys(txdb, keytype = "GENEID")
    tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")

    # importing Salmon count data using tximport
    samples <- read.table('$samples_cond', header = TRUE)
    files <- file.path("salmon", samples$sample, "quant.sf")
    names(files) <- paste0(samples$sample)
    txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

    # DESeq2
    library(DESeq2)
    dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~condition)
    dds <- DESeq(dds)
    res <- results(dds)

    sum <- summary(res)
    capture.output(sum, file = "res_summary.txt")

    """
}
