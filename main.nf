#!/usr/bin/env nextflow

params.original_fastqs = "$baseDir/data/original_fastqs/*"
params.generated_fastqs = "$baseDir/data/generated_fastqs/*"
params.outdir = 'results'

Channel
    .fromPath(params.original_fastqs, checkIfExists: true)
    .map { file -> tuple(file.name, file) }
    .set { original_fastqs }

Channel
    .fromPath(params.generated_fastqs, checkIfExists: true) //checks whether the specified file exists
    .map { file -> tuple(file.name, file) }
    .set { generated_fastqs }


/*
 * Step 0: Update headers of original fastqs to match the generated ones.
 */
process update_headers {
    tag "$name"

    input:
    set val(name), file(fastq) from original_fastqs

    output:
    set val(name), file('*updated_header.fq.gz') into original_fastqs_updated_header

    """
    zcat $fastq | awk '{if (NR%4==1) { split(\$0,a," "); split(a[2],b,":"); print a[1]"/"b[1]} else print \$0 }' > $fastq".updated_header.fq"
    gzip $fastq".updated_header.fq"
    """
}

generated_fastqs.mix(original_fastqs_updated_header)
                .set { to_be_sorted_fastqs }

/*
 * Step 1: sort fastqs
 */
process sort_fastqs {
    tag "$name"

    input:
    set val(name), file(fastq) from to_be_sorted_fastqs

    output:
    set val(name), file('*sorted.fq') into sorted_fastqs

    """
    zcat $fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > $fastq".sorted.fq"
    """
}

/*
 * Step 2: Calculate checksums for the sorted fastqs
 */
process calculate_checksums {
    publishDir params.outdir, mode: 'copy'

    input:
    file fastq from sorted_fastqs.collect()

    output:
    file('checksums.md5')

    """
    md5sum *sorted.fq > checksums.md5
    """
}

