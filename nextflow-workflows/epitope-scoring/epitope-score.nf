#!/usr/bin/env nextflow

params.in = "$baseDir/../workflow-scratchpad/data/alphavirus_proteins/alphavirus_protein_multiseq_test_small.fasta"

alphavirus_fasta = file(params.in)

/*
 * Run Bepipred for B-cell epitope discovery on protein fasta file
 */
process Bepipred {

    input:
    file alphavirus_fasta

    filename =  alphavirus_fasta.Name

    cmd_base = "python /bcell_standalone/predict_antibody_epitope.py -m Bepipred -f /iedb_b_cell_io/${filename} > /iedb_b_cell_io/out.txt"
    cmd = "\"$cmd_base\""
    println cmd

    var full_cmd =
    """
    docker run --rm \
     -v $baseDir/../workflow-scratchpad/data/alphavirus_proteins:/iedb_b_cell_io iedb-b-cell:latest /bin/bash -c \
       $cmd
    """

    println "Running command: $full_cmd"

    output:

    """
    $full_cmd
    """

}

/*
 * Define the workflow
 */
workflow {
    Bepipred(params.in)
}
