#!/usr/bin/env nextflow

params.in = "$baseDir/../workflow-scratchpad/data/alphavirus_proteins/alphavirus_protein_multiseq_test_small.fasta"

alphavirus_fasta = file(params.in)

/*
 * Split a fasta file into multiple files
 */
process Bepipred {

    input:
    file alphavirus_fasta

    println alphavirus_fasta.Name


    script:
    """
    sudo docker run --rm \
     -v $baseDir/../workflow-scratchpad/data/alphavirus_proteins:/iedb_b_cell_io iedb-b-cell:latest /bin/bash -c \
       "python /bcell_standalone/predict_antibody_epitope.py -m Bepipred -f /iedb_b_cell_io/alphavirus_protein_multiseq_test_small.fasta > /iedb_b_cell_io/out.txt"
    """

}

/*
 * Define the workflow
 */
workflow {
    Bepipred(params.in)
}
