#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.protein_file = '/home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/alphavirus_protein_multiseq.fasta'
params.cdhit_similarity_threshold = 0.95

process CDHIT {
    debug true
    label 'CD_HIT'

    input:
    path protein_fasta
    val similarity_threshold

    /*var outfile_name = "cdhit_out_sim_${similarity_threshold_str}.txt"*/

    output:
    path("cdhit_output.${similarity_threshold}"), emit: clstr_fasta
    path("cdhit_output.${similarity_threshold}.clstr"), emit: clstr_file

    script:
    """
    cd-hit \
    -i $protein_fasta \
    -o cdhit_output.$similarity_threshold \
    -c $similarity_threshold
    """
}

process CDHITTOTSV {
    debug true
    publishDir '/home/pathinformatics/epitope_outputs/cdhit_output'

    input:
    path clstr_file
    val similarity_threshold

    output:
    path("cdhit_out_sim_${similarity_threshold}.tsv"), emit: clstr_tsv

    script:
    """
    python3 /cdhit_clstr_to_df.py \
    $clstr_file \
    cdhit_out_sim_${similarity_threshold}.tsv
    """

}

workflow {
    protein_fasta_ch = Channel.fromPath(params.protein_file)
    cdhit_out_ch = CDHIT(protein_fasta_ch, params.cdhit_similarity_threshold)
    CDHITTOTSV(cdhit_out_ch.clstr_file, params.cdhit_similarity_threshold)

}
