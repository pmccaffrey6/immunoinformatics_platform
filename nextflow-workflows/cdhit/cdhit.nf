#!/usr/bin/env nextflow

params.protein_file = '/home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/alphavirus_protein_multiseq.fasta'
protein_fasta = file(params.protein_file)

process CDHIT {
    input:
    path protein_fasta

    var similarity_threshold = 0.95
    var similarity_threshold_str =  "$similarity_threshold".replaceAll(/\./, "")
    var outfile_name = "cdhit_out_sim_${similarity_threshold_str}.txt"

    output:
    path "/home/pathinformatics/example_data_for_nextflow/cdhit_outputs/$outfile_name"

    """
    cd-hit \
    -i /home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/$protein_fasta \
    -o /home/pathinformatics/example_data_for_nextflow/cdhit_outputs/$outfile_name \
    -c $similarity_threshold
    """
}

process CDHITTOTSV {
  input:
  path cdhit_tsv

  println cdhit_tsv

}

workflow {
    CDHIT(params.protein_file)
    CDHITTOTSV(CHDHIT.out)
}
