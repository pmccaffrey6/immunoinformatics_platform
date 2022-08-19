#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.protein_file = '/home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/alphavirus_protein_multiseq.fasta'
protein_fasta = file(params.protein_file)

process CDHIT {
    input:
    path protein_fasta

    var similarity_threshold = 0.95
    var similarity_threshold_str =  "$similarity_threshold".replaceAll(/\./, "")
    var outfile_name = "cdhit_out_sim_${similarity_threshold_str}.txt"

    output:
    """
    cd-hit \
    -i /home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/$protein_fasta \
    -o /home/pathinformatics/example_data_for_nextflow/cdhit_outputs/$outfile_name \
    -c $similarity_threshold
    """
}

process CDHITTOTSV {

  output:
  """
  python3 /cdhit_clstr_to_df.py \
  /home/pathinformatics/example_data_for_nextflow/cdhit_outputs/$outfile_name \
  /home/pathinformatics/example_data_for_nextflow/cdhit_outputs/cdhit_to_tsv
  """

}

workflow {
    CDHIT(params.protein_file)
    CDHITTOTSV(CHDHIT.out)
}
