#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.protein_file = '/home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/alphavirus_protein_multiseq.fasta'
params.cdhit_similarity_threshold = 0.95

params.bepipred = "no"
params.epidope = "no"
params.netmhcpani = "no"
params.netmhcpanii = "no"
params.dc_bcell = "no"

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

process EPIDOPE {
    debug true
    publishDir '/home/pathinformatics/epitope_outputs/'

    input:
    path protein_fasta

    output:
    path("epidope_output"), emit: epidope_output

    script:
    """
    -p 12 \
    -i $protein_fasta \
    -o epidope_output
    """
}

process BEPIPRED {
    debug true

    input:
    path protein_fasta

    output:
    path("bepipred_output"), emit: bepipred_output

    script:
    """
    /bin/bash -c "sed '/^[^>]/s/[B|X|J|U]//g' $protein_fasta > cleaned.fasta \
    && python /bcell_standalone/predict_antibody_epitope.py -m Bepipred -f cleaned.fasta > bepipred_output"
    """
}

process BEPIPREDTOTSV {
    debug true
    publishDir '/home/pathinformatics/epitope_outputs/bepipred_output'

    input:
    path bepipred_output

    output:
    path("bepipred_out.tsv"), emit: bepipred_tsv

    script:
    """
    python3 /bepipred_to_df.py \
    $bepipred_output \
    bepipred_out.tsv
    """
}

process NETMHCPANI {
    debug true
    publishDir '/home/pathinformatics/epitope_outputs/netmhcpan_i_output'

    input:
    path protein_fasta

    output:
    path("netmhcpan_i_out.tsv"), emit: netmhcpan_i_tsv

    script:
    """
    /bin/bash -c "sed '/^[^>]/s/[B|X|J|U]//g' $protein_fasta > cleaned.fasta \
    && /mhc_i/src/predict_binding.py netmhcpan_ba \
    HLA-A*01:01,HLA-A*02:01,HLA-A*03:01,HLA-A*24:02,HLA-A*26:01,HLA-B*07:02,HLA-B*08:01,HLA-B*15:01,HLA-B*27:05,HLA-B*39:01,HLA-B*40:01,HLA-B*58:01 \
    9,9,9,9,9,9,9,9,9,9,9,9 \
    cleaned.fasta > netmhcpan_i_out.tsv"
    """
}

process NETMHCPANII {
    debug true
    publishDir '/home/pathinformatics/epitope_outputs/netmhcpan_ii_output'

    input:
    path protein_fasta

    output:
    path("netmhcpan_ii_out.tsv"), emit: netnhcpan_ii_tsv

    script:
    """
    /bin/bash -c "sed '/^[^>]/s/[B|X|J|U]//g' $protein_fasta > cleaned.fasta \
    && /mhc_ii/mhc_II_binding.py netmhciipan_ba \
    HLA-DRB1*03:01,HLA-DRB1*07:01,HLA-DRB1*15:01,HLA-DRB3*01:01,HLA-DRB3*02:02,HLA-DRB4*01:01,HLA-DRB5*01:01 \
    cleaned.fasta > netmhcpan_ii_out.tsv"
    """
}

workflow {
    protein_fasta_ch = Channel.fromPath(params.protein_file)
    split_fastas_ch = Channel.fromPath('/home/pathinformatics/epitope_outputs/split_fastas/*')

    cdhit_out_ch = CDHIT(protein_fasta_ch, params.cdhit_similarity_threshold)
    CDHITTOTSV(cdhit_out_ch.clstr_file, params.cdhit_similarity_threshold)
    /* B-CELL SCORING */
    if (params.bepipred == "yes") {
        bepipred_out_ch = BEPIPRED(protein_fasta_ch)
        BEPIPREDTOTSV(bepipred_out_ch.bepipred_output)
    }
    if (params.epidope == "yes") {
        EPIDOPE(protein_fasta_ch)
    }
    /* T-CELL SCORING */
    if (params.netmhcpani == "yes") {
        NETMHCPANI(protein_fasta_ch)
    }
    if (params.netmhcpanii == "yes") {
        NETMHCPANII(protein_fasta_ch)
    }

}
