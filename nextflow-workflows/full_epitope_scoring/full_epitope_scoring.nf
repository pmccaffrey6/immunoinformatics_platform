#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*params.protein_file = '/home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/test.fasta'*/
params.protein_file = '/home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/alphavirus_protein_multiseq.fasta'
params.cdhit_similarity_threshold = 0.95
params.alphafold_pdb_folder = '/home/pathinformatics/epitope_outputs/alphafold_predictions/*/*.pdb'

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

process DISCOTOPE {
    debug true
    publishDir '/home/pathinformatics/epitope_outputs/discotope_output'

    input:
    path pdb_file

    output:
    path("discotope_out.txt"), emit: discotope_out

    script:
    """
    /discotope-1.1/discotope -f $pdb_file -chain A > discotope_out.txt
    """
}

process NETMHCPANI {
    debug true
    publishDir '/home/pathinformatics/epitope_outputs/netmhcpan_i_output'

    input:
    path protein_fasta
    val allele

    output:
    path("netmhcpan_i_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_out.tsv"), emit: netmhcpan_i_tsv

    script:
    """
    /bin/bash -c "sed '/^[^>]/s/[B|X|J|U]//g' $protein_fasta > ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta \
    && /mhc_i/src/predict_binding.py netmhcpan_ba \
    $allele 9 \
    ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta > \
    netmhcpan_i_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_out.tsv"
    """
}

process NETMHCPANII {
    debug true
    publishDir '/home/pathinformatics/epitope_outputs/netmhcpan_ii_output'

    input:
    path protein_fasta
    val allele

    output:
    path("netmhcpan_ii_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_out.tsv"), emit: netmhcpan_ii_tsv

    script:
    """
    /bin/bash -c "sed '/^[^>]/s/[B|X|J|U]//g' $protein_fasta > ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta \
    && /mhc_ii/mhc_II_binding.py netmhciipan_ba \
    $allele \
    ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta > \
    netmhcpan_ii_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_out.tsv"
    """
}


workflow {
  protein_fasta_ch = Channel.fromPath(params.protein_file)
  protein_fasta_value_ch = file(params.protein_file)
  alphafold_pdb_ch = Channel.fromPath(params.alphafold_pdb_folder)

  mhc_i_alleles_ch = Channel.from('HLA-A*01:01','HLA-A*02:01','HLA-A*03:01','HLA-A*24:02','HLA-A*26:01','HLA-B*07:02','HLA-B*08:01','HLA-B*15:01>
  mhc_ii_alleles_ch = Channel.from('HLA-DRB1*03:01','HLA-DRB1*07:01','HLA-DRB1*15:01','HLA-DRB3*01:01','HLA-DRB3*02:02','HLA-DRB4*01:01','HLA-DR>

  /*split_fastas_ch = SPLITFASTAS(protein_fasta_ch)*/
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
  if (params.dc_bcell == "yes") {
      discotope_out_ch = DISCOTOPE(alphafold_pdb_ch)
  }
  /* T-CELL SCORING */
  if (params.netmhcpani == "yes") {
      NETMHCPANI(protein_fasta_value_ch, mhc_i_alleles_ch)
  }
  if (params.netmhcpanii == "yes") {
      NETMHCPANII(protein_fasta_value_ch, mhc_ii_alleles_ch)
  }
}
