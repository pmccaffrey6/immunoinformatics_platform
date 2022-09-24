#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*params.protein_file = '/home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/test.fasta'*/
params.protein_file = '/home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/alphavirus_protein_multiseq.fasta'
params.epitope_output_folder = '/home/pathinformatics/epitope_outputs'

params.cdhit_similarity_threshold = 0.95
params.alphafold_pdb_folder = '/home/pathinformatics/epitope_outputs/alphafold_predictions/*/*.pdb'

params.cdhit = "no"
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
    publishDir "${params.epitope_output_folder}/cdhit_output"

    input:
    path clstr_file
    val similarity_threshold
    val clustering_type

    output:
    path("cdhit_out_${clustering_type}_sim_${similarity_threshold}.tsv"), emit: clstr_tsv

    script:
    """
    python3 /cdhit_clstr_to_df.py \
    $clstr_file \
    cdhit_out_${clustering_type}_sim_${similarity_threshold}.tsv
    """

}

process EPIDOPE {
    debug true
    label 'EPIDOPE'
    publishDir "${params.epitope_output_folder}/epidope_output"

    input:
    path protein_fasta

    output:
    path("epidope_output_${protein_fasta}"), emit: epidope_output

    script:
    """
    epidope \
    -p 12 \
    -i $protein_fasta \
    -o epidope_output_$protein_fasta
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
    publishDir "${params.epitope_output_folder}/bepipred_output"

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

process MIGRATEAF2 {
    debug true

    """
    #!/usr/bin/python3
    import glob
    import os
    import shutil

    alphafold_files = glob.glob('/tmp/alphafold/*/ranked_0.pdb')
    print("Migrating AF2 Predictions:")
    print(",".join(alphafold_files))
    for alphafold_file in alphafold_files:
        protname = alphafold_file.split('/')[-2]
        outfile = f"{protname}_ranked_0.pdb"
        alphafold_dest_folder = os.path.join('/home','pathinformatics','epitope_outputs','alphafold_predictions',protname)
        alphafold_dest_file = os.path.join(alphafold_dest_folder,outfile)
        if not os.path.exists(alphafold_dest_folder):
            os.makedirs(alphafold_dest_folder)
        shutil.copy(alphafold_file, alphafold_dest_file)
    """
}

process DISCOTOPE {
    debug true

    input:
    path pdb_file

    output:
    path("${pdb_file}_discotope_out.txt"), emit: discotope_out

    script:
    """
    /discotope-1.1/discotope -f $pdb_file -chain A > ${pdb_file}_discotope_out.txt
    """
}

process DISCOTOPETOTSV {
    debug true
    publishDir "${params.epitope_output_folder}/discotope_output"

    input:
    path discotope_file

    output:
    path("${discotope_file}.tsv"), emit: discotope_out_tsv

    script:
    """
    python3 /discotope_to_tsv.py \
    $discotope_file ${discotope_file}.tsv $discotope_file
    """
}

process NETMHCPANI {
    debug true
    publishDir "${params.epitope_output_folder}/netmhcpan_i_output"

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
    publishDir "${params.epitope_output_folder}/netmhcpan_ii_output"

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

process SPLITFASTAS {
    debug true
    publishDir "${params.epitope_output_folder}/split_fastas"

    input:
    path protein_fasta

    output:
    path ('split_fastas'), emit: split_fastas

    script:
    """
    python3 /split_fastas.py \
    $protein_fasta \
    /home/pathinformatics/epitope_outputs/split_fastas
    """
}

process CONSOLIDATEEPITOPES {
    debug true
    /*publisDir "${params.epitope_output_folder}/consolidated_epitopes"*/

    script:
    """
    #!/usr/bin/python3
    import glob
    import os
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    epitope_paths = glob.glob("${params.epitope_output_folder}/*")
    # Consolidate from Epidope
    epidope_output_files = glob.glob("${params.epitope_output_folder}/epidope_output/epidope_output*/predicted_epitopes.csv")
    epidope_df = pd.concat([pd.read_csv(file, sep='\t') for file in epidope_output_files])[["#Gene_ID","sequence","score"]].drop_duplicates()
    ##[["#Gene_ID","sequence"]]
    epidope_df.rename(columns={"#Gene_ID":"protein_id","score":"epidope_score"}, inplace=True)
    epidope_df["type"] = "epidope"

    # Consolidate from Bepipred
    bepipred_output_files = glob.glob("${params.epitope_output_folder}/bepipred_output/*.tsv")
    bepipred_df = pd.concat([pd.read_csv(file, sep='\t') for file in bepipred_output_files])[["Peptipe","protein_id"]].drop_duplicates()
    ##[["Peptide","protein_id"]]
    bepipred_df.rename(columns={"Peptipe":"sequence"}, inplace=True)
    bepipred_df["type"] = "bepipred"
    bepipred_df["bepipred_score"] = 1

    # Consolidate from Discotope
    discotope_files = glob.glob("${params.epitope_output_folder}/discotope_output/*.tsv")
    discotope_df = pd.concat([pd.read_csv(file, sep='\t') for file in discotope_files])[["mean_score","sequence","protein_name"]].drop_duplicates()
    discotope_df["protein_name"] = [i.split("_ranked")[0] for i in discotope_df["protein_name"].values]
    ##[["sequence","protein_name"]]
    discotope_df.rename(columns={"protein_name":"protein_id"}, inplace=True)
    discotope_df["type"] = "discotope"

    # Consolidate from NetMHCIPan
    netmhcpan_i_files = glob.glob("${params.epitope_output_folder}/netmhcpan_i_output/*.tsv")
    netmhcpan_i_df = pd.concat([pd.read_csv(file, sep='\t') for file in netmhcpan_i_files])
    netmhcpan_i_df["protein_id"] = "prot1234"
    netmhcpan_i_df = netmhcpan_i_df[["peptide","protein_id","allele","ic50","rank"]].drop_duplicates()
    ##[["peptide","proteinid"]]
    netmhcpan_i_df.rename(columns={"peptide":"sequence","ic50":"netmhcpan_i_ic50","rank":"netmhcpan_i_rank"}, inplace=True)
    netmhcpan_i_df["type"] = "netmhcpan_i"

    # Consolidate from NetMHCIIPan
    netmhcpan_ii_files = glob.glob("${params.epitope_output_folder}/netmhcpan_ii_output/*.tsv")
    netmhcpan_ii_df = pd.concat([pd.read_csv(file, sep='\t') for file in netmhcpan_ii_files])
    netmhcpan_ii_df["protein_id"] = "prot1234"
    netmhcpan_ii_df = netmhcpan_ii_df[["peptide","protein_id","allele","ic50","percentile_rank"]].drop_duplicates()
    ##[["peptide","protein_id"]]
    netmhcpan_ii_df.rename(columns={"peptide":"sequence","ic50":"netmhcpan_ii_ic50","percentile_rank":"netmhcpan_ii_rank"}, inplace=True)
    netmhcpan_ii_df["type"] = "netmhcpan_ii"

    # Create consolidated FASTA
    dfs = [epidope_df, bepipred_df, discotope_df, netmhcpan_i_df, netmhcpan_ii_df]
    sequence_records = []
    for df in dfs:
        for idx, row in df[["protein_id","sequence"]].drop_duplicates().iterrows():
            sequence_records.append(SeqRecord(seq=Seq(row["sequence"]), id=row["protein_id"]))

    print(sequence_records)
    all_concat_df_outdir = "${params.epitope_output_folder}/consolidated_outputs"
    if not os.path.exists(all_concat_df_outdir):
        os.makedirs(all_concat_df_outdir)
    with open(f"${params.epitope_output_folder}/consolidated_outputs/consolidated_epitopes.fasta", 'w') as df_fasta_outfile:
        SeqIO.write(sequence_records, df_fasta_outfile, "fasta")

    # Consolidate all dfs together
    all_concat_df = pd.concat([epidope_df, bepipred_df, discotope_df, netmhcpan_i_df, netmhcpan_ii_df])
    print(f"concatenated DF created of shape {all_concat_df.shape}")
    all_concat_df_outdir = "${params.epitope_output_folder}/consolidated_outputs"
    if not os.path.exists(all_concat_df_outdir):
        os.makedirs(all_concat_df_outdir)
    all_concat_df.to_csv("${params.epitope_output_folder}/consolidated_outputs/consolidated_outputs.tsv", sep='\t')
    """
}

workflow {
    protein_fasta_ch = Channel.fromPath(params.protein_file)
    protein_fasta_value_ch = file(params.protein_file)
    alphafold_pdb_ch = Channel.fromPath(params.alphafold_pdb_folder)

    mhc_i_alleles_ch = Channel.from('HLA-A*01:01','HLA-A*02:01','HLA-A*03:01','HLA-A*24:02','HLA-A*26:01','HLA-B*07:02','HLA-B*08:01','HLA-B*15:01','HLA-B*27:05','HLA-B*39:01','HLA-B*40:01','HLA-B*58:01')
    mhc_ii_alleles_ch = Channel.from('HLA-DRB1*03:01','HLA-DRB1*07:01','HLA-DRB1*15:01','HLA-DRB3*01:01','HLA-DRB3*02:02','HLA-DRB4*01:01','HLA-DRB5*01:01')

    /*split_fastas_ch = SPLITFASTAS(protein_fasta_ch)*/
    split_fastas_ch = Channel.fromPath('/home/pathinformatics/epitope_outputs/split_fastas/*')

    /* CDHIT INPUT PROTEINS */
    if (params.cdhit == "yes") {
        cdhit_out_ch = CDHIT(protein_fasta_ch, params.cdhit_similarity_threshold)
        CDHITTOTSV(cdhit_out_ch.clstr_file, params.cdhit_similarity_threshold, "input_proteins")
    }
    /* B-CELL SCORING */
    if (params.bepipred == "yes") {
        bepipred_out_ch = BEPIPRED(protein_fasta_ch)
        BEPIPREDTOTSV(bepipred_out_ch.bepipred_output)
    }
    /*if (params.epidope == "yes") {
        EPIDOPE(protein_fasta_ch)
    }*/
    if (params.epidope == "yes") {
        EPIDOPE(protein_fasta_ch.splitFasta(by: 500, file: true))
    }
    if (params.dc_bcell == "yes") {
        MIGRATEAF2()
        discotope_out_ch = DISCOTOPE(alphafold_pdb_ch)
        DISCOTOPETOTSV(discotope_out_ch)
    }
    /* T-CELL SCORING */
    if (params.netmhcpani == "yes") {
        NETMHCPANI(protein_fasta_value_ch, mhc_i_alleles_ch)
    }
    if (params.netmhcpanii == "yes") {
        NETMHCPANII(protein_fasta_value_ch, mhc_ii_alleles_ch)
    }

    CONSOLIDATEEPITOPES()
}
