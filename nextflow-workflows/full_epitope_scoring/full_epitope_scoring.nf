#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.immunoinformatics_allele_table_path = '/home/pathinformatics/immunoinformatics_platform/host/host-data/allele-frequency-processed-tables/all_allele_frequencies.tsv'
params.tcrpmhc_template_file = '/home/pathinformatics/immunoinformatics_platform/host/host-data/tcrpmhc_templates/TCRpMHC_template_0.fasta'

/*params.protein_file = '/home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/test.fasta'*/
params.protein_file = '/home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/alphavirus_protein_multiseq.fasta'
params.epitope_output_folder = '/home/pathinformatics/epitope_outputs'

params.cdhit_similarity_threshold = 0.95
params.alphafold_pdb_folder = '/home/pathinformatics/epitope_outputs/alphafold_predictions/*/*.pdb'

params.netmhcpan_chunk_size = 500

params.cdhit_input_proteins = "no"
params.bepipred = "no"
params.epidope = "no"
params.netmhcpani = "no"
params.netmhcpanii = "no"
params.dc_bcell = "no"
params.consolidate_epitopes = "no"
params.tcrpmhc = "no"
params.jessev = "no"
params.jessev_top_n = 3
params.include_docking_in_immunogenicity = "no"

params.allele_target_region = "Brazil"

params.netmhcpan_tcr_toprank_threshold = 0.05

params.b_cell_antigen_templates="Q8QZ72,Q8JUX5"

process PROCESSINPUTFASTA {
    debug true

    publishDir "${params.epitope_output_folder}/processed_fastas"

    input:
    path protein_fasta

    output:
    path("${protein_fasta.getBaseName()}_no_descriptors.fasta"), emit: fasta_wo_descriptions

    script:
    """
    #!/usr/bin/python3
    import os
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    sequences = []
    clean_protein_ids = []
    protein_descriptions = []

    for record in SeqIO.parse("${protein_fasta}", "fasta"):
        sequences.append(record.seq)
        clean_protein_ids.append(record.id)
        protein_descriptions.append(record.description)

    clean_records = [SeqRecord(seq=seq, id=id, description="") for seq,id in zip(sequences, clean_protein_ids)]

    with open("${protein_fasta.getBaseName()}_no_descriptors.fasta", "w") as output_handle:
        SeqIO.write(clean_records, output_handle, "fasta")


    if not os.path.exists("${params.epitope_output_folder}/join_tables"):
        os.makedirs("${params.epitope_output_folder}/join_tables")

    input_fasta_table = pd.DataFrame(data={
        "sequence":sequences,
        "protein_id":clean_protein_ids,
        "protein_description":protein_descriptions})

    input_fasta_table.to_csv("${params.epitope_output_folder}/join_tables/input_fasta_table.tsv", sep='\t')
    """

}

process FORMATALLELEFREQUENCIES {
    debug true

    publishDir "${params.epitope_output_folder}/join_tables"

    input:
    val regions

    output:
    path("${regions.replaceAll(/,/,"_")}_regional_allele_frequencies.tsv"), emit: regional_allele_frequencies_tsv

    script:
    """
    #!/usr/bin/python3
    import pandas as pd

    population_terms = "${regions}".split(",")

    allele_df = pd.read_csv("${params.immunoinformatics_allele_table_path}", sep='\t')
    regional_allele_data = pd.concat(allele_df[allele_df["dataset_name"].str.lower().str.contains(popterm.lower())] for popterm in population_terms).drop_duplicates()

    dfs = []

    for allele_type in regional_allele_data[pd.notna(regional_allele_data["allele_type"])]['allele_type'].unique():
        allele_data = regional_allele_data[regional_allele_data["allele_type"]==allele_type]

        allele_data = allele_data[(allele_data["locus"].str.contains(":")) & (pd.notna(allele_data["locus"]))]
        allele_data["locus"] = allele_data["locus"].apply(lambda x: ':'.join(x.split(':')[:2]))

        allele_data["count"] = 1
        allele_denom = len(allele_data)
        allele_grp = allele_data[["locus","count"]].groupby("locus").sum().reset_index()
        allele_grp["prevalence"] = (allele_grp["count"] / allele_denom)

        if allele_type in ["A","B","C"]:
            allele_grp["locus_join"] = ["HLA-"+i.replace("*","") for i in allele_grp["locus"]]
        else:
            allele_grp["locus_join"] = [i.replace("*","_").replace(":","") for i in allele_grp["locus"]]

        dfs.append(allele_grp)

    allele_prevalence = pd.concat(dfs)
    allele_prevalence["region"] = "_".join("${regions}".split(","))

    allele_prevalence.to_csv("${regions.replaceAll(/,/,"_")}_regional_allele_frequencies.tsv", sep='\t')
    """

}

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
    -c $similarity_threshold \
    -l 5 \
    -T 12 \
    -g 1 \
    -aS 0.9 \
    -uS 0.1 \
    -d 200
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

process NETMHCPANIIEDB {
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

process NETMHCPANI {
    debug true
    publishDir "${params.epitope_output_folder}/netmhcpan_i_output"

    input:
    path protein_fasta
    val allele

    output:
    path("netmhcpan_i_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_out.xls"), emit: netmhcpan_i_xls

    script:
    """
    /bin/bash -c "sed '/^[^>]/s/[B|X|J|U]//g' $protein_fasta > ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta \
    && export NETMHCpan=/netMHCpan-4.1/Linux_x86_64/ \
    && export TMPDIR=$HOME \
    && /netMHCpan-4.1/Linux_x86_64/bin/netMHCpan \
    -BA \
    -xls \
    -a $allele \
    -l 9 \
    -t 15 \
    -v \
    -f ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta \
    -xlsfile netmhcpan_i_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_out.xls"
    """
}

process NETMHCPANIIIEDB {
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

process NETMHCPANII {
    debug true
    publishDir "${params.epitope_output_folder}/netmhcpan_ii_output"

    input:
    each protein_fasta
    /*path protein_fasta*/
    val allele

    output:
    path("netmhcpan_ii_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_${protein_fasta.getName()}_out.xls"), emit: netmhcpan_ii_xls

    script:
    """
    /bin/bash -c "sed '/^[^>]/s/[B|X|J|U]//g' $protein_fasta > ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta \
    && export NETMHCIIpan=/netMHCIIpan-4.1 \
    && export TMPDIR=$HOME \
    && /netMHCIIpan-4.1/netMHCIIpan \
    -BA \
    -xls \
    -a $allele \
    -length 15 \
    -filter \
    -rankF 10 \
    -v \
    -f ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta \
    -xlsfile netmhcpan_ii_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_${protein_fasta.getName()}_out.xls"
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

    output:
    val 1

    script:
    """
    #!/usr/bin/python3
    import glob
    import os
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from multiprocessing import Pool

    epitope_paths = glob.glob("${params.epitope_output_folder}/*")
    # Consolidate from Epidope
    print("Consolidating Epidope")
    epidope_output_files = glob.glob("${params.epitope_output_folder}/epidope_output/epidope_output*/predicted_epitopes.csv")
    epidope_df = pd.concat([pd.read_csv(file, sep='\t') for file in epidope_output_files])[["#Gene_ID","sequence","score"]].drop_duplicates()
    ##[["#Gene_ID","sequence"]]
    epidope_df.rename(columns={"#Gene_ID":"protein_id","score":"epidope_score"}, inplace=True)
    epidope_df["type"] = "epidope"
    print(f"Epidope Table Size: {epidope_df.shape}")

    # Consolidate from Bepipred
    print("Consolidating Bepipred")
    bepipred_output_files = glob.glob("${params.epitope_output_folder}/bepipred_output/*.tsv")
    bepipred_df = pd.concat([pd.read_csv(file, sep='\t') for file in bepipred_output_files])[["Peptipe","Length","protein_id"]].drop_duplicates()
    ##[["Peptide","protein_id"]]
    bepipred_df.rename(columns={"Peptipe":"sequence"}, inplace=True)
    bepipred_df["type"] = "bepipred"
    bepipred_df["bepipred_score"] = 1
    # Filter Bepipred predicted epitopes by target size of 5-22 AA
    bepipred_df = bepipred_df[(bepipred_df["Length"]>=5) & (bepipred_df["Length"]<=22)]
    bepipred_df.drop("Length", axis=1, inplace=True)
    print(f"Bepipred Table Size: {bepipred_df.shape}")

    # Consolidate from Discotope
    print("Consolidating Discotope")
    discotope_files = glob.glob("${params.epitope_output_folder}/discotope_output/*.tsv")
    discotope_df = pd.concat([pd.read_csv(file, sep='\t') for file in discotope_files])[["mean_score","sequence","protein_name"]].drop_duplicates()
    discotope_df["protein_name"] = [i.split("_ranked")[0] for i in discotope_df["protein_name"].values]
    ##[["sequence","protein_name"]]
    discotope_df.rename(columns={"protein_name":"protein_id"}, inplace=True)
    discotope_df["type"] = "discotope"
    print(f"Discotope Table Size: {discotope_df.shape}")

    # Consolidate from NetMHCIPan
    TOPRANK = $params.netmhcpan_tcr_toprank_threshold
    toprank = str(TOPRANK).replace(".","_")
    print("Consolidating NetMHCPAN I")
    netmhcpan_i_files = glob.glob("${params.epitope_output_folder}/netmhcpan_i_output/*.xls")
    netmhcpan_i_dfs = []
    for file in netmhcpan_i_files:
        df = pd.read_csv(file, sep='\t', skiprows=1)
        df_for_col = pd.read_csv(file)
        df["allele"] = df_for_col.columns[-1].strip('\t')
        netmhcpan_i_dfs.append(df)
    netmhcpan_i_df = pd.concat(netmhcpan_i_dfs)
    netmhcpan_i_df.rename(columns={"ID":"protein_id"}, inplace=True)
    netmhcpan_i_df = netmhcpan_i_df[["Peptide","protein_id","allele","EL-score","EL_Rank","BA-score","BA_Rank","Ave","NB"]].drop_duplicates()
    ##[["peptide","proteinid"]]
    netmhcpan_i_df.rename(columns={"Peptide":"sequence","EL-score":"netmhcpan_i_el_score","EL_Rank":"netmhcpan_i_el_rank","BA-score":"netmhcpan_i_ba_score","BA_Rank":"netmhcpan_i_ba_rank","Ave":"netmhcpan_i_ave","NB":"netmhcpan_i_nb"}, inplace=True)
    netmhcpan_i_df[["sequence", "protein_id"]].drop_duplicates().to_csv("${params.epitope_output_folder}/join_tables/netmhcpan_i_peptides_protein_ids.tsv", sep='\t')
    netmhcpan_i_df["type"] = "netmhcpan_i"
    netmhcpan_i_df.reset_index().to_feather("${params.epitope_output_folder}/netmhcpan_i_output/all_netmhcpan_i.feather")
    print(f"NetMHCPAN I Table Size: {netmhcpan_i_df.shape}")
    netmhcpan_i_toprank = netmhcpan_i_df[netmhcpan_i_df["netmhcpan_i_ba_rank"]<=TOPRANK]
    netmhcpan_i_toprank.to_csv(f"${params.epitope_output_folder}/netmhcpan_i_output/netmhcpan_i_top_{toprank}.tsv", sep='\t', index=False)
    netmhcpan_i_toprank_records = [SeqRecord(seq=Seq(s), id=s, description="") for s in list(netmhcpan_i_toprank["sequence"].unique())]
    with open(f"${params.epitope_output_folder}/netmhcpan_i_output/netmhcpan_i_top_{toprank}.fasta", "w") as netmhcpan_i_toprank_output_file:
        SeqIO.write(netmhcpan_i_toprank_records, netmhcpan_i_toprank_output_file, "fasta")


    # Consolidate from NetMHCIIPan
    TOPRANK = $params.netmhcpan_tcr_toprank_threshold
    toprank = str(TOPRANK).replace(".","_")
    print("Consolidating NetMHCPAN II")
    netmhcpan_ii_files = glob.glob("${params.epitope_output_folder}/netmhcpan_ii_output/*.xls")
    netmhcpan_ii_dfs = []
    for file in netmhcpan_ii_files:
        df = pd.read_csv(file, sep='\t', skiprows=1)
        df_for_col = pd.read_csv(file)
        df["allele"] = df_for_col.columns[-1].strip('\t')
        netmhcpan_ii_dfs.append(df)
    netmhcpan_ii_df = pd.concat(netmhcpan_ii_dfs)
    netmhcpan_ii_df.rename(columns={"ID":"protein_id"}, inplace=True)
    netmhcpan_ii_df = netmhcpan_ii_df[["Peptide","protein_id","allele","Score","Rank","Score_BA","Rank_BA","Ave","NB"]].drop_duplicates()
    ##[["peptide","protein_id"]]
    netmhcpan_ii_df.rename(columns={"Peptide":"sequence","Score":"netmhcpan_ii_el_score","Rank":"netmhcpan_ii_el_rank","Score_BA":"netmhcpan_ii_ba_score","Rank_BA":"netmhcpan_ii_ba_rank","Ave":"netmhcpan_ii_ave","NB":"netmhcpan_ii_nb"}, inplace=True)
    netmhcpan_ii_df[["sequence", "protein_id"]].drop_duplicates().to_csv("${params.epitope_output_folder}/join_tables/netmhcpan_ii_peptides_protein_ids.tsv", sep='\t')
    netmhcpan_ii_df["type"] = "netmhcpan_ii"
    netmhcpan_ii_df.reset_index().to_feather("${params.epitope_output_folder}/netmhcpan_ii_output/all_netmhcpan_ii.feather")
    print(f"NetMHCPAN II Table Size: {netmhcpan_ii_df.shape}")
    netmhcpan_ii_toprank = netmhcpan_ii_df[netmhcpan_ii_df["netmhcpan_ii_ba_rank"]<=TOPRANK]
    netmhcpan_ii_toprank.to_csv(f"${params.epitope_output_folder}/netmhcpan_ii_output/netmhcpan_ii_top_{toprank}.tsv", sep='\t', index=False)
    netmhcpan_ii_toprank_records = [SeqRecord(seq=Seq(s), id=s, description="") for s in list(netmhcpan_ii_toprank["sequence"].unique())]
    with open(f"${params.epitope_output_folder}/netmhcpan_ii_output/netmhcpan_i_top_{toprank}.fasta", "w") as netmhcpan_ii_toprank_output_file:
        SeqIO.write(netmhcpan_ii_toprank_records, netmhcpan_ii_toprank_output_file, "fasta")

    # Consolidate Outputs
    all_concat_df_outdir = "${params.epitope_output_folder}/consolidated_outputs"
    if not os.path.exists(all_concat_df_outdir):
        os.makedirs(all_concat_df_outdir)



    # Create consolidated FASTA
    dfs = [epidope_df, bepipred_df, discotope_df, netmhcpan_i_df, netmhcpan_ii_df]

    def create_fasta_and_txt_from_df(df):
        sequence_records = []
        df_type = df["type"].values[0]
        df_proc = df[(pd.notna(df["sequence"]))][["protein_id","sequence"]].drop_duplicates()
        txt_lines = list(set(df_proc["sequence"].unique()))

        pd.Series(df_proc["sequence"].unique()).apply(lambda x: sequence_records.append(SeqRecord(seq=Seq(str(x)), id=str(x), description=str(df_type)  )) )

        with open(f"${params.epitope_output_folder}/consolidated_outputs/consolidated_epitopes_{df_type}.fasta", 'w') as df_fasta_outfile:
            SeqIO.write(sequence_records, df_fasta_outfile, "fasta")

        with open(f"${params.epitope_output_folder}/consolidated_outputs/consolidated_epitopes_{df_type}.txt", 'w') as txt_outfile:
            txt_outfile.writelines("\\n".join(txt_lines))

    with Pool(len(dfs)) as p:
        p.map(create_fasta_and_txt_from_df, dfs)

    # Consolidate all dfs together
    for df in [epidope_df, bepipred_df, discotope_df, netmhcpan_i_df, netmhcpan_ii_df]:
        df_type = df["type"].values[0]
        df.reset_index().to_feather(f"${params.epitope_output_folder}/consolidated_outputs/consolidated_outputs_{df_type}.feather")
        df.to_csv(f"${params.epitope_output_folder}/consolidated_outputs/consolidated_outputs_{df_type}.tsv", sep='\t')
    """
}

process GATHEREPITOPEFASTAS {
    debug true

    input:
    val consolidated_output

    output:
    path("all_consolidated_epitopes.fasta"), emit: gathered_epitopes_fasta

    script:
    """
    cat ${params.epitope_output_folder}/consolidated_outputs/consolidated_epitopes*.fasta > \
    all_consolidated_epitopes.fasta
    """

}

process PREPAREDATAFORJESSEV {
    debug true

    publishDir "${params.epitope_output_folder}/jessev_input"

    input:
    path allele_frequencies_table

    output:
    path("jessev_input.csv"), emit: jessev_input_csv

    script:
    """
    #!/usr/bin/python3
    import pandas as pd
    import os

    MAX_NUM_PROTEIN_IDS = 10
    THRESHOLD_VALUE = 05.00
    THRESHOLD_ATTRIBUTE = "netmhcpan_i_ba_rank"

    if not os.path.exists("${params.epitope_output_folder}/jessev_input"):
        os.makedirs("${params.epitope_output_folder}/jessev_input")

    netmhcpan_i_df = pd.read_feather("${params.epitope_output_folder}/netmhcpan_i_output/all_netmhcpan_i.feather")

    allele_freq_df = pd.read_csv("${allele_frequencies_table}", sep='\t')
    netmhcpan_i_df_full = netmhcpan_i_df.merge(allele_freq_df, left_on="allele", right_on="locus_join", how="inner")

    if "${params.include_docking_in_immunogenicity}"=="yes":
        pdb_episa_df = pd.read_csv("${params.epitope_output_folder}/pdb_episa_output/pdb_episa_output_table.tsv", sep='\t').rename(columns={"SEQUENCE":"sequence"})
        netmhcpan_i_df_full = netmhcpan_i_df_full.merge(pdb_episa_df, on="sequence", how="inner")
        print(netmhcpan_i_df_full.head())
        print("shape of epitope table with docking:", netmhcpan_i_df_full.shape)

    netmhcpan_i_df_full["weighted_immunogen"] = netmhcpan_i_df_full["netmhcpan_i_ba_score"] * netmhcpan_i_df_full["prevalence"]
    jess_ev_df = netmhcpan_i_df_full[netmhcpan_i_df_full[THRESHOLD_ATTRIBUTE]<THRESHOLD_VALUE][["sequence", "protein_id", "allele", "weighted_immunogen"]].drop_duplicates().groupby("sequence").agg({"protein_id":";".join, "allele":";".join, "weighted_immunogen":"sum"})
    # Because some of these epitopes appear in MANY proteins, we have to do some deduplication.
    # first we only keep the UNIQUE alleles from the aggregation which makes sense
    # second, we have to set some cap on the number of protein_ids we can carry forward to JessEV
    # and that number is captured in MAX_NUM_PROTEIN_IDS
    jess_ev_df["allele"] = jess_ev_df["allele"].apply(lambda x: ';'.join(list(set(x.split(';')))))
    jess_ev_df["protein_id"] = jess_ev_df["protein_id"].apply(lambda x: ';'.join( sorted(list(set(x.split(';'))))[:MAX_NUM_PROTEIN_IDS]))
    print(f"shape of JessEV input table: {jess_ev_df.shape}")
    jess_ev_df.reset_index().rename(columns={"weighted_immunogen":"immunogen", "sequence":"epitope", "protein_id":"proteins", "allele":"alleles"}).to_csv("jessev_input.csv", index=False)
    """
}

process RUNJESSEV {
    debug true

    input:
    path jessev_input_csv
    val iterations

    script:
    """
    #!/usr/bin/python3
    import pandas as pd
    import subprocess
    import glob
    import os
    import shutil
    import docker
    from docker.types import Mount

    client = docker.from_env()

    NUM_EPITOPES = 3
    MIN_SPACER_LEN = 4

    shutil.copy("${params.epitope_output_folder}/jessev_input/$jessev_input_csv", "${params.epitope_output_folder}/jessev_input/jessev_input_0.csv")

    for iter in range($iterations):
        print("running JessEV iteration:", iter)

        if iter > 0:
            input_csv_pre = f"${params.epitope_output_folder}/jessev_input/jessev_input_{iter-1}.csv"
            df_input_csv_pre = pd.read_csv(input_csv_pre)
            output_csv_pre = f"${params.epitope_output_folder}/jessev_input/jessev_output_{iter-1}.csv"
            df_output_csv_pre = pd.read_csv(output_csv_pre)
            df_input_csv_pre["jessev_used"] = [(i in df_output_csv_pre["vaccine"].values[0]) for i in df_input_csv_pre["epitope"]]
            input_csv_filtered = df_input_csv_pre[df_input_csv_pre["jessev_used"]==False].drop("jessev_used", axis=1)
            input_csv_filtered.to_csv(f"${params.epitope_output_folder}/jessev_input/jessev_input_{iter}.csv", index=False)

        jessev_statement = f"${params.epitope_output_folder}/jessev_input/jessev_input_{iter}.csv ${params.epitope_output_folder}/jessev_input/jessev_output_{iter}.csv"

        client.containers.run("pmccaffrey6/jess_ev:latest",
            command=f"/opt/conda/envs/jessev/bin/python3 /JessEV/design.py -e 3 -s 4 {jessev_statement}",
            auto_remove=True,
            mounts=[Mount('/home/pathinformatics','//home/pathinformatics', type="bind")]
         )

    if not os.path.exists(f"${params.epitope_output_folder}/jessev_output/"):
        os.makedirs(f"${params.epitope_output_folder}/jessev_output/")

    pd.concat([pd.read_csv(i) for i in glob.glob(f"${params.epitope_output_folder}/jessev_input/jessev_output_*.csv")]).to_csv(
        f"${params.epitope_output_folder}/jessev_output/jessev_outputs_top_${iterations}.tsv", sep='\t', index=False)

    """
}

process TCRPMHC {
    debug true

    input:
    path input_fasta
    path template_fasta

    output:
    val "${params.epitope_output_folder}/tcrpmhc_output/${input_fasta}_TCR-pMHC.pdb", emit: tcrpmhc_output

    script:
    """
    echo $input_fasta && \
    echo $template_fasta && \

    cp $template_fasta /home/pathinformatics/epitope_outputs/tcrpmhc_output/template_$input_fasta && \
    cat $input_fasta >> /home/pathinformatics/epitope_outputs/tcrpmhc_output/template_$input_fasta && \
    cat /home/pathinformatics/epitope_outputs/tcrpmhc_output/template_$input_fasta && \

    /opt/conda/envs/TCRpMHCmodels/bin/tcrpmhc_models \
    /home/pathinformatics/epitope_outputs/tcrpmhc_output/template_$input_fasta \
    -n $input_fasta \
    -p $params.epitope_output_folder/tcrpmhc_output
    """
}

process PDBEPISA {
    debug true

    input:
    val tcrpmhc_batch

    output:
    val tcrpmhc_batch

    script:
    """
    mkdir "${file(tcrpmhc_batch.first()).getBaseName()}" && \
    cp ${tcrpmhc_batch.join(" ")} ${file(tcrpmhc_batch.first()).getBaseName()} && \
    /opt/run "${file(tcrpmhc_batch.first()).getBaseName()}" $params.epitope_output_folder/pdb_episa_output
    """
}

process PDBEPISATOTABLE {
    debug true

    input:
    val pdb_files

    script:
    """
    #!/usr/bin/python3
    import pandas as pd
    import os
    import glob
    import xml
    import xml.etree.ElementTree as ET
    from Bio.PDB.Polypeptide import three_to_one

    existing_xml_files = glob.glob(f"${params.epitope_output_folder}/pdb_episa_output/*TCR-pMHC.xml")

    print("received pdb files:")
    print(f"${pdb_files}")

    print("xml files:")
    print(existing_xml_files)

    def parse_pisa_xml(xml_file):
        tree = ET.parse(xml_file)
        root = tree.getroot()
        elems = {i.tag:i for i in list(root)}

        interface_summary = elems['INTERFACESUMMARY']
        structures = list(interface_summary)
        residues = list(elems['RESIDUES'])

        structure_data = {}

        for structure, residue in zip(structures, residues):
            base_data = {i.tag:i for i in list(structure)}
            base_tags = [i.tag for i in list(structure)]

            for i in base_tags:
                if i.startswith("NUMBEROFRED"):
                    numresidues = base_data[i]
                    for i in list(numresidues):
                        base_data[i.tag] = int(i.text)

            sequence = ''.join([ three_to_one(list(i)[0].text.split(':')[-1].strip(' ').split(' ')[0]) for i in list(residue)])
            base_data['SEQUENCE'] = sequence

            base_keys = base_data.keys()
            for base_item in base_keys:
                if base_item.startswith('SOLVENTAREA'):
                    solvent_area = {i.tag:i.text for i in list(base_data[base_item])}
                elif base_item.startswith('SOLVATIONENERGY'):
                    solvation_energy = {i.tag:i.text for i in list(base_data[base_item])}

            for item in solvent_area.keys():
                newkey = f"SOLVENTAREA_{item}"
                base_data[newkey] = float(solvent_area[item])

            for item in solvation_energy.keys():
                newkey = f"SOLVATIONERGY_{item}"
                base_data[newkey] = float(solvation_energy[item])

            structure_data[structure.tag] = {k:v for k,v in base_data.items() if not type(v)==xml.etree.ElementTree.Element}

        df_out = pd.DataFrame(structure_data.values())
        return df_out

    xml_df = pd.concat([parse_pisa_xml(xml_file) for xml_file in existing_xml_files])
    print("XML DF Shape:", xml_df.shape)
    if not os.path.exists(f"${params.epitope_output_folder}/pdb_episa_output"):
        os.makedirs(f"${params.epitope_output_folder}/pdb_episa_output")
    xml_df["SEQUENCE_LENGTH"] = [len(i) for i in xml_df["SEQUENCE"]]
    xml_df[xml_df["SEQUENCE_LENGTH"]<10].to_csv(f"${params.epitope_output_folder}/pdb_episa_output/pdb_episa_output_table.tsv", sep='\t')
    """
}

process GETUNIPROTBYACCESSION {
    debug true

    publishDir "${params.epitope_output_folder}/b_cell_antigen_templates"

    input:
    val uniprot_accession

    output:
    path("${uniprot_accession}.fasta"), emit: b_cell_antigen_fasta

    script:
    """
    #!/usr/bin/python3
    import requests
    import json
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    import os

    uniprot_out = json.loads(requests.get(f"https://www.ebi.ac.uk/proteins/api/proteins/${uniprot_accession}").text)
    fastaout = [SeqRecord(seq=Seq(str(uniprot_out['sequence']['sequence'])), id=uniprot_out['accession'], description="")]

    if not os.path.exists("${params.epitope_output_folder}/b_cell_antigen_templates"):
        os.makedirs("${params.epitope_output_folder}/b_cell_antigen_templates")

    with open(f"${uniprot_accession}.fasta", "w") as uniprot_fasta_out:
        SeqIO.write(fastaout, uniprot_fasta_out, "fasta")
    """
}

process BLASTFROMFILES {
    debug true

    publishDir "${params.epitope_output_folder}/blast_database"

    input:
    val fasta_files
    path query_file

    output:
    path("blast_results.tsv"), emit: blast_results_tsv

    script:
    """
    echo "fasta files ${fasta_files}"
    rm -f ${params.epitope_output_folder}/blast_database/uniprot_all.fasta*
    cat ${fasta_files.join(" ")} >> ${params.epitope_output_folder}/blast_database/uniprot_all.fasta

    makeblastdb -in ${params.epitope_output_folder}/blast_database/uniprot_all.fasta -input_type fasta -dbtype prot -title "b_cell_antigen_db"

    blastp -query ${query_file} \
    -db ${params.epitope_output_folder}/blast_database/uniprot_all.fasta \
    -outfmt 6 -evalue 0.0005 > blast_results.tsv
    """

}

process FILTERPROTEINSBYBLAST {
    debug true

    publishDir "${params.epitope_output_folder}/filtered_blast_results"

    input:
    path blast_results_tsv
    path input_proteins_file

    output:
    path "*.fasta"

    script:
    """
    #!/usr/bin/python3
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    colnames = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
    blast_df = pd.read_csv(f"${blast_results_tsv}", sep='\t', header=None, names=colnames)
    blast_df = blast_df[blast_df['pident']>80.0]
    template_accessions = list(blast_df['sseqid'].unique())
    for template_accession in template_accessions:
        template_blast = blast_df[blast_df['sseqid']==template_accession]
        matched_proteins = list(template_blast['qseqid'].unique())

        filtered_sequences = []

        for record in SeqIO.parse(f"${input_proteins_file}", "fasta"):
            if record.id in matched_proteins:
                filtered_sequences.append(record)

        print(f"number of blast hits for uniprot accession {template_accession}: {len(filtered_sequences)}")

        with open(f"blast_hits_with_{template_accession}.fasta", "w") as output_handle:
            SeqIO.write(filtered_sequences, output_handle, "fasta")
    """
}

process CLUSTALOMEGAMSA {
    debug true

    publishDir "${params.epitope_output_folder}/clustalomega_results"

    input:
    path filtered_fasta

    output:
    path "*.txt"

    script:
    """
    echo filtered_fasta $filtered_fasta && \
    clustalo -v -i $filtered_fasta --outfmt=clu -o ${filtered_fasta.getBaseName()}_aligned.txt
    """
}

workflow {
    protein_fasta_ch = Channel.fromPath(params.protein_file)
    protein_fasta_value_ch = file(params.protein_file)
    protein_fasta_clean_ch = PROCESSINPUTFASTA(protein_fasta_value_ch)
    protein_fasta_clean_ch.view()

    tcrpmhc_templates_value_ch = file(params.tcrpmhc_template_file)

    /*COLLECT B-CELL ANTIGEN TEMPLATES*/
    b_cell_antigen_fastas = GETUNIPROTBYACCESSION(Channel.from(params.b_cell_antigen_templates.split(",")))
    blast_results_tsv_ch = BLASTFROMFILES(b_cell_antigen_fastas.collect(), protein_fasta_value_ch)
    filtered_fasta_ch = FILTERPROTEINSBYBLAST(blast_results_tsv_ch, protein_fasta_value_ch).flatten()
    clustal_omega_output_ch = CLUSTALOMEGAMSA(filtered_fasta_ch)
    clustal_omega_output_ch.collect().view()

    /* CALCULATE ALLELE FREQUENCY TABLES */
    allele_frequencies_table_ch = FORMATALLELEFREQUENCIES(params.allele_target_region)
    allele_frequencies_table_ch.view()

    /*protein_fasta_splits_value_ch = Channel.fromPath(protein_fasta_value_ch).splitFasta(by: params.netmhcpan_chunk_size, file: true).collect()*/
    protein_fasta_splits_value_ch = protein_fasta_clean_ch.splitFasta(by: params.netmhcpan_chunk_size, file: true).collect()

    alphafold_pdb_ch = Channel.fromPath(params.alphafold_pdb_folder)

    mhc_i_alleles_ch = Channel.from('HLA-A01:01','HLA-A02:01','HLA-A03:01','HLA-A24:02','HLA-A26:01','HLA-B07:02','HLA-B08:01','HLA-B15:01','HLA-B27:05','HLA-B39:01','HLA-B40:01','HLA-B58:01')
    mhc_ii_alleles_ch = Channel.from('DRB1_0301','DRB1_0701','DRB1_1501','DRB3_0101','DRB3_0202','DRB4_0101','DRB5_0101')

    /*split_fastas_ch = SPLITFASTAS(protein_fasta_ch)*/
    split_fastas_ch = Channel.fromPath('/home/pathinformatics/epitope_outputs/split_fastas/*')

    /* CDHIT INPUT PROTEINS */
    if (params.cdhit_input_proteins == "yes") {
        cdhit_out_ch = CDHIT(protein_fasta_ch, params.cdhit_similarity_threshold)
        CDHITTOTSV(cdhit_out_ch.clstr_file, params.cdhit_similarity_threshold, "input_proteins")
    }
    /* B-CELL SCORING */
    if (params.bepipred == "yes") {
        bepipred_out_ch = BEPIPRED(protein_fasta_clean_ch)
        BEPIPREDTOTSV(bepipred_out_ch.bepipred_output)
    }
    if (params.epidope == "yes") {
        EPIDOPE(protein_fasta_clean_ch.splitFasta(by: 500, file: true))
    }
    if (params.dc_bcell == "yes") {
        MIGRATEAF2()
        discotope_out_ch = DISCOTOPE(alphafold_pdb_ch)
        DISCOTOPETOTSV(discotope_out_ch)
    }
    /* T-CELL SCORING */
    if (params.netmhcpani == "yes") {
        NETMHCPANI(protein_fasta_clean_ch, mhc_i_alleles_ch)
    }
    if (params.netmhcpanii == "yes") {
        NETMHCPANII(protein_fasta_splits_value_ch, mhc_ii_alleles_ch)
    }

    if (params.consolidate_epitopes == "yes") {
        consolidated_epitopes_output_ch = CONSOLIDATEEPITOPES()
        consolidated_epitopes_fasta_ch = GATHEREPITOPEFASTAS(consolidated_epitopes_output_ch)

        cdhit_epitopes_out_ch = CDHIT(consolidated_epitopes_fasta_ch, params.cdhit_similarity_threshold)
        CDHITTOTSV(cdhit_epitopes_out_ch.clstr_file, params.cdhit_similarity_threshold, "epitopes")
    }
    if (params.tcrpmhc == "yes") {
        protein_records_ch = Channel.fromPath("${params.epitope_output_folder}/netmhcpan_i_output/netmhcpan_i_top*.fasta").splitFasta(by: 1, file:true).take(100)
        tcrpmhc_pdbs = TCRPMHC(protein_records_ch, tcrpmhc_templates_value_ch)
        tcrpmhc_batches = tcrpmhc_pdbs.collate(3)
        tcrpmhc_batches.view()
        pdb_episa_output_ch = PDBEPISA(tcrpmhc_batches)
        PDBEPISATOTABLE(pdb_episa_output_ch.collect())
    }
    if (params.jessev == "yes") {
        jessev_input_file_ch = PREPAREDATAFORJESSEV(allele_frequencies_table_ch)
        jessev_input_file_ch.view()
        RUNJESSEV(jessev_input_file_ch, params.jessev_top_n)
    }
}
