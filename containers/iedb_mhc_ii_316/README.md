# README

## Inputs

fasta file or piped-in sequence

## Outputs

writes tabular-looking data to stdout:

```
allele	seq_num	start	end	length	peptide	consensus_percentile_rank	adjusted_consensus_percentile_rank	comblib_core	comblib_score	comblib_rank	adjusted_comblib_rank	smm_align_core	smm_align_ic50	smm_align_rank	adjusted_smm_align_rank	nn_align_core	nn_align_ic50	nn_align_rank	adjusted_nn_align_rank	sturniolo_core	sturniolo_score	sturniolo_rank	adjusted_sturniolo_rank
HLA-DRB1*03:01	1	8	22	15	EGVSGATWVDLVLEG	59.0	59.00	-	-	-	-	VSGATWVDL	5932.0	59.0	59.00	TWVDLVLEG	3356.9	50.0	50.00	VSGATWVDL	-1.34	88.0	88.00
HLA-DRB1*03:01	1	4	18	15	RDFLEGVSGATWVDL	85.0	85.00	-	-	-	-	FLEGVSGAT	8523.0	69.0	69.00	VSGATWVDL	12707.7	85.0	85.00	FLEGVSGAT	-1.2	87.0	87.00
HLA-DRB1*03:01	1	5	19	15	DFLEGVSGATWVDLV	85.0	85.00	-	-	-	-	FLEGVSGAT	18841.0	85.0	85.00	VSGATWVDL	10634.4	81.0	81.00	FLEGVSGAT	-1.2	87.0	87.00
HLA-DRB1*03:01	1	6	20	15	FLEGVSGATWVDLVL	85.0	85.00	-	-	-	-	FLEGVSGAT	18363.0	85.0	85.00	VSGATWVDL	8154.8	74.0	74.00	FLEGVSGAT	-1.2	87.0	87.00

```

## Command

### Build Container

Note that we just use the epidope:v0.3 container:

`sudo docker build -t iedb-mhc-ii .`

### Run Container

Note that you should pick a host directory to mount as well as a host path in which to save the output:

```
sudo docker run --rm \
 -v <HOST_FOLDER_CONTAINING_MULTIFASTA_INPUT_FILE>:/iedb_io iedb-mhc-ii:latest /bin/bash -c \
 "python /mhc_ii/mhc_II_binding.py <ALGORITHM> <HLA_ALLELE> <SEQUENCE_FASTA_PATH> > <OUTFILE_PATH>"
```

For example:

```
sudo docker run --rm \
 -v $PWD/iedb_mhc_ii_io:/iedb_io iedb-mhc-ii:latest /bin/bash -c \
   "python /mhc_ii/mhc_II_binding.py netmhciipan_ba HLA-DRB1*03:01 test.fasta >  /iedb_io/out.txt"
```
