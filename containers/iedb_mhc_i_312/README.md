# README

## Inputs

fasta file or piped-in sequence

## Outputs

writes tabular-looking data to stdout:

```
allele	seq_num	start	end	length	peptide	core	icore	ic50	rank
HLA-A*02:01	1	6	14	9	TMFEALPHI	TMFEALPHI	TMFEALPHI	5.77	0.04
HLA-A*02:01	1	440	448	9	LMFSTSAYL	LMFSTSAYL	LMFSTSAYL	6.41	0.05
HLA-A*02:01	1	447	455	9	YLVSIFLHL	YLVSIFLHL	YLVSIFLHL	7.59	0.06
HLA-A*02:01	1	137	145	9	TLMSIVSSL	TLMSIVSSL	TLMSIVSSL	10.9	0.1
HLA-A*02:01	1	10	18	9	ALPHIIDEV	ALPHIIDEV	ALPHIIDEV	11.76	0.11
HLA-A*02:01	1	435	443	9	ALMDLLMFS	ALMDLLMFS	ALMDLLMFS	13.04	0.12
HLA-A*02:01	1	45	53	9	ALISFLLLA	ALISFLLLA	ALISFLLLA	24.44	0.21
HLA-A*02:01	1	42	50	9	GIFALISFL	GIFALISFL	GIFALISFL	47.65	0.41
```

## Command

### Build Container

Note that we just use the epidope:v0.3 container:

`sudo docker build -t iedb-mhc-i .`

### Run Container

Note that you should pick a host directory to mount as well as a host path in which to save the output:

```
sudo docker run --rm \
 -v <HOST_FOLDER_CONTAINING_MULTIFASTA_INPUT_FILE>:/iedb_io iedb-mhc-i:latest /bin/bash -c \
 "<ALGORITHM> <HLA_ALLELE> <EPITOPE_LENGTH> <SEQUENCE_FASTA_PATH> > <OUTFILE_PATH>"
```

For example:

```
sudo docker run --rm \
 -v /home/jupyter-pathinformatics/nikos_viral_vaccine:/iedb_io iedb-mhc-i:latest /bin/bash -c \
   "./src/predict_binding.py netmhcpan_ba HLA-A*02:01 9 /mhc_i/examples/input_sequence.fasta >  /iedb_io/out.txt"
```
