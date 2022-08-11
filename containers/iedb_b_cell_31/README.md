# README

## Inputs

fasta file (or SWISSPROT ID)

## Outputs

writes tabular-looking data to stdout:

```
Predicted peptides
No	Start	End	Peptipe	Length
1	20	25	ADVAGH	6
2	38	42	PETLE	5
3	55	63	EMKASEDLK	9
4	81	85	GHHEA	5
5	88	88	K	1
6	90	90	L	1
7	92	97	QSHATK	6
8	119	131	RHPGDFGADAQGA	13
9	150	153	LGYQ	4
Position	Residue	Score	Assignment
1	M	-0.911	.
2	V	-0.439	.
3	L	-0.055	.
4	S	-0.258	.
```

## Command

### Build Container

Note that we just use the epidope:v0.3 container:

`sudo docker build -t iedb-b-cell .`

### Run Container

Note that you should pick a host directory to mount as well as a host path in which to save the output:

```
sudo docker run --rm \
 -v <HOST_FOLDER_CONTAINING_MULTIFASTA_INPUT_FILE>:/iedb_io iedb-b-cell:latest /bin/bash -c \
 "python /bcell_standalone/predict_antibody_epitope.py -m <B_CELL_METHOD_TO_USE> -f <MULTI_FASTA_FILE> > <OUTFILE_PATH>"
```

For example:

```
sudo docker run --rm \
 -v /home/jupyter-pathinformatics/nikos_viral_vaccine:/iedb_io iedb-b-cell:latest /bin/bash -c \
   "python /bcell_standalone/predict_antibody_epitope.py -m Bepipred -f /iedb_io/sequences.txt > iedb_io/out.txt"
```

or with a SWISSPROT ID:

```
sudo docker run --rm \
 -v $PWD/iedb_b_cell_io:/iedb_b_cell_io iedb-b-cell:latest /bin/bash -c \
   "python /bcell_standalone/predict_antibody_epitope.py -m Bepipred -s P02185 > /iedb_b_cell_io/out.txt"
```
