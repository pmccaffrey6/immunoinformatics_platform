# How to use

## Inputs

.PDB file(s)

## Outputs

.txt file that will need to be processed into a .csv or .tsv file, it looks like this:

```
DiscoTope predictions for '1z40.pdb'.
	Looking only at Chain:  A
	Contact Distance = 10.000 Angstroms
	Threshold = -7.700

A	108	ASN	14	-1.459	-8.459
A	109	PRO	11	0.724	-4.776	<=B
A	110	TRP	13	0.804	-5.696	<=B
A	111	THR	12	1.211	-4.789	<=B
A	112	GLU	11	1.331	-4.169	<=B
A	113	TYR	14	0.929	-6.071	<=B
A	114	MET	18	-0.779	-9.779
A	115	ALA	20	-0.444	-10.444
A	116	LYS	21	0.122	-10.378
A	117	TYR	24	-2.172	-14.172
```

## Command

### Build Container

`sudo docker build -t discotope .`

### Run Container

Note that you should pick a host directory to mount as well as a host path in which to save the output:

`sudo docker run --rm -v <HOST_DIRECTORY_WITH_PDB_FILE(S)>:/discotope_in discotope:latest -f /discotope_in/discotope-1.1/test/1z40.pdb -chain A > <PLACE_ON_HOST_TO_SAVE_OUTPUT>`

For example:

`sudo docker run --rm -v $PWD:/discotope_in discotope:latest -f /discotope_in/discotope-1.1/test/1z40.pdb -chain A > /discotope_in/out.txt`
