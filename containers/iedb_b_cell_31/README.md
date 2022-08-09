# README

## Inputs

fasta file (or SWISSPROT ID)

## Outputs



## Command

### Build Container

Note that we just use the epidope:v0.3 container:

`sudo docker build -t iedb-b-cell .`

### Run Container

Note that you should pick a host directory to mount as well as a host path in which to save the output:

```
sudo docker run --rm -t \
 -v <HOST_FOLDER_CONTAINING_MULTIFASTA_INPUT_FILE>:/iedb_io iedb-b-cell:latest \
 -m <B_CELL_METHOD_TO_USE> -f <MULTI_FASTA_FILE>
```

For example:

```
sudo docker run --rm -t \
 -v /home/jupyter-pathinformatics/nikos_viral_vaccine:/iedb_io iedb-b-cell:latest \
   -m Bepipred-2.0 -f /iedb_io/sequences.txt
```
