# README

## Inputs

.PDB file(s)

## Outputs

Multiple files are created:

```
|epidope-|
|        |--[folder of .csv files for each sequence]
|
|plots-|
|      |--[folder of html5 plots for each sequence]
|
|epidope_scores.csv
|predicted_epitopes.csv
```

Importantly, we focus on `predicted_epitopes.csv` and `epidope_scores.csv` as the main outputs

## Command

### Build Container

Note that we just use the epidope:v0.3 container:

`sudo docker pull epidope:v0.3`

or

`singularity pull docker://flomock/epidope:v0.3`

### Run Container

Note that you should pick a host directory to mount as well as a host path in which to save the output:

```
sudo docker run --rm \
  -v <HOST_FOLDER_CONTAINING_MULTIFASTA_INPUT_FILE>:/epidope_io flomock/epidope:v0.3 \
  -i /epidope_io/<MOUNTED_PATH_TO_INPUT_MULTIFASTA_FILE> -o <MOUNTED_PATH_TO_SAVE_EPIDOPE_OUTPUT>
```

For example:

```
sudo docker run --rm \
  -v /home/jupyter-pathinformatics/nikos_viral_vaccine:/epidope_io flomock/epidope:v0.3 \
    -i /epidope_io/source_accessions/alphavirus_proteins/alphavirus_protein_multiseq.fasta -o /epidope_io/epidope_results/
```
