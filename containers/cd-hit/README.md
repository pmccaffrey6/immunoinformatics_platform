# How to Use

## Commands

### Build Container
`sudo docker build -t cd-hit-4-8-1 .`

### Run Container
`sudo docker run --rm -v <HOST_DIRECTORY_WITH_INPUT_FASTA_FILE>:/cdhit_io cd-hit-4-8-1:latest cd-hit -i <MOUNTED_PATH_TO_INPUT_FILE> -o <MOUNTED_PATH_TO_OUTPUT_FILE>`

For example:

`sudo docker run --rm -v $PWD/cdhit_io:/cdhit_io cd-hit-4-8-1:latest cd-hit -i /cdhit_io/alphavirus_protein_multiseq_test_small.fasta -o /cdhit_io/out.txt`
