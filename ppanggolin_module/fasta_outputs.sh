#!/usr/bin/env bash

. "$( dirname "${BASH_SOURCE[0]}" )/filenames.sh"
cd "$pangenome_path"/PPaNGGOLiN_RUN || exit

# Output FASTA files of Gene Data
ppanggolin fasta -p "$base_filename".h5 --output MY_GENES --genes all -f
ppanggolin fasta -p "$base_filename".h5 --output MY_GENES/persistent --genes persistent -f
ppanggolin fasta -p "$base_filename".h5 --output MY_GENES/shell --genes shell -f
ppanggolin fasta -p "$base_filename".h5 --output MY_GENES/cloud --genes cloud -f
ppanggolin fasta -p "$base_filename".h5 --output MY_GENES/core --genes core -f
ppanggolin fasta -p "$base_filename".h5 --output MY_GENES/softcore --genes softcore -f
ppanggolin fasta -p "$base_filename".h5 --output MY_GENES/rgp --genes rgp -f

# Output FASTA files of Protein Families
ppanggolin fasta -p "$base_filename".h5 --output MY_PROT --prot_families all -f
ppanggolin fasta -p "$base_filename".h5 --output MY_PROT/persistent --prot_families persistent -f
ppanggolin fasta -p "$base_filename".h5 --output MY_PROT/shell --prot_families shell -f
ppanggolin fasta -p "$base_filename".h5 --output MY_PROT/cloud --prot_families cloud -f
ppanggolin fasta -p "$base_filename".h5 --output MY_PROT/core --prot_families core -f
ppanggolin fasta -p "$base_filename".h5 --output MY_PROT/softcore --prot_families softcore -f
ppanggolin fasta -p "$base_filename".h5 --output MY_PROT/rgp --prot_families rgp -f
