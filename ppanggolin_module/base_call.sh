#!/usr/bin/env bash

. "$( dirname "${BASH_SOURCE[0]}" )/filenames.sh"
cd "$pangenome_path"/ppan_dataset || exit

# create a clustered .h5 file
ppanggolin annotate --anno organisms.gff.list --fasta organisms.fasta.list --output "$pangenome_path"/PPaNGGOLiN_RUN --basename $base_filename -f --cpu 20
cd "$pangenome_path"/PPaNGGOLiN_RUN || exit
ppanggolin cluster -p "$base_filename".h5 --identity 0.5 --coverage 0.8 -f --cpu 20
#ppanggolin workflow --anno organisms.gff.list -o "$pangenome_path/PPaNGGOLiN_RUN" --basename $base_filename -f --cpu 20|| exit

