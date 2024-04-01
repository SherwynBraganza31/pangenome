#!/usr/bin/env bash

. "$( dirname "${BASH_SOURCE[0]}" )/filenames.sh"
cd "$pangenome_path"/PPaNGGOLiN_RUN || exit

# Generate Partitions, Projections, Gene Families
ppanggolin write -p "$base_filename".h5 --partitions --output partitions -f
ppanggolin write -p "$base_filename".h5 --projection --output projection -f
ppanggolin write -p "$base_filename".h5 --families_tsv --output families-tsv -f