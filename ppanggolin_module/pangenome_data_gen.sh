#!/usr/bin/env bash

. "$( dirname "${BASH_SOURCE[0]}" )/filenames.sh"
cd "$pangenome_path"/PPaNGGOLiN_RUN || exit

# Generate Pangenome Data
ppanggolin graph -p "$base_filename".h5 -f
ppanggolin partition -p "$base_filename".h5 -f
ppanggolin info -p "$base_filename".h5 --content
