#!/usr/bin/env bash

. "$( dirname "${BASH_SOURCE[0]}" )/filenames.sh"
cd "$pangenome_path"/PPaNGGOLiN_RUN || exit

ppanggolin msa -p "$base_filename".h5 --phylo -o "./Phylo" -f
ppanggolin msa -p "$base_filename".h5 -o "./msa" -f