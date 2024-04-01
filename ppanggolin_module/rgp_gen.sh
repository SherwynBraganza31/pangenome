#!/usr/bin/env bash

. "$( dirname "${BASH_SOURCE[0]}" )/filenames.sh"
cd "$pangenome_path"/PPaNGGOLiN_RUN || exit

# Generate Regions of genome plasticity
ppanggolin rgp -p "$base_filename".h5
ppanggolin write -p "$base_filename".h5 --regions --output RGP-Regions -f