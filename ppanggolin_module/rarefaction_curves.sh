#!/usr/bin/env bash

. "$( dirname "${BASH_SOURCE[0]}" )/filenames.sh"
cd "$pangenome_path"/PPaNGGOLiN_RUN || exit

# Generate Regions of genome plasticity


ppanggolin rarefaction -p "$base_filename".h5 --output rarefaction -f