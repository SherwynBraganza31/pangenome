#!/usr/bin/env bash

. "$( dirname "${BASH_SOURCE[0]}" )/filenames.sh"
cd "$pangenome_path"/PPaNGGOLiN_RUN || exit

# Generate Spot Figures
ppanggolin spot -p "$base_filename".h5 -f
ppanggolin write -p "$base_filename".h5 --spots --output spots -f
ppanggolin draw -p "$base_filename".h5 --spots all --output spot-figures -f