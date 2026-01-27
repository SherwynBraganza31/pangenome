#!/usr/bin/env bash

. "$( dirname "${BASH_SOURCE[0]}" )/filenames.sh"
cd "$pangenome_path"/PPaNGGOLiN_RUN || exit

# Generate Curves and Matrix Files
ppanggolin draw -p "$base_filename".h5 --ucurve --tile_plot --output "$base_filename"_withcloud -f
ppanggolin draw -p "$base_filename".h5 --tile_plot --nocloud --output "$base_filename"_nocloud -f
ppanggolin write -p "$base_filename".h5 --gexf --csv --output matrix -f
ppanggolin write -p "$base_filename".h5 --light_gexf --output light_gexf -f
