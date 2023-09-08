#!/bin/bash

DIRNAME=$(dirname "$0")
echo "The script is located in: $DIRNAME"
cd "$DIRNAME"
for file in *minimization-out.mae; do
    base_name="${file%.*}"; # This strips the extension from the file name
    /usr/local/schrodinger/2022-2/utilities/structconvert "$file" "$base_name.sdf";
done
