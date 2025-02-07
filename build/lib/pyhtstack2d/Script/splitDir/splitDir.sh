#!/bin/bash

src_dir="path/to/source_directory"
dest_dir_base="path/to/destination_directory_base"
max_subdirs=50
counter=0
dir_index=1

# Get all files in the source directory
files=("$src_dir"/*)
files_only=()
for item in "${files[@]}"; do
    if [ -f "$item" ]; then
        files_only+=("$item")
    fi
done

for subdir in "$src_dir"/*; do
    if [ -d "$subdir" ]; then
        if [ $counter -eq 0 ]; then
            dest_dir="${dest_dir_base}_${dir_index}"
            mkdir -p "$dest_dir"
            ((dir_index++))

            # Copy files to destination directory
            for file in "${files_only[@]}"; do
                cp "$file" "$dest_dir"
            done
        fi

        mv "$subdir" "$dest_dir"
        ((counter++))

        if [ $counter -eq $max_subdirs ]; then
            counter=0
        fi
    fi
done