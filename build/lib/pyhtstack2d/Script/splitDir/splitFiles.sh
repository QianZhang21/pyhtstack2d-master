#!/bin/bash

src_dir="path/to/source_directory"
dest_dir_base="path/to/destination_directory_base"
max_subdirs=35
counter=0
dir_index=1


for subdir in "$src_dir"/*; do
    if [ -f "$subdir" ]; then
        if [ $counter -eq 0 ]; then
            dest_dir="${dest_dir_base}_${dir_index}"
            mkdir -p "$dest_dir"
            ((dir_index++))

        fi

        mv "$subdir" "$dest_dir"
        ((counter++))

        if [ $counter -eq $max_subdirs ]; then
            counter=0
        fi
    fi
done
