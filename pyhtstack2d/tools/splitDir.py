import os
import shutil

# Bash script equivalent to the Python function below
""" 
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
"""


def split_directories(src_dir, dest_dir_base, max_subdirs):
    subdirs = [d for d in os.listdir(src_dir) if os.path.isdir(os.path.join(src_dir, d))]
    files = [f for f in os.listdir(src_dir) if os.path.isfile(os.path.join(src_dir, f))]

    # Create destination directories and move subdirectories
    for i in range(0, len(subdirs), max_subdirs):
        dest_dir = f"{dest_dir_base}_{i // max_subdirs + 1}"
        os.makedirs(dest_dir, exist_ok=True)

        # Copy files to destination directory
        for file in files:
            shutil.copy(os.path.join(src_dir, file), os.path.join(dest_dir, file))

        for subdir in subdirs[i:i + max_subdirs]:
            shutil.move(os.path.join(src_dir, subdir), os.path.join(dest_dir, subdir))


if __name__ == "__main__":
    """
    When there are a large number of calculations, use this script to split subdirectories and then perform calculations.
    The split_directories function takes a source directory, a destination directory base, and a maximum number of subdirectories as arguments.
    """

    src_directory = "path/to/source_directory"
    dest_directory_base = "path/to/destination_directory_base"
    max_subdirectories = 50

    split_directories(src_directory, dest_directory_base, max_subdirectories)


