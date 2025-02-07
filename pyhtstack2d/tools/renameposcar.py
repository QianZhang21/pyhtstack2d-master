import os


def swap_poscar_filename(folder_path):
    if not os.path.exists(folder_path):
        print("The folder does not exist.")
        return

    files = os.listdir(folder_path)
    for filename in files:
        # Check if the file name starts with "POSCAR-"
        if filename.startswith("POSCAR-"):
            # Get the new name of the file
            new_name = filename.replace("POSCAR-", "", 1) + "-POSCAR"

            # Get the full path of the file
            old_file = os.path.join(folder_path, filename)
            new_file = os.path.join(folder_path, new_name)

            # Rename the file
            os.rename(old_file, new_file)
            print(f"File '{filename}' renamed to '{new_name}'.")
