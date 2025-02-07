from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure
import os


def get_all_files(folder_path):
    file_paths = []
    for dirpath, dirnames, filenames in os.walk(folder_path):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            file_paths.append(file_path)
    return file_paths


def get_subfolders(folder_path):
    subfolders = []
    for item in os.listdir(folder_path):
        item_path = os.path.join(folder_path, item)
        if os.path.isdir(item_path):
            subfolders.append(item_path)
    return subfolders


def delete_empty_subfolders(folder_path):
    for dirpath, dirnames, filenames in os.walk(folder_path, topdown=False):
        if not filenames and not dirnames:
            try:
                # print(f"Deleting empty folder: {dirpath}")
                os.rmdir(dirpath)
            except OSError as e:
                print(f"Error: {e.strerror} - {dirpath}")

    for dirpath, dirnames, filenames in os.walk(folder_path, topdown=False):
        if not filenames and not dirnames:
            try:
                # print(f"Deleting empty folder: {dirpath}")
                os.rmdir(dirpath)
            except OSError as e:
                print(f"Error: {e.strerror} - {dirpath}")


def remove_duplicates(posdir, multist=True, threshold=0.01, deleteempty=True):
    """
    Remove duplicates from a list of structures.

    posdir: str
        Path to the directory containing the POSCAR files.
    multist: bool, optional, default=True
        Flag indicating whether the directory contains multiple different systems.
        If True, the function will handle multiple systems; if False, it assumes
        there is only one system in the directory.
    threshold: float, optional, default=0.01
        Tolerance level for comparing the structures. The function considers
        two structures as duplicates if their distance is below this threshold.
    deleteempty: bool, optional, default=True
        Flag indicating whether to delete empty directories after removing duplicates.
        If True, any empty directories will be deleted.
    """
    matcher = StructureMatcher()
    if multist:
         posdirlist = get_subfolders(posdir)
         for psd in posdirlist:
             remove_duplicates(psd, multist=False, deleteempty=False)
         delete_empty_subfolders(posdir)
    else:
        similar_pos = []
        stfiles = get_all_files(posdir)
        for i in range(len(stfiles)):
            for j in range(i+1, len(stfiles)):
                stfile1 = stfiles[i]
                stfile2 = stfiles[j]
                st1 = Structure.from_file(stfile1)
                st2 = Structure.from_file(stfile2)
                are_similar = matcher.fit(st1, st2)
                if are_similar:
                    score = matcher.get_rms_dist(st1, st2)
                    if score[1] < threshold:
                        similar_pos.append(stfile2)
        for rpos in set(similar_pos):
            os.remove(rpos)
        if deleteempty:
            delete_empty_subfolders(posdir)
