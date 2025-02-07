import os
"""
Changing elements in POSCAR file.
"""


def modify_elements_in_poscar(file_path, new_elements, save_path=None):
    # read the original POSCAR
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Modify elements
    elements_line_index = 5
    lines[elements_line_index] = '     '.join(new_elements) + '\n'
    elem_num = lines[elements_line_index + 1].split()
    assert len(elem_num) == len(new_elements), "The number of elements is not equal to the number of new elements."
    # Write new POSCAR
    fomula_w = ""
    for i, num in enumerate(elem_num):
        fomula_w += new_elements[i]+str(int(num))
    if save_path:
        savefile = os.path.join(save_path, fomula_w + '-POSCAR')
    else:
        savefile = fomula_w + '-POSCAR'
    with open(savefile, 'w') as file:
        file.writelines(lines)

# # Example
# if __name__ == '__main__':
#     file_path = "POSCAR"
#     new_elements = ['Ag', 'V', 'P', 'S']
#     save_path = "POSCAR_dir"
#     modify_elements_in_poscar(file_path, new_elements, save_path)
