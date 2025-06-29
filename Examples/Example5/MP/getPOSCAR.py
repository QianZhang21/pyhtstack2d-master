from mp_api.client import MPRester
from pymatgen.core.periodic_table import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import os
import json

my_api_key = "Your_Materials_Project_API_Key_Here"  # Replace with your actual Materials Project API key

if not os.path.exists("POSCAR_dir"):
    os.mkdir("POSCAR_dir")

metals = [el for el in Element if Element(el).is_metal]  # Get all metals from the periodic table
lanthanides = {57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
               89, 90, 91, 92, 93, 94}  # Lanthanides and actinides
metals = [el for el in metals if Element(el).Z not in lanthanides]  # Exclude lanthanides and actinides
dict_metal = {}  # Create a dictionary to store metal information
with MPRester(api_key=my_api_key) as mpr:  # Initialize MPRester with your API key
    for metal in metals:  # Iterate through each metal
        ls_mat_ids = mpr.get_material_ids(chemsys_formula=metal.symbol)  # Get material IDs for the metal
        for mat_id in ls_mat_ids:  # Iterate through each material ID
            # Search for materials with specific criteria
            mat_summary = mpr.summary.search(material_ids=mat_id, num_elements=(1, 1), num_sites=(1, 6),
                                             theoretical=False, total_magnetization=(0, 0.5), energy_above_hull=(0, 0),
                                             fields=["material_id", "structure", "nsites", "is_stable", "symmetry",
                                                     "band_gap", "formula_pretty", "energy_above_hull", "is_magnetic",
                                                     "nelements", "chemsys"])

            for mat in mat_summary:  # Iterate through each material summary
                structure = mat.structure  # Get the structure of the material
                sga = SpacegroupAnalyzer(structure)  # Analyze the symmetry of the structure
                primitive = sga.get_primitive_standard_structure()  # Get the primitive standard structure
                filename = f"{mat.formula_pretty}_{mat.material_id}-POSCAR"  # Create a filename based on the formula and material ID
                primitive.to(filename=os.path.join("POSCAR_dir", filename))  # Save the primitive structure to a POSCAR file
                print(f"Saved: {filename}")

                dict_metal[mat.material_id] = {}
                dict_metal[mat.material_id]['element'] = mat.chemsys
                dict_metal[mat.material_id]['nelements'] = mat.nelements
                dict_metal[mat.material_id]["nsites"] = mat.nsites
                dict_metal[mat.material_id]["symmetry"] = mat.symmetry.symbol
                dict_metal[mat.material_id]["is_stable"] = mat.is_stable
                dict_metal[mat.material_id]["is_magnetic"] = mat.is_magnetic
                dict_metal[mat.material_id]["energy_above_hull"] = mat.energy_above_hull
                dict_metal[mat.material_id]["band_gap"] = mat.band_gap

with open("nomagnetic_metal.json", "w") as f:
    json.dump(dict_metal, f, indent=2)
