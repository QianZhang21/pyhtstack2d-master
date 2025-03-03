from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp import Poscar
import os

if not os.path.exists("POSCAR_dir"):
    os.mkdir("POSCAR_dir")

# Your Materials Project API key
API_KEY = "your_api_key_here"


def download_and_save_poscars(mpr, material_ids):
    """
    Download structures from Materials Project and save them as POSCAR files.
    """
    for material_id in material_ids:
        try:
            structure = mpr.get_structure_by_material_id(material_id)
            poscar = Poscar(structure)
            filename = os.path.join("POSCAR_dir", f"{material_id}-POSCAR")
            poscar.write_file(filename)
            print(f"✅ {material_id} saved as {filename}")
        except Exception as e:
            print(f"❌ Error downloading {material_id}: {e}")


def query_and_save_poscars(mpr, criteria):
    """
    Query materials from Materials Project based on criteria and save them as POSCAR files.
    """
    properties = ["material_id", "pretty_formula", "structure"]
    try:
        materials = mpr.query(criteria=criteria, properties=properties)
        print(f"🔍 Found {len(materials)} materials matching the criteria.")

        for material in materials:
            material_id = material["material_id"]
            structure = material["structure"]
            poscar = Poscar(structure)
            filename = os.path.join("POSCAR_dir", f"{material_id}-POSCAR")
            poscar.write_file(filename)
            print(f"✅ {material['pretty_formula']} ({material_id}) saved as {filename}")

    except Exception as e:
        print(f"❌ Error querying materials: {e}")


# Connect to Materials Project API
with MPRester(API_KEY) as mpr:
    # =================================
    # Manually specified material IDs
    material_ids = ["mp-5304", "mp-5229"]
    download_and_save_poscars(mpr, material_ids)
    # =================================

    # =================================
    # Query materials based on criteria
    criteria = {
        "elements": {"$all": ["Si", "O"]},  # Must contain Si and O
        "band_gap": {"$gte": 1.0, "$lte": 2.0},  # Band gap: 1.0 - 2.0 eV
        "density": {"$gte": 2},  # Density > 2 g/cm³
        "spacegroup.symbol": "Fm-3m"  # Space group Fm-3m
    }
    query_and_save_poscars(mpr, criteria)
    # =================================

