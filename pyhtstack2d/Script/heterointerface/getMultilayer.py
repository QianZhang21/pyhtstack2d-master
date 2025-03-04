import os
import shutil
from pyhtstack2d.buildbilayer.stackBilayer import Bilayer

# =============== Specify the Input Parameters ==============
# Intralayer distances (atomic spacing within each monolayer)
d_intra = [1.93161, 1.91135]
# Number of layers in each multilayer structure (nlayer + nlayer will form the heterointerface)
nlayer = 6
# Monolayer structure files for the upper and bottom layers
monoup = "LaAlO3-POSCAR.vasp"  # "mp-5304-POSCAR"
monobottom = "SrTiO3-POSCAR.vasp"  # "mp-5229-POSCAR"
# Interlayer distance for the final heterointerface
d_inter = 1.92358
# ==========================================================

# Generate multilayer structures for both the upper and bottom layers
for i, lstr in enumerate(["up", "bottom"]):  # Loop through both layers
    di = d_intra[i]  # Assign the corresponding intralayer distance
    mono = monoup if i == 0 else monobottom  # Select the correct monolayer file
    for ni in range(nlayer-1):  # Iteratively build up multilayers
        # For the first layer, use the monolayer structure; for subsequent layers, use the growing multilayer structure
        mono_multilayer = mono if ni == 0 else f"{lstr}_multilayer-POSCAR"
        # Ensure the working directory is clean before generating new structures
        if os.path.exists("BiPOSCAR_dir"):
            shutil.rmtree("BiPOSCAR_dir")
        # Generate the multilayer structure with the specified intralayer distance
        bigen = Bilayer(mono, mono_multilayer, d_inter=di, skip_xy_rev=True, lv=nlayer * 5)
        bigen.WritePOSCAR()
        # Define the path for the newly generated POSCAR file
        biposcar = os.path.join(bigen.savepath, bigen.formula_w, "AA", "cord1", "POSCAR")
        # Copy the generated POSCAR file to serve as input for the next layer
        shutil.copy(biposcar, f"{lstr}_multilayer-POSCAR")

# Construct the final heterointerface structure by stacking the two multilayers
bigen = Bilayer("up_multilayer-POSCAR", "bottom_multilayer-POSCAR", d_inter=d_inter, skip_xy_rev=True, lv=nlayer * 10)
bigen.WritePOSCAR()
# Retrieve the final heterointerface POSCAR file
biposcar = os.path.join(bigen.savepath, bigen.formula_w, "AA", "cord1", "POSCAR")
# Save the heterointerface POSCAR file as the final output
shutil.copy(biposcar, "heterointerface-POSCAR")
# Cleanup: Remove temporary directories and files to keep the workspace clean
shutil.rmtree("BiPOSCAR_dir")  # Remove temporary directory
os.remove("up_multilayer-POSCAR")  # Delete the intermediate multilayer file for the upper layer
os.remove("bottom_multilayer-POSCAR")  # Delete the intermediate multilayer file for the bottom layer

print("Heterointerface structure generation complete.")

