# import re
#
# mid_list = []
#
# with open("entropy.txt", "r") as f:
#     for line in f:
#         line_info = line.split(",")
#         material_id = line_info[0]
#         entropy = float(line_info[1])
#         mismatch = re.findall(r"\d+\.\d+|\d+", (material_id.split("/")[0]).split("_")[1])[0]
#         if float(mismatch) >= 3:
#             continue
#
#         if abs(entropy) > 0.001:
#             mid_list.append(material_id)
#
#
# with open("entropymid.txt", "w", newline="\n") as f:
#     for material_id in mid_list:
#         f.write(material_id + "\n")


from pyhtstack2d.tools.genInput import GetMagEntropy

GetMagEntropy(multilevel=4, magpath="opt/").separateEntropy(recalheader="header.sh")
