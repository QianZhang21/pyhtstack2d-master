import re
import math
from statistics import mode

mag_dict = {}
mag_list = []
nomag_list = []
count = 0

with open("mag.txt", "r") as f:
    for line in f:
        line_info = line.split(",")
        material_id = line_info[0]
        mid = material_id.split("/")[0]
        mag = float(line_info[1])
        mismatch = re.findall(r"\d+\.\d+|\d+", (material_id.split("/")[0]).split("_")[1])[0]
        if float(mismatch) >= 3:
            count += 1
            continue

        if mid not in mag_dict:
            mag_dict[mid] = {"mid": [], "mag": []}

        mag_dict[mid]["mid"].append(material_id)
        mag_dict[mid]["mag"].append(mag)

        if abs(mag) > 0.1:
            mag_list.append(material_id)
        else:
            nomag_list.append(material_id)

print(count, len(mag_list), len(nomag_list))
rc_list = []
for mid in mag_dict:
    # print(mid, mag_dict[mid]["mag"])
    mid_mag_list = [abs(num) for num in mag_dict[mid]["mag"]]
    # most_common_mag = round(max(mid_mag_list))
    mid_mag_list_round = [round(num) for num in mid_mag_list]
    most_common_mag = mode(mid_mag_list_round)
    non_mode_indexes = [index for index, value in enumerate(mid_mag_list_round) if value != most_common_mag]
    for index in non_mode_indexes:
        rc_list.append(mag_dict[mid]["mid"][index])

with open("magrc.txt", "w", newline="\n") as f:
    for material_id in rc_list:
        mid = material_id.split("/")[0]
        print(material_id, mag_dict[mid]["mag"])
        f.write(material_id + "\n")

