import re

mag_list = []
mag_mid_list = []
nomag_list = []
nomag_mid_list = []
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
        if abs(mag) > 0.1:
            mag_list.append(material_id)
            if mid not in mag_mid_list:
                mag_mid_list.append(mid)
        else:
            nomag_list.append(material_id)
            if mid not in nomag_mid_list:
                nomag_mid_list.append(mid)

print(count, len(mag_list), len(nomag_list))
# print(list(set(mag_mid_list) & set(nomag_mid_list)))
rc = list(set(mag_mid_list) & set(nomag_mid_list))
rc_list = []
for rc_ in nomag_list:
    rc_head = rc_.split("/")[0]
    if rc_head in rc:
        rc_list.append(rc_)

with open("maglist.txt", "w", newline="\n") as f:
    for material_id in mag_list:
        f.write(material_id + "\n")

with open("nomaglist.txt", "w", newline="\n") as f:
    for material_id in nomag_list:
        f.write(material_id + "\n")

with open("magrc.txt", "w", newline="\n") as f:
    for material_id in rc_list:
        f.write(material_id + "\n")

