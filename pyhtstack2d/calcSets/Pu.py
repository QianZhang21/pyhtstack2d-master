import itertools

# Common valence states of elements
import os
import json

uinfopath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'U_info.json')
with open(uinfopath, "r") as f:
    Udict = json.load(f)


def INCARPu(elemlist, incar=None, u_setting=None, magatomlist=None, openmixparam=True, mixparam=None):
    """
    Add U setting to the INCAR file.

    elemlist: list of elements in the structure
    incar: path to the INCAR file
    u_setting: dictionary of U settings for specific elements
    magatomlist: list of magnetic atoms, if not None, only magnetic atoms will be considered for U setting
    openmixparam: if True, open the mixing parameters in the INCAR file
    mixparam: mixing parameters
    incarbasic: path to the basic INCAR file, if not None, the basic INCAR file will be copied to the new INCAR file
    """
    LDAU_info = {'LDAU': '.TRUE.', 'LDAUTYPE': 2, "LDAUL": [], "LDAUU": [], "LDAUJ": []}
    mixparam_default = {'LMAXMIX': 4, 'AMIX': 0.2, 'BMIX': 0.0001, 'AMIX_MAG': 0.8, 'BMIX_MAG': 0.0001}
    if mixparam:
        mixparam_default.update(mixparam)
    if openmixparam:
        LDAU_info.update(mixparam_default)

    for elem_ in elemlist:
        if magatomlist is not None and elem_ not in magatomlist:
            LDAU_info["LDAUL"].append(0)
            LDAU_info["LDAUU"].append(0)
            LDAU_info["LDAUJ"].append(0)
        elif u_setting is not None and elem_ in u_setting:
            LDAU_info["LDAUL"].append(u_setting[elem_][0])
            LDAU_info["LDAUU"].append(u_setting[elem_][1])
            LDAU_info["LDAUJ"].append(u_setting[elem_][2])
        elif elem_ in Udict["period_info"]:
            LDAU_info["LDAUL"].append(Udict["U_info"][Udict["period_info"][elem_]][0])
            LDAU_info["LDAUU"].append(Udict["U_info"][Udict["period_info"][elem_]][1])
            LDAU_info["LDAUJ"].append(Udict["U_info"][Udict["period_info"][elem_]][2])
        else:
            LDAU_info["LDAUL"].append(Udict["U_info"]["other"][0])
            LDAU_info["LDAUU"].append(Udict["U_info"]["other"][1])
            LDAU_info["LDAUJ"].append(Udict["U_info"]["other"][2])

    if incar is not None:
        assert os.path.exists(incar), "INCAR file not found!"
        with open(incar, 'a+', newline='\n') as wf:
            wf.write("LDAU = .TRUE.\n")
            wf.write("LDAUTYPE = 2\n")
            if openmixparam:
                for key, value in mixparam_default.items():
                    wf.write(key + " = " + str(value) + "\n")
            wf.write("LDAUL = ")
            for i in range(len(elemlist)):
                wf.write(str(LDAU_info["LDAUL"][i]) + " ")
            wf.write("\n")
            wf.write("LDAUU = ")
            for i in range(len(elemlist)):
                wf.write(str(LDAU_info["LDAUU"][i]) + " ")
            wf.write("\n")
            wf.write("LDAUJ = ")
            for i in range(len(elemlist)):
                wf.write(str(LDAU_info["LDAUJ"][i]) + " ")
            wf.write("\n")

    return LDAU_info
