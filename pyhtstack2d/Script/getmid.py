def getmidfunc(multilevel=None, scf=""):
    if multilevel is None or multilevel == 2:
        multilevel = None
    else:
        assert isinstance(multilevel, int) and multilevel > 2, "multilevel should be an integer greater than 2."
        multilevel = multilevel - 1

    getsh = "getmaterial.sh"
    materialfile = "material.txt"
    if scf != "" and scf.endswith("/"):
        scfdir = scf
    elif scf != "":
        scfdir = scf + "/"
    else:
        scfdir = ""
    with open(getsh, "w", newline="\n") as f:
        f.write("#!/bin/bash\n\n")
        f.write(f"rm -f {materialfile}\n\n")
        if multilevel is None:
            f.write("for file in $(ls -d */); do\n")
            f.write(
                f"    if [ -f $file/{scfdir}OUTCAR ] && [ `grep -c \"reached\" $file/{scfdir}OUTCAR` -ne 0 ] && [ `grep -c \"Voluntary\" $file/{scfdir}OUTCAR` -ne 0 ]; then\n")
            f.write(f"        echo $file >> {materialfile}\n")
            f.write("    fi\n")
            f.write("done\n")
        else:
            file_path = "$file"
            cd_com = "../" * multilevel
            f.write("for file in $(ls -d */); do\n")
            f.write("    cd $file\n")
            for level_i in range(multilevel - 1):
                file_path += f"$file_{level_i}"
                f.write("    " * level_i + f"    for file_{level_i} in $(ls -d */)\n")
                f.write("    " * level_i + "    do\n")
                f.write("    " * level_i + f"        cd $file_{level_i}\n")

            f.write("    " * (
                        multilevel - 2) + f"        if [ -f {scfdir}OUTCAR ] && [ `grep -c \"reached\" {scfdir}OUTCAR` -ne 0 ] && [ `grep -c \"Voluntary\" {scfdir}OUTCAR` -ne 0 ]; then\n")
            f.write("    " * (
                    multilevel - 2) + f"                echo {file_path} >> {cd_com}{materialfile}\n")
            f.write("    " * (multilevel - 2) + "        fi\n")
            for level_i in range(multilevel - 1):
                f.write("    " * (multilevel - 2 - level_i) + "        cd ..\n")
                f.write("    " * (multilevel - 2 - level_i) + "    done\n")
            f.write("    cd ..\n")
            f.write("done\n")
