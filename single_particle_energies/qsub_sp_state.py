import os
import time

edf = input("Choose EDF (UNE0~2,SKM*,SKP,SLY4,SV-min):")
for isotope in ["100Sn","132Sn","16O","208Pb","40Ca","48Ca","56Ni","78Ni","90Zr","z126_n184"]:
    os.chdir(isotope)
    os.system("rm *.hel thoout* mass*.msub.* \n")
    f1 = open("hfbtho_NAMELIST.dat","r")
    l1 = f1.readlines()
    out_str = ""
    for line in l1:
        if "&" not in line: out_str += " "
        ss = line.split()
        if "functional" in line:
            ss[2] = "'" + edf + "',"
        for word in ss:
            out_str += word + " "
        out_str += "\n"
    out_file = open("hfbtho_NAMELIST.dat","w")
    out_file.write(out_str)
    out_file.close(); f1.close()
    os.system("cp ~/hfbthoV300/src/hfbtho/single_qp_energy/hfbtho_main ./ \n")
    os.system("qsub mass_table.msub \n")
    if False:
     for isospin in ["neutron","proton"]:
         os.chdir(isospin)
         time.sleep(1)
         lsLoc = os.listdir("./")
         print (lsLoc)
         for ff in lsLoc:
             print (ff)
             os.chdir(ff)
             os.system("rm *.hel thoout* mass*.msub.* \n")
             f1 = open("hfbtho_NAMELIST.dat","r")
             l1 = f1.readlines()
             out_str = ""
             for line in l1:
                 if "&" not in line: out_str += " "
                 ss = line.split()
                 if "functional" in line:
                     ss[2] = "'" + edf + "',"
                 for word in ss:
                     out_str += word + " "
                 out_str += "\n"
             out_file = open("hfbtho_NAMELIST.dat","w")
             out_file.write(out_str)
             out_file.close(); f1.close()
             os.system("cp ~/hfbthoV300/src/hfbtho/single_qp_energy/hfbtho_main ./ \n")
             os.system("qsub mass_table.msub \n")
             os.chdir("..")
         os.chdir("..")
    os.chdir("..")
