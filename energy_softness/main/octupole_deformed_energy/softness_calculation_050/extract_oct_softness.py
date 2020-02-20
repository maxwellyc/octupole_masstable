import numpy as np

def softness(edf):
    f_b2b3 = open("beta2_min_inputs/"+edf+"_octupole_energy_050.dat","r")
    f_soft = open("reduced/"+edf+"_reduced_softness_050.dat","r")
    
    l0 = f_b2b3.readlines();l1 = f_soft.readlines()
    f_b2b3.close(); f_soft.close()
    
    d_b2b3,d_soft = {},{}

    for ind,line in enumerate(l0):
        if ind:
            ss = line.split()
            # BE_gs is BE of true ground state in octupole masstable
            # BE_b2 is the BE of previous control (b2 only, hence _b2) group,
            # with reflection sym. turned on.
            Z,N,BE_gs,BE_b2 = int(ss[0]),int(ss[1]),float(ss[4]),float(ss[5])
            # b2 of gs, b3 of gs, b2 of control(b2 only, hence _b2)
            b2_gs,b3_gs,b2_b2 = float(ss[7]),float(ss[6]),float(ss[8])
            d_b2b3[(Z,N)] = (BE_gs,BE_b2,b2_gs,b3_gs,b2_b2)
        
    for ind,line in enumerate(l1):
        if ind:
            ss = line.split()
            if ss[-1] != "YES": print (f"Unconverged: {edf}\t{line}")
            Z,N,BE,b2,b3,conv = int(ss[0]),int(ss[1]),float(ss[3]),float(ss[4]),float(ss[5]),ss[-1]
            # if 0.35 set also has this nuclei
            d_soft[(Z,N)] = (BE,b2,b3,conv)

    # write octupole softness data to file
    if True:
        f_out = open("softness_data/"+edf+"_octupole_energy_050.dat","w")
        outputStr = "Z".ljust(7)+"N".ljust(7)+"A".ljust(7)+"E_soft(MeV)".ljust(20)+ \
        "BE_b3=0 (MeV)".ljust(20)+"BE_b3=0.01 (MeV)".ljust(20)+"beta2".ljust(15)+ \
        "beta3 (0.01)".ljust(15)+"delta_beta2".ljust(15)+"CONV".ljust(9)+"\n"
        for k in d_soft:
            Z,N = k[0],k[1]
            outputStr += str(k[0]).ljust(7)+str(k[1]).ljust(7)+str(k[0]+k[1]).ljust(7)+\
            str(round(d_soft[k][0]-d_b2b3[k][1],6)).ljust(20)+\
            str(d_b2b3[k][1]).ljust(20)+str(d_soft[k][0]).ljust(20)+\
            str(d_soft[k][1]).ljust(15)+str(d_soft[k][2]).ljust(15)+\
            str(round(d_soft[k][1]- d_b2b3[k][2],6)).ljust(15) + d_soft[k][3]+"\n"

        f_out.write(outputStr)

EDF = ["UNEDF0","UNEDF1","UNEDF2","SLY4","SV-MIN","SKMS","SKP"]
for e in EDF:
    softness(e)
        
        


