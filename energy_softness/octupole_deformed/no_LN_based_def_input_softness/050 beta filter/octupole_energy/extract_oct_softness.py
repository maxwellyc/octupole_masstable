import numpy as np

def softness(edf):
    f_gs = open("octupole_gs_reduced/"+edf+"_reduced_beta2_3_050.dat","r")
    f_control = open("control_gs_extract/reduced/Q20/"+edf+\
    "_control_groundstate_reduced_050.dat","r")
    f_035 = open("035_reduced/Q20/"+edf+"_control_groundstate_reduced_Q20.dat","r")
    
    l0 = f_gs.readlines()
    l1 = f_control.readlines(); l2 = f_035.readlines()
    f_gs.close(); f_control.close(); f_035.close()
    
    d_gs,d_control,d_soft = {},{},{}

    for ind,line in enumerate(l1):
        if ind:
            ss = line.split()
            if ss[-1] != "YES": print (f"Check error: {line}")
            Z,N,BE,beta2 = int(ss[0]),int(ss[1]),float(ss[3]),float(ss[4])
            d_control[(Z,N)] = (BE,beta2)
        
    for ind,line in enumerate(l2):
        if ind:
            ss = line.split()
            if ss[-1] != "YES": print (f"Check error: {line}")
            Z,N,BE,beta2 = int(ss[0]),int(ss[1]),float(ss[3]),float(ss[4])
            # if 0.35 set also has this nuclei
            if (Z,N) in d_control and BE < d_control[(Z,N)][0]:
                d_control[(Z,N)] = (BE,beta2)
            else:
                d_control[(Z,N)] = (BE,beta2)
    
        
    count = 0
    for line in l0:
        if not count:
            count+=1
            continue
        ss = line.split()
        if ss[-1] != "YES": print (f"Check error: {line}")
        Z,N,BE,beta2,beta3 = int(ss[0]),int(ss[1]),float(ss[3]),float(ss[4]),float(ss[5])
        d_gs[(Z,N)] = (BE,beta2,beta3)

    # write octupole softness data to file
    if True:
        f_out = open("octupole_energy/"+edf+"_octupole_energy_050.dat","w")
        outputStr = "Z".ljust(7)+"N".ljust(7)+"A".ljust(7)+"E_oct(MeV)".ljust(20)+ \
        "BE_gs(MeV)".ljust(20)+"BE_control(MeV)".ljust(20)+"beta3_gs".ljust(15)+ \
        "beta2_gs".ljust(15)+"Q20_control (b)".ljust(15)+"\n"
        for k in d_control:
            Z,N = k[0],k[1]
            outputStr += str(k[0]).ljust(7)+str(k[1]).ljust(7)+str(k[0]+k[1]).ljust(7)+\
            str(round(d_gs[k][0]-d_control[k][0],6)).ljust(20)+\
            str(d_gs[k][0]).ljust(20)+str(d_control[k][0]).ljust(20)+\
            str(d_gs[k][2]).ljust(15)+str(d_gs[k][1]).ljust(15)+str(d_control[k][1]).ljust(15)+"\n"

        f_out.write(outputStr)

EDF = ["UNEDF0","UNEDF1","UNEDF2","SLY4","SV-MIN","SKMS","SKP"]
for e in EDF:
    softness(e)
        
        


