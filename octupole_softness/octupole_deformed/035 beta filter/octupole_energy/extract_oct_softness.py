import numpy as np

def softness():
    edf = input("SELECT::EDF(default SLY4 -OR- SKMS, SKP, SV-MIN, UNEDF0~2):") or "SLY4"
    
    f_gs = open("octupole_gs_reduced/"+edf+"_gs_reduced.dat","r")
    f_control = open("control_gs_extract/reduced/"+edf+"_control_groundstate_reduced.dat","r")
    
    l0 = f_gs.readlines(); l1 = f_control.readlines()
    f_gs.close(); f_control.close()
    
    d_gs,d_control,d_soft = {},{},{}
    count = 0
    for line in l1:
        if not count:
            count+=1
            continue
        ss = line.split()
        if ss[-1] != "YES": print (f"Check error: {line}")
        Z,N,BE,beta2 = int(ss[0]),int(ss[1]),float(ss[3]),float(ss[4])
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
        f_out = open("results/"+edf+"_octupole_energy.dat","w")
        outputStr = "Z".ljust(7)+"N".ljust(7)+"A".ljust(7)+"E_oct(MeV)".ljust(20)+ \
        "BE_gs(MeV)".ljust(20)+"BE_control(MeV)".ljust(20)+"beta3_gs".ljust(15)+ \
        "beta2_gs".ljust(15)+"beta2_control".ljust(15)+"delta_beta2".ljust(15)+"\n"
        for k in d_control:
            Z,N = k[0],k[1]
            outputStr += str(k[0]).ljust(7)+str(k[1]).ljust(7)+str(k[0]+k[1]).ljust(7)+\
            str(round(d_gs[k][0]-d_control[k][0],6)).ljust(20)+\
            str(d_gs[k][0]).ljust(20)+str(d_control[k][0]).ljust(20)+\
            str(d_gs[k][2]).ljust(15)+str(d_gs[k][1]).ljust(15)+str(d_control[k][1]).ljust(15)+\
            str(round(d_gs[k][1]-d_control[k][1],6)).ljust(15)+"\n"
            if d_gs[k][0]-d_control[k][0] > 0 and (180<=N<=210 or 120<=N<=140):
                print (Z,N,d_gs[k][0]-d_control[k][0])

        f_out.write(outputStr)

softness()
        
        


