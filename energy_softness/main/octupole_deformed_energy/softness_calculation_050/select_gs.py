import numpy as np
import math

def read_file(f_dir,f_name):
    f1 = open(f_dir+f_name, "r")
    l1 = f1.readlines()
    return l1

def select_gs(edf,cutoff=40):
    #edf = input("SELECT::EDF(default SLY4 -OR- SKMS, SKP, SV-MIN, UNEDF0~2):")
    if len(edf) == 0 : edf = "SLY4"
    f_dir = "raw/"
    f_name = edf +"_octupole_table_combined.dat"
    output_f1 = "new_data/" + edf + "_all_gs_data_beta2_3_0"+str(cutoff)+".dat"
    l1 = read_file(f_dir,f_name)
    visitedNuc = []
    outputFile = open("./"+output_f1,"w")
    outputStr = ""; i = 0; cc=0
    for line in l1:
        if i == 0:
            i += 1
            outputStr += line
        elif i > 0:
            ss = line.split()
            Z = ss[0];N = ss[1]; A = int(ss[2]); c_conv = ss[20]
            beta2 = float(ss[6]); Q30 = float(ss[12])
            r0=1.2*math.pow(A,1/3);r03=r0**3;c3=4*math.pi/3;c3=4*math.pi/3
            beta3=round((Q30/A/r03)*1000*c3,6)
            if (Z,N) not in visitedNuc and c_conv == "YES" and abs(beta3)<=cutoff/100 and abs(beta2)<=cutoff/100:
                if int(Z)>= 100 and abs(beta3) < 0.00001 and c_conv == "NO":
                    print (Z,N,beta3)
                visitedNuc.append((Z,N))
                outputStr += line
                cc+=1
    print (f"Total nuclei: {cc}, {edf}")

    outputFile.write(outputStr)
    outputFile.close()

def reduce_var(edf,cutoff=40):
    #edf = input("REDUCE::EDF(default SLY4 -OR- SKMS, SKP, SV-MIN, UNEDF0~2):")
    if len(edf) == 0 : edf = "SLY4"
    input_f1 = "raw/"+edf+"_softness_050.dat"
    output_f2 = edf+"_reduced_softness_050.dat"
    outputFile = open("reduced/"+output_f2,"w")
    l2 = read_file("",input_f1)
    i = 0
    uc = 0
    for line in l2:
        if i == 0:
            i += 1
            default_label = line
            outputStr = "Z".ljust(7)+"N".ljust(7)+"A".ljust(7)+"Binding_Energy_(MeV)".ljust(24)+"beta2_total".ljust(15)\
            +"beta3_total".ljust(15)+"fileID".ljust(10)+"Convergence".ljust(15)+"\n"
        elif i > 0:
            ss = line.split()
            #Z      N      A       Binding_Energy_(MeV)    Quad_Def_Beta2_P     Quad_Def_Beta2_N     Quad_Def_Beta2_total   Quad_Moment_Q2_P_(fm^2)    Quad_Moment_Q2_N_(fm^2)    Quad_Moment_Q2_total_(fm^2)     Octupole_Moment_Q3_P_(fm^3)     Octupole_Moment_Q3_N_(fm^3)     Octupole_Moment_Q3_total_(fm^3)      Pairing_gap_P_(MeV)     Pairing_gap_N_(MeV)   RMS_radius_P_(fm)    RMS_radius_N_(fm)    RMS_radius_total_(fm)    Charge_Radius_(fm)     File_ID
            Z=int(ss[0]);N=int(ss[1]);A=int(ss[2])
            f_id=ss[19]; convergence=ss[20]
            bind_E=ss[3];pairGap_P=ss[13];pairGap_N=ss[14]
            beta2_P=ss[4];beta2_N=ss[5];beta2_T=ss[6];
            Q2_P=ss[7];Q2_N=ss[8];Q2_T=ss[9];
            Q3_P=float(ss[10]);Q3_N=float(ss[11]);Q3_T=float(ss[12])
            r_P=ss[15];r_N=ss[16];r_T=ss[17];r_charge=ss[18]
            r0=1.2*math.pow(A,1/3);r03=r0**3;c3=4*math.pi/3
            beta3_P=round((Q3_P/Z/r03)*1000*c3,6)
            beta3_N=round((Q3_N/N/r03)*1000*c3,6)
            beta3_T=round((Q3_T/A/r03)*1000*c3,6)
            if convergence != "YES": uc+=1; print (Z,N,"unconverged",uc)
            outputStr += str(Z).ljust(7)+str(N).ljust(7)+str(A).ljust(7)+str(bind_E).ljust(24)+str(beta2_T).ljust(15)+str(beta3_T).ljust(15)+f_id.ljust(10)+convergence.ljust(15)+"\n"
    outputFile.write(outputStr)
    outputFile.close()

def compare_old():
    edf = input("EDF (defaults to SLY4):")
    if len(edf) == 0 : edf = "SLY4"
    input_f1 = edf+"_binding_beta2_3_filtered" + ".dat"
    #input_f2 = edf + "_even-even_massexplorer.dat"
    input_f2 = edf + "_even-even_nuclei-no_drip_lines.dat"
    #input_f3 = "AME2016/AME2016.dat"
    output_f1 = "compare_"+edf+".dat"
    outputFile = open(output_f1,"w")
    l1 = read_file("new_data/",input_f1)   #change this to new_data/ when comparing with octupole table
    l2 = read_file("old_data/",input_f2)
    visited1=[]; commonNuc=[]
    new_data={}
    i=0;j=0;count=0;compared=0
#   # For comparing Erik's old data with massexplorer
#    for line in l1:
#        if i==0: i+=1
#        elif i>0:
#            ss=line.split()
#            bindE=float(ss[4]);beta2=ss[9];beta3="0";f_id="old"
#            Z=int(ss[1]);N=int(ss[2]);A=int(ss[3])
#            visited1.append((Z,N))
#            new_data[(Z,N)]=[bindE,beta2,beta3,f_id]
    # For extracting data from my octupole table.
    for line in l1:
        if i==0: i+=1
        elif i>0:
            ss=line.split()
            Z=int(ss[0]);N=int(ss[1]);A=int(ss[2])
            visited1.append((Z,N))
            bindE=float(ss[3]);beta2=ss[4];beta3=ss[5];f_id=ss[6]
            new_data[(Z,N)]=[bindE,beta2,beta3,f_id]
    outputStr = "Z".ljust(7)+"N".ljust(7)+"A".ljust(7)+"Energy_new(MeV)".ljust(20)+"Energy_old(MeV)".ljust(20)+"Energy_diff(MeV)".ljust(20)+"beta2_new".ljust(13)+"beta2_old".ljust(13)+"beta3_new".ljust(13)+"fileID".ljust(10)+"\n"
    print (outputStr)
    for line in l2:
        if j==0: j+=1
        elif j>0:
            ss=line.split()
            Z=int(ss[1]);N=int(ss[2]);A=int(ss[3])
            if (Z,N) in visited1:
                compared+=1
                commonNuc.append((Z,N))
                bindE_old=float(ss[4]);beta2_old=ss[9]
                bindE_new=new_data[(Z,N)][0]
                diff = round(bindE_new - bindE_old,6)
                beta2_new=new_data[(Z,N)][1]
                beta3_new=new_data[(Z,N)][2]
                diff_beta2 = float(beta2_new) - float(beta2_old)
                f_id=new_data[(Z,N)][3]
                outputStr += str(Z).ljust(7)+str(N).ljust(7)+str(A).ljust(7)+str(bindE_new).ljust(20)+str(bindE_old).ljust(20)+str(diff).ljust(20)+beta2_new.ljust(13)\
            +beta2_old.ljust(13)+beta3_new.ljust(13)+f_id.ljust(10)+"\n"
                if Z>=106 and abs(float(beta3_new))<0.00001:#and abs(diff_beta2)<0.0001:
#                if abs(float(beta3_new))>0.5 or float(beta2_new) > 0.5:
                    count+=1
                    print (str(Z).ljust(7)+str(N).ljust(7)+str(A).ljust(7)+str(bindE_new).ljust(20)+str(bindE_old).ljust(20)+str(diff).ljust(20)+beta2_new.ljust(13)\
            +beta2_old.ljust(13)+beta3_new.ljust(13)+f_id.ljust(10)+str(count).ljust(10))
    print (f"Compared {compared} nucleus\t{edf}")
    outputFile.write(outputStr)
    outputFile.close()

#for edf in ["SLY4","SKMS","SKP","SV-MIN","UNEDF0","UNEDF1","UNEDF2"]:
#    select_gs(edf,100)
for edf in ["SLY4","SKMS","SKP","SV-MIN","UNEDF0","UNEDF1","UNEDF2"]:
    reduce_var(edf,100)
#compare_old()
