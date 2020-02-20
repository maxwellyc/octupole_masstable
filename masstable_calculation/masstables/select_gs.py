import numpy as np
import math

def read_file(f_dir,f_name):
    f1 = open(f_dir+f_name, "r")
    l1 = f1.readlines()
    return l1

def moment_to_beta(Q2, Q3, A):
    c2 = math.sqrt(5 * math.pi) / 3
    c3 = 4 * math.pi / 3  # ~=4.18879
    r0 = 1.2 * A**(1.0/3.0)
    b2 = Q2 / (A * r0**2 / c2 / 100)
    b3 = Q3 / (A * r0**3 / c3 / 1000)
    return round(b2,3),round(b3,3)

def select_gs(edf,cutoff=50):
    #edf = input("SELECT::EDF(default SLY4 -OR- SKMS, SKP, SV-MIN, UNEDF0~2):")
    if len(edf) == 0 : edf = "SLY4"
    f_dir = "EDF_extracted_raw/LN_raw/"
    f_name = "HFBTHOv300_"+edf+"_All_Data_20_shells_LN_deformation-masstable.dat"
    output_f1 = "new_data/LN_all_data/0"+str(cutoff)+"/"+ edf + "_all_gs_data_beta2_3_0"+str(cutoff)+".dat"
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
            Q20, Q30 = float(ss[9]), float(ss[12])
            beta2, beta3 = moment_to_beta(Q20, Q30, A)
            if (Z,N) not in visitedNuc and c_conv == "YES":
                if abs(beta3)<=cutoff/100 and abs(beta2)<=cutoff/100:
                    visitedNuc.append((Z,N))
                    outputStr += line
                    cc+=1
    print (f"Total nuclei: {cc}, {edf}")

    outputFile.write(outputStr)
    outputFile.close()

def reduce_var(edf,cutoff=50):
    outputStr = ""
    #edf = input("REDUCE::EDF(default SLY4 -OR- SKMS, SKP, SV-MIN, UNEDF0~2):")
    if len(edf) == 0 : edf = "SLY4"
    input_f1 = "new_data/LN_all_data/0"+str(cutoff)+"/"+edf+"_all_gs_data_beta2_3_0"+str(cutoff)+".dat"
    output_f2 = edf+"_reduced_beta2_3_0"+str(cutoff)+".dat"
    outputFile = open("new_data/LN_reduced/050/"+output_f2,"w")
    l2 = read_file("",input_f1)
    i = 0
    uc = 0
    for line in l2:
        if i == 0:
            i += 1
            default_label = line
            outputStr = "Z".ljust(7)+"N".ljust(7)+"A".ljust(7)+"Binding_Energy_(MeV)".ljust(24)+"beta2".ljust(15)\
            +"beta3".ljust(15)+"fileID".ljust(22)+"Convergence".ljust(15)+"\n"
        elif i > 0:
            ss = line.split()
            #Z      N      A       Binding_Energy_(MeV)    Quad_Def_Beta2_P     Quad_Def_Beta2_N     Quad_Def_Beta2_total   Quad_Moment_Q2_P_(fm^2)    Quad_Moment_Q2_N_(fm^2)    Quad_Moment_Q2_total_(fm^2)     Octupole_Moment_Q3_P_(fm^3)     Octupole_Moment_Q3_N_(fm^3)     Octupole_Moment_Q3_total_(fm^3)      Pairing_gap_P_(MeV)     Pairing_gap_N_(MeV)   RMS_radius_P_(fm)    RMS_radius_N_(fm)    RMS_radius_total_(fm)    Charge_Radius_(fm)     File_ID
            Z=int(ss[0]);N=int(ss[1]);A=int(ss[2])
            f_id=ss[19]; convergence=ss[20]
            bind_E=ss[3]
            # pairGap_P=ss[13];pairGap_N=ss[14]
            # beta2_P=ss[4];beta2_N=ss[5];beta2_T=ss[6];
            # r_P=ss[15];r_N=ss[16];r_T=ss[17];r_charge=ss[18]
            Q20, Q30 = float(ss[9]), float(ss[12])
            beta2, beta3 = moment_to_beta(Q20, Q30, A)
            if convergence != "YES": uc+=1; print (Z,N,"unconverged",uc)
            outputStr += str(Z).ljust(7)+str(N).ljust(7)+str(A).ljust(7)+str(bind_E).ljust(24)+str(beta2).ljust(15)+str(beta3).ljust(15)+f_id.ljust(22)+convergence.ljust(15)+"\n"
    outputFile.write(outputStr)
    outputFile.close()

for edf in ["SLY4","SV-MIN","UNEDF0","UNEDF1","UNEDF2"]:
    select_gs(edf,35)
for edf in ["SLY4","SV-MIN","UNEDF0","UNEDF1","UNEDF2"]:
    reduce_var(edf,35)
#compare_old()
