import os
import sys
import numpy as np
import math
import itertools as it

def beta_to_moment(beta2, beta3, A):
    c2 = math.sqrt(5 * math.pi) / 3
    c3 = 4 * math.pi / 3  # ~=4.18879
    r0 = 1.2 * A**(1.0/3.0)
    Q2 = beta2 * A * r0**2 / c2 / 100
    Q3 = beta3 * A * r0**3 / c3 / 1000
    return Q2,Q3



def create_masstable_out(z_list):
    # already calculated beta = 0.35 set, select only non-calculated
    f1 = open("beta3_001-050_count.dat", "r")
    f2 = open("UNEDF0_reduced_beta2_3_050.dat","r")
    l1 = f1.readlines(); l2 = f2.readlines()
    EDF_list, no_oct, oct = [],{},{}
    f1.close(); f2.close()
    # beta=0.35 set
    for line in l1:
        ss = line.split()
        if not EDF_list:
            EDF_list = ss[2:-1]; print (EDF_list)
            continue
        Z,N = int(ss[0]),int(ss[1])
        if Z in z_list:
            for ind,edf in enumerate(EDF_list):
                if int(ss[2+ind]):
                    if edf not in oct:
                        oct[edf] = {(Z,N)}
                    else:
                        oct[edf].add((Z,N))
    #print (oct["UNEDF0"])
    for ind,line in enumerate(l2):
        if ind:
            ss = line.split()
            Z,N = int(ss[0]),int(ss[1])
            if Z in z_list:
                for e in EDF_list:
                    if (Z,N) not in oct[e]:
                        no_oct[e] = no_oct.setdefault(e,[]) + [(Z,N)]


    # for each edf, find nuclei in column of their own.
    beta2_min, beta2_max, beta2_step, beta3 = -0.5, 0.5, 0.05, 0.0
    for edf in EDF_list:
        count,output_content = 0,""
        f = open(f"hfbtho_MASSTABLE_{edf}.dat","w")
        print (edf,len(no_oct[edf]))
        for z,n in no_oct[edf]:
            A = z+n
            # range for beta2 calculations
            for b2 in np.arange(beta2_min, beta2_max+0.00001,beta2_step):
                Q2,Q3 = beta_to_moment(b2,beta3,A)
                Q2,Q3,b2 = round(Q2,10), round(Q3,10), round(b2,2)
                output_content += str(z).rjust(6) + str(n).rjust(6) + str(Q2).rjust(18) +\
                                  str(Q3).rjust(18) + str(b2).rjust(12) + str(beta3).rjust(12) +" \n"
                count += 1

        output_head = " Nrows \n" + "   " + str(count) + " \n" + "    Z     N          Q20                 Q30           beta2       beta3\n"
        output_content = output_head + output_content
        f.write(output_content)
        f.close()
        print (f"{edf}, Row count: {count}")

# Zr, Xe, Ba, Ce, Hg; 
# Rn, Ra, Th, U, Pu.
z_list = [40,54,56,58,80] + list(range(86,102,2))
create_masstable_out(z_list)
