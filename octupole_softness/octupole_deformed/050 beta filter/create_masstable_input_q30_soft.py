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



def create_masstable_out(b_min=1,b_max=50):
    # already calculated beta = 0.35 set, select only non-calculated
    f0_input = open("beta3_001-035_count.dat", "r")
    f_input = open("beta3_00"+str(b_min)+"-0"+str(b_max)+"_count.dat", "r")
    lines0 = f0_input.readlines()
    lines = f_input.readlines()
    EDF_list, oct_dict, oct_dict0 = [],{},{}
    
    # beta=0.35 set
    for line in lines0:
        ss = line.split()
        if not EDF_list:
            EDF_list = ss[2:-1]; print (EDF_list)
            continue
        Z,N = int(ss[0]),int(ss[1])
        for ind,edf in enumerate(EDF_list):
            if int(ss[2+ind]):
                oct_dict0[edf] = oct_dict0.setdefault(edf,[]) + [(Z,N)]


    EDF_list = [];
    # cross holds nuclei that're both in set 0.35 and set 0.5, we need only
    # do the +/- 0.35~0.5 part for these.
    oct_dict_cross = {}
    for line in lines:
        ss = line.split()
        if not EDF_list:
            EDF_list = ss[2:-1]; print (EDF_list)
            continue
        Z,N = int(ss[0]),int(ss[1])
        for ind,edf in enumerate(EDF_list):
            if int(ss[2+ind]) and (Z,N) not in oct_dict0[edf]:
                oct_dict[edf] = oct_dict.setdefault(edf,[]) + [(Z,N)]
            if int(ss[2+ind]) and (Z,N) in oct_dict0[edf]:
                oct_dict_cross[edf] = oct_dict_cross.setdefault(edf,[]) + [(Z,N)]

    # for each edf, find nuclei in column of their own.
    beta2_min, beta2_max, beta2_step = -0.4, 0.4, 0.05
    for edf in EDF_list:
        count,output_content = 0,""
        f = open(f"hfbtho_MASSTABLE_{edf}.dat","w")
        print (edf,len(oct_dict[edf]))
        for z,n in oct_dict[edf]:
            A = z+n
            # range for beta2 calculations
            for k in np.arange(-b_max/100, b_max/100+0.00001,0.05):
                Q2,Q3 = beta_to_moment(k,0.0,A)
                kk,Q2,ll,Q3 = round(k,2), round(Q2,10), 0.0, 0.0
                output_content += str(z).rjust(6) + str(n).rjust(6) + str(Q2).rjust(18) +\
                                  str(Q3).rjust(18) + str(kk).rjust(12) + str(ll).rjust(12) +" \n"
                count += 1
        for z,n in oct_dict_cross[edf]:
            A = z+n
            for k in it.chain(np.arange(-b_max/100,-0.449,0.05),
            np.arange(0.45,b_max/100+0.00001,0.05)):
                Q2,Q3 = beta_to_moment(k,0.0,A)
                kk,Q2,ll,Q3 = round(k,2), round(Q2,10), 0.0, 0.0
                output_content += str(z).rjust(6) + str(n).rjust(6) + str(Q2).rjust(18) +\
                                  str(Q3).rjust(18) + str(kk).rjust(12) + str(ll).rjust(12) +" \n"
                count += 1

        output_head = " Nrows \n" + "   " + str(count) + " \n" + "    Z     N          Q20                 Q30           beta2       beta3\n"
        output_content = output_head + output_content
        f.write(output_content)
        f.close()
        print (f"{edf}, Row count: {count}")

create_masstable_out()
