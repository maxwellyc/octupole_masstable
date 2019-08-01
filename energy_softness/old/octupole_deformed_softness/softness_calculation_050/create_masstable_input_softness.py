import os
import sys
import numpy as np
import math

def beta_to_moment(beta2, beta3, A):
    c2 = math.sqrt(5 * math.pi) / 3
    c3 = 4 * math.pi / 3  # ~=4.18879
    r0 = 1.2 * A**(1.0/3.0)
    Q2 = beta2 * A * r0**2 / c2 / 100
    Q3 = beta3 * A * r0**3 / c3 / 1000
    return Q2,Q3



def create_masstable_out(edf):
    b3_in = 0.03
    f_input = open("beta2_min_inputs/"+edf+"_octupole_energy_050.dat","r")
    l = f_input.readlines()
    f = open("hfbtho_MASSTABLE_"+edf+".dat","w")
    output_content,count = "",0
    for ind,line in enumerate(l):
        if ind:
            ss = line.split()
            Z,N,b2 = int(ss[0]),int(ss[1]),float(ss[-2])
            Q2,Q3 = beta_to_moment(b2,b3_in,Z+N)
            Q2 = round(Q2,10);Q3 = round(Q3,10)
            output_content += str(Z).rjust(6) + str(N).rjust(6) + str(Q2).rjust(18) +\
                              str(Q3).rjust(18) + str(b2).rjust(12) + str(b3_in).rjust(12) +" \n"
            count += 1
    output_head = " Nrows \n" + "   " + str(count) + " \n" + "    Z     N          Q20                 Q30           beta2       beta3\n"
    output_content = output_head + output_content
    f.write(output_content)
    f.close()
    print (f"{edf}, Row count: {count}")

EDF = ["UNEDF0","UNEDF1","UNEDF2","SLY4","SV-MIN","SKMS","SKP"]
for e in EDF:
    create_masstable_out(e)
