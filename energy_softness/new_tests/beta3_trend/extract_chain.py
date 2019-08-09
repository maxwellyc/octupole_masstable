import math
import numpy as np
def isNum(x):
    try:
        float(x)
        return True
    except:
        return False

def beta3_to_moment(Q2, b2, b3, A):
    c2 = math.sqrt(5 * math.pi) / 3
    c3 = 4 * math.pi / 3  # ~=4.18879
    # Use Q2 and b2 directly from calculation to figure out value of rms radius
    # So we can derive a more precise Q3 value to constrain on.
    # Unless there's cases where input b2 is too small, division of zero occurs, we turn
    # back to r0 = 1.2 * A(1/3)
    try:
        r0 = math.sqrt(Q2 / b2 / A * 100 * c2)
    except:
        r0 = 1.2 * A**(1.0/3.0)
    Q3 = round( b3 * A * r0**3 / c3 / 1000,6)
    return Q3

def softness_inputs(edf,z_ext = 56):
    l = open(edf+"_octupole_energy_050.dat",'r').readlines()[1:]
    outstr = ""
    outfile = open("hfbtho_MASSTABLE_"+edf+".dat","w")
    count = 0
    for line in l:
        ss = line.split()
        z, n, b2, Q2= int(ss[0]), int(ss[1]), round(float(ss[-2]),6), round(float(ss[-1]),6)
        for b3 in np.arange(0,0.021,0.001):
            count += 1
            Q3 = beta3_to_moment(Q2, b2, b3, z+n)
            b3 = round(b3,6)
            outstr += str(z).rjust(6) + str(n).rjust(6) + str(Q2).rjust(18) +\
                              str(Q3).rjust(18) + str(b2).rjust(18) + str(b3).rjust(18) +" \n"
    outstr = " Nrows \n" + "   " + str(count) + " \n" + "    Z     N           " +\
    "Q20               Q30               beta2                 beta3\n" + outstr
    outfile.write(outstr)
    outfile.close()

for edf in ["SLY4"]:
    softness_inputs(edf,56)
