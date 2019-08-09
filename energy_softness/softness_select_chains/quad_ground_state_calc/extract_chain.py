import math
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

#print (beta3_to_moment(7.421899, 0.166117, 0.01, 144 ))

def extract(edf, z_ext = 56):
    outfile = open("s2_inputs/"+edf + "_z" + str(z_ext) + "_chain.dat",'w')
    # octupole nuclei
    l1 = open("s1_inputs/"+edf+"_oct.dat",'r').readlines()
    # no octupole nuclei, the files' format are different
    l2 = open("s1_inputs/"+edf+"_no_oct.dat",'r').readlines()
    outstr = "Z      N      quad_gs_beta2\n"
    temp = {}
    for line in l1:
        ss = line.split()
        if not ss or not isNum(ss[0]): continue
        if int(ss[0]) == z_ext:
            # neutron number vs Z,N,beta2,Q20 will sort on neutron numbers later
            temp[int(ss[1])] = ss[0] + "      "+ss[1] + "      " + ss[-2] + "      " + ss[-1] + "\n"
    for line in l2:
        ss = line.split()
        if not ss or not isNum(ss[0]): continue
        if int(ss[0]) == z_ext:
            if int(ss[1]) in temp: print (edf, "repeat", int(ss[1]))
            # neutron number vs Z,N,beta2,Q20 will sort on neutron numbers later
            temp[int(ss[1])] = ss[0] + "      " + ss[1] + "      " + ss[4] + "      " + ss[5] + "\n"

    for N in sorted(temp.keys()):
        outstr += temp[N]
    outfile.write(outstr)
    outfile.close()

def softness_inputs(edf,z_ext = 56, b3 = 0.01):
    l = open("s2_inputs/" + edf + "_z" + str(z_ext) + "_chain.dat",'r').readlines()[1:]
    outstr = " Nrows \n" + "   " + str(len(l)) + " \n" + "    Z     N           " +\
    "Q20               Q30               beta2                 beta3\n"
    outfile = open("s3_inputs/hfbtho_MASSTABLE_"+edf+".dat","w")
    for line in l:
        ss = line.split()
        z, n, b2, Q2= int(ss[0]), int(ss[1]), round(float(ss[2]),6), round(float(ss[3]),6)
        Q3 = beta3_to_moment(Q2, b2, b3, z+n)
        outstr += str(z).rjust(6) + str(n).rjust(6) + str(Q2).rjust(18) +\
                          str(Q3).rjust(18) + str(b2).rjust(18) + str(b3).rjust(18) +" \n"
    outfile.write(outstr)
    outfile.close()

for edf in ["SLY4","SV-MIN","UNEDF0","UNEDF1","UNEDF2"]:
    extract(edf,56)
    softness_inputs(edf,56,0.01)
