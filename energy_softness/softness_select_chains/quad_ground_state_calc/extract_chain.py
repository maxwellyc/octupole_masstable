import math
def isNum(x):
    try:
        float(x)
        return True
    except:
        return False

def beta_to_moment(Q2, b2, beta3, A):
    c2 = math.sqrt(5 * math.pi) / 3
    c3 = 4 * math.pi / 3  # ~=4.18879
#    r0 = 1.2 * A**(1.0/3.0)
    # Use Q2 and b2 directly from calculation to figure out value of rms radius
    # So we can derive a more precise Q3 value to constrain on.
    r0 = math.sqrt(Q2 / b2 / A * 100 * c2)
#    Q2 = round(beta2 * A * r0**2 / c2 / 100,6)
    b2_out = round(100* Q2 * c2 / A / (r0**2),3)
    Q3 = round(beta3 * A * r0**3 / c3 / 1000,6)
    return beta2_out,Q3

print (beta_to_moment(7.421899, , 0.01, 144 ))

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
            # neutron number vs line. we can sort on neutron numbers later
            temp[int(ss[1])] = ss[0]+"      "+ss[1]+"      "+ss[-1] + "\n"
    for line in l2:
        ss = line.split()
        if not ss or not isNum(ss[0]): continue
        if int(ss[0]) == z_ext:
            if int(ss[1]) in temp: print (edf, "repeat", int(ss[1]))
            # neutron number vs line. we can sort on neutron numbers later
            temp[int(ss[1])] = ss[0]+"      "+ss[1]+"      "+ss[5] + "\n"

    for N in sorted(temp.keys()):
        outstr += temp[N]
    outfile.write(outstr)
    outfile.close()

def softness_inputs(edf,z_ext = 56):
    l = open("s2_inputs/" + edf + "_z" + str(z_ext) + "_chain.dat",'r').readlines()[1:]
    outstr = " Nrows \n" + "   " + str(len(l)) + " \n" + "    Z     N           " +\
    "Q20               Q30               beta2                 beta3\n"
    outfile = open("s3_inputs/hfbtho_MASSTABLE_"+edf+".dat","w")
    for line in l:
        ss = line.split()
        z, n, Q2, b3 = int(ss[0]), int(ss[1]), round(float(ss[2]),6), 0.001
        b2, Q3 = beta_to_moment(Q2, b3, z+n)
        outstr += str(z).rjust(6) + str(n).rjust(6) + str(Q2).rjust(18) +\
                          str(Q3).rjust(18) + str(b2).rjust(18) + str(b3).rjust(18) +" \n"
    outfile.write(outstr)
    outfile.close()

for edf in ["SLY4","SV-MIN","UNEDF0","UNEDF1","UNEDF2"]:
    extract(edf,56)
    softness_inputs(edf,56)
