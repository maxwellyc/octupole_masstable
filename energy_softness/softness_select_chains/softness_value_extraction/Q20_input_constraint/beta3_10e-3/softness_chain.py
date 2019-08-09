# read-in softness base, there'll be two parts
# baseline can be achieved for octupole deformed and non-oct nuclei with:
# 1. octupole nuclei's quad. minimum has been calculated in the octupole energy
# series calculations.
# 2. non-oct quad. are the recently calculated ones that are used to initialize
# softness calculations.

# need to make sure beta2 are exactly the same between the two subtracted states
# so we'll need two columns to compare before and after beta2 values.

# We'll need 3 input files. beta3 = 0.01 file; baseline file for oct and non-oct.
# Also we need to find out non-converged states and see what can be done there.

# Let's first gather the latter two files while waiting for beta3 = 0.01 calculations
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def extract_t1(lines):
    dct = {}
    for line in lines:
        ss = line.split()
        Z,N,BE,b2 = int(ss[0]),int(ss[1]),float(ss[3]),float(ss[4])
        dct[(Z,N)] = (BE,b2)
    return dct

def extract_t2(lines):
    dct = {}
    for line in lines:
        ss = line.split()
        Z,N,BE,b2 = int(ss[0]),int(ss[1]),float(ss[5]),float(ss[8])
        dct[(Z,N)] = (BE,b2)
    return dct

def softness_plot(edf,zz=56):
    plt.clf()
    # beta3 = 0.01 binding energy
    l0 = open("data0/"+edf+"_small_beta3_Ba_chain_reduced.dat",'r').readlines()[1:]
    # quadrupole ground state binding energy for non-octupole deformed nuclei
    l1 = open("data1/"+edf+"_no_oct_quad_gs_reduced.dat",'r').readlines()[1:]
    # quadrupole ground state binding energy for octupole deformed nuclei
    l2 = open("data2/"+edf+"_octupole_energy_050.dat",'r').readlines()[1:]

    d0 = extract_t1(l0)
    d1 = extract_t1(l1)
    d2 = extract_t2(l2)
    soft = {}
    for Z,N in d1:
        if Z == zz and (Z,N) in d0:
            soft[(Z,N)] = -d1[(Z,N)][0] + d0[(Z,N)][0]
            #print (Z,N, d1[(Z,N)][0],  d0[(Z,N)][0], "no-oct")
    for Z,N in d2:
        if Z == zz and (Z,N) in d0:
            soft[(Z,N)] = -d2[(Z,N)][0] + d0[(Z,N)][0]
            #print (Z,N, d2[(Z,N)][0],  d0[(Z,N)][0], "oct")
    xx, yy = [], []
    xx_oct, yy_oct = [], []
    for Z,N in sorted(soft.keys()):
        # dripline cut
        if 58 <= N <= 126:
            xx.append(Z+N)
            yy.append(soft[(Z,N)])
            if (Z,N) in d2:
                xx_oct.append(Z+N)
                yy_oct.append(soft[(Z,N)])

    print (edf,xx_oct,yy_oct)
    plt.xticks(range(114,187,6))
    plt.ylim([-0.03,0.11])
    plt.grid(axis='x',alpha=0.5)
    plt.axhline(alpha=0.5,lw=0.5,ls=':')
    plt.plot(xx,yy,'bo-',label='not octupole-deformed',ms=5,zorder=1)
    plt.scatter(xx_oct,yy_oct,c='r',s=30,zorder=2,label='octupole-deformed')
    plt.legend()
    plt.title("Ba chain "+edf+"  beta3 = 0.001")
    plt.xlabel("A")
    plt.ylabel(r"$BE_{\beta3=0.001} - BE_{\beta3=0}$"+" (MeV)")
    plt.savefig("softness_plots/"+edf+"_softness.pdf",bbox_inches='tight')


for edf in ["SLY4","SV-MIN","UNEDF0","UNEDF1","UNEDF2"]:
    softness_plot(edf,56)
