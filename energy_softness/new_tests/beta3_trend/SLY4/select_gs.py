import numpy as np
import math
import matplotlib.pyplot as plt

def Q30_b3(Q2, b2, Q3, A):
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
    #Q3 = round( b3 * A * r0**3 / c3 / 1000,6)
    b3 = round( 1000 * c3 * Q3 / A / r0**3 ,6)
    return b3

def reduce():
    BE_trend = {}
    Q20_values = {}
    lines = open("HFBTHOv300_SLY4_All_Data_20_shells_z56to56-masstable.dat",'r').readlines()[1:]
    outputFile = open("SLY4_octupole_Ba_trends_reduced.dat","w")
    outputStr = "Z".ljust(7)+"N".ljust(7)+"A".ljust(7)+"Binding_Energy_(MeV)".ljust(24)+"beta2_total".ljust(15)\
    +"beta3_total".ljust(15)+"fileID".ljust(10)+"Convergence".ljust(15)+"\n"
    for line in lines:
        ss = line.split()
        # Z      N      A
        # Binding_Energy_(MeV)    Quad_Def_Beta2_P     Quad_Def_Beta2_N     Quad_Def_Beta2_total
        # Quad_Moment_Q2_P_(fm^2)    Quad_Moment_Q2_N_(fm^2)    Quad_Moment_Q2_total_(fm^2)
        # Octupole_Moment_Q3_P_(fm^3)     Octupole_Moment_Q3_N_(fm^3)     Octupole_Moment_Q3_total_(fm^3)
        # Pairing_gap_P_(MeV)     Pairing_gap_N_(MeV)   RMS_radius_P_(fm)
        # RMS_radius_N_(fm)    RMS_radius_total_(fm)    Charge_Radius_(fm)     File_ID
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
        if convergence != "YES":
            print (Z,N,"unconverged, file# ",ss[-2])
            continue
        if N not in BE_trend:
            BE_trend[N] = []
        b3 = Q30_b3(float(Q2_T), float(beta2_T), Q3_T, N+Z)
        BE_trend[N].append( ( round(b3,6), round(float(bind_E),6) ) )
        Q20_values[N] = round(float(Q2_T),6)
        outputStr += str(Z).ljust(7)+str(N).ljust(7)+str(A).ljust(7)+str(bind_E).ljust(24)+\
        str(beta2_T).ljust(15)+str(beta3_T).ljust(15)+f_id.ljust(10)+convergence.ljust(15)+"\n"
    outputFile.write(outputStr)
    outputFile.close()

    for N in sorted(BE_trend.keys()):
        # print (sorted(BE_trend[N]))
        # print ("++")
        fig, ax = plt.subplots(figsize=(15,5))
        xx,yy = zip(*sorted(BE_trend[N]))
        xx, yy = list(xx), list(yy)
        plt.clf()
        plt.plot(xx,yy,'ro-',label=str(56+N)+"Ba, Q20="+str(Q20_values[N]))
        plt.legend()
        plt.xticks(np.arange(0,0.021,0.001))
        plt.grid(alpha=0.3)
        plt.xlabel(r"$\beta_3$")
        plt.ylabel("BE (MeV)")
        plt.title(str(56+N)+"Ba")
        plt.savefig( str(N+56)+"Ba_BE_vs_Q30.pdf",bbox_inches='tight')


reduce()
