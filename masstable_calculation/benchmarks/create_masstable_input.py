import os
import sys
import numpy as np
import math

# beta3 = 2424.068 * Q30 / (proton_number + neutron_number)**2

def default_dripline(zz):
    """ default dripline, for a given Z, returns (N_min,N_max) in a tuple
    """
    # dictionary, key is proton, value is (min_N,max_N) tuple
    # first create min and max, then zip together.
    drip_min={2:2, 4:2, 6:2, 8:2, 10:2, 12:2, 14:2, 16:4, 18:6, 20:8, 22:12, 24:12, 26:14, 28:14, 30:18,
    32:22, 34:24, 36:26, 38:28, 40:30, 42:32, 44:34, 46:36, 48:38, 50:40, 52:48, 54:50, 56:52, 58:54, 60:56,
    62:60, 64:62, 66:66, 68:68, 70:72, 72:74, 74:76, 76:80, 78:84, 80:86, 82:90, 84:96, 86:100, 88:104, 90:108,
    92:112, 94:116, 96:122, 98:124, 100:126, 102:128, 104:134, 106:140, 108:144, 110:148, 112:152, 114:154,
    116:158, 118:162, 120:166}
    drip_max={2:6, 4:20, 6:22, 8:26, 10:38, 12:40, 14:40, 16:48, 18:56, 20:62, 22:64, 24:66, 26:68,
    28:78, 30:84,32:88, 34:92, 36:96, 38:100, 40:102, 42:110, 44:116, 46:130, 48:132, 50:132, 52:132,
    54:132, 56:138, 58:150, 60:152,62:158, 64:162, 66:188, 68:190, 70:190, 72:190, 74:190, 76:190,
    78:190, 80:190, 82:192, 84:208, 86:214, 88:220, 90:226,92:238, 94:262, 96:262, 98:262, 100:264,
    102:264, 104:264, 106:264, 108:272, 110:280, 112:286, 114:296, 116:300, 118:300, 120:304}
    drip_N = (drip_min[(zz)], drip_max[(zz)])
    return drip_N

def nuclei_range_input():
    # Lower limit of neutron number
    Z_min = input("Lower limit of Z (2 ~ 130): ")
    try:
        Z_min = int(Z_min)
        if not (2 <= Z_min <= 130):
                raise ValueError("Please enter positive integer between 2~300 !")
    except ValueError:
        raise ValueError("Please enter positive integer between 2 ~ 130 !")
    # Upper limit of proton number
    Z_max = input("Upper limit of Z (2 ~ 130): ")
    try:
        Z_max = int(Z_max)
        if not (2 <= Z_max <= 130):
                raise ValueError("Please enter positive integer between 2 ~ 130 !")
    except ValueError:
        raise ValueError("Please enter positive integer between 2 ~ 120 !")
    # Lower limit of neutron number
    N_min = input("Lower limit of N (2 ~ 300, hit return for default dripline): ")
    if len(N_min) > 0:
        try:
            N_min = int(N_min)
            if not (2 <= N_min <= 300):
                raise ValueError("Please enter positive integer between 2 ~ 300 !")
        except ValueError:
            raise ValueError("Please enter positive integer between 2 ~ 300 !")
    elif len(N_min) == 0:
        N_min = -1
    # Upper limit of neutron number
    N_max = input("Upper limit of N (2 ~ 300, hit return for default dripline): ")
    if len(N_max) > 0:
        try:
            N_max = int(N_max)
            if not (2 <= N_max <= 300):
                raise ValueError("Please enter positive integer between 2 ~ 300 !")
        except ValueError:
            raise ValueError("Please enter positive integer between 2 ~ 300 !")
    elif len(N_max) == 0:
        N_max = -1
    
    return Z_min,Z_max,N_min,N_max
    # Need a list for each Z the lower and upper limit for N

def beta_range_input():
    # Set beta2 range and stepsize
    beta2_min = input("Lower limit of beta2 (-0.5 ~ 0.5,  hit return for default -0.35): ")
    if len(beta2_min) > 0:
        try:
            beta2_min = float(beta2_min)
            if not (-0.5 <= beta2_min <= 0.5):
                raise ValueError("Please enter number between -0.5 ~ 0.5 !")
        except ValueError:
            raise ValueError("Please enter positive number!")
    elif len(beta2_min) == 0:
        beta2_min = -0.35
    beta2_max = input("Upper limit of beta2 (-0.5 ~ 0.5,  hit return for default 0.35): ")
    if len(beta2_max) > 0:
        try:
            beta2_max = float(beta2_max)
            if not (-0.5 <= beta2_max <= 0.5):
                raise ValueError("Please enter number between -0.5 ~ 0.5 !")
            if beta2_min > beta2_max:
                raise  ValueError("beta2_max is smaller than beta2_min !")
        except ValueError:
            raise ValueError("Please enter positive number!")
    elif len(beta2_max) == 0:
        beta2_max = 0.35

    beta2_step = input("Step size of beta2 (hit return for default 0.05): ")
    if len(beta2_step) > 0:
        try:
            beta2_step = float(beta2_step)
        except ValueError:
            raise ValueError("Please enter positive number!")
    elif len(beta2_step) == 0:
        beta2_step = 0.05
    
    # Set beta3 range and stepsize
    beta3_min = input("Lower limit of beta3 (0 ~ 0.4,  hit return for default 0.0): ")
    if len(beta3_min) > 0:
        try:
            beta3_min = float(beta3_min)
            if not (0 <= beta3_min <= 0.4):
                raise ValueError("Please enter positive number between 0 ~ 0.4 !")
        except ValueError:
            raise ValueError("Please enter positive number between 0 ~ 0.4 !")
    elif len(beta3_min) == 0:
        beta3_min = 0.0

    beta3_max = input("Upper limit of beta3 (0 ~ 0.4, hit return for default 0.4): ")
    if len(beta3_max) > 0:
        try:
            beta3_max = float(beta3_max)
            if not (0 <= beta3_max <= 0.4):
                raise ValueError("Please enter positive number between 0 ~ 0.4 !")
            elif beta3_min > beta3_max:
                raise  ValueError("beta3_max is smaller than beta3_min !")
        except ValueError:
            raise ValueError("Please enter positive number between 0 ~ 0.4 !")
    elif len(beta3_max) == 0:
        beta3_max = 0.4

    beta3_step = input("Step size of beta3 (hit return for default 0.1): ")
    if len(beta3_step) > 0:
        try:
            beta3_step = float(beta3_step)
        except ValueError:
            raise ValueError("Please enter positive number!")
    elif len(beta3_step) == 0:
        beta3_step = 0.1
    return beta2_min, beta2_max, beta2_step,beta3_min, beta3_max, beta3_step

def beta_to_moment(beta2, beta3, A):
    c2 = math.sqrt(5 * math.pi) / 3
    c3 = 4 * math.pi / 3  # ~=4.18879
    r0 = 1.2 * A**(1.0/3.0)
    Q2 = beta2 * A * r0**2 / c2 / 100
    Q3 = beta3 * A * r0**3 / c3 / 1000
    return Q2,Q3



def create_masstable_out(fname="hfbtho_MASSTABLE.dat"):
    f = open(fname,"w")
    output_content = ""
    Z_min,Z_max,N_min,N_max = nuclei_range_input()
    beta2_min, beta2_max, beta2_step,beta3_min, beta3_max, beta3_step = beta_range_input()
    count = 0
    print("Z\t N_min N_max")
    for i in range(Z_min,Z_max+1,2):
        if N_min == -1:
            n0 = default_dripline(i)[0]
        elif N_min > 0:
            n0 = N_min
        if N_max == -1:
            n1 = default_dripline(i)[1]
        elif N_max > 0:
            n1 = N_max
        print (i,"\t",n0,n1)
        for j in range(n0,n1+1,2):
            A = i+j
            for l in np.arange(beta3_min, beta3_max+0.000001,beta3_step):
                for k in np.arange(beta2_min, beta2_max+0.000001,beta2_step):
                    Q2,Q3 = beta_to_moment(k,l,A)
                    kk = round(k,2); ll = round(l,3); Q2 = round(Q2,10); Q3 = round(Q3,10)
                    output_content += str(i).rjust(6) + str(j).rjust(6) + str(Q2).rjust(18) +\
                                      str(Q3).rjust(18) + str(kk).rjust(12) + str(ll).rjust(12) +" \n"
                    count += 1
    output_head = " Nrows \n" + "   " + str(count) + " \n" + "    Z     N          Q20                 Q30           beta2       beta3\n"
    output_content = output_head + output_content
    f.write(output_content)
    f.close()
    print (f"Row count: {count}")

create_masstable_out()


