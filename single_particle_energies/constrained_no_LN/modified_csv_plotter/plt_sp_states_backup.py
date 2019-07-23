import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec

# For Adobe illustrator text
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def readcsv(ff,d,d1,d2,d3):
    lines = ff.readlines()
    lines.pop(0)
    for line in lines:
        ss = line.split(',')
        if "Ground" in line:
            edf,label,Z,N,A = ss[0],ss[1],int(ss[2]),int(ss[3]),int(ss[4])
            if edf not in d2: d2.append(edf)
            if label not in d1:
                d1[(label,"neutron")],d1[(label,"proton")] = [],[]
            continue
        isospin = ss[5]
        spOrbit = ss[7]
        p_h = ss[6]
        d1[(label,isospin)].append(spOrbit)
#        print (ss[8].split('"'))
#        print (ss[12].split('"'))
        nj = int(ss[8].split('"')[1])
        nl = int(ss[12].split('"')[0])
        spEnergy = round(float(ss[14]),6)
        d[(edf,label,isospin,spOrbit)] = spEnergy
        if spOrbit not in d3:
            d3[spOrbit] = (nl,nj)


def sp_orbits():
    EDF_list = ["UNEDF0","UNEDF1","UNEDF2","SV-MIN","SLY4","SKMS","SKP"]
    # keys of sp_dict are edf, element label, neutron/proton
    # values are sp energy, orbit's label, l,j quantum number
    sp_dict = {}; elem_spo = {}; edf_xtick = []; orbit_lj = {}; E3 = {}
    marker_type_list = ["o","*","s","v","<",">","^","1","+","x","D","p","H","P"]
    marker_color_list = ['b', 'g', 'r', 'c', 'm', 'y','k','b', 'g', 'r', 'c', 'm', 'y', 'k']#"C7","C8","C9","xkcd:gold","xkcd:lavender","xkcd:orange","xkcd:navy"]
    for i in EDF_list:
        f0 = open(i+"_constrained_no_LN.csv","r")
        readcsv(f0,sp_dict,elem_spo,edf_xtick,orbit_lj)
        f0.close()
        for k1 in orbit_lj:
            nl1,nj1 = orbit_lj[k1]
            for k2 in orbit_lj:
                nl2,nj2 = orbit_lj[k2]
                if abs(nj1-nj2) == 6 and abs(nl1-nl2)==3:
                    if (k1,k2) not in E3:
                        E3[(k1,k2)] = 1
    # k is list of element labels, each constitute 2 plots, neutron & proton
    # each big loop is one plot
    nuclear_shell = ["1s_(1/2)","1p_(3/2)","1p_(1/2)","1d_(5/2)","2s_(1/2)","1d_(3/2)","1f_(7/2)",
    "2p_(3/2)","1f_(5/2)","2p_(1/2)","1g_(9/2)","1g_(7/2)","2d_(5/2)","2d_(3/2)","3s_(1/2)",
    "1h_(11/2)","1h_(9/2)","2f_(7/2)","2f_(5/2)","3p_(3/2)","3p_(1/2)","1i_(13/2)","2g_(9/2)",
    "3d_(5/2)","1i_(11/2)","2g_(7/2)","4s_(1/2)","3d_(3/2)","1j_(15/2)","2h_(11/2)","1j_(13/2)",
    "3f_(7/2)","2h_(9/2)","4p_(3/2)","3f_(5/2)","1k_(17/2)"]
    for i,n in enumerate(nuclear_shell):
        n = n.split("/")
        nx = n[0].split("(")
        nuclear_shell[i] = nx[0]+nx[1]+"_"+n[1][0]
    if True:
        for k,isospin in elem_spo:
            orbits = elem_spo[(k,isospin)]
            # same orbit use same marker shape and color
            visited = set()
            for i in orbits:
                for j in orbits:
                    if (i,j) in E3 and (j,i) not in visited and (i,j) not in visited:
                        visited.add((i,j))
                        print (k,isospin,i,j)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title(f"{k} - {isospin}")
            ax.set_ylabel("Single Particle Energy(MeV)")
            ind = 0
            for spo in reversed(nuclear_shell):
                if spo in orbits:
                    s_spo = spo.split("_")
                    l_spo = s_spo[0]+"$_{"+s_spo[1]+"/"+s_spo[2]+"}$"
                    mt = marker_type_list[ind]
                    mc = marker_color_list[ind]
                    x_links = []; y_links = []
                    x0 = np.array([0.5,1.2])
                    for e in edf_xtick:
                        sp_e = sp_dict[(e,k,isospin,spo)]
                        y0 = np.array([sp_e,sp_e])
                        if e == "UNE0":
                            ax.plot(x0,y0,marker=mt,linestyle='-',
                            color=mc,linewidth=1.0,markersize=5,alpha=0.7,label=l_spo)
                        else:
                            ax.plot(x0,y0,marker=mt,linestyle='-',
                            color=mc,linewidth=1.0,markersize=5,alpha=0.7)
                        x_links.append(x0[0]);x_links.append(x0[1])
                        y_links.append(sp_e); y_links.append(sp_e)
                        x0 += 1
                    ind+=1
                    # plot links between EDFs, not really necessary
        #            ax.plot(np.array(x_links),np.array(y_links),ls="--",lw=0.5,color=mc)
            ax.set_xlim([0,9.7])
            ax.set_xticks([0.82,1.82,2.82,3.82,4.82,5.82,6.82])
            ax.set_xticklabels(edf_xtick)
            ax.tick_params(axis='x',length=0)
            ax.legend(numpoints=2)
            plt.savefig("plots/sp_orbits/"+k+"_"+isospin+".pdf",format="pdf")
            plt.clf()

def octupole_split():
    EDF_list = ["UNEDF0","UNEDF1","UNEDF2","SV-MIN","SLY4","SKMS","SKP"]
    #edf_list = ["UNE0","UNE1","UNE2","SV-min","SLy4","SKM*","SKP"]
    isotope_list = ["16O","40Ca","48Ca","56Ni","78Ni","90Zr","100Sn","132Sn","208Pb","z126_n184"]
    isospin_list = ["neutron","proton"]
    # In order of nucleon number: 184, 126, 82, 50, 28
    j3_pairs = [("1k_17_2","2h_11_2"),("1j_15_2","2g_9_2"),("1i_13_2","2f_7_2"),
    ("1h_11_2","2d_5_2"),("1g_9_2","2p_3_2")]
    # keys of sp_dict are edf, element label, neutron/proton
    # values are sp energy, orbit's label, l,j quantum number
    sp_dict = {}; elem_spo = {}; edf_xtick = []; orbit_lj = {}; j3_split = {}
    for i in EDF_list:
        f0 = open(i+"_constrained_no_LN.csv","r")
        readcsv(f0,sp_dict,elem_spo,edf_xtick,orbit_lj)
        f0.close()
    for i in edf_xtick:
        for j in isotope_list:
            for k in isospin_list:
                for a,b in j3_pairs:
                    # a,b order as listed in j3_pairs, we choose only particle states
                    # which will be listed in the front of the list, if we found the pair
                    # we break
                    if (i,j,k,a) and (i,j,k,b) in sp_dict:
                        # edf, isospin, isotope, (delta_j3 pairs)
                        j3_split[(i,k,j,(a,b))] = round(sp_dict[(i,j,k,a)]- sp_dict[(i,j,k,b)],6)
                        break
    # Regions of interest:
    # -1: 34, -2: 56, -3: 88, -4: 134, -5: 184+
    # Proton  -1: 56Ni, -2: 100Sn, -3: 132Pb, -4: 208Pb
    # Neutron -1: 56Ni, -2: 100Sn, -3: 132Sn, -4: 208Pb, -5: z126_n184
    # (N,Z): (-1,-1); (-2,-1); (-2,-2); (-3,-2); (-4,-2); (-4,-3); (-5,-3);
    # Each region should generate a plot
    region_labels = ["80Zr","100Zr","112Ba","146Ba","190Ba","224Ra","286Th"]
    # in order of above region, in order of neutron / proton for each pair
    # a pair goes into the same plot
    oct_region = [[("56Ni",(j3_pairs[-1])),("56Ni",(j3_pairs[-1]))], #80Zr,z40,n40
    [("100Sn",(j3_pairs[-2])),("56Ni",(j3_pairs[-1]))],  # 100Zr,z40,n60
    [("100Sn",(j3_pairs[-2])),("100Sn",(j3_pairs[-2]))], # 112Ba,z56,n56
    [("132Sn",(j3_pairs[-3])),("132Sn",(j3_pairs[-2]))], # 146Ba,z56,n90
    [("208Pb",(j3_pairs[-4])),("132Sn",(j3_pairs[-2]))], # 190Ba,z56,n134
    [("208Pb",(j3_pairs[-4])),("208Pb",(j3_pairs[-3]))], # 224Ra,z88,n136
    [("z126_n184",(j3_pairs[-5])),("208Pb",(j3_pairs[-3]))]] # 286Th,z90,n196
    # plot or not
    if False:
        for ind1,(a,b) in enumerate(oct_region):
            x0, y1, y2 =  [], [], []
            n_orbits1 = a[1][0].split("_")
            n_orbits2 = a[1][1].split("_")
            n_orb_latex1 = n_orbits1[0]+"$_{"+n_orbits1[1]+"/"+n_orbits1[2]+"}$"
            n_orb_latex2 = n_orbits2[0]+"$_{"+n_orbits2[1]+"/"+n_orbits2[2]+"}$"
            p_orbits1 = b[1][0].split("_")
            p_orbits2 = b[1][1].split("_")
            p_orb_latex1 = p_orbits1[0]+"$_{"+p_orbits1[1]+"/"+p_orbits1[2]+"}$"
            p_orb_latex2 = p_orbits2[0]+"$_{"+p_orbits2[1]+"/"+p_orbits2[2]+"}$"
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for ind, k in enumerate(edf_xtick):
                # a from neutron pairs, b from proton pairs
                l = (k,"neutron",)+a; m = (k,"proton",)+b
                # neutron
                y1.append(j3_split[l])
                # proton
                y2.append(j3_split[m])
                # x-axis is functionals
                x0.append([ind+1])
            x0 = np.array(x0); y1 = np.array(y1); y2 = np.array(y2)
            plt.plot(x0,y1,'b-',marker='s',label="neutron "+
            n_orb_latex1+"-"+n_orb_latex2)
            plt.plot(x0,y2,'r-',marker='s',label="proton   "+
            p_orb_latex1+"-"+p_orb_latex2)
            plt.title(region_labels[ind1]+" Region")
            plt.ylabel(r"E$_{\Delta{l,j}= 3}$ (MeV)",fontsize=14)
            ax.set_ylim([-0.5,8])
            ax.set_xticks(x0)
            ax.set_xticklabels(edf_xtick)
            ml = MultipleLocator(5)
            ax.tick_params(axis='x',length=0)
            plt.axes().yaxis.set_minor_locator(ml)
            plt.legend()
            plt.savefig("plots/splitting/"+region_labels[ind1]+
            "_region_splitting.pdf",format="pdf")
            plt.clf()

#    print(j3_split)
    print (j3_split[("UNE0","neutron","100Sn",("1h_11_2","2d_5_2"))])
    print (j3_split[("UNE0","proton","100Sn",("1h_11_2","2d_5_2"))])

    print (j3_split[("UNE0","neutron","132Sn",("1i_13_2","2f_7_2"))])
    print (j3_split[("UNE0","proton","132Sn",("1h_11_2","2d_5_2"))])


octupole_split()
