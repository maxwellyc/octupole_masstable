import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
# For Adobe illustrator text
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


def plot_softness():
#
    edf = ["UNEDF0", "UNEDF1","UNEDF2","SV-MIN","SLY4"]#, "SKMS", "SKP"]
    colors  = ['r','b','k','y','m','c','g']
    markers = ['s','*','o','^','d','+','x']


    f,l,E_oct,f_gs,gs_def,l_gs = {},{},{},{},{},{}
    for e in edf:
        f[e] = open("octupole_energy/"+e+"_octupole_energy_050.dat","r")
        f_gs[e] = open("octupole_gs_reduced/"+e+"_reduced_beta2_3_050.dat","r")
        l[e] = f[e].readlines()
        l_gs[e] = f_gs[e].readlines()
        count = 0; E_oct[e] = {}
        for line in l[e]:
            if not count: count+=1;continue
            ss=line.split()
            Z,N,E,b2,b3 = int(ss[0]),int(ss[1]),float(ss[3]),float(ss[7]),float(ss[6])
            E_oct[e][(Z,N)] = (E,b2,b3)

        count = 0; gs_def[e] = {}
        for line in l_gs[e]:
            if not count: count+=1;continue
            ss=line.split()
            Z,N,b2,b3 = int(ss[0]),int(ss[1]),float(ss[4]),float(ss[5])
            gs_def[e][(Z,N)] = (b2,b3)
        f[e].close()
        f_gs[e].close()

    # plot isotope, 3 plots for one chain. Each edf gets one chain
    z_plot = 94; plt_title = "Pu"
    plt_list = [(z_plot,plt_title)]
    #[(56,"Ba"),(86,"Rn"),(88,"Ra"),(90,"Th"),(92,"U"),
    #(94,"Pu"),(96,"Cm"),(98,"Cf"),(100,"Fm"),(102,"No"),(104,"Rf"),(106,"Sg"),
    #(108,"Hs"),(110,"Ds"),(112,"Cn"),(114,"Fl"),(116,"Lv"),(118,"Og"),
    #(120,"z120"),(58,"Ce"),(54,"Xe"),(80,"Hg")]
    x_range = {}
    for z_plot,plt_title in plt_list:
        x_range[z_plot] = (1000,0)
        for e in edf:
            for z,n in E_oct[e].keys():
                if z == z_plot:
                    x_range[z_plot] = (min(x_range[z_plot][0],n),max(x_range[z_plot][1],n))

    for z_plot,plt_title in plt_list:
        fig, axes = plt.subplots(3,1, sharex = True, sharey = False)#, figsize = (10,30))
        axes[0].set_title(plt_title)
        for j,e in enumerate(edf):
            x_n, y_Eoct,y_b2, y_b3 = [] , [], [], []

            # one-time use
            n1, n2 = 126,142
            for N in range(n1,n2+1,2):
            # original code:
            #for N in range(x_range[z_plot][0]-4,x_range[z_plot][1]+1,2):

                x_n.append(N+z_plot)
                if (z_plot,N) in E_oct[e]:
                    y_Eoct.append(E_oct[e][(z_plot,N)][0])
                    y_b2.append(E_oct[e][(z_plot,N)][1])
                    y_b3.append(E_oct[e][(z_plot,N)][2])
                else:
                    y_Eoct.append(0)
                    y_b2.append(gs_def[e][(z_plot,N)][0])
                    y_b3.append(gs_def[e][(z_plot,N)][1])
            axes[0].plot(x_n,y_b2,c=colors[j],linestyle='-',linewidth=0.5,marker=markers[j],
            markersize=2,label=e,alpha=0.7)
            axes[1].plot(x_n,y_b3,c=colors[j],linestyle='-',linewidth=0.5,marker=markers[j],
            markersize=2,label=e,alpha=0.7)
            axes[2].plot(x_n,y_Eoct,c=colors[j],linestyle='-',linewidth=0.5,marker=markers[j],
            markersize=2,label=e,alpha=0.7)
        axes[0].set_ylabel(r"$\beta_2$")
        axes[1].set_ylabel(r"$\beta_3$")
        axes[2].set_ylabel(r"$E_{oct} (MeV)$")
        axes[2].set_xlabel("A")

        # one-time use
        plt.xticks(range(n1+z_plot,n2+z_plot+1,2),fontsize=12)
        # original code:
        #axes[2].xaxis.set_minor_locator(AutoMinorLocator(3))
        #plt.xticks(range(x_range[z_plot][0],x_range[z_plot][1]+1,6),fontsize=6)

        plt.legend(loc=4,fontsize=8)#prop={'size':6}
        plt.savefig("plots/Z"+str(z_plot)+"_"+plt_title+".pdf",format="pdf")
        plt.clf()
plot_softness()
