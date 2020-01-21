#Takes in a data file eliminates all lines with energy values above the ground state
#and gives as output a file whose energies are only ground state values.  These files
#are then used as input by erik.py
import matplotlib.pyplot as plt
import numpy as np

def ctFile(b_min=1,b_max=50):

 if b_max != 0:
     f0 = open("Data/UNEDF0_reduced_beta2_3_0"+str(b_max)+".dat","r")
     f1 = open("Data/UNEDF1_reduced_beta2_3_0"+str(b_max)+".dat","r")
     f2 = open("Data/UNEDF2_reduced_beta2_3_0"+str(b_max)+".dat","r")
     f3 = open("Data/SV-MIN_reduced_beta2_3_0"+str(b_max)+".dat","r")
     f4 = open("Data/SLY4_reduced_beta2_3_0"+str(b_max)+".dat","r")
     f5 = open("Data/SKMS_reduced_beta2_3_0"+str(b_max)+".dat","r")
     f6 = open("Data/SKP_reduced_beta2_3_0"+str(b_max)+".dat","r")
 else:
     f0 = open("Data/UNEDF0_reduced_no_filter.dat","r")
     f1 = open("Data/UNEDF1_reduced_no_filter.dat","r")
     f2 = open("Data/UNEDF2_reduced_no_filter.dat","r")
     f3 = open("Data/SV-MIN_reduced_no_filter.dat","r")
     f4 = open("Data/SLY4_reduced_no_filter.dat","r")
     f5 = open("Data/SKMS_reduced_no_filter.dat","r")
     f6 = open("Data/SKP_reduced_no_filter.dat","r")
 
 lines0 = f0.readlines()
 lines1 = f1.readlines()
 lines2 = f2.readlines()
 lines3 = f3.readlines()
 lines4 = f4.readlines()
 lines5 = f5.readlines()
 lines6 = f6.readlines()
 
 linew = '';ctf = []
 ct0,ct1,ct2,ct3,ct4,ct5,ct6 = {},{},{},{},{},{},{}
 
# UNEDF0
 for line in lines0:

    ss = line.split()
    try:                                    
     N=int(float(ss[1])+0.00001)      #Number of Neutrons
     Z=int(float(ss[0])+0.00001)      #Number of Protons
     beta3=float(float(ss[5]))               #beta3 of a particular nucleus

     if(abs(beta3) >= b_min/100) :
       ct0[(N,Z)] = 1
     else:
       ct0[(N,Z)] = 0
       continue
    except (ValueError, IndexError):
       continue                      #N,Z, or, Efn are not numbers
 f0.close()


# UNEDF1
 for line in lines1:

    ss = line.split()
    try:                                    
     N=int(float(ss[1])+0.00001)      #Number of Neutrons
     Z=int(float(ss[0])+0.00001)      #Number of Protons
     beta3=float(float(ss[5]))               #beta3 of a particular nucleus

     if(abs(beta3) >= b_min/100) :
       ct1[(N,Z)] = 1
     else:
       ct1[(N,Z)] = 0
       continue
    except (ValueError, IndexError):
       continue                      #N,Z, or, Efn are not numbers
 f1.close()


# UNEDF2
 for line in lines2:

    ss = line.split()
    try:                                    
     N=int(float(ss[1])+0.00001)      #Number of Neutrons
     Z=int(float(ss[0])+0.00001)      #Number of Protons
     beta3=float(float(ss[5]))               #beta3 of a particular nucleus

     if(abs(beta3) >= b_min/100) :      #Looks for the smallest Efn
       ct2[(N,Z)] = 1
     else:
       ct2[(N,Z)] = 0
       continue
    except (ValueError, IndexError):
       continue                      #N,Z, or, Efn are not numbers
 f2.close()

# SV-min
 for line in lines3:

    ss = line.split()
    try:                                    
     N=int(float(ss[1])+0.00001)      #Number of Neutrons
     Z=int(float(ss[0])+0.00001)      #Number of Protons
     beta3=float(float(ss[5]))               #beta3 of a particular nucleus

     if(abs(beta3) >= b_min/100) :
       ct3[(N,Z)] = 1
     else:
       ct3[(N,Z)] = 0
       continue
    except (ValueError, IndexError):
       continue                      #N,Z, or, Efn are not numbers
 f3.close()

# SLy4
 for line in lines4:

    ss = line.split()
    try:                                    
     N=int(float(ss[1])+0.00001)      #Number of Neutrons
     Z=int(float(ss[0])+0.00001)      #Number of Protons
     beta3=float(float(ss[5]))               #beta3 of a particular nucleus

     if(abs(beta3) >= b_min/100) :
       ct4[(N,Z)] = 1
     else:
       ct4[(N,Z)] = 0
       continue
    except (ValueError, IndexError):
       continue                      
 f4.close()

# SKM*
 for line in lines5:

    ss = line.split()
    try:                                    
     N=int(float(ss[1])+0.00001)      #Number of Neutrons
     Z=int(float(ss[0])+0.00001)      #Number of Protons
     beta3=float(float(ss[5]))               #beta3 of a particular nucleus

     if(abs(beta3) >= b_min/100) :
       ct5[(N,Z)] = 1
     else:
       ct5[(N,Z)] = 0
       continue
    except (ValueError, IndexError):
       continue                     
 f5.close()

# SKP
 for line in lines6:

    ss = line.split()
    try:
     N=int(float(ss[1])+0.00001)      #Number of Neutrons
     Z=int(float(ss[0])+0.00001)      #Number of Protons
     beta3=float(float(ss[5]))               #beta3 of a particular nucleus

     if(abs(beta3) >= b_min/100) :
       ct6[(N,Z)] = 1
     else:
       ct6[(N,Z)] = 0
       continue
    except (ValueError, IndexError):
       continue
 f6.close()

 target = open("beta3_00"+str(b_min)+"-0"+str(b_max)+"_count.dat", "w")
 title = 'Z'.rjust(3)+'    '+'N'.rjust(3)+'    '+'UNEDF0'.rjust(6)+'    '+'UNEDF1'.rjust(6)+'    '+'UNEDF2'.rjust(6)+'    '+'SV_MIN'.rjust(6)+'    '+'SLY4'.rjust(6)+'    '+'SKMS'.rjust(6)+'    '+'SKP'.rjust(6)+'    '+'Total'.rjust(6)
 target.write(title)
 #print title
 target.write("\n")
 tt = 0
 linect = 0
 
 n1,n2,n3,n4,n5,n6,n7 = [],[],[],[],[],[],[]
 z1,z2,z3,z4,z5,z6,z7 = [],[],[],[],[],[],[]
 
 for Z in range(0,122):
   
   for N in range (0,302):
     if (N,Z) in ct0 or (N,Z) in ct1 or (N,Z) in ct2 or (N,Z) in ct3 or (N,Z) in ct4 or (N,Z) in ct5 or (N,Z) in ct6:
       tt = 0
       linew = ''
       linect += 1
       linew = str(Z).rjust(3)+'    '+str(N).rjust(3)+'    '
       #count UNEDF0
       if (N,Z) in ct0:
          linew += str(ct0[(N,Z)]).rjust(6) + '    '
          tt += ct0[(N,Z)]
       else:
          linew += '0'.rjust(6) + '    '
       #count UNEDF1
       if (N,Z) in ct1:
          linew += str(ct1[(N,Z)]).rjust(6) + '    '
          tt += ct1[(N,Z)]
       else:
          linew += '0'.rjust(6) + '    '
       #count UNEDF2
       if (N,Z) in ct2:
          linew += str(ct2[(N,Z)]).rjust(6) + '    '
          tt += ct2[(N,Z)]
       else:
          linew += '0'.rjust(6) + '    '
       #count SKMS
       if (N,Z) in ct3:
          linew += str(ct3[(N,Z)]).rjust(6) + '    '
          tt += ct3[(N,Z)]
       else:
          linew += '0'.rjust(6) + '    '
       #count SLY4
       if (N,Z) in ct4:
          linew += str(ct4[(N,Z)]).rjust(6) + '    '
          tt += ct4[(N,Z)]
       else:
          linew += '0'.rjust(6) + '    '
       #count SV-MIN
       if (N,Z) in ct5:
          linew += str(ct5[(N,Z)]).rjust(6) + '    '
          tt += ct5[(N,Z)]
       else:
          linew += '0'.rjust(6) + '    '

       if (N,Z) in ct6:
          linew += str(ct6[(N,Z)]).rjust(6) + '    '
          tt += ct6[(N,Z)]
       else:
          linew += '0'.rjust(6) + '    '
       
       if (tt!=0):
         #count TOTAL
         linew += str(tt).rjust(6)
         #print linew+'\t'+str(int(linect))
         target.write(linew)
         target.write('\n')
       if (tt == 1):
         n1.append(N)
         z1.append(Z)
       
       if (tt == 2):
         n2.append(N)
         z2.append(Z)

       if (tt == 3):
         n3.append(N)
         z3.append(Z)
         
       if (tt == 4):
         n4.append(N)
         z4.append(Z)

       if (tt == 5):
         n5.append(N)
         z5.append(Z)
         
       if (tt == 6):
         n6.append(N)
         z6.append(Z)

       if (tt == 7):
         n7.append(N)
         z7.append(Z)
       
       
     else:
       tt = 0
       continue
 fig, ax = plt.subplots()
 ax.plot(n1,z1,'o',c='0.9',label = 1, ms = 1.5)
 ax.plot(n2,z2,'o',c='0.8',label = 2, ms = 1.5)
 ax.plot(n3,z3,'o',c='0.7',label = 3, ms = 1.5)
 ax.plot(n4,z4,'o',c='0.6',label = 4, ms = 1.5)
 ax.plot(n5,z5,'o',c='0.5',label = 5, ms = 1.5)
 ax.plot(n6,z6,'o',c='0.2',label = 6, ms = 1.5)
 ax.plot(n7,z7,'D',c='b',label = 7, ms = 1)
 ax.set_aspect('equal')
 plt.xlabel('Neutron Number')
 plt.ylabel('Proton Number')
 plt.title("Octupole "+str(b_max/100)+r"$\geq\mid{\beta_3}\mid\geq$0.01 Counts")#+"\n"+"UNEDF0~2,SV-min,SLy4,SkM*,SkP")
 plt.xlim([0,302])
 plt.ylim([0,122])
 plt.xticks(np.arange(0, 302,20))
 plt.yticks(np.arange(20, 122,20))
 plt.grid(alpha=0.5,linestyle=':',linewidth=0.5)
 plt.legend(numpoints = 1,loc = 4)
 plt.savefig("mult_beta3_0"+str(b_max)+".pdf",format='pdf')

# count0=0
#
# for k in ct6:
#     if ct6[k]:
#        count0+=1
# print (count0)

ctFile(b_min=1,b_max=75)




 


 



