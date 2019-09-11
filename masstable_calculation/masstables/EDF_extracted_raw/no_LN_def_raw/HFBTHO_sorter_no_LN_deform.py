#This code organizes the individual output files for HFBTHOv300 into two separate files:
# 1.  HFBTHOv300_"functional-name"_All_Data.dat (containing the data from every constrained calculation)
# 2.  HFBTHOv300_"functional-name"_Ground_State_Data.dat (containing the ground state data for each nucleus from all of the constrained calculations)
#
#The input files are assumed to have the following form:  thoout_000_000001.dat
#where "000" is the Team ID and "000001" is the File ID
#
#To run the code, type the following into the terminal:
#
#python SortHFBTHOv300.py "Name_of_functional" "Lowest_Team_ID_Number" "Highest_Team_ID_Number" "Lowest_File_ID_Number" "Highest_File_ID_Number"

import math
import os
import sys
import decimal
import re  #Real expressions
import subprocess  #Shell commands

def ElementName(Z):
  #=======================
  #Labels for each element
  #=======================
  if (Z == 1):
    elementlabel = 'H'
  elif (Z == 2):
    elementlabel = 'He'
  elif (Z == 3):
    elementlabel = 'Li'
  elif (Z == 4):
    elementlabel = 'Be'
  elif (Z == 5):
    elementlabel = 'B'
  elif (Z == 6):
    elementlabel = 'C'
  elif (Z == 7):
    elementlabel = 'N'
  elif (Z == 8):
    elementlabel = 'O'
  elif (Z == 9):
    elementlabel = 'F'
  elif (Z == 10):
    elementlabel = 'Ne'
  elif (Z == 11):
    elementlabel = 'Na'
  elif (Z == 12):
    elementlabel = 'Mg'
  elif (Z == 13):
    elementlabel = 'Al'
  elif (Z == 14):
    elementlabel = 'Si'
  elif (Z == 15):
    elementlabel = 'P'
  elif (Z == 16):
    elementlabel = 'S'
  elif (Z == 17):
    elementlabel = 'Cl'
  elif (Z == 18):
    elementlabel = 'Ar'
  elif (Z == 19):
    elementlabel = 'K'
  elif (Z == 20):
    elementlabel = 'Ca'
  elif (Z == 21):
    elementlabel = 'Sc'
  elif (Z == 22):
    elementlabel = 'Ti'
  elif (Z == 23):
    elementlabel = 'V'
  elif (Z == 24):
    elementlabel = 'Cr'
  elif (Z == 25):
    elementlabel = 'Mn'
  elif (Z == 26):
    elementlabel = 'Fe'
  elif (Z == 27):
    elementlabel = 'Co'
  elif (Z == 28):
    elementlabel = 'Ni'
  elif (Z == 29):
    elementlabel = 'Cu'
  elif (Z == 30):
    elementlabel = 'Zn'
  elif (Z == 31):
    elementlabel = 'Ga'
  elif (Z == 32):
    elementlabel = 'Ge'
  elif (Z == 33):
    elementlabel = 'As'
  elif (Z == 34):
    elementlabel = 'Se'
  elif (Z == 35):
    elementlabel = 'Br'
  elif (Z == 36):
    elementlabel = 'Kr'
  elif (Z == 37):
    elementlabel = 'Rb'
  elif (Z == 38):
    elementlabel = 'Sr'
  elif (Z == 39):
    elementlabel = 'Y'
  elif (Z == 40):
    elementlabel = 'Zr'
  elif (Z == 41):
    elementlabel = 'Nb'
  elif (Z == 42):
    elementlabel = 'Mo'
  elif (Z == 43):
    elementlabel = 'Tc'
  elif (Z == 44):
    elementlabel = 'Ru'
  elif (Z == 45):
    elementlabel = 'Rh'
  elif (Z == 46):
    elementlabel = 'Pd'
  elif (Z == 47):
    elementlabel = 'Ag'
  elif (Z == 48):
    elementlabel = 'Cd'
  elif (Z == 49):
    elementlabel = 'In'
  elif (Z == 50):
    elementlabel = 'Sn'
  elif (Z == 51):
    elementlabel = 'Sb'
  elif (Z == 52):
    elementlabel = 'Te'
  elif (Z == 53):
    elementlabel = 'I'
  elif (Z == 54):
    elementlabel = 'Xe'
  elif (Z == 55):
    elementlabel = 'Cs'
  elif (Z == 56):
    elementlabel = 'Ba'
  elif (Z == 57):
    elementlabel = 'La'
  elif (Z == 58):
    elementlabel = 'Ce'
  elif (Z == 59):
    elementlabel = 'Pr'
  elif (Z == 60):
    elementlabel = 'Nd'
  elif (Z == 61):
    elementlabel = 'Pm'
  elif (Z == 62):
    elementlabel = 'Sm'
  elif (Z == 63):
    elementlabel = 'Eu'
  elif (Z == 64):
    elementlabel = 'Gd'
  elif (Z == 65):
    elementlabel = 'Tb'
  elif (Z == 66):
    elementlabel = 'Dy'
  elif (Z == 67):
    elementlabel = 'Ho'
  elif (Z == 68):
    elementlabel = 'Er'
  elif (Z == 69):
    elementlabel = 'Tm'
  elif (Z == 70):
    elementlabel = 'Yb'
  elif (Z == 71):
    elementlabel = 'Lu'
  elif (Z == 72):
    elementlabel = 'Hf'
  elif (Z == 73):
    elementlabel = 'Ta'
  elif (Z == 74):
    elementlabel = 'W'
  elif (Z == 75):
    elementlabel = 'Re'
  elif (Z == 76):
    elementlabel = 'Os'
  elif (Z == 77):
    elementlabel = 'Ir'
  elif (Z == 78):
    elementlabel = 'Pt'
  elif (Z == 79):
    elementlabel = 'Au'
  elif (Z == 80):
    elementlabel = 'Hg'
  elif (Z == 81):
    elementlabel = 'Tl'
  elif (Z == 82):
    elementlabel = 'Pb'
  elif (Z == 83):
    elementlabel = 'Bi'
  elif (Z == 84):
    elementlabel = 'Po'
  elif (Z == 85):
    elementlabel = 'At'
  elif (Z == 86):
    elementlabel = 'Rn'
  elif (Z == 87):
    elementlabel = 'Fr'
  elif (Z == 88):
    elementlabel = 'Ra'
  elif (Z == 89):
    elementlabel = 'Ac'
  elif (Z == 90):
    elementlabel = 'Th'
  elif (Z == 91):
    elementlabel = 'Pa'
  elif (Z == 92):
    elementlabel = 'U'
  elif (Z == 93):
    elementlabel = 'Np'
  elif (Z == 94):
    elementlabel = 'Pu'
  elif (Z == 95):
    elementlabel = 'Am'
  elif (Z == 96):
    elementlabel = 'Cm'
  elif (Z == 97):
    elementlabel = 'Bk'
  elif (Z == 98):
    elementlabel = 'Cf'
  elif (Z == 99):
    elementlabel = 'Es'
  elif (Z == 100):
    elementlabel = 'Fm'
  elif (Z == 101):
    elementlabel = 'Md'
  elif (Z == 102):
    elementlabel = 'No'
  elif (Z == 103):
    elementlabel = 'Lr'
  elif (Z == 104):
    elementlabel = 'Rf'
  elif (Z == 105):
    elementlabel = 'Db'
  elif (Z == 106):
    elementlabel = 'Sg'
  elif (Z == 107):
    elementlabel = 'Bh'
  elif (Z == 108):
    elementlabel = 'Hs'
  elif (Z == 109):
    elementlabel = 'Mt'
  elif (Z == 110):
    elementlabel = 'Ds'
  elif (Z == 111):
    elementlabel = 'Rg'
  elif (Z == 112):
    elementlabel = 'Cn'
  elif (Z == 113):
    elementlabel = 'Nh'
  elif (Z == 114):
    elementlabel = 'Fl'
  elif (Z == 115):
    elementlabel = 'Mc'
  elif (Z == 116):
    elementlabel = 'Lv'
  elif (Z == 117):
    elementlabel = 'Ts'
  elif (Z == 118):
    elementlabel = 'Og'
  elif (Z == 119):
    elementlabel = 'Uue'
  elif (Z == 120):
    elementlabel = 'Ubn'
  #========================
  #Outputs the element name
  #========================
  return elementlabel
#=====================================================================================================================
#This function takes all of the HFBTHOv300 output files from 'masstable' mode and puts the relevant data into one file
#=====================================================================================================================
def Read_HFBTHO_Masstable(thoout,bl,Output_Dict, Incomp_No, Incomp_Other):
    #==============================================================================================================
    #Goes through every available "thoout" output file, extracts useful data, and puts it onto the list Output_List
    #==============================================================================================================
    lines = open(thoout,'r').readlines()
    tho_name = thoout.split(".dat")[0]
    tho_name = tho_name.split("_")[-1]
    file_ID = bl + "-" + tho_name

    # Initialize variables
    convergence = "YES"; cpu_count=0
    N, Z, BE = 9999,9999,9999
    pairing_gap_N, pairing_gap_P = 0,0
    rms_radius_N, rms_radius_P, rms_radius_T = 0,0,0
    charge_radius, quad_def_beta2_N, quad_def_beta2_P, quad_def_beta2_T = 0,0,0,0
    quad_moment_Q2_N, quad_moment_Q2_P, quad_moment_Q2_T = 0,0,0
    oct_moment_Q3_N, oct_moment_Q3_P, oct_moment_Q3_T = 0,0,0

    for line in lines:
        if "iterations limit interrupt after1001" in line: convergence = "NO"
        if "CPU" in line: cpu_count += 1
        ss = line.split()
        try:
            #-----------------------------------------
            #Identifies the proton and neutron numbers
            #-----------------------------------------
            if (ss[0] == "Requested"):
                N = int(float(ss[2]) + 0.0001)
                Z = int(float(ss[3]) + 0.0001)
            #---------------------------
            #Identifies the pairing gaps
            #---------------------------
            elif((ss[0] == "delta(n,p),") and (ss[1] == "pwi")):
                pairing_gap_N = float(ss[3])     #Neutron pairing gap
                pairing_gap_P = float(ss[4])     #Proton pairing gap
            #------------------------
            #Identifies the rms-radii
            #------------------------
            elif ((ss[0] == "rms-radius") and (ss[1] == "..........")):
                rms_radius_N = float(ss[2])      #Neutron rms radius
                rms_radius_P = float(ss[3])      #Proton rms radius
                rms_radius_T = float(ss[4])      #Total rms radius
            #----------------------------
            #Identifies the charge radius
            #----------------------------
            elif ((ss[0] == "charge-radius,") and (ss[1] == "r0")):
                charge_radius = float(ss[3])     #Charge radius
            #------------------------------------------------
            #Identifies the quadrupole deformation parameters
            #------------------------------------------------
            elif((ss[0] == "deformation") and (ss[1] == "beta")):
                quad_def_beta2_N = float(ss[3])  #Neutron quadrupole deformation parameter
                quad_def_beta2_P = float(ss[4])  #Proton quadrupole deformation parameter
                quad_def_beta2_T = float(ss[5])  #Total quadrupole deformation parameter
            #---------------------------------
            #Identifies the quadrupole moments
            #---------------------------------
            elif((ss[0] == "quadrupole") and (ss[1] == "moment[b]")) and not quad_moment_Q2_T:
              # This gathers no LN deformation, these are needed for constraint calculation benchmark
              # Current HFBTHO has constraint on no LN deformation.
              quad_moment_Q2_N = float(ss[2])  #Neutron quadrupole moment
              quad_moment_Q2_P = float(ss[3])  #Proton quadrupole moment
              quad_moment_Q2_T = float(ss[4])  #Total quadrupole moment
            #-------------------------------
            #Identifies the octupole moments
            #-------------------------------
            elif((ss[0] == "octupole") and (ss[1] == "moment")) and not oct_moment_Q3_T:
              # This gathers no LN deformation, these are needed for constraint calculation benchmark
              # Current HFBTHO has constraint on no LN deformation.
              oct_moment_Q3_N = float(ss[3])  #Neutron octupole moment
              oct_moment_Q3_P = float(ss[4])  #Proton octupole moment
              oct_moment_Q3_T = float(ss[5])  #Total octupole moment
            #-----------------------------
            #Identifies the binding energy
            #-----------------------------
            elif ((ss[0] == 'tEnergy:') and (ss[1] == 'ehfb(qp)+LN')):
                BE = float(ss[2])                #Binding Energy
            #---------------------------------------------------------
            #No useful pieces of information, moves onto the next line
            #---------------------------------------------------------
            else:
                continue
        except IndexError:
            continue
    if Z > 200: return
    if cpu_count != 2: convergence = "***"
    if (Z,N) not in Output_Dict: Output_Dict[(Z,N)] = []
    Output_Dict[(Z,N)].append((Z,N,BE,quad_def_beta2_P,quad_def_beta2_N,quad_def_beta2_T,quad_moment_Q2_P,quad_moment_Q2_N,quad_moment_Q2_T,oct_moment_Q3_P,oct_moment_Q3_N,oct_moment_Q3_T,rms_radius_P,
                      rms_radius_N,rms_radius_T,charge_radius,pairing_gap_N,pairing_gap_P,file_ID,convergence))
    if convergence == "NO":
        Incomp_No.append((Z,N,file_ID,"No convergence"))
    if convergence == "***":
        Incomp_Other.append((Z,N,file_ID,"No convergence other"))
    return
#===========
#User Inputs
#===========
EDFs = ['SLY4', 'SV-MIN', 'UNEDF0', 'UNEDF1', 'UNEDF2']  # 'SKMS', 'SKP', 'UNEDF1-SO'
number_of_shells = 20

for functional in EDFs:
    # Locate block directories
    os.system("shopt -s extglob\n"+"rm HFBTHOv300_"+functional+"*.dat")
    os.chdir(functional)
    block_ls = os.listdir()
    blocks = []
    for bl in block_ls:
        if 'block' in bl and "." not in bl:
            blocks.append(bl)

    Output_Dict = {}  #Dict for output data
    Incomp_No, Incomp_Other = [], []

    #----------------------------------------------------------
    #Writes and properly formats the titles for the output file
    #----------------------------------------------------------
    all_data_str = '{:6} {:6} {:9} {:23} {:20} {:22} {:20} {:26} {:30} {:27} {:31} {:34} {:34} {:23} {:21} {:20} {:22} {:22} {:13} {:20} {:6} \n'.format(
                            'Z', 'N', 'A', 'Binding_Energy_(MeV)', 'Quad_Def_Beta2_P', 'Quad_Def_Beta2_N', 'Quad_Def_Beta2_total', 'Quad_Moment_Q2_P_(fm^2)', 'Quad_Moment_Q2_N_(fm^2)',
                            'Quad_Moment_Q2_total_(fm^2)', 'Octupole_Moment_Q3_P_(fm^3)', 'Octupole_Moment_Q3_N_(fm^3)', 'Octupole_Moment_Q3_total_(fm^3)', 'Pairing_gap_P_(MeV)',
                            'Pairing_gap_N_(MeV)', 'RMS_radius_P_(fm)', 'RMS_radius_N_(fm)', 'RMS_radius_total_(fm)', 'Charge_Radius_(fm)', 'File_ID',"Converged")

    for bl in blocks:
        os.chdir(bl)
        tho_ls = os.listdir()
        tho_list = []
        for fn in tho_ls:
            if "thoout" in fn and ".dat" in fn:
                tho_list.append(fn)
        print (functional,"\t",bl, "\tFile Count: ", len(tho_list))
        for ind,thoout in enumerate(tho_list):
            if not (ind+1) % 1000 or ind+1 == len(tho_list): print (ind+1,"/",len(tho_list))
            Read_HFBTHO_Masstable(thoout,bl,Output_Dict, Incomp_No, Incomp_Other)
        os.chdir("..")

    # All data of a single EDF should be stored in Output_Dict at this point, now we sort in order of Z,N,BE
    for key in sorted(Output_Dict):
        print (key)
        nuc_all = Output_Dict[key]
        # Sort on binding energy
        for entry in sorted(nuc_all, key = lambda x:x[2]):

            Z, N, BE = entry[0], entry[1], entry[2]
            file_ID, convergence = entry[18], entry[19]
            quad_def_beta2_P, quad_def_beta2_N, quad_def_beta2_T = entry[3],entry[4],entry[5]
            quad_moment_Q2_P, quad_moment_Q2_N, quad_moment_Q2_T = entry[6],entry[7],entry[8]
            oct_moment_Q3_P,  oct_moment_Q3_N,  oct_moment_Q3_T  = entry[9],entry[10],entry[11]
            rms_radius_P,     rms_radius_N,     rms_radius_T     = entry[12],entry[13],entry[14]
            charge_radius,    pairing_gap_N,    pairing_gap_P    = entry[15],entry[16],entry[17]

            all_data_str += '{:6} {:6} {:9} {:23} {:20} {:22} {:20} {:26} {:30} {:27} {:31} {:34} {:34} {:23} {:21} {:20} {:22} {:22} {:13} {:20} {:6}\n'.format(
                str(Z), str(N), str(Z+N), str(BE).rjust(13, ), str(quad_def_beta2_P).rjust(10, ), str(quad_def_beta2_N).rjust(10, ), str(quad_def_beta2_T).rjust(10, ),
                str(quad_moment_Q2_P).rjust(12, ), str(quad_moment_Q2_N).rjust(12, ), str(quad_moment_Q2_T).rjust(12, ), str(oct_moment_Q3_P).rjust(12, ),
                str(oct_moment_Q3_N).rjust(12, ), str(oct_moment_Q3_T).rjust(12, ), str(pairing_gap_P).rjust(10, ), str(pairing_gap_N).rjust(10, ),
                str(rms_radius_P).rjust(10, ), str(rms_radius_N).rjust(10, ), str(rms_radius_T).rjust(10, ), str(charge_radius).rjust(10, ), str(file_ID).rjust(12, ), str(convergence).rjust(6,))

    os.chdir("..")
    Data_File_Out = "HFBTHOv300_"+functional+"_All_Data_"+str(number_of_shells)+"_shells_no_LN_deformation-masstable.dat"  #Output file for Read_HFBTHO_Masstable_Output
    all_data_output = open(Data_File_Out, "w")    #Output file for all data
    all_data_output.write(all_data_str)
    all_data_output.close()
    print ("Incomplete:\n")
    for inp in Incomp_No:
        print (inp[0],"\t",inp[1],"\t",inp[2])
    print ("Incomplete Other:\n")
    for inp in Incomp_Other:
        print (inp[0],"\t",inp[1],"\t",inp[2])
