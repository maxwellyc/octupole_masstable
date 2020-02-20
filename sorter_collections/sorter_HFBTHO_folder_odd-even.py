# This code organizes the individual output files for HFBTHOv300 into one single file, for each nucleus the outputs are sorted by
# ascending binding energy, meaning the first occurance of each nucleus is from the ground state output file.
# WARNING:: Please check last column of each row for whether this row is from converged calculation. This convergence check is
# only valid for even-even nuclei kick-off mode calculations, for all other calculations this may not be accurate. 09/18/2019

# HFBTHOv300_"functional-name"_All_Data.dat (containing the data from every constrained calculation)

# The input files are assumed to have the following form:  thoout_000001.dat

# For HFBTHOv300 odd-even nuclei non-blocking kickoff calculation, sometimes (unknown pattern as of this date) the calculation
# doesn't terminate and will repeat itself again and again.

# Maxwell Cao 09/18/2019

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
               #  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,  9 , 10 ,
  elementList = [
                 'H','He','Li','Be', 'B', 'C', 'N', 'O', 'F','Ne',\
                'Na','Mg','Al','Si', 'P', 'S','Cl','Ar', 'K','Ca',\
                'Sc','Ti', 'V','Cr','Mn','Fe','Co','Ni','Cu','Zn',\
                'Ga','Ge','As','Se','Br','Kr','Rb','Sr', 'Y','Zr',\
                'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',\
                'Sb','Te', 'I','Xe','Cs','Ba','La','Ce','Pr','Nd',\
                'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',\
                'Lu','Hf','Ta', 'W','Re','Os','Ir','Pt','Au','Hg',\
                'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',\
                'Pa', 'U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',\
                'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds',\
                'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og','Uue','Ubn'
                ]
  elementLabel = {ind+1:label for ind,label in enumerate(elementList)}
  return elementLabel[Z]
#=====================================================================================================================
#This function takes all of the HFBTHOv300 output files from 'masstable' mode and puts the relevant data into one file
#=====================================================================================================================
def Read_HFBTHO_Masstable(thoout,Output_Dict, Incomp_No):
    #==============================================================================================================
    #Goes through every available "thoout" output file, extracts useful data, and puts it onto the list Output_List
    #==============================================================================================================
    lines = open(thoout,encoding="ISO-8859-1").readlines()
    tho_name = thoout.split(".dat")[0]
    tho_name = tho_name.split("_")[-1]
    file_ID = tho_name

    # Initialize variables
    N, Z, BE = 9999,9999,9999
    pairing_gap_N, pairing_gap_P = 0,0
    rms_radius_N, rms_radius_P, rms_radius_T = 0,0,0
    charge_radius, quad_def_beta2_N, quad_def_beta2_P, quad_def_beta2_T = 0,0,0,0
    quad_moment_Q2_N, quad_moment_Q2_P, quad_moment_Q2_T = 0,0,0
    oct_moment_Q3_N, oct_moment_Q3_P, oct_moment_Q3_T = 0,0,0

    # Need to refresh count when encountering "HARTREE-FOCK-BOGOLIUBOV" because this
    # means a new calculation was initialized
    # In each new calculation, the convergence is reached when one "CPU" keyword is seen
    # followed by "iteration converged"
    convergence = "NO"; cpu_count = 0
    iter_converged = 0

    for line in lines:
        if "HARTREE-FOCK-BOGOLIUBOV" in line:
            cpu_count = 0
        elif "CPU" in line:
            cpu_count = 1
        if cpu_count == 1 and "converged" in line: convergence = "YES"
        ss = line.split()
        try:
            #-----------------------------------------
            #Identifies the proton and neutron numbers
            #-----------------------------------------
            if (ss[0] == "Requested"):
                N, Z = int(float(ss[2]) + 0.0001), int(float(ss[3]) + 0.0001)
            #---------------------------
            #Identifies the pairing gaps
            #---------------------------
            elif((ss[0] == "delta(n,p),") and (ss[1] == "pwi")):
                pairing_gap_N, pairing_gap_P = float(ss[3]), float(ss[4])
            #------------------------
            #Identifies the rms-radii
            #------------------------
            elif ((ss[0] == "rms-radius") and (ss[1] == "..........")):
                rms_radius_N, rms_radius_P, rms_radius_T = float(ss[2]), float(ss[3]), float(ss[4])
            #----------------------------
            #Identifies the charge radius
            #----------------------------
            elif ((ss[0] == "charge-radius,") and (ss[1] == "r0")):
                charge_radius = float(ss[3])     #Charge radius
            #------------------------------------------------
            #Identifies the quadrupole deformation parameters
            #------------------------------------------------
            elif((ss[0] == "deformation") and (ss[1] == "beta")):
                quad_def_beta2_N, quad_def_beta2_P, quad_def_beta2_T = float(ss[3]), float(ss[4]), float(ss[5])
            #---------------------------------
            #Identifies the quadrupole moments
            #---------------------------------
            elif((ss[0] == "quadrupole") and (ss[1] == "moment[b]")) and not quad_moment_Q2_T:
              # This gathers no LN deformation, these are needed for constraint calculation benchmark
              # Current HFBTHO has constraint on no LN deformation.
              quad_moment_Q2_N, quad_moment_Q2_P, quad_moment_Q2_T  = float(ss[2]), float(ss[3]), float(ss[4])
            #-------------------------------
            #Identifies the octupole moments
            #-------------------------------
            elif((ss[0] == "octupole") and (ss[1] == "moment")) and not oct_moment_Q3_T:
              # This gathers no LN deformation, these are needed for constraint calculation benchmark
              # Current HFBTHO has constraint on no LN deformation.
              oct_moment_Q3_N, oct_moment_Q3_P, oct_moment_Q3_T = float(ss[3]), float(ss[4]), float(ss[5])
            #-----------------------------
            #Identifies the binding energy
            #-----------------------------
            elif ((ss[0] == 'tEnergy:') and (ss[1] == 'ehfb(qp)+LN')):
                BE = float(ss[2])                #Binding Energy
            #---------------------------------------------------------
            #No useful pieces of information, moves onto the next line
            #---------------------------------------------------------
            elif "Blocking candidates" in line and convergence == "YES":
                break
            else:
                continue
        except IndexError:
            continue
    if Z > 200: return
    if (Z,N) not in Output_Dict: Output_Dict[(Z,N)] = []
    Output_Dict[(Z,N)].append((Z,N,BE,quad_def_beta2_P,quad_def_beta2_N,quad_def_beta2_T,quad_moment_Q2_P,quad_moment_Q2_N,quad_moment_Q2_T,oct_moment_Q3_P,oct_moment_Q3_N,oct_moment_Q3_T,rms_radius_P,
                      rms_radius_N,rms_radius_T,charge_radius,pairing_gap_N,pairing_gap_P,file_ID,convergence))
    if convergence == "NO":
        Incomp_No.append((Z,N,file_ID,"No convergence"))
    return

def read_entry(entry):
    Z, N, BE = entry[0], entry[1], entry[2]
    file_ID, convergence = entry[18], entry[19]
    quad_def_beta2_P, quad_def_beta2_N, quad_def_beta2_T = entry[3],entry[4],entry[5]
    quad_moment_Q2_P, quad_moment_Q2_N, quad_moment_Q2_T = entry[6],entry[7],entry[8]
    oct_moment_Q3_P,  oct_moment_Q3_N,  oct_moment_Q3_T  = entry[9],entry[10],entry[11]
    rms_radius_P,     rms_radius_N,     rms_radius_T     = entry[12],entry[13],entry[14]
    charge_radius,    pairing_gap_N,    pairing_gap_P    = entry[15],entry[16],entry[17]
    return '{:6} {:6} {:9} {:23} {:20} {:22} {:20} {:26} {:30} {:27} {:31} {:34} {:34} {:23} {:21} {:20} {:22} {:22} {:13} {:20} {:6}\n'.format(
        str(Z), str(N), str(Z+N), str(BE).rjust(13, ), str(quad_def_beta2_P).rjust(10, ), str(quad_def_beta2_N).rjust(10, ), str(quad_def_beta2_T).rjust(10, ),
        str(quad_moment_Q2_P).rjust(12, ), str(quad_moment_Q2_N).rjust(12, ), str(quad_moment_Q2_T).rjust(12, ), str(oct_moment_Q3_P).rjust(12, ),
        str(oct_moment_Q3_N).rjust(12, ), str(oct_moment_Q3_T).rjust(12, ), str(pairing_gap_P).rjust(10, ), str(pairing_gap_N).rjust(10, ),
        str(rms_radius_P).rjust(10, ), str(rms_radius_N).rjust(10, ), str(rms_radius_T).rjust(10, ), str(charge_radius).rjust(10, ), str(file_ID).rjust(12, ), str(convergence).rjust(6,))

#===========
#User Inputs
#===========
EDFs = ['SLY4', 'SV-MIN', 'UNEDF0', 'UNEDF1', 'UNEDF2']  # 'SKMS', 'SKP', 'UNEDF1-SO'
number_of_shells = 20

for functional in EDFs:
    # Locate block directories
    os.system("shopt -s extglob\n"+"rm HFBTHOv300_"+functional+"*.dat")
    os.chdir(functional)
    Output_Dict, Incomp_No = {}, []  # Dict for output data, incomplete list

    #----------------------------------------------------------
    #Writes and properly formats the titles for the output file
    #----------------------------------------------------------
    all_data_str = '{:6} {:6} {:9} {:23} {:20} {:22} {:20} {:26} {:30} {:27} {:31} {:34} {:34} {:23} {:21} {:20} {:22} {:22} {:13} {:20} {:6} \n'.format(
                            'Z', 'N', 'A', 'Binding_Energy_(MeV)', 'Quad_Def_Beta2_P', 'Quad_Def_Beta2_N', 'Quad_Def_Beta2_total', 'Quad_Moment_Q2_P_(fm^2)', 'Quad_Moment_Q2_N_(fm^2)',
                            'Quad_Moment_Q2_total_(fm^2)', 'Octupole_Moment_Q3_P_(fm^3)', 'Octupole_Moment_Q3_N_(fm^3)', 'Octupole_Moment_Q3_total_(fm^3)', 'Pairing_gap_P_(MeV)',
                            'Pairing_gap_N_(MeV)', 'RMS_radius_P_(fm)', 'RMS_radius_N_(fm)', 'RMS_radius_total_(fm)', 'Charge_Radius_(fm)', 'File_ID',"Converged")


    tho_ls = os.listdir()
    tho_list = []
    for fn in tho_ls:
        if "thoout" in fn and ".dat" in fn:
            tho_list.append(fn)
    for ind,thoout in enumerate(tho_list):
        Read_HFBTHO_Masstable(thoout,Output_Dict, Incomp_No)
    os.chdir("..")

    # All data of a single EDF should be stored in Output_Dict at this point, now we sort in order of Z,N,BE
    for key in sorted(Output_Dict):
        # Sort Z, N value of nucleus
        nuc_all = Output_Dict[key]
        # Split into converged and unconverged and then sort on binding energy:
        nuc_conv, nuc_unconv = [], []
        for entry in nuc_all:
            if entry[-1] == "YES":
                nuc_conv.append(entry)
            elif entry[-1] == "NO":
                nuc_unconv.append(entry)
        nuc_conv.sort(key = lambda x:x[2])
        nuc_unconv.sort(key = lambda x:x[2])
        for entry in nuc_conv:
            all_data_str += read_entry(entry)
        for entry in nuc_unconv:
            all_data_str += read_entry(entry)

    Data_File_Out = "HFBTHOv300_"+functional+"_All_Data_"+str(number_of_shells)+"_shells_no_LN_deformation-masstable.dat"  #Output file for Read_HFBTHO_Masstable_Output
    all_data_output = open(Data_File_Out, "w")    #Output file for all data
    all_data_output.write(all_data_str)
    all_data_output.close()
    print ("\n",functional,"\tFile Count: ", len(tho_list),', unconverged: ',len(Incomp_No))
    print ("BE_converged - BE_unconverged (MeV):", round(nuc_conv[0][2] - nuc_unconv[0][2],6))
    print ("Converged ground state file ID: ",nuc_conv[0][-2])
    print ("========================================")

    # for inp in Incomp_No:
    #     print (inp[0],"\t",inp[1],"\t",inp[2],"\t unconverged\t",functional)
