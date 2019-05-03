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

#==========================================================
#This function reads data from the file hfbtho_NAMELIST.dat
#==========================================================
def read_hfbtho_NAMELIST():
  Input_File = open("hfbtho_NAMELIST.dat")    #File to be read: hfbtho_NAMELIST.dat
  lines = Input_File.readlines()              #Reads all the lines of the input file
  for line in lines:
    ss = line.split()
    try:
      if (ss[0] == 'number_of_shells'):
        number_of_shells = re.sub("[^\d\.]","",ss[2])  #Regular expressions needed for formatting; number_of_shells number has a comma at the end of it
    except IndexError:
      continue
  #------------------------------------
  #Gives the number of shells as output
  #------------------------------------
  return number_of_shells
#===========================================================
#This function reads data from the file hfbtho_MASSTABLE.dat
#===========================================================
def read_hfbtho_MASSTABLE():
  Input_File = open("hfbtho_MASSTABLE.dat")   #File to be read: hfbtho_MASSTABLE.dat
  lines = Input_File.readlines()              #Reads all of the lines of the input file
  i=0    #Iteration variable for each line of the input file
  for line in lines:
    ss = line.split()   #Splits each element of the line separated by empty space
    #-------------------------------------------------------------------------------
    #Extracts data for Highest_File_ID, lower_proton_number, and upper_proton_number
    #-------------------------------------------------------------------------------
    if (i == 1):  #Second line of the code: gives HighestFileID
      Highest_File_ID = int(ss[0])
      i=i+1
      continue
    elif (i >= 3):                       #First line with numbers; gives lower_proton_number
      if (i == 3):
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        #Note: the numbers here must be defined as integers otherwise the code has trouble going from strings of length 2 (like 98) to length 3 (like 100) and gives incorrect results
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        lower_proton_number = int(ss[0])      #Lowest proton number for this particular mass table
        i=i+1
        continue
      else:
        if (i == (Highest_File_ID-1)):
          upper_proton_number = int(ss[0])
          exit
        else:
          i=i+1
          continue
    else:
      i=i+1
      continue
  #----------------------------------------------------------------------
  #Gives the lower and upper proton numbers and highest file ID as output
  #----------------------------------------------------------------------
  return (str(lower_proton_number), str(upper_proton_number), Highest_File_ID)
#=====================================================
#This function reads data from the file hfbtho_PES.dat
#=====================================================
def read_hfbtho_PES():
  Input_File = open("hfbtho_PES.dat")
  lines = Input_File.readlines()              #Reads all of the lines of the input file
  i=0    #Iteration variable for each line of the input file
  for line in lines:
    ss = line.split()   #Splits each element of the line separated by empty space
    #-------------------------------------------------------------------------------
    #Extracts data for Highest_File_ID, lower_proton_number, and upper_proton_number
    #-------------------------------------------------------------------------------
    if (i == 1):  #Second line of the code: gives HighestFileID
      Highest_File_ID = int(ss[0])
      i=i+1
      continue
    elif (i >= 5):
      if (i == 5):
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        #Note: the numbers here must be defined as integers otherwise the code has trouble going from strings of length 2 (like 98) to length 3 (like 100) and gives incorrect results
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        lower_proton_number = int(ss[0])      #Lowest proton number for this particular mass table
        i=i+1
        continue
      else:
        if (i == (Highest_File_ID-1)):
          upper_proton_number = int(ss[0])
          exit
        else:
          i=i+1
          continue
    else:
      i=i+1
      continue
  #----------------------------------------------------------------------
  #Gives the lower and upper proton numbers and highest file ID as output
  #----------------------------------------------------------------------
  return (str(lower_proton_number), str(upper_proton_number), Highest_File_ID)
#=======================================================================================
#This function takes in a proton number and outputs the abbreviation of its element name
#=======================================================================================
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
#====================================================================================================================
#This function takes all of the HFBTHOv300 output files from 'dripline' mode and puts the relevant data into one file
#====================================================================================================================
def Read_HFBTHO_Dripline_Output(functional,Data_File_In,Data_File_Out,LowestTeamID,TeamID,LowestFileID,FileID):
  lowest_team_ID = int(LowestTeamID)
  current_team_ID = int(TeamID)
  lowest_file_ID = int(LowestFileID)
  current_file_ID = int(FileID)
  #=========================================================================
  #Deletes the current HFBTHOv300d_All_Data file to be replaced by a new one
  #=========================================================================
  if ((lowest_team_ID == current_team_ID) and (lowest_file_ID == current_file_ID) and (os.path.isfile("HFBTHOv300_"+functional+"_All_Data.dat"))):
    os.remove("HFBTHOv300_"+functional+"_All_Data.dat")
  f1 = open(str(Data_File_In))
  all_data_output = open(str(Data_File_Out), "a")
  lines = f1.readlines()
  BEdic={}        #Dictionary for binding energies
  Output_List=[]  #List for output data
  for line in lines:
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
      elif((ss[0] == "quadrupole") and (ss[1] == "moment[b]")):
        quad_moment_Q2_N = float(ss[2])  #Neutron quadrupole moment
        quad_moment_Q2_P = float(ss[3])  #Proton quadrupole moment
        quad_moment_Q2_T = float(ss[4])  #Total quadrupole moment
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
  f1.close()
  #==================================
  #Transfers the data above to a list
  #==================================
  Output_List.append((Z,N,BE,quad_def_beta2_P,quad_def_beta2_N,quad_def_beta2_T,quad_moment_Q2_P,quad_moment_Q2_N,quad_moment_Q2_T,rms_radius_P,rms_radius_N,rms_radius_T,charge_radius,pairing_gap_N,pairing_gap_P))
  Output_List.sort()
  #==========================================================
  #Writes and properly formats the titles for the output file
  #==========================================================
  if ((current_team_ID == lowest_team_ID) and (current_file_ID == lowest_file_ID)):
    all_data_output.write( '{:6} {:6} {:7} {:23} {:20} {:20} {:22} {:26} {:26} {:31} {:23} {:21} {:20} {:20} {:24} {:22} \n'.format('Z', 'N', 'A', 'Binding_Energy_(MeV)', 'Quad_Def_Beta2_P', 'Quad_Def_Beta2_N', 'Quad_Def_Beta2_total', 'Quad_Moment_Q2_P_(fm^2)', 'Quad_Moment_Q2_N_(fm^2)', 'Quad_Moment_Q2_total_(fm^2)', 'Pairing_gap_P_(MeV)', 'Pairing_gap_N_(MeV)', 'RMS_radius_P_(fm)', 'RMS_radius_N_(fm)', 'RMS_radius_total_(fm)', 'Charge_Radius_(fm)'))
  #=================================
  #Writes the list to an output file
  #=================================
  for item in Output_List:
    #---------------------------------
    #Proton, neutron, and mass numbers
    #---------------------------------
    Z = Output_List[0][0]
    N = Output_List[0][1]
    A = int(Z) + int(N)
    #--------------
    #Binding energy
    #--------------
    BE = Output_List[0][2]
    #-----------------------
    #Quadrupole deformations
    #-----------------------
    quad_def_beta2_P = Output_List[0][3]
    quad_def_beta2_N = Output_List[0][4]
    quad_def_beta2_T = Output_List[0][5]
    #------------------
    #Quadrupole moments
    #------------------
    quad_moment_Q2_P = Output_List[0][6]
    quad_moment_Q2_N = Output_List[0][7]
    quad_moment_Q2_T = Output_List[0][8]
    #---------
    #RMS radii
    #---------
    rms_radius_P = Output_List[0][9]
    rms_radius_N = Output_List[0][10]
    rms_radius_T = Output_List[0][11]
    #-------------
    #Charge radius
    #-------------
    charge_radius = Output_List[0][12]
    #------------
    #Pairing gaps
    #------------
    proton_pairing_gap = Output_List[0][13]
    neutron_pairing_gap = Output_List[0][14]
    #---------------------------------
    #Writes the data to an output file
    #---------------------------------
    all_data_output.write( '{:6} {:6} {:9} {:23} {:20} {:22} {:20} {:26} {:30} {:27} {:23} {:21} {:20} {:22} {:22} {:17} \n'.format(str(Z), str(N), str(A), str(BE).rjust(13, ), str(quad_def_beta2_P).rjust(10, ), str(quad_def_beta2_N).rjust(10, ), str(quad_def_beta2_T).rjust(10, ), str(quad_moment_Q2_P).rjust(12, ), str(quad_moment_Q2_N).rjust(12, ), str(quad_moment_Q2_T).rjust(12, ), str(proton_pairing_gap).rjust(10, ), str(neutron_pairing_gap).rjust(10, ), str(rms_radius_P).rjust(10, ), str(rms_radius_N).rjust(10, ), str(rms_radius_T).rjust(10, ), str(charge_radius).rjust(10, )))
  all_data_output.close()
#=====================================================================================================================
#This function takes all of the HFBTHOv300 output files from 'masstable' mode and puts the relevant data into one file
#=====================================================================================================================
def Read_HFBTHO_Masstable_and_PES_Output(functional,Data_File_Out,Lowest_File_ID,Highest_File_ID):
  lowest_file_ID = int(Lowest_File_ID)
  highest_file_ID = int(Highest_File_ID)
  all_data_output = open(str(Data_File_Out), "a")    #Output file for all data
  Output_List=[]  #List for output data
  #==============================================================================================================
  #Goes through every available "thoout" output file, extracts useful data, and puts it onto the list Output_List 
  #==============================================================================================================
  for i in range(lowest_file_ID,highest_file_ID+1):
    #-------------------------------------------------------------------------
    #Deletes the current HFBTHOv300d_All_Data file to be replaced by a new one
    #-------------------------------------------------------------------------
    if ((i == lowest_file_ID) and (os.path.isfile("HFBTHOv300_"+functional+"_All_Data.dat"))):
      os.remove("HFBTHOv300_"+functional+"_All_Data.dat")
    try:
      FileID = str(i).rjust(6, '0')   #Proper formatting for the number in each file
      f1 = open("thoout_"+FileID+".dat")    #Open input file to be read
      lines = f1.readlines()
      convergence = "YES"; cpu_count=0
      #-----------------------------------
      #Reading the lines of the input file
      #-----------------------------------
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
          elif((ss[0] == "quadrupole") and (ss[1] == "moment[b]")):
            quad_moment_Q2_N = float(ss[2])  #Neutron quadrupole moment
            quad_moment_Q2_P = float(ss[3])  #Proton quadrupole moment
            quad_moment_Q2_T = float(ss[4])  #Total quadrupole moment
          #-------------------------------
          #Identifies the octupole moments
          #-------------------------------
          elif((ss[0] == "octupole") and (ss[1] == "moment")):
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
      if cpu_count != 2: convergence = "***"
      f1.close()
      current_file_ID=i    #File ID of this calculation (useful for troubleshooting)
      #----------------------------------
      #Transfers the data above to a list
      #----------------------------------
      Output_List.append((Z,N,BE,quad_def_beta2_P,quad_def_beta2_N,quad_def_beta2_T,quad_moment_Q2_P,quad_moment_Q2_N,quad_moment_Q2_T,oct_moment_Q3_P,oct_moment_Q3_N,oct_moment_Q3_T,rms_radius_P,
                          rms_radius_N,rms_radius_T,charge_radius,pairing_gap_N,pairing_gap_P,current_file_ID,convergence))
      #----------------------------------------------------------
      #Writes and properly formats the titles for the output file
      #----------------------------------------------------------
      if (i == lowest_file_ID):   #Ensures the column labels are only printed once on the output
        all_data_output.write( '{:6} {:6} {:7} {:23} {:20} {:20} {:22} {:26} {:26} {:31} {:31} {:31} {:36} {:23} {:21} {:20} {:20} {:24} {:22} {:14} {:20} \n'.format(
                                'Z', 'N', 'A', 'Binding_Energy_(MeV)', 'Quad_Def_Beta2_P', 'Quad_Def_Beta2_N', 'Quad_Def_Beta2_total', 'Quad_Moment_Q2_P_(fm^2)', 'Quad_Moment_Q2_N_(fm^2)',
                                'Quad_Moment_Q2_total_(fm^2)', 'Octupole_Moment_Q3_P_(fm^3)', 'Octupole_Moment_Q3_N_(fm^3)', 'Octupole_Moment_Q3_total_(fm^3)', 'Pairing_gap_P_(MeV)',
                                'Pairing_gap_N_(MeV)', 'RMS_radius_P_(fm)', 'RMS_radius_N_(fm)', 'RMS_radius_total_(fm)', 'Charge_Radius_(fm)', 'File_ID',"Converged"))
      #---------------------
      #Looks for the file ID
      #---------------------
      if cpu_count == 2: print (i)
      else: print(i,"Incomplete")
      i=i+1
    except IOError:
      continue
  #=====================================================================
  #Sorts the data so all calculations with the same Z and N are together
  #=====================================================================
  Output_List.sort()
  i=0  #Reinitializes the iterative variable
  #---------------------------------
  #Writes the list to an output file
  #---------------------------------
  for item in Output_List:
    #---------------------------------
    #Proton, neutron, and mass numbers
    #---------------------------------
    Z = Output_List[i][0]
    N = Output_List[i][1]
    A = int(Z) + int(N)
    #--------------
    #Binding energy
    #--------------
    BE = Output_List[i][2]
    #-----------------------
    #Quadrupole deformations
    #-----------------------
    quad_def_beta2_P = Output_List[i][3]
    quad_def_beta2_N = Output_List[i][4]
    quad_def_beta2_T = Output_List[i][5]
    #------------------
    #Quadrupole moments
    #------------------
    quad_moment_Q2_P = Output_List[i][6]
    quad_moment_Q2_N = Output_List[i][7]
    quad_moment_Q2_T = Output_List[i][8]
    #----------------
    #Octupole moments
    #----------------
    oct_moment_Q3_P = Output_List[i][9]
    oct_moment_Q3_N = Output_List[i][10]
    oct_moment_Q3_T = Output_List[i][11]
    #---------
    #RMS radii
    #---------
    rms_radius_P = Output_List[i][12]
    rms_radius_N = Output_List[i][13]
    rms_radius_T = Output_List[i][14]
    #-------------
    #Charge radius
    #-------------
    charge_radius = Output_List[i][15]
    #------------
    #Pairing gaps
    #------------
    proton_pairing_gap = Output_List[i][16]
    neutron_pairing_gap = Output_List[i][17]
    #-------
    #File ID
    #-------
    current_file_ID = Output_List[i][18]
    #-------
    #Convergence check
    #-------
    file_convergence = Output_List[i][19]
    #---------------------------------
    #Writes the data to an output file
    #---------------------------------
    all_data_output.write( '{:6} {:6} {:9} {:23} {:20} {:22} {:20} {:26} {:30} {:27} {:31} {:34} {:34} {:23} {:21} {:20} {:22} {:22} {:13} {:20} {:6}\n'.format(
                              str(Z), str(N), str(A), str(BE).rjust(13, ), str(quad_def_beta2_P).rjust(10, ), str(quad_def_beta2_N).rjust(10, ), str(quad_def_beta2_T).rjust(10, ),
                              str(quad_moment_Q2_P).rjust(12, ), str(quad_moment_Q2_N).rjust(12, ), str(quad_moment_Q2_T).rjust(12, ), str(oct_moment_Q3_P).rjust(12, ),
                              str(oct_moment_Q3_N).rjust(12, ), str(oct_moment_Q3_T).rjust(12, ), str(proton_pairing_gap).rjust(10, ), str(neutron_pairing_gap).rjust(10, ),
                              str(rms_radius_P).rjust(10, ), str(rms_radius_N).rjust(10, ), str(rms_radius_T).rjust(10, ), str(charge_radius).rjust(10, ), str(current_file_ID).rjust(10, ), str(file_convergence).rjust(6,)))
    i=i+1
  #======================
  #Closes the output file
  #======================
  all_data_output.close()
#
#
#
#=====================================================================================
#Function which finds the oblate, prolate, and spherical minima from the All_Data file
#=====================================================================================
def All_Minima(Total_Input_File,Minima_Output_File):
  Input_File=open(str(Total_Input_File))
  Output_File=open(str(Minima_Output_File), "w")
  lines=Input_File.readlines()
  BE_Oblate={}               #Dictionary for binding energies of oblate constraint calculations
  BE_Spherical={}            #Dictionary for the binding energy of the spherical constraint calculation
  BE_Prolate={}              #Dictionary for the binding energies of prolate constraint calculations
  Oblate_Minimum_Line={}     #Dictionary for the lines of data for oblate minima
  Spherical_Minimum_Line={}  #Dictionary for the lines of data for the spherical minimum
  Prolate_Minimum_Line={}    #Dictionary for the lines of data for prolate minima
  All_Minima_Output=[]       #List of data for all the minima of a nucleus
  for line in lines:
    ss=line.split()          #Splits each element of the line separated by a space
    try:
      Z=int(ss[0])                    #Proton number
      N=int(ss[1])                    #Neutron number
      BE=float(ss[3])                 #Binding energy
      Q2=round(float(ss[9]))          #Quadrupole deformation
      #==============
      #Oblate minimum
      #==============
      if (Q2 < -0.02):
        if ((N,Z) in BE_Oblate):
          if (BE < BE_Oblate[(N,Z)]):
            BE_Oblate[(N,Z)]=BE
            Oblate_Minimum_Line[(N,Z)]=line
          else:
            continue
        else:
          BE_Oblate[(N,Z)]=BE
          Oblate_Minimum_Line[(N,Z)]=line
      #=================
      #Spherical minimum
      #=================
      elif (-0.02 <= Q2 <= 0.02):
        if ((N,Z) in BE_Spherical):
          if (BE < BE_Spherical[(N,Z)]):
            BE_Spherical[(N,Z)]=BE
            Spherical_Minimum_Line[(N,Z)]=line
          else:
            continue
        else:
          BE_Spherical[(N,Z)]=BE
          Spherical_Minimum_Line[(N,Z)]=line
      #===============
      #Prolate minimum
      #===============
      elif (Q2 > 0.02):
        if ((N,Z) in BE_Prolate):
          if (BE < BE_Prolate[(N,Z)]):
            BE_Prolate[(N,Z)]=BE
            Prolate_Minimum_Line[(N,Z)]=line
          else:
            continue
        else:
          BE_Prolate[(N,Z)]=BE
          Prolate_Minimum_Line[(N,Z)]=line
    except ValueError:
      continue
  Input_File.close()
  #===========================================================
  #Transfers the data from the dictionaries into a single list
  #===========================================================
  for key in BE_Oblate.keys():
    N=key[0]; Z=key[1]
    All_Minima_Output.append((Z,N,Oblate_Minimum_Line[(N,Z)]))
    All_Minima_Output.append((Z,N,Spherical_Minimum_Line[(N,Z)]))
    All_Minima_Output.append((Z,N,Prolate_Minimum_Line[(N,Z)]))
  All_Minima_Output.sort()
  #----------------------------------------------------------
  #Writes and properly formats the titles for the output file
  #----------------------------------------------------------
  Output_File.write( '{:6} {:6} {:7} {:23} {:20} {:20} {:22} {:26} {:26} {:31} {:31} {:31} {:36} {:23} {:21} {:20} {:20} {:24} {:22} {:22} \n'.format(
                          'Z', 'N', 'A', 'Binding_Energy_(MeV)', 'Quad_Def_Beta2_P', 'Quad_Def_Beta2_N', 'Quad_Def_Beta2_total', 'Quad_Moment_Q2_P_(fm^2)', 'Quad_Moment_Q2_N_(fm^2)',
                          'Quad_Moment_Q2_total_(fm^2)', 'Octupole_Moment_Q3_P_(fm^3)', 'Octupole_Moment_Q3_N_(fm^3)', 'Octupole_Moment_Q3_total_(fm^3)','Pairing_gap_P_(MeV)', 'Pairing_gap_N_(MeV)',
                          'RMS_radius_P_(fm)', 'RMS_radius_N_(fm)', 'RMS_radius_total_(fm)', 'Charge_Radius_(fm)', 'File_ID'))  
  #=================================================
  #Transfers the data on the list to the output file
  #=================================================
  for item in All_Minima_Output:
    sp=item[2].split()
    #---------------------------------
    #Proton, neutron, and mass numbers
    #---------------------------------
    Z = sp[0]
    N = sp[1]
    A = sp[2]
    #--------------
    #Binding energy
    #--------------
    BE = sp[3]
    #-----------------------
    #Quadrupole deformations
    #-----------------------
    quad_def_beta2_P = sp[4]
    quad_def_beta2_N = sp[5]
    quad_def_beta2_T = sp[6]
    #------------------
    #Quadrupole moments
    #------------------
    quad_moment_Q2_P = sp[7]
    quad_moment_Q2_N = sp[8]
    quad_moment_Q2_T = sp[9]
    #----------------
    #Octupole moments
    #----------------
    oct_moment_Q3_P = sp[10]
    oct_moment_Q3_N = sp[11]
    oct_moment_Q3_T = sp[12]
    #------------
    #Pairing gaps
    #------------
    proton_pairing_gap = sp[13]
    neutron_pairing_gap = sp[14]
    #---------
    #RMS radii
    #---------
    rms_radius_P = sp[15]
    rms_radius_N = sp[16]
    rms_radius_T = sp[17]
    #-------------
    #Charge radius
    #-------------
    charge_radius = sp[18]
    #-------
    #File ID
    #-------
    current_file_ID = sp[19]
    #---------------------------------
    #Writes the data to an output file
    #---------------------------------
    Output_File.write( '{:6} {:6} {:9} {:23} {:20} {:22} {:20} {:26} {:30} {:27} {:31} {:34} {:34} {:23} {:21} {:20} {:22} {:22} {:17} {:20} \n'.format(
                              str(Z), str(N), str(A), str(BE).rjust(13, ), str(quad_def_beta2_P).rjust(10, ), str(quad_def_beta2_N).rjust(10, ), str(quad_def_beta2_T).rjust(10, ),
                              str(quad_moment_Q2_P).rjust(12, ), str(quad_moment_Q2_N).rjust(12, ), str(quad_moment_Q2_T).rjust(12, ), str(oct_moment_Q3_P).rjust(12, ),
                              str(oct_moment_Q3_N).rjust(12, ), str(oct_moment_Q3_T).rjust(12, ),str(proton_pairing_gap).rjust(10, ), str(neutron_pairing_gap).rjust(10, ),
                              str(rms_radius_P).rjust(10, ), str(rms_radius_N).rjust(10, ), str(rms_radius_T).rjust(10, ), str(charge_radius).rjust(10, ), str(current_file_ID).rjust(10, )))
  Output_File.close()
#
#
#
#============================================================================================
#Function which finds the properties of each nucleus corresponding to the ground state energy
#============================================================================================
def HFBTHOv300_Ground_State_Mass_Table(Total_Input_File,GS_Output_File,az,an,ae,LNrp,LNrn,PPG,NPG):
  f1 = open(str(Total_Input_File))                    #Defines input file
  ground_state_data = open(str(GS_Output_File), "w")  #Defines output file
  lines = f1.readlines()                              #Reads the input file
  BEdic={}                                            #Dictionary of binding energy values
  DATA_Line={}                                        #Dictionary for the lines of data on the mass table file
  neutron_skin_data={}                                #Dictionary for calculated neutron skin values
  ProtonPairingGap={}                                 #Dictionary of proton pairing gap values
  NeutronPairingGap={}                                #Dictionary of neutron pairing gap values
  EvenMassTableDATA=[]                                #List for the drip line sorted mass table (even-even nuclei)
  Mass_Table_Data=[]                                  #List of output data
  i=0                                                 #Iteration variable (for current nucleus)
  j=-1                                                #Iteration variable (for previous nucleus)
  k=1                                                 #Iteration variable (for next nucleus)
  pdripreached=0                                      #Counter to see if the list has hit the proton drip line yet
  protondrip=0		                              #Current N value of the 2 proton drip line
  #===========================================================
  #Writes the values of the input file to various dictionaries
  #===========================================================
  for line in lines:
    ss = line.split()
    try:
      Z=int(float(ss[int(az)])+0.0001)         #Proton number
      N=int(float(ss[int(an)])+0.0001)         #Neutron number
      BE=float(ss[int(ae)])                    #Binding energy
      rms_proton_radius=float(ss[int(LNrp)])   #Root-mean-square proton radius
      rms_neutron_radius=float(ss[int(LNrn)])  #Root-mean-square neutron radius
      ProtonPair=float(ss[int(PPG)])           #Proton pairing gap
      NeutronPair=float(ss[int(NPG)])          #Neutron pairing gap
      #---------------------------
      #Calculates the neutron skin
      #---------------------------
      neutron_skin = decimal.Decimal(rms_neutron_radius - rms_proton_radius).quantize(decimal.Decimal('.000001'), rounding=decimal.ROUND_UP)
      #======================================================================================
      #Looks for the smallest binding energy between spherical, prolate, and oblate solutions
      #(Only works if binding energy values are negative)
      #======================================================================================
      if ((N,Z) in BEdic):
        if(BE < BEdic[(N,Z)]):
          neutron_skin_data[(N,Z)]=neutron_skin
          ProtonPairingGap[(N,Z)]=ProtonPair
          NeutronPairingGap[(N,Z)]=NeutronPair       
          BEdic[(N,Z)]=BE
          DATA_Line[(N,Z)]=line
        else:
          continue
      else:
        neutron_skin_data[(N,Z)]=neutron_skin
        ProtonPairingGap[(N,Z)]=ProtonPair
        NeutronPairingGap[(N,Z)]=NeutronPair
        BEdic[(N,Z)]=BE             #The binding energy is added to the dictionary BEdic
        DATA_Line[(N,Z)]=line       #The current line is added to the dictionary DATALine
    except ValueError:
      continue                      #N,Z, or, BE are not numbers    
  f1.close()
  #===========================================================
  #Data for even-even nuclei added to the list Mass_Table_Data
  #===========================================================
  for key in BEdic.keys():
    N=key[0]; Z=key[1]
    #-----------------------------------
    #Adds most of the nuclei to the list
    #-----------------------------------
    if (((N,Z) in BEdic) and ((N-2,Z) in BEdic) and ((N,Z-2) in BEdic)):  #Makes sure data exists for calculation
      #--------------------------------------------------
      #Calculates S_p, S_{2p}, S_n, S_{2n}, and Q_{alpha}
      #--------------------------------------------------
      BE1 = BEdic[(N,Z)]                     #BE(N,Z)
      BE2 = BEdic[(N,Z-2)]                   #BE(N,Z-2)
      PPG1 = ProtonPairingGap[(N,Z)]         #ProtonPairingGap(N,Z)
      PPG2 = ProtonPairingGap[(N,Z-2)]       #ProtonPairingGap(N,Z-2)
      BEavg = (BE1 + BE2)/2.0
      PPGavg = (PPG1 + PPG2)/2.0
      BEoddZ = BEavg + PPGavg                #BE(N,Z-1)
      s1p = decimal.Decimal(BEoddZ - BEdic[(N,Z)]).quantize(decimal.Decimal('.000001'), rounding=decimal.ROUND_UP)
      s2p = decimal.Decimal(BEdic[(N,Z-2)] - BEdic[(N,Z)]).quantize(decimal.Decimal('.000001'), rounding=decimal.ROUND_UP)
      BE3 = BEdic[(N-2,Z)]                   #BE(N-2,Z)
      NPG1 = NeutronPairingGap[(N,Z)]        #NeutronPairingGap(N,Z)
      NPG2 = NeutronPairingGap[(N-2,Z)]      #NeutronPairingGap(N-2,Z)
      BEavgN = (BE1 + BE3)/2.0
      NPGavg = (NPG1 + NPG2)/2.0
      BEoddN = BEavgN + NPGavg               #BE(N-1,Z)
      s1n = decimal.Decimal(BEoddN - BEdic[(N,Z)]).quantize(decimal.Decimal('.000001'), rounding=decimal.ROUND_UP)
      s2n = decimal.Decimal(BEdic[(N-2,Z)] - BEdic[(N,Z)]).quantize(decimal.Decimal('.000001'), rounding=decimal.ROUND_UP)
      Qalpha = decimal.Decimal(28.3 + BEdic[(N,Z)] - BEdic[(N-2,Z-2)]).quantize(decimal.Decimal('.000001'), rounding=decimal.ROUND_UP)
      #-------------------------------------
      #Neutron skin values calculated before
      #-------------------------------------  
      neutron_skin = neutron_skin_data[(N,Z)]
      #--------------------------------
      #Adds all of the data to the list
      #--------------------------------
      Mass_Table_Data.append((Z,N,s2p,s2n,neutron_skin,Qalpha,s1p,s1n,DATA_Line[(N,Z)]))
    #--------------------------------------------------------------------------------------------
    #Ensures that proton and neutron indexing doesn't prevent certain nuclei from making the list
    #--------------------------------------------------------------------------------------------
    elif (((N,Z-2) not in BEdic) and ((N-2,Z) not in BEdic)):
      #----------------------------------------------------------
      #Not enough data for S_p, S_n, S_{2p}, S_{2n}, or Q_{alpha}
      #----------------------------------------------------------
      s1p = "No_Data"
      s1n = "No_Data"
      s2p = "No_Data"
      s2n = "No_Data"
      Qalpha = "No_Data"
      #-------------------------------------
      #Neutron skin values calculated before
      #-------------------------------------
      neutron_skin = neutron_skin_data[(N,Z)]
      #--------------------------------
      #Adds all of the data to the list
      #--------------------------------
      Mass_Table_Data.append((Z,N,s2p,s2n,neutron_skin,Qalpha,s1p,s1n,DATA_Line[(N,Z)]))
    #---------------------------------------------------------------------------------
    #Ensures that neutron indexing doesn't prevent certain nuclei from making the list
    #---------------------------------------------------------------------------------
    elif ((N-2,Z) not in BEdic):
      #--------------------------------------------
      #Not enough data for S_n, S_{2n} or Q_{alpha}
      #--------------------------------------------
      s1n = "No_Data"
      s2n = "No_Data"
      Qalpha = "No_Data"
      #--------------------------
      #Calculating S_p and S_{2p}
      #--------------------------
      BE1 = BEdic[(N,Z)]                   #BE(N,Z)
      BE2 = BEdic[(N,Z-2)]                 #BE(N,Z-2)
      PPG1 = ProtonPairingGap[(N,Z)]       #ProtonPairingGap(N,Z)
      PPG2 = ProtonPairingGap[(N,Z-2)]     #ProtonPairingGap(N,Z-2)
      BEavg = (BE1 + BE2)/2.0
      PPGavg = (PPG1 + PPG2)/2.0
      BEoddZ = BEavg + PPGavg              #BE(N,Z-1)
      s1p = decimal.Decimal(BEoddZ - BEdic[(N,Z)]).quantize(decimal.Decimal('.000001'), rounding=decimal.ROUND_UP)
      s2p = decimal.Decimal(BEdic[(N,Z-2)] - BEdic[(N,Z)]).quantize(decimal.Decimal('.000001'), rounding=decimal.ROUND_UP)
      #-------------------------------------
      #Neutron skin values calculated before
      #-------------------------------------     
      neutron_skin = neutron_skin_data[(N,Z)]
      #--------------------------------
      #Adds all of the data to the list
      #--------------------------------
      Mass_Table_Data.append((Z,N,s2p,s2n,neutron_skin,Qalpha,s1p,s1n,DATA_Line[(N,Z)]))
    #--------------------------------------------------------------------------------
    #Ensures that proton indexing doesn't prevent certain nuclei from making the list
    #--------------------------------------------------------------------------------
    elif ((N,Z-2) not in BEdic):
      #--------------------------------------------
      #Not enough data for S_p, S_{2p} or Q_{alpha}
      #--------------------------------------------
      s1p = "No_Data"
      s2p = "No_Data"
      Qalpha = "No_Data"
      #--------------------------
      #Calculating S_n and S_{2n}
      #--------------------------
      BE1 = BEdic[(N,Z)]                    #BE(N,Z)
      BE2 = BEdic[(N-2,Z)]                  #BE(N-2,Z)
      NPG1 = NeutronPairingGap[(N,Z)]       #NeutronPairingGap(N,Z)
      NPG2 = NeutronPairingGap[(N-2,Z)]     #NeutronPairingGap(N-2,Z)
      BEavg = (BE1 + BE2)/2.0
      NPGavg = (NPG1 + NPG2)/2.0
      BEoddN = BEavg + NPGavg               #BE(N-1,Z)
      s1n = decimal.Decimal(BEoddN - BEdic[(N,Z)]).quantize(decimal.Decimal('.000001'), rounding=decimal.ROUND_UP)
      s2n = decimal.Decimal(BEdic[(N-2,Z)] - BEdic[(N,Z)]).quantize(decimal.Decimal('.000001'), rounding=decimal.ROUND_UP)
      #-------------------------------------
      #Neutron skin values calculated before
      #-------------------------------------     
      neutron_skin = neutron_skin_data[(N,Z)]
      #--------------------------------
      #Adds all of the data to the list
      #--------------------------------
      Mass_Table_Data.append((Z,N,s2p,s2n,neutron_skin,Qalpha,s1p,s1n,DATA_Line[(N,Z)]))
    else:
      continue
  #---------------------------------------------------------------------------------
  #Sorts the data on the list to be in order of increasing proton and neutron number
  #---------------------------------------------------------------------------------
  Mass_Table_Data.sort()
  #==========================================================
  #Writes and properly formats the titles for the output file
  #==========================================================
  ground_state_data.write(' {:10} {:6} {:6} {:7} {:23} {:20} {:20} {:22} {:26} {:26} {:31} {:31} {:31} {:36} {:23} {:21} {:20} {:20} {:24} {:22} {:14} {:16} {:13} {:16} {:21}  \n'.format('Symbol', 'Z',
                          'N', 'A', 'Binding_Energy_(MeV)', 'Quad_Def_Beta2_P', 'Quad_Def_Beta2_N', 'Quad_Def_Beta2_total', 'Quad_Moment_Q2_P_(fm^2)', 'Quad_Moment_Q2_N_(fm^2)',
                          'Quad_Moment_Q2_total_(fm^2)', 'Octupole_Moment_Q3_P_(fm^3)', 'Octupole_Moment_Q3_N_(fm^3)', 'Octupole_Moment_Q3_total_(fm^3)', 'Pairing_gap_P_(MeV)', 'Pairing_gap_N_(MeV)',
                          'RMS_radius_P_(fm)', 'RMS_radius_N_(fm)', 'RMS_radius_total_(fm)', 'Charge_Radius_(fm)', 'S_p_(MeV)', 'S_{2p}_(MeV)', 'S_n_(MeV)', 'S_{2n}_(MeV)', 'Q_{alpha}_(MeV)'))
  #--------------------------------------------------------
  #Initialize 2-proton drip line (location at Z=98 for now)
  #--------------------------------------------------------
  protondrip=130  
  #===========================================
  #Writes the even-even list to an output file
  #===========================================
  for item in Mass_Table_Data:
    #---------------------------------
    #Proton, neutron, and mass numbers
    #---------------------------------
    Z = Mass_Table_Data[i][0]
    N = Mass_Table_Data[i][1]
    A = int(Z) + int(N)
    #------------
    #Element name
    #------------
    elementlabel = ElementName(Z) 
    #-----------------------------------------------
    #Breaks up each line of data from the input file
    #-----------------------------------------------
    Modline = item[8]
    sp = Modline.split()
    #--------------
    #Binding energy
    #--------------
    BE = sp[3]
    #-----------------------
    #Quadrupole deformations
    #-----------------------
    quad_def_beta2_P = sp[4]
    quad_def_beta2_N = sp[5]
    quad_def_beta2_T = sp[6]
    #------------------
    #Quadrupole moments
    #------------------
    quad_moment_Q2_P = sp[7]
    quad_moment_Q2_N = sp[8]
    quad_moment_Q2_T = sp[9]
    #----------------
    #Octupole moments
    #----------------
    oct_moment_Q3_P = sp[10]
    oct_moment_Q3_N = sp[11]
    oct_moment_Q3_T = sp[12]
    #------------
    #Pairing gaps
    #------------
    proton_pairing_gap = sp[13]
    neutron_pairing_gap = sp[14]
    #---------
    #RMS radii
    #---------
    rms_radius_P = sp[15]
    rms_radius_N = sp[16]
    rms_radius_T = sp[17]
    #-------------
    #Charge radius
    #-------------
    charge_radius = sp[18]
    #----------------------------------------------------------------------
    #Testing purposes only; outputs all data and does no drip line analysis
    #----------------------------------------------------------------------
    #ground_state_data.write('   {:8} {:6} {:6} {:9} {:23} {:20} {:22} {:20} {:26} {:30} {:27} {:31} {:34} {:34} {:23} {:21} {:20} {:22} {:22} {:17} {:15} {:15} {:14} {:18} {:21} \n'.format(
    #                             elementlabel, str(Z), str(N), str(A), str(BE).rjust(13, ), str(quad_def_beta2_P).rjust(10, ), str(quad_def_beta2_N).rjust(10, ), str(quad_def_beta2_T).rjust(10, ),
    #                             str(quad_moment_Q2_P).rjust(12, ), str(quad_moment_Q2_N).rjust(12, ), str(quad_moment_Q2_T).rjust(12, ), str(oct_moment_Q3_P).rjust(12, ),
    #                             str(oct_moment_Q3_N).rjust(12, ), str(oct_moment_Q3_T).rjust(12, ), str(proton_pairing_gap).rjust(10, ), str(neutron_pairing_gap).rjust(10, ),
    #                             str(rms_radius_P).rjust(10, ), str(rms_radius_N).rjust(10, ), str(rms_radius_T).rjust(10, ), str(charge_radius).rjust(10, ),
    #                             str(Mass_Table_Data[i][6]).rjust(10, ), str(Mass_Table_Data[i][2]).rjust(10, ), str(Mass_Table_Data[i][7]).rjust(10, ), str(Mass_Table_Data[i][3]).rjust(10, ),
    #                             str(Mass_Table_Data[i][5]).rjust(10, )))
    #i=i+1
    #j=j+1
    #k=k+1    
    #----------------------------
    #Adds data to the output file
    #----------------------------
    try:
      #----------------------------------
      #Adds Helium isotopes to the output
      #----------------------------------
      if((Mass_Table_Data[i][0]==2) and (Mass_Table_Data[i][1]<=6)):
        protondrip=2
        #---------------------------------
        #Writes the data to an output file
        #---------------------------------
        ground_state_data.write('   {:8} {:6} {:6} {:9} {:23} {:20} {:22} {:20} {:26} {:30} {:27} {:23} {:21} {:20} {:22} {:22} {:17} {:15} {:15} {:14} {:18} {:21} \n'.format(
                                 elementlabel, str(Z), str(N), str(A), str(BE).rjust(13, ), str(quad_def_beta2_P).rjust(10, ), str(quad_def_beta2_N).rjust(10, ), str(quad_def_beta2_T).rjust(10, ),
                                 str(quad_moment_Q2_P).rjust(12, ), str(quad_moment_Q2_N).rjust(12, ), str(quad_moment_Q2_T).rjust(12, ), str(proton_pairing_gap).rjust(10, ),
                                 str(neutron_pairing_gap).rjust(10, ), str(rms_radius_P).rjust(10, ), str(rms_radius_N).rjust(10, ), str(rms_radius_T).rjust(10, ), str(charge_radius).rjust(10, ),
                                 str(Mass_Table_Data[i][6]).rjust(10, ), str(Mass_Table_Data[i][2]).rjust(10, ), str(Mass_Table_Data[i][7]).rjust(10, ), str(Mass_Table_Data[i][3]).rjust(10, ),
                                 str(Mass_Table_Data[i][5]).rjust(10, )))
        i=i+1
        j=j+1
        k=k+1
      #===========================================================================================================================================
      #Finds the case for the proton drip line (-s2p to +s2p) or proceeds if proton drip line has already been reached for this particular element
      #===========================================================================================================================================
      elif(((Mass_Table_Data[i][2] > 0) and (Mass_Table_Data[j][2] < 0)) or (pdripreached==1)):
        #---------------------------------------------------------------------
        #Removes all nuclei below Z=100 (since we don't have those values yet)
        #---------------------------------------------------------------------
        if (Z < 100):
          i=i+1
          j=j+1
          k=k+1
          continue
        #=============================================
        #Removes early entries of the proton drip line
        #=============================================
        if((Mass_Table_Data[i][2] > 0) and (Mass_Table_Data[j][2] < 0)):
          if(N < protondrip):   #Pre-entrant proton drip line nucleus
            i=i+1
            j=j+1
            k=k+1
            continue
          else:
            protondrip = N
        pdripreached=1  #Proton drip line has been found and passed
        #=====================================
        #Tests for cases where s2p>0 and s2n>0
        #=====================================
        if((Mass_Table_Data[i][2] > 0) and (Mass_Table_Data[i][3] > 0)):
          #---------------------------------
          #Writes the data to an output file
          #---------------------------------
          ground_state_data.write('   {:8} {:6} {:6} {:9} {:23} {:20} {:22} {:20} {:26} {:30} {:27} {:31} {:34} {:34} {:23} {:21} {:20} {:22} {:22} {:17} {:15} {:15} {:14} {:18} {:21} \n'.format(
                                       elementlabel, str(Z), str(N), str(A), str(BE).rjust(13, ), str(quad_def_beta2_P).rjust(10, ), str(quad_def_beta2_N).rjust(10, ), str(quad_def_beta2_T).rjust(10, ),
                                       str(quad_moment_Q2_P).rjust(12, ), str(quad_moment_Q2_N).rjust(12, ), str(quad_moment_Q2_T).rjust(12, ), str(oct_moment_Q3_P).rjust(12, ),
                                       str(oct_moment_Q3_N).rjust(12, ), str(oct_moment_Q3_T).rjust(12, ), str(proton_pairing_gap).rjust(10, ), str(neutron_pairing_gap).rjust(10, ),
                                       str(rms_radius_P).rjust(10, ), str(rms_radius_N).rjust(10, ), str(rms_radius_T).rjust(10, ), str(charge_radius).rjust(10, ),
                                       str(Mass_Table_Data[i][6]).rjust(10, ), str(Mass_Table_Data[i][2]).rjust(10, ), str(Mass_Table_Data[i][7]).rjust(10, ), str(Mass_Table_Data[i][3]).rjust(10, ),
                                       str(Mass_Table_Data[i][5]).rjust(10, )))
          i=i+1
          j=j+1
          k=k+1
        #------------------------------------------------------------------
        #Neutron drip line reached; search continues for reentrant isotopes
        #------------------------------------------------------------------
        elif((Mass_Table_Data[i][0]==Mass_Table_Data[k][0]) or (Mass_Table_Data[i][0]==120)):
          i=i+1
          j=j+1
          k=k+1
        #----------------------------------------------------------
        #All neutron drip lines passed; move on to the next element
        #----------------------------------------------------------
        else:
          pdripreached=0
          i=i+1
          j=j+1
          k=k+1
      else:
        i=i+1
        j=j+1
        k=k+1
    except IndexError:
      continue
  #----------------------
  #Closes the output file
  #----------------------
  ground_state_data.close()
#===========
#User Inputs
#===========
if __name__ == "__main__":    #If this file is imported into another, this will prevent this file from running and only use the function called
  EDFs = ['SKMS', 'SKP', 'SLY4', 'SV-MIN', 'UNEDF0', 'UNEDF1', 'UNEDF2', 'UNEDF1-SO']
  type_of_calculation = "masstable"#input('Enter the type of mode used for this calculation ("masstable", "PES", "minima", or "ground"): ')    #Type of calculation
  functional = "SLY4"
  #LowestFileID = input('Enter the smallest FileID number (thoout_"FileID".dat): ')  #Smallest numerial label for the input file:  thoout_"FileID".dat
  #=====================================
  #Input files are from "masstable" mode
  #=====================================
  if (type_of_calculation == 'masstable'):
    os.chdir("./")
    os.system("shopt -s extglob\n"+"rm HFBTHO*.dat")
    number_of_shells = read_hfbtho_NAMELIST()  #Reads from file hfbtho_NAMELIST.dat to get the number of shells used for this calculation
    Lowest_File_ID = 1  #The first thoout file to be read
    (lower_proton_number, upper_proton_number, Highest_File_ID) = read_hfbtho_MASSTABLE()  #Reads from file hfbtho_MASSTABLE.dat to get the lower and upper proton numbers and the highest file ID
    Data_File_Out = "HFBTHOv300_"+functional+"_All_Data_"+number_of_shells+"_shells_z"+lower_proton_number+"to"+upper_proton_number+"-masstable.dat"  #Output file for Read_HFBTHO_Masstable_Output
    Read_HFBTHO_Masstable_and_PES_Output(functional,Data_File_Out,Lowest_File_ID,Highest_File_ID)
  #====================================
  #Input files are from "dripline" mode
  #====================================
  elif (type_of_calculation == 'dripline'):
    LowestTeamID = input('Enter the smallest TeamID number (thoout_"TeamID"_000001.dat): ')  #Smallest numerical label for the team part:  thoout_"TeamID"_000001.dat
    HighestTeamID = input('Enter the largest TeamID number (thouut_"TeamID"_000001.dat): ')  #Largest numerical label for the team part:   thoout_"TeamID"_000001.dat
    for i in range(int(LowestTeamID),int(HighestTeamID)+1):
      try:
        TeamID = str(i).rjust(3, '0')     #Proper formatting for the team number of each file
        for j in range(int(LowestFileID),int(HighestFileID)+1):
          try:
            FileID = str(j).rjust(6, '0')   #Proper formatting for the number in each file
            #---------------------------------------------------------
            #Parameters for the function Read_HFBTHO_Output_File above
            #---------------------------------------------------------
            Data_File_In = "thoout_"+TeamID+"_"+FileID+".dat"
            Data_File_Out = "HFBTHOv300_"+functional+"_All_Data.dat"
            Read_HFBTHO_Dripline_Output(functional,Data_File_In,Data_File_Out,LowestTeamID,TeamID,LowestFileID,FileID)
            #---------------------
            #Looks for the file ID
            #---------------------
            j=j+1
          except IOError:
            continue
        #--------------------------
        #Looks for the next team ID
        #--------------------------
        i=i+1
      except IOError:
        continue
  #===============================
  #Input files are from "PES" mode
  #===============================
  elif (type_of_calculation == 'PES'):
    os.chdir(functional+"/")
    #os.chdir("Maxwell_Test/")
    number_of_shells = read_hfbtho_NAMELIST()  #Reads from file hfbtho_NAMELIST.dat to get the number of shells used for this calculation
    Lowest_File_ID = 1  #The first thoout file to be read
    (lower_proton_number, upper_proton_number, Highest_File_ID) = read_hfbtho_PES()  #Reads from file hfbtho_PES.dat to get the lower and upper proton numbers and highest file ID
    Data_File_Out = "HFBTHOv300_"+functional+"_All_Data_"+number_of_shells+"_shells_z"+lower_proton_number+"to"+upper_proton_number+"-PES.dat"  #Output file for Read_HFBTHO_Masstable_Output
    Read_HFBTHO_Masstable_and_PES_Output(functional,Data_File_Out,Lowest_File_ID,Highest_File_ID)
  #===================================================================================
  #Input file is an "All_Data" file, we want the oblate, prolate, and spherical minima
  #===================================================================================
  elif (type_of_calculation == 'minima'):
    os.chdir(functional+"/")
    number_of_shells = read_hfbtho_NAMELIST()  #Reads from file hfbtho_NAMELIST.dat to get the number of shells used for this calculation
    Total_Input_File = "HFBTHOv300_"+functional+"_All_Data_"+number_of_shells+"_shells_z98to120-PES.dat"
    Minima_Output_File = "HFBTHOv300_"+functional+"_All_Minima_"+number_of_shells+"_shells_z98to120-PES.dat"
    All_Minima(Total_Input_File,Minima_Output_File)  
  #========================================
  #Input file is the output file from above
  #========================================
  elif (type_of_calculation == 'ground'):
    #os.chdir(functional+"/")
    number_of_shells = read_hfbtho_NAMELIST()  #Reads from file hfbtho_NAMELIST.dat to get the number of shells used for this calculation
    Total_Input_File = "HFBTHOv300_"+functional+"_All_Data_"+number_of_shells+"_shells_z98to120-masstable-norestartfile-kickoff-UNEDF1_pairing.dat"
    #GS_Output_File = "HFBTHOv300_"+functional+"_Ground_State_Data_"+number_of_shells+"_shells_z98to120-masstable-norestartfile-kickoff.dat"
    GS_Output_File = "HFBTHOv300_"+functional+"_Ground_State_Data_"+number_of_shells+"_shells_z98to120-UNEDF1_pairing-Drip_Line.dat"
    HFBTHOv300_Ground_State_Mass_Table(Total_Input_File,GS_Output_File,0,1,3,15,16,13,14)
  #===================================
  #Error in user input, the code exits
  #===================================
  else:
    sys.exit('Type of calculation is not one of the acceptable options, please try again.')
