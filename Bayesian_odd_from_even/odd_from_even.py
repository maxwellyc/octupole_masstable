import itertools as it
import math

def isNum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def read_file( file_loc ):
    f1 = open(file_loc)
    l1 = f1.readlines()
    S1n = {}; S2n = {}; BE = {}; gap_z = {}; gap_n = {}; Sym = {}
    S1p = {}; S2p = {}; Qa = {}
    # Input file WITHOUT first column as element name
    for line in l1:
        print (line)
        try:
            ss = line.split()
            Z = int(ss[0]); N = int(ss[1]); A = int(ss[2])
            if isNum(ss[3]):  BE[(Z,N)] = round( float(ss[3]),6)   #BE
            if isNum(ss[13]): gap_z[(Z,N)] = float(ss[13])  #proton pairing gap
            if isNum(ss[14]): gap_n[(Z,N)] = float(ss[14])  #neutron pairing gao
        except (ValueError, IndexError): continue
    # Input file WITH first column as element name
#    for line in l1:
#        try:
#            ss = line.split()
#            Z = int(ss[1]); N = int(ss[2]); A = int(ss[3])
#            Sym[(Z,N)] = ss[0]
#            if isNum(ss[4]):  BE[(Z,N)] = round( float(ss[4]),6)   #BE
#            if isNum(ss[13]): gap_z[(Z,N)] = float(ss[13])  #proton pairing gap
#            if isNum(ss[14]): gap_n[(Z,N)] = float(ss[14])  #neutron pairing gao
#            if isNum(ss[20]): S1p[(Z,N)] = round( float(ss[20]),6)  #S1p
#            if isNum(ss[21]): S2p[(Z,N)] = round( float(ss[21]),6)  #S2p
#            if isNum(ss[22]): S1n[(Z,N)] = round( float(ss[22]),6)  #S1n
#            if isNum(ss[23]): S2n[(Z,N)] = round( float(ss[23]),6)  #S2n
#            if isNum(ss[24]): Qa[(Z,N)]  = round( float(ss[24]),6)  #S2n
#        except (ValueError, IndexError): continue

    # proton pairing gap of odd:
    for Z in range(2,121,2):
        for N in range(2,301,2):
            if (Z,N) in gap_z and (Z,N-2) in gap_z:
                gap_z[(Z,N-1)] = 0.5 * ( gap_z[(Z,N)] + gap_z[(Z,N-2)] )
    
    # BE of odd nuclei
    for Z in range(2,121,2):
        for N in range(2,301,2):
            #evenZ-evenN data stored in previous stage, not need to calculate
            #evenZ-oddN
            if (Z,N) in BE and (Z,N-2) in BE and (Z,N) in gap_n and (Z,N-2) in gap_n:
                BE[(Z,N-1)] = round(0.5 * (   BE[(Z,N)] + BE[(Z,N-2)] + gap_n[(Z,N)] + gap_n[(Z,N-2)]    ),6)
            #oddZ-evenN
            if (Z,N) in BE and (Z-2,N) in BE and (Z,N) in gap_z and (Z-2,N) in gap_z:
                BE[(Z-1,N)] = round(0.5 * (   BE[(Z,N)] + BE[(Z-2,N)] + gap_z[(Z,N)] + gap_z[(Z-2,N)]   ),6)
    #oddZ-oddN
    for Z in range(2,121,2):
        for N in range(2,301,2):
            if (Z,N-1) in BE and (Z-2,N-1) in BE and (Z,N-1) in gap_z and (Z-2,N-1) in gap_z:
                BE[(Z-1,N-1)] = round(0.5 * (   BE[(Z,N-1)] + BE[(Z-2,N-1)] + gap_z[(Z,N-1)] + gap_z[(Z-2,N-1)]   ),6)


    # S1/2n of odd nuclei
    for Z in range(2,121):
        for N in range(2,301):
            if True: #Z%2 or N%2:
                if (Z,N) in BE and (Z,N-1) in BE:
                    S1n[(Z,N)] = -round( BE[(Z,N)] - BE[(Z,N-1)] , 6)
                if (Z,N) in BE and (Z,N-2) in BE:
                    S2n[(Z,N)] = -round( BE[(Z,N)] - BE[(Z,N-2)] , 6)
    # S1/p of odd nuclei
    for Z in range(2,121):
        for N in range(2,301):
            if True:#Z%2 or N%2:
                if (Z,N) in BE and (Z-1,N) in BE:
                    S1p[(Z,N)] = -round( BE[(Z,N)] - BE[(Z-1,N)] , 6)
                if (Z,N) in BE and (Z-2,N) in BE:
                    S2p[(Z,N)] = -round( BE[(Z,N)] - BE[(Z-2,N)] , 6)

    # Q alpha of odd nuclei:
    for Z in range(2,121):
        for N in range(2,301):
            if True:#Z%2 or N%2:
                if (Z,N) in BE and (Z-2,N-2) in BE:
                    Qa[(Z,N)] = round( 28.3 + BE[(Z,N)] - BE[(Z-2,N-2)] , 6)


    edfName = file_loc.split("_")[0]
    outputFile = open("./"+edfName+"_octupole_separation_filtered.dat", "w")
    outputStr = "  Symbol     Z      N      A      Binding_Energy_(MeV)    S_p_(MeV)      S_{2p}_(MeV)     S_n_(MeV)     S_{2n}_(MeV)     Q_{alpha}_(MeV) \n"
    for Z in range(2,121):
        for N in range(2,301):
            if (Z,N) in BE:
                outputStr = outputStr + "    " + "elem".ljust(9) + str(Z).ljust(7) + str(N).ljust(7) + str(Z+N).ljust(10) + \
                            str(BE[(Z,N)]).ljust(22)
                if (Z,N) in S1p: outputStr = outputStr + str(S1p[(Z,N)]).ljust(15)
                else: outputStr  = outputStr + "No_Data".ljust(15)

                if (Z,N) in S2p: outputStr = outputStr + str(S2p[(Z,N)]).ljust(17)
                else: outputStr  = outputStr + "No_Data".ljust(17)

                if (Z,N) in S1n: outputStr = outputStr + str(S1n[(Z,N)]).ljust(14)
                else: outputStr  = outputStr + "No_Data".ljust(14)

                if (Z,N) in S2n: outputStr = outputStr + str(S2n[(Z,N)]).ljust(17)
                else: outputStr  = outputStr + "No_Data".ljust(17)

                if (Z,N) in Qa:  outputStr = outputStr + str(Qa[(Z,N)]).ljust(19)
                else: outputStr  = outputStr + "No_Data".ljust(19)
                outputStr = outputStr + "\n"
            else: continue
    outputFile.write(outputStr)


#read_file("UNEDF2_even-even_nuclei-no_drip_lines.dat")
#read_file("UNEDF0_even-even_nuclei-no_drip_lines.dat")
#read_file("UNEDF1_even-even_nuclei-no_drip_lines.dat")
#read_file("SKMS_even-even_nuclei-no_drip_lines.dat")
#read_file("SKP_even-even_nuclei-no_drip_lines.dat")
#read_file("SLY4_even-even_nuclei-no_drip_lines.dat")
#read_file("SV-MIN_even-even_nuclei-no_drip_lines.dat")
#read_file("UNEDF2_gs_HFBTHOv300.dat")

read_file("../code/sorted_masstable_output/new_data/SLY4_groundstate_octupole_filtered.dat")
