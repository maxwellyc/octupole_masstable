import numpy as np
import re
import itertools as it
import math
import os
from operator import itemgetter
from functools import reduce

def read_file(folder_loc,file_name):
    f1 = open(folder_loc+file_name, "r")
    l1 = f1.readlines()
    return l1

class run_stats:
    def __init__(self, num_threads, tot_time, iter_steps, tlog):
        self.num_threads = num_threads
        self.tot_time = tot_time
        self.iter_steps = iter_steps
        self.tlog = tlog

class nuc_stats:
    def __init__(self, edf_name, ehfb_ln, nucleus, beta2_in, beta3_in, q20_in, q30_in, q20_out, q30_out, tlog, interr):
        self.edf_name = edf_name
        self.ehfb_ln = ehfb_ln
        self.nucleus = nucleus
        self.beta2_in = beta2_in
        self.beta3_in = beta3_in
        self.q20_in = q20_in
        self.q30_in = q30_in
        self.q20_out = q20_out
        self.q30_out = q30_out
        self.tlog = tlog
        self.interr = interr

def find_stats(folder_loc="./"):
    lsLoc = os.listdir(folder_loc)
    num_threads ={}; edf_name ={}; tot_time ={}; iter_steps ={}; tlog = []
    ehfb_ln = {}; nucleus = {}; beta2_in = {}; beta3_in = {}
    q20_in = {}; q30_in = {}; q20_out = {}; q30_out = {}; interr = []
    for ff in lsLoc:
        try:
            if "thoout_" in ff:
                f_index = re.findall('\d+',ff)
                f_index = int(f_index[0])
                tlog.append(f_index)
                edf_skip = 0; time_skip = 0; deform_skip = 0
                l1 = read_file(folder_loc,ff)
                num_threads[(f_index)] = 1  # defaults to 1 if no OpenMP enabled
                print (len(tlog))
                for line in l1:
                    if "Multi-threading" in line:
                        read_num = re.findall('\d+',line)
                        num_threads[(f_index)] = int(read_num[0])
                        continue
                    if "Nucleus: " in line:
                        read_num = re.findall('\d+',line)
                        nucleus[(f_index)] = (int(read_num[2]),int(read_num[1])) #Z, N of nucleus
                        continue
                    if "functional" in line and edf_skip != 1:
                        edf_name[(f_index)] = (line.split())[0]
                        edf_skip = 1
                        continue
                    if "b2_ws" in line:
                        read_num = line.split()
                        beta2_in[(f_index)] = round(float(read_num[-1]),3)  # beta2 input
                        continue
                    if "b3_ws" in line:
                        read_num = line.split()
                        beta3_in[(f_index)] = round(float(read_num[-1]),3)  # beta2 input
                        continue
                    if "Initial constraints Q20, Q30" in line:
                        read_num = line.split()
                        q20_in[(f_index)] = round(float(read_num[-2]),6)
                        q30_in[(f_index)] = round(float(read_num[-1]),6)
                    if "Total CPU time" in line and time_skip == 0:
                        time_skip = 1
                        continue
                    if "Total CPU time" in line and time_skip == 1:
                        read_num = line.split()
                        if "*****" in line:
                            tot_time[(f_index)] = 51    # This temporarily solves the time = "********" issue
                        else:
                            tot_time[(f_index)]  = round(float(read_num[-2])*60,3)    #total cpu time in seconds
                        continue
                    if "after" in line and time_skip == 1:
                        read_num = re.findall('\d+',line)
                        iter_steps[(f_index)] = int(read_num[0])
                    # IF YOU WANT TO EXCLUDE CALCULATION WHERE THERE'S
                    # LOW ACCURACY ITERATIONS, UNCOMMENT BELOW
#                    if "Low accuracy" in line:
#                        iter_steps[(f_index)] = 1001
#                        interr.append(f_index)
#                        break
                    if "interrupt" in line and time_skip == 1:
                        interr.append(f_index)
                        continue
                    if "tEnergy" in line and "LN" in line:
                        read_num = line.split()
                        ehfb_ln[(f_index)] = round(float(read_num[-1]),6)
                        continue
                    if "With Lipkin-Nogami Corrections" in line:
                        deform_skip = 1
                        continue
                    if deform_skip == 1:
                        if "quadrupole moment" in line:
                            read_num = line.split()
                            q20_out[(f_index)] = float(read_num[-1])
                            continue
                        if "octupole moment" in line:
                            read_num = line.split()
                            q30_out[(f_index)] = float(read_num[-1])
                            deform_skip = 0
                            continue
        except:
            continue
    rs1 = run_stats(num_threads, tot_time, iter_steps, tlog)
    ns1 = nuc_stats(edf_name, ehfb_ln, nucleus, beta2_in, beta3_in, q20_in, q30_in, q20_out, q30_out, tlog, interr)
    return rs1,ns1

def sort_ehfbln( ns, rs, hfb_accuracy, grid_s2, grid_s3):
    # For each nuclei, find lowest hfb_ln energy
    nuc_set = set((ns.nucleus).values())
    nuc_row = {}; gs_ehfb = {}; gs_row = {} #key is row #
    out_str = "(Z,N)".rjust(12) + "E_HFB_LN(MeV)".rjust(18) + "beta2_in".rjust(15) + "beta3_in".rjust(15) +\
                "Q20_out".rjust(15) + "Q30_out".rjust(15) + "file#".rjust(8) + "iter#".rjust(10) + "\n"
    for (z,n) in sorted(nuc_set):
        nuc_row[(z,n)] = []
        gs_ehfb[(z,n)] = 0
        gs_row[(z,n)] = []
        # this loop groups the same nucleus from difference rows together, and find its g.s. energy from within
        # add a filter of q20 and q30 steps to check if subgroup of grid suffice
        for rows in ns.nucleus.keys():
            beta2_check = 100 * ns.beta2_in[(rows)]
            beta3_check = 100 * ns.beta3_in[(rows)]
            beta2_check = beta2_check % 10      # beta*_check == 0 select only 0.1 step sized grid points
            beta3_check = beta3_check % 10
            if grid_s2 == 1: beta2_check = 0
            if grid_s3 == 1: beta3_check = 0
            if ns.nucleus[(rows)] == (z,n) and beta2_check == 0 and beta3_check == 0:
                nuc_row[(z,n)].append(rows)
        # find actual minimum g.s.
        for f_ind in nuc_row[(z,n)]:
            # add the statement "and not (f_ind in ns.interr)" to keep only converged results
            if f_ind in ns.ehfb_ln and not (f_ind in ns.interr):
                if ns.ehfb_ln[(f_ind)] < gs_ehfb[(z,n)]:
                    gs_ehfb[(z,n)] = ns.ehfb_ln[(f_ind)]
        # round g.s. energy and include all that's within certain decimal accuracy defined by var: hfb_accuracy
        for f_ind in nuc_row[(z,n)]:
            # add the statement "and not (f_ind in ns.interr)" to keep only converged results
            if f_ind in ns.ehfb_ln and not (f_ind in ns.interr):
                if round(ns.ehfb_ln[(f_ind)],hfb_accuracy) <= round(gs_ehfb[(z,n)],hfb_accuracy):
                    gs_row[(z,n)].append(f_ind)

        f_iter_min = 2000
        for f_ind in sorted(gs_row[(z,n)]):
            if f_ind in ns.interr:
                print (f"Warning: Ground state energy of file {f_ind} did not converge")
            out_str += str(ns.nucleus[(f_ind)]).rjust(12) + str(ns.ehfb_ln[(f_ind)]).rjust(18) +\
                str(ns.beta2_in[(f_ind)]).rjust(15) + str(ns.beta3_in[(f_ind)]).rjust(15) + str(ns.q20_out[(f_ind)]).rjust(15) +\
                str(ns.q30_out[(f_ind)]).rjust(15) + str(f_ind).rjust(8) + str(rs.iter_steps[(f_ind)]).rjust(10) + "\n"
            if int(rs.iter_steps[(f_ind)]) <= f_iter_min:
                f_iter_min = rs.iter_steps[(f_ind)]
        out_str += f"*************** Minimum iteration to achieve lowest g.s. energy (rounded to {hfb_accuracy} decimals) is: {f_iter_min}\n"
        print ((z,n), f" min_iter({hfb_accuracy} decimals)#:\t",f_iter_min)
    out_str = out_str + "======================================================"*2 + "\n"
    return out_str


def gen_stats(rs):
    tot_time = rs.tot_time
    iters = rs.iter_steps
    tlist = rs.tlog
    t_sum = 0; it_sum = 0
    file_count = 0
    for key in tlist:
        # some files don't print final cpu time for reasons unknown
        if key in tot_time and key in iters:
            t_sum += tot_time[(key)]
            it_sum += iters[(key)]
            file_count += 1
    n_thread = (rs.num_threads).values()
    avg_thread = reduce(lambda x, y: x+y, n_thread)/len(n_thread)
    return t_sum, it_sum, file_count, avg_thread

def median_stat(rs):
    tot_time = (rs.tot_time).values()
    iters = (rs.iter_steps).values()
    sorted_time = sorted(tot_time)
    sorted_iter = sorted(iters)
    med_time = sorted_time[int(len(tot_time)/2)]
    med_iter = sorted_iter[int(len(iters)/2)]
    max_time = sorted_time[len(tot_time)-1]
    return med_time, med_iter, max_time

def task_times(a,b,rs,avg_threads):
    f_str = f_loc(a,b)
    lsLoc = os.listdir(f_str)
    task_row = {}; task_time = {}; task_rows_count = {}; max_rows = []; max_rows_it = []
    tot_time = rs.tot_time; iter_steps = rs.iter_steps
    for ff in lsLoc:
        if "mass_table.msub.o" in ff:
            l1 = read_file(f_str,ff)
            break
    for line in l1:
        if "finished row" in line:
            nums = re.findall('\d+',line)
            task_id = int(nums[0])
            row_id = int(nums[1])
            if not (task_id in task_row):
                task_row[(task_id)] = []
                task_row[(task_id)].append(row_id)
            elif task_id in task_row:
                task_row[(task_id)].append(row_id)
    for task_id in task_row.keys():
        task_time[(task_id)] = 0; task_rows_count[(task_id)] = 0
        for row_id in task_row[(task_id)]:
            if row_id in tot_time:
                task_time[(task_id)] = task_time[(task_id)] + tot_time[(row_id)]
                task_rows_count[(task_id)] += 1
        task_time[(task_id)] = round(task_time[(task_id)]/60,2)
    max_task_time = 0
    for key in task_time:
        if float(task_time[(key)]) >= max_task_time:
            max_task_time = float(task_time[(key)])
            max_task_key = key
    max_task_row_done = len(task_row[(max_task_key)])
    avg_time_task_row_max = round(max_task_time / max_task_row_done * avg_threads,2)
    for row in task_row[(max_task_key)]:
        max_rows.append(row)
        if row in iter_steps:
            max_rows_it.append(iter_steps[(row)])
    max_it_avg = int(reduce(lambda a,b:a+b, max_rows_it)/len(max_rows_it))
    print ("*****")
    print (f"Maximum cpu time used by task {max_task_key}, total rows calculated: {len(task_row[(max_task_key)])}, total cpu time spent {round(avg_threads*max_task_time/60,2)} hours, actual time spent {round(max_task_time/60,2)} hours")
    print (f"Rows calculated by task {max_task_key} are:\n{max_rows}")
    print ("Corresponding iterations used in calculating above rows")
    print (max_rows_it)
    print (f"Average iterations per row for this task is: {max_it_avg}")
    print ("*****")
    print (f"OpenMP threads: {avg_threads}, total tasks: {len(task_row.keys())}, total cpu asked: {avg_threads*len(task_row.keys())}")
    print (f"Maximum average cpu time for tasks per row:\t{avg_time_task_row_max} minutes")


def stat_output(a,b,hfb_accuracy = 3, grid_s2 = 1, grid_s3 = 1):
    f_str = "./" + f_loc(a,b)
    #print (f_str)
    rs,ns = find_stats(f_str)
    time_sum, iter_sum, file_count, avg_threads = gen_stats(rs)
    # To get results of medians and max time consumption, use the 3 lines below
    #med_t, med_it, max_t = median_stat(rs)
    #med_t = round(med_t * avg_threads / 60, 2)
    #max_t = round(max_t * avg_threads / 60, 2)
    # average time per iteration, in seconds
    avg_t_iter = round(time_sum / iter_sum,2)
    # average time used per row (the hfbtho_MASSTABLE.dat nrow), in minutes
    avg_t_row = round( (time_sum / file_count / 60),2)
    # average cpu time per iteration, in seconds
    avg_cput_iter = round(avg_t_iter * avg_threads,2)
    # average cpu time per row, in minutes
    avg_cput_row = round((avg_t_row * avg_threads),2)
    # average iteration used
    avg_iter = int(iter_sum / file_count)
    out_str = "Folder:\t" +f_str +"\n" + "File count:\t" + str(int(file_count))+"\n"
    out_str = out_str + "Average time per iteration:\t" + str(avg_t_iter) + "\ts\n" +\
                        "Average time per row:\t\t"+str(avg_t_row) + "\tmin\n" +\
                        "Average cpu time per iter:\t" + str(avg_cput_iter) + "\ts***\n" +\
                        "Average cpu time per row:\t" + str(avg_cput_row) + "\tmin***\n" +\
                        "Average iterations per row:\t" + str(avg_iter) + "\n" #+\
#                        "Median cpu time per row:\t" + str(med_t) + "\tmin \n" +\
#                        "Maximum cpu time per row:\t" + str(max_t) + "\tmin \n" +\
#                        "Median iterations per row:\t" + str(med_it)
    # out_str contains runtime stats
    #print (out_str)
    task_times(a,b,rs,avg_threads)
    # For hfb calculation output, uncomment the 2 lines below
    # For no grid point filters, use sort_ehfbln(ns,rs,1,1); to filter q2/q3 for 0.1 step only, use any other numbers for last two inputs of sort_ehfbln()
    save_str = f_str + "\n" + sort_ehfbln(ns,rs,hfb_accuracy,grid_s2,grid_s3)
    print ("\n===============================================================")
    del rs; del ns
    return save_str


# As of 02/10/19, raw data has been moved to Document directory, it's taking too much space of google drive and was not needed anyways


def f_loc(a,b):
    edf = ["UNEDF0","UNEDF1","UNEDF2","SLY4","SKMS","SKP","SV-MIN"]
    return edf[a]+"/"+"block"+b+"/"

out_f = open("UNEDF0_block1a.txt","a")
wrt_str = stat_output(0,"1a",hfb_accuracy = 6, grid_s2 = 1, grid_s3 = 1)
out_f.write(wrt_str)

#os.system("shopt -s extglob\n"+"cd output_stats \n rm *\n")
#for a in [1]:
#    for b in range(0,2):
#        f_str = f_loc(a,b)
#        f_ss = ""
#        for ss in f_str.split("/"):
#            f_ss += ss + "_"
#        f_ss = (f_ss.split("test_"))[1]
#        if "_grid_large" in f_ss:
#            f_ss = (f_ss.split("_grid_large"))[0]
#        for i in range(0,2):
#            if f_ss[-1] == "_": f_ss = f_ss[:-1]
#        print (f_str)
#        out_f = open("output_stats/"+f_ss+".txt","a")
#        #a,b relates to file location, check f_loc(), 3rd number is for hfb_accuracy filter
#        wrt_str = stat_output(a,b,4,1,1)
#        out_f.write(wrt_str)


#def f_loc(test_num, part):
#    f_str = "raw_data/"
#    if test_num == 0:
#        f_str += "grid_test_Fm/"
#        if part == 0:
#            f_str += "230-240/"
#        elif part == 1:
#            f_str += "288-302/"
#    elif test_num == 1:
#        f_str += "heavy_test_Th/"
#        if part == 0:
#            f_str += "220-236/grid_large/"
#        elif part == 1:
#            f_str += "280-294/grid_large/"
#    elif test_num == 2:
#        f_str += "light_test/"
#        if part == 0:
#            f_str += "Ca/grid_large/"
#        elif part == 1:
#            f_str += "O/grid_large/"
#    elif test_num == 3:
#        f_str += "medium_test/"
#        if part == 0:
#            f_str += "Dy/grid_large/"
#        elif part == 1:
#            f_str += "Gd/grid_large/"
#    elif test_num == 4:
#        f_str += "Maxwell_Output/"
#    elif test_num == 5:
#        f_str += "epsi_test/"
#        if part == 0:
#            f_str += "Fm/"
#        elif part == 1:
#            f_str += "O/"
#    elif test_num == 6:
#        f_str += "EDFs_test/"
#        if part == 0:
#            f_str += "SKMS/"
#        elif part == 1:
#            f_str += "SKP/"
#        elif part == 2:
#            f_str += "SLY4/"
#        elif part == 3:
#            f_str += "SV-min/"
#        elif part == 4:
#            f_str += "UNE0/"
#        elif part == 5:
#            f_str += "UNE1/"
#        elif part == 6:
#            f_str += "UNE2/"
#    elif test_num == 7:
#        f_str += "10e0_test/"
#        if part == 0:
#            f_str += "Fm_288_302/"
#        elif part == 1:
#            f_str += "O/"
#    elif test_num == 8:
#        f_str += "10e-2_test/"
#        if part == 0:
#            f_str += "Fm_288_302/"
#        elif part == 1:
#            f_str += "O/"
#    elif test_num == 9:
#        f_str += "10e-3_test/"
#        if part == 0:
#            f_str += "Fm_288_302/"
#        elif part == 1:
#            f_str += "O/"
#    return f_str
