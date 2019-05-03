import numpy as np
from functools import reduce
import math
import datetime as dt
import os
import time
import pexpect
import calendar as cld

def read_file(f_loc,file_name):
    f1 = open(f_loc+file_name, "r")
    l1 = f1.readlines()
    # Time tick when the file was updated
    ticks_mfile = os.path.getmtime(f_loc+file_name)
    t_file = time.ctime(ticks_mfile)
    # Time tick when the job started, read from output
    start_time = time.strptime(l1[0][:28],"%a %b %d %H:%M:%S %Z %Y")
    start_epoch = cld.timegm(start_time)
    # Current time tick
    time_now = time.strptime(str(dt.datetime.now())[:-7],"%Y-%m-%d %H:%M:%S")
    now_epoch = cld.timegm(time_now)
    ticks_diff = start_epoch - now_epoch
    print (f"Job started: {l1[0][0:28]}\t Last update: {t_file}")
    return l1,ticks_diff

def scp_file(edf,num_block,job_id,f_loc,f_head,f_suff):
    pwf = open("pw").readlines()[0]
    # Used for masstable calculation that uses block directories
    #hpc_path = "scp HPCC:~/hfbthoTest/src/hfbtho/Octupole_Table/"+edf+"/block"
    #hpc_path = hpc_path + num_block + f_suff + "/" + f_head + str(job_id) + " "
    #hpc_path = hpc_path + f_loc
    hpc_path = "scp HPCC:~/hfbthoTest/src/hfbtho/Octupole_Table/"+\
    "non_octupole_deformed_softness/"
    hpc_path += edf+"/mass_table.msub.o"+str(job_id) + " " + f_loc
    child = pexpect.spawn(hpc_path)
    child.expect("d:")
    child.sendline(pwf)
    child.expect(pexpect.EOF, timeout=10)

def read_time(edf,job_id,num_block='',tags=0,f_update=True,f_loc="hpcc_job_output/", f_head="mass_table.msub.o",show=False,f_suff=""):
    print (f"BLOCK-{num_block+'-'+str(tags)+'-'+edf}:")
    if f_update:
        scp_file(edf,num_block,job_id,f_loc,f_head,f_suff)
    l1,ticks_diff = read_file(f_loc,f_head+str(job_id))
    task_row_time = {}; recv_row = []; fin_row = []; unfin_row = []; task_recv_row = {}
    task_t_accum = {};
    task_id = 0; task_time = 0.1; time_sum = 0; row_count = 0; exit_count = 0
    max_task_id = 0; exit_rate = 0; run_time = 0; max_row = 0; max_time = 0
    est_t_run = []; est_t_left = []; task_ing = set()
    max_accum = 0; min_time = 10000; min_time_row = 0; max_time_row = 0
    for line in l1:
        ss = line.split()
        if "Total row number" in line:  max_row = int(ss[-1])
        if "received row" in line:
            task_id = int(ss[1])
            row_id = int(ss[4])
            if task_id > max_task_id: max_task_id = task_id
            if task_id not in task_recv_row.keys() :
                task_recv_row[(task_id)] = []
            else:
                task_recv_row[(task_id)].append(row_id)
            recv_row.append(row_id)
            if row_id <= max_row:   task_ing.add(task_id)
        if "minutes" in line:
            # In newer version of dynamic, the time is the index = -2 number
            # In the old version, time is stuck with the string "minutes"
            task_id = int(ss[1])
            if task_id > max_task_id:   max_task_id = task_id
            # initialize list in dictionaries
            if task_id not in task_row_time.keys(): task_row_time[(task_id)] = []
            if len(ss[-1]) > 9:
            # in newer version, do ss[-2] without slicing
                task_time = float((ss[-1])[0:-7])
            elif len(ss[-1]) < 9:
                task_time = float(ss[-2])
            if task_time > max_time:
                max_time = round(task_time,2)
                max_time_row = int(ss[4])
            if task_time < min_time:
                min_time = round(task_time,2)
                min_time_row = int(ss[4])
            # task_id vs list of (row_id, time)
            task_row_time[(task_id)].append([int(ss[4]),round(task_time,3)])
            # task_id vs list of rows completed
            time_sum += task_time
            row_count += 1  # number of rows finished
            task_ing.discard(task_id)
        if "exit signal" in line:   exit_count += 1
        if "finished row" in line:  fin_row.append(int(ss[4]))
        if "RunTime" in line:
            ss = line.split("=")
            run_time_str = ss[1]
            days = 0
            if "-" in run_time_str:
                days = int(run_time_str.split("-")[0])
                run_time_str = run_time_str.split("-")[1]
            tt1 = run_time_str.split(":")
            run_time = int(tt1[0])*60 + int(tt1[1]) + 1 + days*24*60
    fin_row.sort(); recv_row.sort()
    if exit_count > 0: max_row = recv_row[-1]-1
    up_bound = max_row
    if max_row in fin_row:
        for i in range(len(fin_row)-1):
            this_row = fin_row[len(fin_row)-i-2]
            last_row = fin_row[len(fin_row)-i-1]
            if (last_row - this_row) == 1 :
                up_bound = this_row-1
            else:
                break
    unfin_row_set = set(recv_row) - set(fin_row)
    unfin_row_set.discard(max_row+1)
    unfin_row_set.discard(-1)
    unfin_row = sorted(unfin_row_set)
    exit_rate = round(exit_count / max_task_id *100,2)
    avg_time_per_row = round(time_sum/row_count,2)
    unfin_row_count = len(unfin_row)
    avg_run_time = round(run_time*max_task_id/row_count,2)
    margin = round(avg_run_time / avg_time_per_row * 100,2)
    if len(unfin_row) == 0:
        unfin_row.append("None")
        unfin_row_count = 0
    if exit_rate == 100:
        up_bound = "None"
    # Find potential error and frozen MPI tasks:
    # 1. Find total time accumuated in each tasks
    # 2. Find the longest time accumulated, that gives us an estimate of how much time was elapsed
    # 3. Use the longest time accumulated, and the maximum time for a row, ideally this max time is due to iteration limit being hit, so more or less should be the max needed.
    # 4. Sum each task's accumulated time with max_time, if smaller than max_accum, then it's almost certain there's a crash. Test on block-2, 6 samples
    min_row_comp = 100; max_row_comp = 0
    for key in task_row_time.keys():
        row_comp = len(task_row_time[(key)])
        if min_row_comp > row_comp:
            min_row_comp = row_comp
            min_row_task = key
        if max_row_comp < row_comp:
            max_row_comp = row_comp
            max_row_task = key
    for i in range(1,max_task_id+1):
        try:
            task_row_time[(i)] = np.array(task_row_time[(i)])
            task_t_accum[(i)] = round(reduce((lambda a,b : a+b), (task_row_time[(i)])[:,1]),2)
            if max_accum < task_t_accum[(i)]: max_accum = task_t_accum[(i)]
        except:
            continue
    crash_count=0; warning_rows=set()
    for i in range(1,max_task_id+1):
        try:
            if task_t_accum[(i)] + max(max_time,120) < max_accum:
                cache_row = list(set(task_recv_row[(i)]) - set(task_row_time[(i)][:,0]))[0]
                froze_time = round(max_accum - task_t_accum[(i)],2)
                crash_count+=1
                warning_rows.add(i)
                print (f"\n### WARNING ###: Task {i} potential failure. Check row: {cache_row}, time frozen: {froze_time}min \t{crash_count}\n")
        except:
            continue
    c_margin = max_accum/(row_count/max_task_id)/avg_time_per_row
    print (f"Job ID: {job_id}, ntasks={max_task_id+1}, total rows: {max_row}")
    try:
        print (f"Min: Task {min_row_task} completed {min_row_comp} rows, currently working on row {task_recv_row[(min_row_task)][-1]}")
    except:
        pass
    try:
        print (f"Max: Task {max_row_task} completed {max_row_comp} rows, currently working on row {task_recv_row[(max_row_task)][-1]}")
    except:
        pass
    print ("Slowest row time:",max_time,"minutes, row:",max_time_row)
    print ("Fastest row time:",min_time,"minutes, row:",min_time_row)
    if run_time != 0 :
        print (f"Actual avg walltime: {avg_run_time} minutes, total job time: {str(dt.timedelta(minutes=run_time))[:-3]} hrs, extra time margin: {margin}%")
        print ("*********Complete*********")
    else:
        comp_rate = round(row_count / (max_row-recv_row[0]+1) *100,2)
        t_cache = int(avg_time_per_row*((max_row-row_count-recv_row[0]+1)/max_task_id)*c_margin) + 60
#    if exit_rate > 0:
        t_cache2 = 0
        for key,value in task_t_accum.items():
            # you want to check if the task_id is still running
            if key in task_ing:
                task_time_left = max_time - (max_accum - value)
                if task_time_left > t_cache2: t_cache2 = int(task_time_left)
        worst_end_time = str(dt.datetime.now()+dt.timedelta(minutes=(t_cache2+max_accum))+dt.timedelta(seconds=ticks_diff))[:-10]
        max_accum = int(max_accum)
        if t_cache2 > t_cache:
            t_cache = t_cache2; flag = "worst case time left"
        else: flag = "est. time left"
        est_end_time = str(dt.datetime.now()+dt.timedelta(minutes=(t_cache+max_accum))+dt.timedelta(seconds=ticks_diff))[:-10]
        est_end_time = max(worst_end_time,est_end_time)
        print (f"Completion rate: {comp_rate}%, {flag}: {str(dt.timedelta(minutes=t_cache))[:-3]}, est. time elapsed: {str(dt.timedelta(minutes=max_accum))[:-3]}\nEstimated end time:\t{est_end_time} EST, est. time margin: {round(c_margin*100,2)}%")
        print ("*********IN PROGRESS*********")
    print (f"Average {avg_time_per_row} minutes per row, finished rows: {row_count}, rows left: {max_row-row_count-recv_row[0]+1}, exit rate: {exit_rate}%\nSmallest row in process: {unfin_row[0]}, re-run upper bound: {up_bound}, rows in process: {unfin_row_count}")
    print ("========================================================================================")
    # check how many rows did each task do, debug purpose.
    # Tasks tend to crash when doing their 82nd row for unknown reasons 1/22/19
    for key in task_recv_row.keys():
        len_dict = len(task_recv_row[(key)])
        if len_dict > 200:
            print (key, len_dict)
    largest_unfin=0
    if show:
        all_unfin = []
        for n in range(1,max_row+1):
            if n not in fin_row: all_unfin.append(n)
        all_unfin.sort()
        for i in range(1,len(all_unfin)):
            if all_unfin[-i] - all_unfin[-i-1] == 1:
                largest_unfin = all_unfin[-i-1]
            else: break
        print (all_unfin)#(unfin_row_set)
        print (len(all_unfin))#len(unfin_row_set))
        print (largest_unfin)


################################################################################################
#read_time("SLY4",8108656,"1a")  #halted due to interaction check
#read_time("SLY4",8108889,"1b")  #halted due to interaction check
#read_time("SKMS",8112000,"8")  #halted due to interaction check
################################################################################################

################################################################################################
#active tasks:



################################################################################################


################################################################################################
#completed:
# non-Octupole deformed pure quadrupole minimum
#read_time("UNEDF0",16463204,show=1)
#read_time("UNEDF1",16463210,show=1)
#read_time("UNEDF2",16463216,show=1)
#read_time("SV-MIN",16463238,show=0)
#read_time("SLY4",16463250,show=1)
#read_time("SKMS",16464687,show=0)
#read_time("SKP",16464657,show=0)



# Octupole softness no LN deformation input
#read_time("UNEDF0",16277700,show=0)
#read_time("UNEDF1",16277660,show=0)
#read_time("UNEDF2",16277721,show=0)
#read_time("SV-MIN",16277804,show=0)
#read_time("SLY4",16277764,show=0)
#read_time("SKMS",16277843,show=0)
#read_time("SKP",16278179,show=0)

# Octupole softness
#read_time("UNEDF0",15065552,show=1)
#read_time("UNEDF1",15065565,show=1)
#read_time("UNEDF2",15065497,show=1)
#read_time("SLY4",15065572,show=1)
#read_time("SKMS",15065543,show=1)
#read_time("SV-MIN",15065576,show=1)
#read_time("SKP",15065583,show=1)



# Octupole energy follow up
#read_time("UNEDF0",14417394,f_update=False)
#read_time("UNEDF1",14405172) # 3 unfin [3329, 3636, 3732]
#read_time("UNEDF2",14447413,show=1)# 15 unfin [5352, 5674, 5829, 6036, 6173, 6266, 6342, 6359, 6444, 6452, 6519, 6617, 6633, 6676, 6803]
#read_time("SLY4",14418141,show=1) # 9 unfin [5417, 5596, 5597, 5756, 5769, 5867, 5912, 5954, 5997]
#read_time("SV-MIN",14417799)

#read_time("SKP",14417989)
#read_time("SKMS",14496759)

# Octupole energy 0.50 set
#read_time("UNEDF0",14988092,show=1)
#read_time("UNEDF1",14988168,show=1)
#read_time("UNEDF2",14988159,show=1)
#read_time("SLY4",14988015,show=1)
#read_time("SKMS",14988001,show=1)
#read_time("SV-MIN",14988005,show=1)
#read_time("SKP",14988081,show=1)



# tests
#read_time("UNEDF0",8290915,"Oxy",f_update=False)
#read_time("UNEDF0",8314518,"hbzero_20")
#read_time("SLY4",8339402 ,"1",f_update=False)


#SKP
#read_time("SKP",9848807,"1a",1,f_update=False)
#read_time("SKP",10070271,"1a",2,f_update=False)
#read_time("SKP",9848817,"1b",1,f_update=False)
#read_time("SKP",10136256,"1b",2,f_update=False)
#read_time("SKP",9896566,"2",1,f_update=False)
#read_time("SKP",9919747,"2",2,f_update=False)
#read_time("SKP",9919720,"3",f_update=False)
#read_time("SKP",10012223,"4",f_update=False)
#read_time("SKP",10079727,"5",f_update=False)
#read_time("SKP",10038828,"6",f_update=False) #18504
#read_time("SKP",10152253,"6",2,f_update=False) #406 need finished=336
#read_time("SKP",10152488,"6",3,f_suff="/batch3",f_update=False) #365 ok
#read_time("SKP",10152151,"7",2,f_update=False) #502 ok
#read_time("SKP",10152392,"7",3,f_suff="/batch3",f_update=False) #592 ok
#read_time("SKP",10038829,"7",f_update=False) #19906
#read_time("SKP",10079728,"8",1,f_update=False)
#read_time("SKP",10169844,"8",2,f_update=False) #15599-19575
#read_time("SKP",9848774,"9",1,f_update=False)
#read_time("SKP",9896638,"10",1,f_update=False)
#read_time("SKP",10199156,"10",2)

#SV-min
#read_time("SV-min",9672265,"1a",1,f_update=False)
#read_time("SV-min",9848380,"1a",2,f_update=False)
#read_time("SV-min",9672283,"1b",1,f_update=False)
#read_time("SV-min",9961338,"1b",2,f_update=False)
#read_time("SV-min",9763416,"2",f_update=False)
#read_time("SV-min",9754102,"3",1,f_update=False)
#read_time("SV-min",9848552,"3",2,f_update=False)
#read_time("SV-min",9863128,"4",f_update=False)
#read_time("SV-min",9806916,"5",1,f_update=False)
#read_time("SV-min",9896258,"5",2,f_update=False)
#read_time("SV-min",9863162,"6",f_update=False)
#read_time("SV-min",9863173,"7",1,f_update=False)
#read_time("SV-min",9919731,"7",2,f_update=False)
#read_time("SV-min",10070479,"7",3,f_update=False)
#read_time("SV-min",10082145,"7",4,f_update=False)
#read_time("SV-min",9691889,"8",f_update=False)
#read_time("SV-min",9671970,"9",f_update=False)
#read_time("SV-min",9691858,"10",f_update=False)


#SKMS
#read_time("SKMS",9046455,"1a",1,f_update=False)
#read_time("SKMS",9153594,"1a",2,f_update=False)
#read_time("SKMS",9046475,"1b",f_update=False)
#read_time("SKMS",8074073,"2",f_update=False)
#read_time("SKMS",9235637,"3",f_update=False)
#read_time("SKMS",9235712,"4",f_update=False)
#read_time("SKMS",9269031,"5",1,f_update=False)
#read_time("SKMS",9327777,"5",2,f_update=False)

#read_time("SKMS",9301070,"6",1,f_update=False)
#read_time("SKMS",9301079,"7",1,f_update=False)
#read_time("SKMS",9634804,"7",2,f_update=False)
#read_time("SKMS",9331471,"8",f_update=False)
#read_time("SKMS",9235771,"9",f_update=False)
#read_time("SKMS",9331491,"10",f_update=False)

#UNEDF0
#read_time("UNEDF0",9179643,"1a",f_update=False)
#read_time("UNEDF0",9218783,"1b",f_update=False)
#read_time("UNEDF0",8757009,"2",f_update=False)
#read_time("UNEDF0",8835374,"3",1,f_update=False)
#read_time("UNEDF0",9009293,"3",2,f_update=False)
#read_time("UNEDF0",9035857,"3",3,f_update=False)
#read_time("UNEDF0",8835420,"4",f_update=False)
#read_time("UNEDF0",9235601,"5",1,f_update=False)
#read_time("UNEDF0",9292987,"5",2,f_update=False)
#read_time("UNEDF0",8994372,"6",f_update=False)
#read_time("UNEDF0",9068429,"7",1,f_update=False)
#read_time("UNEDF0",9152844,"7",2,f_update=False)
#read_time("UNEDF0",9100068,"8",f_update=False)
#read_time("UNEDF0",8741671,"9",1,f_update=False,show_unfin=True)
#read_time("UNEDF0",8774591,"9",2,f_update=False)
#read_time("UNEDF0",9169676,"10")


#UNEDF1
#read_time("UNEDF1",8131261,"1a",f_update=False)
#read_time("UNEDF1",8261576,"1b",f_update=False)
#read_time("UNEDF1",8219632,"2a",1,f_update=False)
#read_time("UNEDF1",8261717,"2b",1,f_update=False)
#read_time("UNEDF1",8280804,"3",1,f_update=False)
#read_time("UNEDF1",8498217,"4",1,f_update=False)
#read_time("UNEDF1",8537518,"4",2,f_update=False)
#read_time("UNEDF1",8539881,"5a",f_update=False)
#read_time("UNEDF1",8164020,"5b",1,f_update=False)
#read_time("UNEDF1",8700645,"6",1,f_update=False)
#read_time("UNEDF1",8573049,"7",1,f_update=False)
#read_time("UNEDF1",8570325,"8a",1,f_update=False)
#read_time("UNEDF1",8280900,"8b",f_update=False)
#read_time("UNEDF1",8121542,"9",1,f_update=False)
#read_time("UNEDF1",8164470,"9",2,f_update=False)
#read_time("UNEDF1",8276555,"10",f_update=False)
#read_time("UNEDF1",10733473,"1",f_update=0)
#read_time("UNEDF1",10733504,"2",f_update=0,show=1)

#UNEDF2
#read_time("UNEDF2",7439690,"1",f_update=False)
#read_time("UNEDF2",6964196,"2","ptg",f_update=False)
#read_time("UNEDF2",7245783,"3a",f_update=False)
#read_time("UNEDF2",7492108,"3b",f_update=False)
#read_time("UNEDF2",7526383,"4a",1,f_update=False)
#read_time("UNEDF2",7582562,"4a",2,f_update=False)
#read_time("UNEDF2",7602218,"4a",3,f_update=False)
#read_time("UNEDF2",7247935,"4b",f_update=False)  #re-run partial
#read_time("UNEDF2",7248400,"5",f_update=False)
#read_time("UNEDF2",6964112,"6","ptg",f_update=False)
#read_time("UNEDF2",7526364,"7",1,f_update=False)
#read_time("UNEDF2",7582583,"7",2,f_update=False)
#read_time("UNEDF2",7589318,"7",3,f_update=False)
#read_time("UNEDF2",7565767,"8",f_update=False)
#read_time("UNEDF2",7589713,"9",f_update=False)
#read_time("UNEDF2",7646402,"10",f_update=False)
#read_time("UNEDF2",7609175,"11",f_update=False)
#read_time("UNEDF2",7631105,"12",f_update=False)

#SLY4
#read_time("SLY4",8108656,"1a",1,f_update=False)
#read_time("SLY4",8131388,"1a",2,f_update=False)
#read_time("SLY4",8108889,"1b",1,f_update=False)
#read_time("SLY4",8131413,"1b",2,f_update=False)
#read_time("SLY4",7820004,"2",f_update=False)
#read_time("SLY4",7934255,"3",f_update=False)
#read_time("SLY4",8004466,"4",f_update=False)
#read_time("SLY4",7820120,"5",1,f_update=False)
#read_time("SLY4",8064502,"5",2,f_update=False)
#read_time("SLY4",8004491,"6",1,f_update=False)
#read_time("SLY4",8014572,"7",1,f_update=False)
#read_time("SLY4",8039078,"8",f_update=False)
#read_time("SLY4",8039077,"9",f_update=False)
#read_time("SLY4",8064690,"10",f_update=False)
#

################################################################################################

