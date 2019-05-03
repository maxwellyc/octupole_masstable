import numpy as np

def read_file(folder_loc,file_name):
    f1 = open(folder_loc+file_name, "r")
    l1 = f1.readlines()
    return l1

def read_incomplete(job_id,folder_loc=""):
    l1 = read_file(folder_loc,str(job_id))
    row_inc = []
    for line in l1:
        if "thoout_" in line:
            # In newer version of dynamic, the time is the index = -2 number
            # In the old version, time is stuck with the string "minutes"
            ss = line.split("thoout_")
            row_str = ss[-1]
            row_inc.append(int(row_str[:6]))
    row_inc.sort()
    print ("incomplete row count:",len(row_inc))
    print ("smallest incomplete row #: ",row_inc[0])


read_incomplete("incomplete.txt")


