def main(cut):
    EDF = ["UNEDF0","UNEDF1","UNEDF2","SLY4","SV-MIN","SKMS","SKP"]
    for e in EDF:
        f = open("new_data/reduced/no_filter/"+e+"_reduced_no_filter.dat","r")
        l = f.readlines()
        out = open("abnormal_deformation/"+e+"_def_beyond_0"+str(int(cut*100))+".dat","w")
        out_str = ""
        for ind,line in enumerate(l):
            if not ind:
                out_str += line
            if ind:
                ss=line.split()
                Z,N,b2,b3 = int(ss[0]),int(ss[1]),float(ss[4]),float(ss[5])
                if abs(b2)>cut or abs(b3)>cut:
                    out_str += line
        out.write(out_str)
        f.close();out.close()

main(0.4)
