import numpy as np

f = open("CFM_LEAP_mod.txt","r")
lines = []
for line in f:
    lines.append(line)

vals = []
for l in lines:
    spl = l.split()
    vals.append(spl)

f.close()

f = open("CFM_LEAP_civil.txt", "w")

print(vals)

for i in range(len(vals[0])):
    #f.write(str(vals[0][i])+ "   "+str(vals[1][i])+"   "+str(vals[2][i])+"   "+str(vals[3][i])+"\n")
    f.write('%-50s %-20s %-30s %-30s\n'%(str(vals[0][i]),str(vals[1][i]),str(vals[2][i]),str(vals[3][i])))