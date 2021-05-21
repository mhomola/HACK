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
print(len(vals[3]))

for i in range(len(vals[0])):
    f.write(str(vals[0][i])+ "   "+str(vals[1][i])+"   "+str(vals[2][i])+"\n")