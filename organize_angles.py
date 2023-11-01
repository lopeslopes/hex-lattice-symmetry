import numpy as np

with open("angles.dat") as f:
    data_a = f.readlines()

angles_a = []
for d in data_a:
    angles_a.append(float(d))

angles_a = sorted(angles_a)

for a in angles_a:
    a_dg = a*(180/np.pi)
    if (a_dg < 1.2): print(a_dg, a)
