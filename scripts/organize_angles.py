import numpy as np

with open("angles.dat") as f:
    data_a = f.readlines()

angles_a = []
for d in data_a:
    angles_a.append(float(d))

angles_a = sorted(angles_a)

treated = set()
for a in angles_a:
    a_dg = a*(180/np.pi)
    if ((a_dg < 1.1) and (a_dg > 1.09)):
        treated.add(a)
treated2 = sorted(treated)

for ang in treated2:
    print(ang, ang*180/np.pi)
