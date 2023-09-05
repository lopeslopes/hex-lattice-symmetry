import numpy as np

with open("anglesA.dat") as f:
    data_a = f.readlines()

with open("anglesB.dat") as f:
    data_b = f.readlines()

angles_a = []
for d in data_a:
    angles_a.append(float(d))

angles_b = []
for d in data_b:
    angles_b.append(float(d))

angles_a = sorted(angles_a)
angles_b = sorted(angles_b)

f = open("treated_anglesA.dat", "w")
for a in angles_a:
    if (a <= np.pi/3):
        f.write(f"{a}\n")

f = open("treated_anglesB.dat", "w")
for b in angles_b:
    if (b <= np.pi/3):
        f.write(f"{b}\n")
