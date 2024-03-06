import numpy as np
import matplotlib.pyplot as plt


with open("latticeAA.dat") as f:
    dados_aa = f.readlines()
    aa = []
    for d in dados_aa:
        aux = d.strip().split(";")
        aa.append([float(aux[0]), float(aux[1])])
    aa = np.array(aa)

with open("latticeAB.dat") as f:
    dados_ab = f.readlines()
    ab = []
    for d in dados_ab:
        aux = d.strip().split(";")
        ab.append([float(aux[0]), float(aux[1])])
    ab = np.array(ab)

with open("latticeBA.dat") as f:
    dados_ba = f.readlines()
    ba = []
    for d in dados_ba:
        aux = d.strip().split(";")
        ba.append([float(aux[0]), float(aux[1])])
    ba = np.array(ba)

with open("latticeBB.dat") as f:
    dados_bb = f.readlines()
    bb = []
    for d in dados_bb:
        aux = d.strip().split(";")
        bb.append([float(aux[0]), float(aux[1])])
    bb = np.array(bb)

try:
    with open("latticeA1.dat") as f:
        dados_a1 = f.readlines()
        a1 = []
        for d in dados_a1:
            aux = d.strip().split(";")
            a1.append([float(aux[0]), float(aux[1])])
        a1 = np.array(a1)
    
    with open("latticeA2.dat") as f:
        dados_a2 = f.readlines()
        a2 = []
        for d in dados_a2:
            aux = d.strip().split(";")
            a2.append([float(aux[0]), float(aux[1])])
        a2 = np.array(a2)
    
    with open("latticeB1.dat") as f:
        dados_b1 = f.readlines()
        b1 = []
        for d in dados_b1:
            aux = d.strip().split(";")
            b1.append([float(aux[0]), float(aux[1])])
        b1 = np.array(b1)
    
    with open("latticeB2.dat") as f:
        dados_b2 = f.readlines()
        b2 = []
        for d in dados_b2:
            aux = d.strip().split(";")
            b2.append([float(aux[0]), float(aux[1])])
        b2 = np.array(b2)
except:
    print("No files found for lattices A1, A2, B1 or B2")


diff = ab[35]
origin = aa[44]
print(f"norm of depicted vector: {np.linalg.norm(diff)}")

# creating new lattices adding the diff
aa2 = aa + diff
ab2 = ab + diff
ba2 = ba + diff
bb2 = bb + diff

# plots
ax1 = plt.subplot(111)
ax1.scatter(aa[:,0],aa[:,1], color='blue')
ax1.scatter(ab[:,0],ab[:,1], color='orange')
# ax1.scatter(ba[:,0],ba[:,1], color='orange')
# ax1.scatter(bb[:,0],bb[:,1], color='blue')
ax1.set_aspect("equal")

ax2 = plt.subplot(111)
ax2.quiver(origin[0], origin[1], diff[0], diff[1], angles='xy', scale_units='xy', scale=1)

# ax3 = plt.subplot(111)
# ax3.scatter(aa2[:,0],aa2[:,1])
# ax3.scatter(ab2[:,0],ab2[:,1])
# ax3.scatter(ba2[:,0],ba2[:,1])
# ax3.scatter(bb2[:,0],bb2[:,1])
# ax3.set_aspect("equal")

plt.show()