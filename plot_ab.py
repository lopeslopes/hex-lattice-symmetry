import matplotlib.pyplot as plt
import numpy as np

with open(file="lattice1t.dat") as f:
    data1 = f.readlines()

points1 = []
for line in data1:
    aux = line.split(";")
    points1.append([float(aux[0]), float(aux[1]), float(aux[2])])
points1 = np.array(points1)

with open(file="lattice2t.dat") as f:
    data2 = f.readlines()

points2 = []
for line in data2:
    aux = line.split(";")
    points2.append([float(aux[0]), float(aux[1]), float(aux[2])])
points2 = np.array(points2)

with open(file="lattice3t.dat") as f:
    data3 = f.readlines()

points3 = []
for line in data3:
    aux = line.split(";")
    points3.append([float(aux[0]), float(aux[1]), float(aux[2])])
points3 = np.array(points3)

with open(file="lattice4t.dat") as f:
    data4 = f.readlines()

points4 = []
for line in data4:
    aux = line.split(";")
    points4.append([float(aux[0]), float(aux[1]), float(aux[2])])
points4 = np.array(points4)

# 2D PLOT (FASTER)
ax1 = plt.subplot(111)
graf1 = ax1.scatter(points1[:,0], points1[:,1], s=10)

ax2 = plt.subplot(111)
graf2 = ax2.scatter(points2[:,0], points2[:,1], s=10)

ax3 = plt.subplot(111)
graf3 = ax3.scatter(points3[:,0], points3[:,1], s=10)

ax4 = plt.subplot(111)
graf4 = ax4.scatter(points4[:,0], points4[:,1], s=10)

ax1.set_aspect(1)
ax2.set_aspect(1)
ax3.set_aspect(1)
ax4.set_aspect(1)

plt.show()
