import matplotlib.pyplot as plt
import numpy as np

with open(file="lattice1A.dat") as f:
    data1A = f.readlines()

points1A = []
for line in data1A:
    aux = line.split(";")
    points1A.append([float(aux[0]), float(aux[1]), float(aux[2])])
points1A = np.array(points1A)

with open(file="lattice1B.dat") as f:
    data1B = f.readlines()

points1B = []
for line in data1B:
    aux = line.split(";")
    points1B.append([float(aux[0]), float(aux[1]), float(aux[2])])
points1B = np.array(points1B)

with open(file="lattice2A.dat") as f:
    data2A = f.readlines()

points2A = []
for line in data2A:
    aux = line.split(";")
    points2A.append([float(aux[0]), float(aux[1]), float(aux[2])])
points2A = np.array(points2A)

with open(file="lattice2B.dat") as f:
    data2B = f.readlines()

points2B = []
for line in data2B:
    aux = line.split(";")
    points2B.append([float(aux[0]), float(aux[1]), float(aux[2])])
points2B = np.array(points2B)

with open(file="lattice03.dat") as f:
    data3 = f.readlines()

points3 = []
for line in data3:
    aux = line.split(";")
    points3.append([float(aux[0]), float(aux[1]), float(aux[2])])
points3 = np.array(points3)

# 2D PLOT (FASTER)
ax1A = plt.subplot(111)
graf1A = ax1A.scatter(points1A[:,0], points1A[:,1], s=10)

ax2A = plt.subplot(111)
graf2A = ax2A.scatter(points2A[:,0], points2A[:,1], s=10)

ax1B = plt.subplot(111)
graf1B = ax1B.scatter(points1B[:,0], points1B[:,1], s=10)

ax2B = plt.subplot(111)
graf2B = ax2B.scatter(points2B[:,0], points2B[:,1], s=10)

ax3 = plt.subplot(111)
graf3 = ax3.scatter(points3[:,0], points3[:,1], s=10)

ax1A.set_aspect(1)
ax2A.set_aspect(1)
ax1B.set_aspect(1)
ax2B.set_aspect(1)
ax3.set_aspect(1)

# # 3D PLOT (TOO SLOW)
# ax1 = plt.subplot(111, projection='3d')
# graf1 = ax1.scatter(points1[:,0], points1[:,1], points1[:,2], s=1)
# ax2 = plt.subplot(111, projection='3d')
# graf2 = ax2.scatter(points2[:,0], points2[:,1], points2[:,2], s=1)
# ax3 = plt.subplot(111, projection='3d')
# graf3 = ax3.scatter(points3[:,0], points3[:,1], points3[:,2], s=1)
# 
# ax1.set_box_aspect((np.ptp(points1[:,0]), np.ptp(points1[:,1]), np.ptp(points1[:,2])))
# ax2.set_box_aspect((np.ptp(points2[:,0]), np.ptp(points2[:,1]), np.ptp(points2[:,2])))
# ax1.azim = 90
# ax1.elev = 90

plt.show()
