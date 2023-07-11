import matplotlib.pyplot as plt
import numpy as np

with open(file="lattice01.dat") as f:
    data1 = f.readlines()

points1 = []
for line in data1:
    aux = line.split(";")
    points1.append([float(aux[0]), float(aux[1]), float(aux[2])])
points1 = np.array(points1)

with open(file="lattice02.dat") as f:
    data2 = f.readlines()

points2 = []
for line in data2:
    aux = line.split(";")
    points2.append([float(aux[0]), float(aux[1]), float(aux[2])])
points2 = np.array(points2)

with open(file="lattice03.dat") as f:
    data3 = f.readlines()

points3 = []
for line in data3:
    aux = line.split(";")
    points3.append([float(aux[0]), float(aux[1]), float(aux[2])])
points3 = np.array(points3)


ax1 = plt.subplot(111)
graf1 = ax1.scatter(points1[:,0], points1[:,1], s=5)
ax2 = plt.subplot(111)
graf2 = ax2.scatter(points2[:,0], points2[:,1], s=5)
ax3 = plt.subplot(111)
graf3 = ax3.scatter(points3[:,0], points3[:,1], s=5)

ax1.set_aspect(1)
ax2.set_aspect(1)
ax3.set_aspect(1)

ax1.set_xlim(-70, 70)
ax1.set_ylim(0, 140)

## 3D PLOT (TOO SLOW)
# ax1 = plt.subplot(111, projection='3d')
# graf1 = ax1.scatter(points1[:,0], points1[:,1], points1[:,2], s=1)
# ax2 = plt.subplot(111, projection='3d')
# graf2 = ax2.scatter(points2[:,0], points2[:,1], points2[:,2], s=1)
# ax3 = plt.subplot(111, projection='3d')
# graf3 = ax3.scatter(points3[:,0], points3[:,1], points3[:,2], s=1)

# ax1.set_box_aspect((1,1,1/10))
# ax2.set_box_aspect((1,1,1/10))
# ax3.set_box_aspect((1,1,1/10))

plt.show()
