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


ax1 = plt.subplot(111, projection='3d')
graf1 = ax1.scatter(points1[:,0], points1[:,1], points1[:,2])
ax2 = plt.subplot(111, projection='3d')
graf2 = ax2.scatter(points2[:,0], points2[:,1], points2[:,2])

ax1.set_box_aspect((1,1,1/10))
ax2.set_box_aspect((1,1,1/10))

plt.show()
