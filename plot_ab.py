import matplotlib.pyplot as plt
import numpy as np

with open(file="lattice01.dat") as f:
    data1 = f.readlines()

pointsA = []
pointsB = []

for i, line in enumerate(data1):
    aux = line.split(";")
    ponto = [float(aux[0]), float(aux[1]), float(aux[2])]
    if (i % 2 == 0):
        pointsA.append(ponto)
    else:
        pointsB.append(ponto)

with open(file="lattice02.dat") as f:
    data2 = f.readlines()

for i, line in enumerate(data2):
    aux = line.split(";")
    ponto = [float(aux[0]), float(aux[1]), float(aux[2])]
    if (i % 2 == 0):
        pointsA.append(ponto)
    else:
        pointsB.append(ponto)

pointsA = np.array(pointsA)
pointsB = np.array(pointsB)

with open(file="lattice03.dat") as f:
    data3 = f.readlines()

points3 = []
for line in data3:
    aux = line.split(";")
    points3.append([float(aux[0]), float(aux[1]), float(aux[2])])
points3 = np.array(points3)

# 2D PLOT (FASTER)
axA = plt.subplot(111)
grafA = axA.scatter(pointsA[:,0], pointsA[:,1], s=10)

axB = plt.subplot(111)
grafB = axB.scatter(pointsB[:,0], pointsB[:,1], s=10)

ax3 = plt.subplot(111)
graf3 = ax3.scatter(points3[:,0], points3[:,1], s=10)

axA.set_aspect(1)
axB.set_aspect(1)
ax3.set_aspect(1)

plt.show()
