import matplotlib.pyplot as plt
import numpy as np

with open(file="latticeA1.dat") as f:
    dataA1 = f.readlines()

pointsA1 = []
for line in dataA1:
    aux = line.split(";")
    pointsA1.append([float(aux[0]), float(aux[1]), float(aux[2])])
pointsA1 = np.array(pointsA1)

with open(file="latticeB1.dat") as f:
    dataB1 = f.readlines()

pointsB1 = []
for line in dataB1:
    aux = line.split(";")
    pointsB1.append([float(aux[0]), float(aux[1]), float(aux[2])])
pointsB1 = np.array(pointsB1)

with open(file="latticeA2.dat") as f:
    dataA2 = f.readlines()

pointsA2 = []
for line in dataA2:
    aux = line.split(";")
    pointsA2.append([float(aux[0]), float(aux[1]), float(aux[2])])
pointsA2 = np.array(pointsA2)

with open(file="latticeB2.dat") as f:
    dataB2 = f.readlines()

pointsB2 = []
for line in dataB2:
    aux = line.split(";")
    pointsB2.append([float(aux[0]), float(aux[1]), float(aux[2])])
pointsB2 = np.array(pointsB2)

# 2D PLOT (FASTER)
axA1 = plt.subplot(111)
grafA1 = axA1.scatter(pointsA1[:,0], pointsA1[:,1], s=30)

axB1 = plt.subplot(111)
grafB1 = axB1.scatter(pointsB1[:,0], pointsB1[:,1], s=30)

axA2 = plt.subplot(111)
grafA2 = axA2.scatter(pointsA2[:,0], pointsA2[:,1], s=10)

axB2 = plt.subplot(111)
grafB2 = axB2.scatter(pointsB2[:,0], pointsB2[:,1], s=10)

axA1.set_aspect(1)
axB1.set_aspect(1)
axA2.set_aspect(1)
axB2.set_aspect(1)

plt.legend(['A1', 'B1', 'A2', 'B2'])

plt.show()
