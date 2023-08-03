import matplotlib.pyplot as plt
import numpy as np

with open(file="latticeA1.dat") as f:
    dataA1 = f.readlines()

pointsA1 = []

for i, line in enumerate(dataA1):
    aux = line.split(";")
    ponto = [float(aux[0]), float(aux[1]), float(aux[2])]
    pointsA1.append(ponto)

pointsB1 = []

with open(file="latticeB1.dat") as f:
    dataB1 = f.readlines()

for i, line in enumerate(dataB1):
    aux = line.split(";")
    ponto = [float(aux[0]), float(aux[1]), float(aux[2])]
    pointsB1.append(ponto)

pointsA1 = np.array(pointsA1)
pointsB1 = np.array(pointsB1)

# 2D PLOT (FASTER)
axA = plt.subplot(111)
grafA = axA.scatter(pointsA1[:,0], pointsA1[:,1], s=10)

axB = plt.subplot(111)
grafB = axB.scatter(pointsB1[:,0], pointsB1[:,1], s=10)

axA.set_aspect(1)
axB.set_aspect(1)

plt.show()
