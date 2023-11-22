import matplotlib.pyplot as plt
import numpy as np

with open(file="latticeA1.dat") as f:
    dataA1 = f.readlines()

pointsA1 = []
for line in dataA1:
    aux = line.split(";")
    pointsA1.append([float(aux[0]), float(aux[1])])
pointsA1 = np.array(pointsA1)

with open(file="latticeB1.dat") as f:
    dataB1 = f.readlines()

pointsB1 = []
for line in dataB1:
    aux = line.split(";")
    pointsB1.append([float(aux[0]), float(aux[1])])
pointsB1 = np.array(pointsB1)

with open(file="latticeA2.dat") as f:
    dataA2 = f.readlines()

pointsA2 = []
for line in dataA2:
    aux = line.split(";")
    pointsA2.append([float(aux[0]), float(aux[1])])
pointsA2 = np.array(pointsA2)

with open(file="latticeB2.dat") as f:
    dataB2 = f.readlines()

pointsB2 = []
for line in dataB2:
    aux = line.split(";")
    pointsB2.append([float(aux[0]), float(aux[1])])
pointsB2 = np.array(pointsB2)

# 2D PLOT (FASTER)
axA1 = plt.subplot(111)
grafA1 = axA1.scatter(pointsA1[:,0], pointsA1[:,1], s=30, color="blue")

axB1 = plt.subplot(111)
grafB1 = axB1.scatter(pointsB1[:,0], pointsB1[:,1], s=30, color="blue")

axA2 = plt.subplot(111)
grafA2 = axA2.scatter(pointsA2[:,0], pointsA2[:,1], s=30, color="orange")

axB2 = plt.subplot(111)
grafB2 = axB2.scatter(pointsB2[:,0], pointsB2[:,1], s=30, color="orange")

try:
    with open(file="latticeAA.dat") as f:
        dataAA = f.readlines()

    pointsAA = []
    for line in dataAA:
        aux = line.split(";")
        pointsAA.append([float(aux[0]), float(aux[1])])
    pointsAA = np.array(pointsAA)
    axAA = plt.subplot(111)
    grafAA = axAA.scatter(pointsAA[:,0], pointsAA[:,1], s=20, color="red")
except:
    print("No AA points found")

try:
    with open(file="latticeAB.dat") as f:
        dataAB = f.readlines()

    pointsAB = []
    for line in dataAB:
        aux = line.split(";")
        pointsAB.append([float(aux[0]), float(aux[1])])
    pointsAB = np.array(pointsAB)
    axAB = plt.subplot(111)
    grafAB = axAB.scatter(pointsAB[:,0], pointsAB[:,1], s=20, color="green")
except:
    print("No AB points found")

try:
    with open(file="latticeBA.dat") as f:
        dataBA = f.readlines()

    pointsBA = []
    for line in dataBA:
        aux = line.split(";")
        pointsBA.append([float(aux[0]), float(aux[1])])
    pointsBA = np.array(pointsBA)
    axBA = plt.subplot(111)
    grafBA = axBA.scatter(pointsBA[:,0], pointsBA[:,1], s=20, color="purple")
except:
    print("No BA points found")

try:
    with open(file="latticeBB.dat") as f:
        dataBB = f.readlines()

    pointsBB = []
    for line in dataBB:
        aux = line.split(";")
        pointsBB.append([float(aux[0]), float(aux[1])])
    pointsBB = np.array(pointsBB)
    axBB = plt.subplot(111)
    grafBB = axBB.scatter(pointsBB[:,0], pointsBB[:,1], s=20, color="yellow")
except:
    print("No BB points found")

axA1.set_aspect(1)
axB1.set_aspect(1)
axA2.set_aspect(1)
axB2.set_aspect(1)

plt.legend(["AA", "AB", "BA", "BB"], loc="upper right")
plt.show()
