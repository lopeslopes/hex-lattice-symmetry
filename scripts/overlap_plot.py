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

# 2D PLOT (FASTER)
axA1 = plt.subplot(111)
grafA1 = axA1.scatter(pointsA1[:,0], pointsA1[:,1], s=30, color="blue")

axB1 = plt.subplot(111)
grafB1 = axB1.scatter(pointsB1[:,0], pointsB1[:,1], s=30, color="blue")

try:
    with open(file="latticeOA.dat") as f:
        dataOVA = f.readlines()

    pointsOVA = []
    for line in dataOVA:
        aux = line.split(";")
        pointsOVA.append([float(aux[0]), float(aux[1]), float(aux[2])])
    pointsOVA = np.array(pointsOVA)
    axOVA = plt.subplot(111)
    grafOVA = axOVA.scatter(pointsOVA[:,0], pointsOVA[:,1], s=20, color="green")
except:
    print("No overlap close enough")


try:
    with open(file="latticeOB.dat") as f:
        dataOVB = f.readlines()

    pointsOVB = []
    for line in dataOVB:
        aux = line.split(";")
        pointsOVB.append([float(aux[0]), float(aux[1]), float(aux[2])])
    pointsOVB = np.array(pointsOVB)
    axOVB = plt.subplot(111)
    grafOVB = axOVB.scatter(pointsOVB[:,0], pointsOVB[:,1], s=20, color="red")
except:
    print("No overlap close enough")

axA1.set_aspect(1)
axB1.set_aspect(1)

plt.legend(["A1", "B1", "Ov_A", "Ov_B"])
plt.show()
