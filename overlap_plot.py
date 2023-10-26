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

# with open(file="latticeA2.dat") as f:
#     dataA2 = f.readlines()
# 
# pointsA2 = []
# for line in dataA2:
#     aux = line.split(";")
#     pointsA2.append([float(aux[0]), float(aux[1]), float(aux[2])])
# pointsA2 = np.array(pointsA2)
# 
# with open(file="latticeB2.dat") as f:
#     dataB2 = f.readlines()
# 
# pointsB2 = []
# for line in dataB2:
#     aux = line.split(";")
#     pointsB2.append([float(aux[0]), float(aux[1]), float(aux[2])])
# pointsB2 = np.array(pointsB2)

# 2D PLOT (FASTER)
axA1 = plt.subplot(111)
grafA1 = axA1.scatter(pointsA1[:,0], pointsA1[:,1], s=30, color="blue")

axB1 = plt.subplot(111)
grafB1 = axB1.scatter(pointsB1[:,0], pointsB1[:,1], s=30, color="orange")

# axA2 = plt.subplot(111)
# grafA2 = axA2.scatter(pointsA2[:,0], pointsA2[:,1], s=30, color="red")
#
# axB2 = plt.subplot(111)
# grafB2 = axB2.scatter(pointsB2[:,0], pointsB2[:,1], s=30, color="green")

try:
    with open(file="latticeAA.dat") as f:
        dataAA = f.readlines()

    pointsAA = []
    for line in dataAA:
        aux = line.split(";")
        pointsAA.append([float(aux[0]), float(aux[1]), float(aux[2])])
    pointsAA = np.array(pointsAA)
    axAA = plt.subplot(111)
    grafAA = axAA.scatter(pointsAA[:,0], pointsAA[:,1], s=20, color="red")
except:
    print("No AA or BB points found")

try:
    with open(file="latticeAB.dat") as f:
        dataAB = f.readlines()

    pointsAB = []
    for line in dataAB:
        aux = line.split(";")
        pointsAB.append([float(aux[0]), float(aux[1]), float(aux[2])])
    pointsAB = np.array(pointsAB)
    axAB = plt.subplot(111)
    grafAB = axAB.scatter(pointsAB[:,0], pointsAB[:,1], s=20, color="green")
except:
    print("No AB or BA points found")

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
# axA2.set_aspect(1)
# axB2.set_aspect(1)

plt.show()
