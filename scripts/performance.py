import matplotlib.pyplot as plt

times = [    28,    167,   1216,    5003,   32720]
sizes = [100000, 200000, 500000, 1000000, 2500000]

ax1 = plt.subplot(111)
graf1 = ax1.scatter(sizes, times)
graf2 = ax1.plot(sizes, times)

plt.show()
