import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("result.txt")
x_data = data[:, 0]
y1_data = data[:, 1]
y2_data = data[:, 2]
y3_data = data[:, 3]
y4_data = data[:, 4]
y5_data = data[:, 5]
y6_data = data[:, 6]


y1_data = -y1_data
y2_data = -y2_data
y3_data = -y3_data  
y4_data = -y4_data
y5_data = -y5_data
y6_data = -y6_data

y1_data = y1_data - y1_data[0]
y2_data = y2_data - y2_data[0]
y3_data = y3_data - y3_data[0]
y4_data = y4_data - y4_data[0]
y5_data = y5_data - y5_data[0]
y6_data = y6_data - y6_data[0]





fig, axs = plt.subplots(3, 2)

axs[0, 0].plot(x_data, y1_data, color='blue')
axs[0, 1].plot(x_data, y2_data, color='blue')
axs[1, 0].plot(x_data, y3_data, color='blue')
axs[1, 1].plot(x_data, y4_data, color='blue')
axs[2, 0].plot(x_data, y5_data, color='blue')
axs[2, 1].plot(x_data, y6_data, color='blue')

plt.show()