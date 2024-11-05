import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("result.txt")
x_data = data[:, 0]
y_data = data[:, 1]
plt.plot(x_data, y_data, color='blue')
plt.show()