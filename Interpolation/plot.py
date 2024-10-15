import numpy as np
import matplotlib.pyplot as plt
with open('result.txt', 'r') as file:
    lines = file.readlines()
with open('result_new.txt', 'w') as file:
    file.writelines(lines[1:])



# a, b = 0, 30
a, b = -5, 5
# a,b = -10, 10


def f(x):
    # return 1.0/(1.0 + x*x)
    # return np.cos(x)
    return np.abs(x)
    # return pow(x,10) + 5.0*pow(x,5) - 2.0*pow(x,2) + 7.0*x -15.0





data = np.loadtxt("result_new.txt")
x_data = data[:, 0]
y_polinomial = data[:, 2]
y_lagrange = data[:, 4]

x_range = np.linspace(a, b, 1000)
y_range = f(x_range)

fig, axs = plt.subplots(4, figsize=(10, 12))
fig.suptitle('1/(1 + x^2)')

axs[0].plot(x_range, y_range,  color='blue')
axs[0].scatter(x_data, y_polinomial, color='red')
axs[0].set_xlabel('x')
axs[0].set_ylabel('y')
axs[0].set_title('P(x_k)')

axs[1].scatter(x_data, data[:,3],  color='green')
# axs[1].scatter(x_data, y_data, color='red')
axs[1].set_xlabel('x')
axs[1].set_ylabel('y')
axs[1].set_title('|y(x_k) - P(x_k)| ')

axs[2].plot(x_range, y_range,  color='red')
axs[2].scatter(x_data,y_lagrange, color='purple')
axs[2].set_xlabel('x')
axs[2].set_ylabel('y')
axs[2].set_title('L(x_k)')

axs[3].scatter(x_data, data[:,5],  color='black')
# axs[3].scatter(x_data, y_data, color='red')
axs[3].set_xlabel('x')
axs[3].set_ylabel('y')
axs[3].set_title('|y(x_k) - L(x_k)|')

plt.tight_layout()
plt.show()