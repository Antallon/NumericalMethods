import numpy as np
import matplotlib.pyplot as plt
with open('result.txt', 'r') as file:
    lines = file.readlines()
with open('result_new.txt', 'w') as file:
    file.writelines(lines[1:])



# a, b = 0, 30
a, b = -3, 3
# a,b = -10, 10


def f(x):
    return 1.0/(1.0 + 25*x*x)
    # return np.cos(x)
     #return np.abs(x)
    # return pow(x,10) + 5.0*pow(x,5) - 2.0*pow(x,2) + 7.0*x -15.0





data = np.loadtxt("result_new.txt")
x_data = data[:, 0]
y_polinomial = data[:, 2]
y_lagrange = data[:, 4]
y_mnrp = data[:, 6]
y_vp = data[:, 8]

x_range = np.linspace(a, b, 1000)
y_range = f(x_range)


fig, axs = plt.subplots(4,2, figsize=(12, 12))
fig.suptitle('1/(1 + x^2)')

axs[0,0].plot(x_range, y_range,  color='blue')
axs[0,0].scatter(x_data, y_polinomial, color='red')
axs[0,0].set_xlabel('x')
axs[0,0].set_ylabel('y')
axs[0,0].set_title('P(x_k)')

axs[0,1].scatter(x_data, data[:,3],  color='green')
# axs[1].scatter(x_data, y_data, color='red')
axs[0,1].set_xlabel('x')
axs[0,1].set_ylabel('y')
axs[0,1].set_title('|y(x_k) - P(x_k)| ')

axs[1,0].plot(x_range, y_range,  color='red')
axs[1,0].scatter(x_data,y_lagrange, color='purple')
axs[1,0].set_xlabel('x')
axs[1,0].set_ylabel('y')
axs[1,0].set_title('L(x_k)')

axs[1,1].scatter(x_data, data[:,5],  color='red')
# axs[3].scatter(x_data, y_data, color='red')
axs[1,1].set_xlabel('x')
axs[1,1].set_ylabel('y')
axs[1,1].set_title('|y(x_k) - L(x_k)|')

axs[2,0].plot(x_range, y_range,  color='red')
axs[2,0].scatter(x_data,y_mnrp, color='green')
axs[2,0].set_xlabel('x')
axs[2,0].set_ylabel('y')
axs[2,0].set_title('MNRP(x_k)')

axs[2,1].scatter(x_data, data[:,7],  color='green')
# axs[3].scatter(x_data, y_data, color='red')
axs[2,1].set_xlabel('x')
axs[2,1].set_ylabel('y')
axs[2,1].set_title('|y(x_k) - MNRP(x_k)|')

axs[3,0].plot(x_range, y_range,  color='red')
axs[3,0].scatter(x_data,y_vp, color='brown')
axs[3,0].set_xlabel('x')
axs[3,0].set_ylabel('y')
axs[3,0].set_title('VP(x_k)')

axs[3,1].scatter(x_data, data[:,9],  color='blue')
# axs[3].scatter(x_data, y_data, color='red')
axs[3,1].set_xlabel('x')
axs[3,1].set_ylabel('y')
axs[3,1].set_title('|y(x_k) - VP(x_k)|')



plt.tight_layout()
plt.show()