import numpy as np
import matplotlib.pyplot as plt

def f(x):
    # return x * (1 - x**2) * np.sin( np.pi * x)
    # return x * (1 - x) * np.log(1 + 10*x)
    # return np.cos(2* np.pi * x)/2 + np.cos(np.pi*x)
    # return np.sin( np.pi * x)
    return x**2 * (1 - x**2) * np.sin(np.pi * x)
    # return np.cos( np.pi * x) * np.cos( np.pi * x)
    # return x**2 * (1 - x)**2 #+ np.cos(np.pi * x)
    # return np.cos( 2*np.pi * x)

data = np.loadtxt("result.txt")
x_data = data[:, 0]
y_data = data[:, 1]

x_range = np.linspace(0, 1, 1000)
y_range = f(x_range)

plt.figure(figsize=(10, 6))
plt.plot(x_range, y_range,  color='blue')
plt.scatter(x_data, y_data, color='red')
plt.xlabel('x')
plt.ylabel('y')
plt.title('cos^2(pi x)')
plt.show()
