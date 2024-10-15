import numpy as np
import matplotlib.pyplot as plt
h = 1.0 / ( 25 - 0.5 )
x = np.linspace(-h/2, 1 + h/2, 100)  
y = np.linspace(-h/2, 1 + h/2, 100)  
x, y = np.meshgrid(x, y)     

data = np.loadtxt('result.txt')
x_data = data[:, 0]
y_data = data[:, 1]
z_data = data[:, 2]

def f(x, y):
    # return np.cos(2 * np.pi * x) * np.cos(2 * np.pi * y)
    return np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y)
    # return np.cos(2 * np.pi*x) * np.sin(2 * np.pi*y)
    # return np.cos(2 * np.pi * x) * np.cos(2 * np.pi * y)*np.cos(2 * np.pi * x) * np.cos(2 * np.pi * y)

    # return np.cos(2 * np.pi * x) * np.cos(2 * np.pi * y)*np.cos(2 * np.pi * x) * np.cos(2 * np.pi * y)
    # return x*x*(1-x)*(1-x)*y*y*(1-y)*(1-y)
    # return x*x*(1-x*x)*y*y*(1-y*y) * np.sin(np.pi * x) * np.sin(np.pi * y)
    # return x*y 



z = f(x, y)
z_real = f(x_data, y_data)
z_diff = (z_real - z_data)

fig = plt.figure(figsize=(12, 8))
fig.suptitle("sin(2pi x) sin(2pi y)")

ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(x, y, z, cmap='viridis', alpha=0.5)
ax1.scatter(x_data, y_data, z_data, color='r', s=50, label='Data Points')
ax1.set_xlim(-0.1, 1.1)
ax1.set_ylim(-0.1, 1.1)
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')
ax1.set_title("График")



ax2 = fig.add_subplot(122, projection='3d')
ax2.scatter(x_data, y_data, z_diff, color='r', s=50, label='Data Points')
ax2.set_zlim(-1e-8,-1e-5)
ax2.set_title("Разность")


plt.show()