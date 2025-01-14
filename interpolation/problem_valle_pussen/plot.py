import numpy as np
import matplotlib.pyplot as plt


# def f(x):
#     return 1.0/(1.0 + 25*x*x)

# x_range = np.linspace(a, b, 1000)
# y_range = f(x_range)
data = np.loadtxt("result.txt")
x_data = data[:, 0]

y_fun = data[:, 1]
y_poly = data[:, 2]
diff = data[:, 3]
max_error = diff.max()
print(max_error)




_, axs = plt.subplots(2,1, figsize=(10, 10))

axs[0].plot(x_data, y_fun,  color='red')
axs[0].scatter(x_data,y_poly, color='green')

# axs[1].scatter(x_data, diff,  color='green')


plt.show()