import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('result2.txt')

fig, axs = plt.subplots(1, figsize=(6, 6))

axs.plot(data[:, 0], data[:, 1], color='g', label = '$Residual$', marker = '+', markersize = 3, alpha = 0.8, linewidth = 0.8, markeredgecolor = 'red')
axs.plot(data[:, 0], data[:, 2], color = 'blue', label = 'Theoretical residual', marker = '+', markersize = 2, alpha = 0.7, markeredgecolor = 'purple')
axs.set_xlabel('iter')
axs.set_ylabel('Resid')
axs.set_title('Погрешность')
axs.grid(alpha=0.2, linestyle = '--')
axs.legend()

fig.tight_layout()
plt.show()