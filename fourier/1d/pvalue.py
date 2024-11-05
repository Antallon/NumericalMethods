import numpy as np
import matplotlib.pyplot as plt


N3 = 100
err3 = 7.97245e-07
N4 = 500
err4 = 6.30205e-09
N5 = 1000
err5 = 7.86585e-10
N1 = 10
err1 = 0.000909355
N2 = 50
err2 = 6.47377e-06



p1 = np.log(err1/err2)/np.log(N2/N1)
p2 = np.log(err2/err3)/np.log(N3/N2)



plt.figure(figsize=(10, 6))
plt.scatter([-np.log(1/N1), -np.log(1/N2), -np.log(1/N3), -np.log(1/N4), -np.log(1/N5)], [-np.log(err1), -np.log(err2), -np.log(err3), -np.log(err4),-np.log(err5)],  color='blue')
plt.plot([-np.log(1/N1), -np.log(1/N2), -np.log(1/N3), -np.log(1/N4), -np.log(1/N5)], [-np.log(err1), -np.log(err2), -np.log(err3), -np.log(err4),-np.log(err5)],  color='blue')

plt.xlabel('-np.log(1/N)')
plt.ylabel('-np.log(err)')
