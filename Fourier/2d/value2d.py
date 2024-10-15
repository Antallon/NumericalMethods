import numpy as np
import matplotlib.pyplot as plt

N1 = 10
err1 = 0.000200179

N2 = 15
err2 = 5.64538e-05

N3 = 25
err3 = 1.17379e-05

N4 = 50
err4 = 1.42428e-06

N5 = 100
err5 = 1.75412e-07





p1 = np.log(err1/err2)/np.log(N2/N1)
p2 = np.log(err2/err3)/np.log(N3/N2)
p3 = np.log(err3/err4)/np.log(N4/N3)
p4 = np.log(err4/err5)/np.log(N5/N4)
print(p1,p2,p3,p4)



plt.figure(figsize=(10, 6))
plt.scatter([-np.log(1/N1), -np.log(1/N2), -np.log(1/N3), -np.log(1/N4), -np.log(1/N5)], [-np.log(err1), -np.log(err2), -np.log(err3), -np.log(err4),-np.log(err5)],  color='blue')
plt.plot([-np.log(1/N1), -np.log(1/N2), -np.log(1/N3), -np.log(1/N4), -np.log(1/N5)], [-np.log(err1), -np.log(err2), -np.log(err3), -np.log(err4),-np.log(err5)],  color='blue')

plt.xlabel('-np.log(1/N)')
plt.ylabel('-np.log(err)')
plt.show()