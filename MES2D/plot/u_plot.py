import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


data = np.loadtxt('../outdir/C.dat')
x = data[:, 0]
y = data[:, 1]
u = data[:, 2]


plt.figure(figsize=(10, 8))


contour = plt.tricontourf(x, y, u, levels=100, cmap='inferno')

plt.colorbar(contour, label='Rozwiązanie U')
plt.xlabel('Oś X')
plt.ylabel('Oś Y')
plt.title(f'Rozwiązanie MES (Zagęszczone)')
plt.axis('equal')
plt.tight_layout()
plt.show()