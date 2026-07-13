import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

# Wczytanie danych
data = np.loadtxt('../outdir/C.dat')
x = data[:, 0]
y = data[:, 1]
u = data[:, 2]

# Tworzenie wykresu
plt.figure(figsize=(10, 8))

# tricontourf tworzy idealnie gładkie izolinie z chmury punktów
# levels=100 oznacza wysoką płynność przejść kolorów
contour = plt.tricontourf(x, y, u, levels=100, cmap='inferno')

plt.colorbar(contour, label='Rozwiązanie U')
plt.xlabel('Oś X')
plt.ylabel('Oś Y')
plt.title(f'Rozwiązanie MES (Zagęszczone)')
plt.axis('equal') # Ważne, żeby koło nie było jajkiem!
plt.tight_layout()
plt.show()