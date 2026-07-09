import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

print("Wczytywanie danych z ../outdir/C.dat ...")

try:

    x, y, u = np.loadtxt('../outdir/C.dat', delimiter=',', unpack=True)
except OSError:
    print("Błąd: Nie znaleziono pliku ../outdir/C.dat! Sprawdź ścieżkę.")
    exit()


triang = mtri.Triangulation(x, y)


plt.figure(figsize=(8, 8))
contour = plt.tricontourf(triang, u, levels=50, cmap='inferno')


plt.colorbar(contour, label='Wartość u (np. potencjał/temperatura)')
plt.title('Rozwiązanie MES (Równanie Laplace\'a)')
plt.xlabel('Oś X')
plt.ylabel('Oś Y')

plt.axis('equal') 

plt.tight_layout()
plt.show()