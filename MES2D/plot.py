import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import os

# Ścieżki do plików wygenerowanych przez C++
points_file = './outdir/pointsList.dat'
triangles_file = './outdir/trianglesPointsList.dat'

# Sprawdzenie czy pliki istnieją
if not os.path.exists(points_file) or not os.path.exists(triangles_file):
    print("Błąd: Nie znaleziono plików. Upewnij się, że program C++ wykonał się poprawnie.")
    exit()

# 1. Wczytanie danych
# points: tablica 2D, gdzie każda kolumna to [x, y]
points = np.loadtxt(points_file)
x = points[:, 0]
y = points[:, 1]

# triangles: tablica 2D z indeksami węzłów [p1, p2, p3]
# Wymuszamy typ całkowity (int), bo to są indeksy
triangles = np.loadtxt(triangles_file, dtype=int)

# 2. Utworzenie obiektu triangulacji
# Łączy on współrzędne węzłów z macierzą topologii elementów
triang = mtri.Triangulation(x, y, triangles)

# 3. Rysowanie wykresu
plt.figure(figsize=(10, 8))

# Rysowanie krawędzi trójkątów (siatka)
plt.triplot(triang, 'k-', lw=1.0, alpha=0.7, label='Krawędzie elementów')

# Rysowanie węzłów (punkty)
plt.plot(x, y, 'ro', markersize=3, label='Węzły (punkty)')

# Opcjonalnie: dodanie numeracji węzłów na wykresie (przydatne przy małych siatkach do debugowania)
# for i, (px, py) in enumerate(zip(x, y)):
#     plt.text(px, py, str(i), color='blue', fontsize=8, ha='bottom', va='bottom')

# Ustawienia wizualne
plt.title('Siatka trójkątna MES')
plt.xlabel('Oś X')
plt.ylabel('Oś Y')
plt.gca().set_aspect('equal') # Zapewnia brak zniekształceń proporcji
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()

# Wyświetlenie wyniku
plt.show()