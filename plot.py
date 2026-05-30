import numpy as np
import matplotlib.pyplot as plt

# ==================== WCZYTANIE DANYCH ====================
# Wczytanie siatki (pozycji węzłów)
mesh_data = np.loadtxt("./outdir_stationary/mesh.dat", skiprows=1)
x_nodes = mesh_data[:, 1]
y_nodes_zero = np.zeros_like(x_nodes)

# Wczytanie rozwiązania (kolumna 0: x, kolumna 1: u(x))
# Plik u.dat zawiera gęste punkty wygenerowane przez m_dxSave
sol_data = np.loadtxt("./outdir_stationary/u.dat")
x_sol = sol_data[:, 0]
u_sol = sol_data[:, 1]


# ==================== WYKRES 1: GEOMETRIA SIATKI ====================
plt.figure(1, figsize=(12, 2))

# Rysowanie ciągłej linii pręta i węzłów
plt.plot([x_nodes[0], x_nodes[-1]], [0, 0], color='blue', linestyle='-', zorder=1)
plt.scatter(x_nodes, y_nodes_zero, color='red', s=30, label='Węzły siatki', zorder=2)

plt.title("Rozkład węzłów na pręcie")
plt.xlabel("Pozycja x")
plt.yticks([])  # Ukrywamy oś Y
plt.grid(True, axis='x', linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig("./outdir_stationary/nodes.png")


# ==================== WYKRES 2: ROZWIĄZANIE U(X) ====================
plt.figure(2, figsize=(10, 6))

# 1. Rysujemy ciągły wykres funkcji u(x) z pliku u.dat
plt.plot(x_sol, u_sol, color='darkblue', linewidth=2, label='Rozwiązanie MES $u(x)$', zorder=1)

# 2. Rysujemy "kropki" węzłów na osi OX (tam gdzie u=0)
plt.scatter(x_nodes, y_nodes_zero, color='red', s=25, marker='o', label='Położenie węzłów ($u=0$)', zorder=2)

# Opcjonalnie: Jeśli chcesz zobaczyć wartości u(x) dokładnie W WĘZŁACH, 
# możesz odkomentować poniższą linijkę (wykorzystuje funkcję R_U_1D z MES):
# plt.scatter(x_nodes, [u_sol[np.abs(x_sol - x).argmin()] for x in x_nodes], color='green', s=40, zorder=3, label='Wartości w węzłach')

# Ustawienia drugiego wykresu
plt.title("Wykres rozwiązania MES $u(x)$", fontsize=14)
plt.xlabel("Pozycja x", fontsize=12)
plt.ylabel("Wartość rozwiązania $u$", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='best')

# Dodatkowe wyróżnienie osi OX na drugim wykresie
plt.axhline(0, color='black', linewidth=0.8, linestyle='-', alpha=0.5)

plt.tight_layout()

# Pokazanie obu okien na raz
plt.savefig("./outdir_stationary/u(x).png")