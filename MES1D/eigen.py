import numpy as np
import matplotlib.pyplot as plt

# ==================== WCZYTANIE DANYCH ====================
mode=196

# Wczytanie rozwiązania (kolumna 0: x, kolumna 1: u(x))
# Plik u.dat zawiera gęste punkty wygenerowane przez m_dxSave
sol_data = np.loadtxt(f"./outdir_eigen/psi_{mode}.dat")
x_sol = sol_data[:, 0]
psi_sol = sol_data[:, 1]




# ==================== WYKRES 2: ROZWIĄZANIE U(X) ====================
plt.figure(2, figsize=(10, 6))

# 1. Rysujemy ciągły wykres funkcji u(x) z pliku u.dat
plt.plot(x_sol, psi_sol, color='darkblue', linewidth=2, label='Rozwiązanie MES $psi(x)$', zorder=1)



plt.title(f"Wykres rozwiązania MES $psi_{mode}(x)$", fontsize=14)
plt.xlabel("Pozycja x", fontsize=12)
plt.ylabel("Wartość rozwiązania $psi(x)$", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='best')

# Dodatkowe wyróżnienie osi OX na drugim wykresie
plt.axhline(0, color='black', linewidth=0.8, linestyle='-', alpha=0.5)

plt.tight_layout()

# Pokazanie obu okien na raz
plt.show()