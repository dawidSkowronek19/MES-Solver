import numpy as np
import matplotlib.pyplot as plt

# ==================== Load Data ====================
mode=3


sol_data = np.loadtxt(f"./outdir/psi_{mode}.dat")
x_sol = sol_data[:, 0]
psi_sol = sol_data[:, 1]




# ==================== Plot Solution ====================
plt.figure(2, figsize=(10, 6))

plt.plot(x_sol, psi_sol, color='darkblue', linewidth=2, label='Rozwiązanie MES $psi(x)$', zorder=1)



plt.title(f"Wykres rozwiązania MES $psi_{mode}(x)$", fontsize=14)
plt.xlabel("Pozycja x", fontsize=12)
plt.ylabel("Wartość rozwiązania $psi(x)$", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='best')


plt.axhline(0, color='black', linewidth=0.8, linestyle='-', alpha=0.5)

plt.tight_layout()


plt.show()