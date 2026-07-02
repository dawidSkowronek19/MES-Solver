import numpy as np
import matplotlib.pyplot as plt


mesh_data = np.loadtxt("./outdir_stationary/mesh.dat", skiprows=1)
x_nodes = mesh_data[:, 1]
y_nodes_zero = np.zeros_like(x_nodes)

sol_data = np.loadtxt("./outdir_stationary/u0.000000.dat")
x_sol = sol_data[:, 0]
u_sol = sol_data[:, 1]


plt.figure(1, figsize=(12, 2))


plt.plot([x_nodes[0], x_nodes[-1]], [0, 0], color='blue', linestyle='-', zorder=1)
plt.scatter(x_nodes, y_nodes_zero, color='red', s=30, label='Węzły siatki', zorder=2)

plt.title("Rozkład węzłów na pręcie")
plt.xlabel("Pozycja x")
plt.yticks([])
plt.grid(True, axis='x', linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig("./outdir_stationary/nodes.png")


plt.figure(2, figsize=(10, 6))


plt.plot(x_sol, u_sol, color='darkblue', linewidth=2, label='$u(x)$', zorder=1)

plt.scatter(x_nodes, y_nodes_zero, color='red', s=25, marker='o', label='nodes', zorder=2)

plt.xlabel("x", fontsize=12)
plt.ylabel("$u$", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(loc='best')

plt.axhline(0, color='black', linewidth=0.8, linestyle='-', alpha=0.5)

plt.tight_layout()
plt.savefig("./u(x).png")