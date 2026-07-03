import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.style.use('dark_background')

def sh_funNb(file_path):
    dir = Path(file_path)
    if not dir.exists() or not dir.is_dir():
        print("Error")
        return 0
    number = sum(1 for element in dir.iterdir()
                 if element.is_file() and "phi_" in element.name)
    
    return number//2

for idx in range(sh_funNb("../outdir/")):
    data = np.loadtxt(f"../outdir/phi_{idx}.dat")
    ksi = data[:,0]
    eta = data[:,1]
    phi = data[:,2]

    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_trisurf(ksi, eta, phi, cmap='viridis', edgecolor='none', antialiased=True)
    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, pad=0.15)
    cbar.set_label(f"\u03D5_{idx} (\u03BE, \u03B7)")
    ax.set_axis_off()

    ax.set_title("Shape function in elemental triangle")
    plt.show()