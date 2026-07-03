import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib import cm
import matplotlib.tri as mtri

plt.style.use('dark_background')
num_frames=180

def sh_funNb(file_path):
    dir = Path(file_path)
    if not dir.exists() or not dir.is_dir():
        print("Error")
        return 0
    number = sum(1 for element in dir.iterdir()
                 if element.is_file() and "phi_" in element.name)
    
    return number//2



for idx in range(sh_funNb("../outdir/")):
    data_phi = np.loadtxt(f"../outdir/phi_{idx}.dat")
    data_dphi = np.loadtxt(f"../outdir/dphi_{idx}.dat")
    ksi = data_phi[:,0]
    eta = data_phi[:,1]
    phi = data_phi[:,2]
    dphi_dksi=data_dphi[:,2]
    dphi_deta=data_dphi[:,3]

    triang = mtri.Triangulation(ksi, eta)
    x0, y0 = ksi[triang.triangles[:, 0]], eta[triang.triangles[:, 0]]
    x1, y1 = ksi[triang.triangles[:, 1]], eta[triang.triangles[:, 1]]
    x2, y2 = ksi[triang.triangles[:, 2]], eta[triang.triangles[:, 2]]


    L1_sq = (x0 - x1)**2 + (y0 - y1)**2
    L2_sq = (x1 - x2)**2 + (y1 - y2)**2
    L3_sq = (x2 - x0)**2 + (y2 - y0)**2

    max_edge_sq = np.max([L1_sq, L2_sq, L3_sq], axis=0)

    mask = max_edge_sq > 0.05**2
    triang.set_mask(mask)

    fig = plt.figure(figsize=(12,8))
    fig.suptitle("Shape Functions in reference element")
    


    plot_data = [
    (phi, r"$\phi$", r"$\phi(\xi, \eta)$"),
    (dphi_dksi, r"$\partial_{\xi}\phi$", r"$\partial_{\xi}\phi(\xi, \eta)$"),
    (dphi_deta, r"$\partial_{\eta}\phi$", r"$\partial_{\eta}\phi(\xi, \eta)$")
    ]

    axs =[]
    surfs=[]

    for i, (data, title, label) in enumerate(plot_data):
        ax = fig.add_subplot(1,3,1+i, projection='3d')
        ax.set_title(title, fontsize=12)

        surf = ax.plot_trisurf(triang, data, cmap='viridis', edgecolor='none', antialiased=True)
        cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, pad=0.1)

        ax.set_xlim(-1,1)
        ax.set_ylim(-1,1)
        ax.set_zlim(np.min(data), np.max(data))

        ax.set_axis_off()
        axs.append(ax)
        surfs.append(surf)
    plt.tight_layout(rect=[0,0.03, 1, 0.95])

    def update(frame):
        angle = frame*(360/num_frames)
        for ax in axs:
            ax.view_init(elev=30, azim=angle)
        
        return surfs

    anim = FuncAnimation(fig, update, frames=num_frames, interval=70, blit=False)
    anim.save(f"../graphs/phi_{idx}.gif")
    plt.close()

