import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from matplotlib.animation import FuncAnimation


files = sorted(glob.glob('../outdir/u/sol_*.dat'), 
               key=lambda x: int(os.path.basename(x).split('_')[1].split('.')[0]))


data_first = np.loadtxt(files[0])
vmin, vmax = data_first[:, 2].min(), data_first[:, 2].max()

fig, ax = plt.subplots(figsize=(10, 8))

def update(frame_file):
    ax.clear()
    data = np.loadtxt(frame_file)
    x, y, u = data[:, 0], data[:, 1], data[:, 2]
    

    contour = ax.tricontourf(x, y, u, levels=100, cmap='inferno', vmin=vmin, vmax=vmax)
    
    iter_nb = os.path.basename(frame_file).split('_')[1].split('.')[0]
    ax.set_title(f'Rozwiązanie MES - Iteracja: {iter_nb}')
    ax.axis('equal')
    return contour


ani = FuncAnimation(fig, update, frames=files, interval=100, repeat=False)

plt.show()


ani.save('symulacja_mes.mp4', writer='ffmpeg')