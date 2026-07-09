import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import os

points_file = '../outdir/pointsList.dat'
triangles_file = '../outdir/trianglesPointsList.dat'




points = np.loadtxt(points_file)
x = points[:, 0]
y = points[:, 1]


triangles = np.loadtxt(triangles_file, dtype=int)


triang = mtri.Triangulation(x, y, triangles)


plt.figure(figsize=(10, 8))


plt.triplot(triang, 'k-', lw=1.0, alpha=0.7, label='edges')


plt.plot(x, y, 'ro', markersize=3, label='nodes')


plt.title('Triangular mesh')
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()

plt.savefig("../graphs/triangular_circle.png")