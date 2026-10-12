"""
Example 4: build your own mesh from raw coordinates and connectivity.

A tapered bracket: 0.4 m tall at the wall, 0.1 m tall at the free end.
Nodes are generated row by row; each element lists its 4 nodes
COUNTER-CLOCKWISE (bottom-left, bottom-right, top-right, top-left).

Also shows: plane strain, nodes_where() with a custom rule, a sloped load.
"""

import os

import numpy as np

from _setup import OUTPUT

from engine import FEAModel, Material, Mesh
from visualizer import plot_results

L, H_ROOT, H_TIP = 1.0, 0.4, 0.1
NX, NY = 40, 8

nodes = []
for j in range(NY + 1):                       # rows bottom -> top
    for i in range(NX + 1):                   # columns left -> right
        x = L * i / NX
        h = H_ROOT + (H_TIP - H_ROOT) * x / L
        y = -h / 2 + h * j / NY
        nodes.append([x, y])

elements = []
for j in range(NY):
    for i in range(NX):
        n0 = j * (NX + 1) + i                 # bottom-left
        elements.append([n0, n0 + 1, n0 + NX + 2, n0 + NX + 1])   # CCW

mesh = Mesh(nodes, elements)                  # raises MeshTopologyError if anything is clockwise
print(mesh)

concrete = Material(E=30e9, nu=0.2, thickness=1.0, plane="strain", name="concrete")
model = FEAModel(mesh, concrete)

wall = mesh.nodes_where(lambda x, y: np.isclose(x, 0.0))
model.fix(wall, "xy")
tip = mesh.closest_node(L, 0.0)
model.add_nodal_force(tip, fx=20e3, fy=-50e3)  # pushes right and down

r = model.solve()
print(f"tip displacement: ux = {r.displacements[tip, 0]:+.4e} m, uy = {r.displacements[tip, 1]:+.4e} m")
print(f"max Von Mises   : {r.max_von_mises / 1e6:.3f} MPa")
plot_results(model, r, filename=os.path.join(OUTPUT, "04_custom_mesh.png"),
             title="Tapered bracket (plane strain)", stress_unit="Pa")
print("picture saved to output/04_custom_mesh.png")
