"""
Example 2: a beam resting on two supports with pressure on its top face.

        p (uniform pressure, N/m^2) pushing down on the top edge
     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
     |============================================|
     ^ pin (x and y held)                roller ^ (y held only)

Shows: add_edge_traction (distributed load), pin vs roller supports,
reactions, and comparison with 5 q L^4 / (384 E I).
"""

import os

import numpy as np

from _setup import OUTPUT

from engine import FEAModel, Material, Mesh
from visualizer import plot_results, plot_stress_components

L, H, T = 3.0, 0.15, 0.05
E, NU = 200e9, 0.3
PRESSURE = 2.0e5                              # Pa on the top face

mesh = Mesh.structured_rectangle(L, H, nx=240, ny=12, origin=(0.0, -H / 2))
model = FEAModel(mesh, Material(E, NU, thickness=T, yield_strength=355e6))

left = mesh.closest_node(0.0, 0.0)            # supports on the neutral axis
right = mesh.closest_node(L, 0.0)
model.fix(left, "xy")                         # pin
model.fix(right, "y")                         # roller

top = mesh.nodes_at(y=H / 2)
top = top[np.argsort(mesh.nodes[top, 0])]     # traction chains must be ordered along the edge
model.add_edge_traction(top, (0.0, -PRESSURE))

r = model.solve()

q = PRESSURE * T                              # load per unit length, N/m
I = T * H**3 / 12
mid = mesh.closest_node(L / 2, 0.0)
print(f"mid-span deflection FEA  : {-r.displacements[mid, 1]:.6e} m")
print(f"mid-span deflection beam : {5 * q * L**4 / (384 * E * I):.6e} m  (5 q L^4 / 384 E I)")
print(f"support reactions (y)    : left {r.reactions[2 * left + 1]:.1f} N, "
      f"right {r.reactions[2 * right + 1]:.1f} N, expected {q * L / 2:.1f} N each")
print(f"max Von Mises            : {r.max_von_mises / 1e6:.2f} MPa")

plot_results(model, r, filename=os.path.join(OUTPUT, "02_simply_supported.png"),
             title="Simply supported beam, uniform pressure")
plot_stress_components(model, r, filename=os.path.join(OUTPUT, "02_simply_supported_components.png"))
print("pictures saved to output/02_simply_supported*.png")
