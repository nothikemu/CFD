"""
Example 3: a plate with a round hole pulled from both sides (stress concentration).

Only a quarter is modelled; symmetry is enforced with rollers:
    x = 0 edge : u_x = 0        y = 0 edge : u_y = 0
The curved, graded mesh comes from Mesh.plate_with_hole_quarter.
"""

import os

import numpy as np

from _setup import OUTPUT

from engine import FEAModel, Material, Mesh
from visualizer import plot_results

A, W, T = 0.01, 0.10, 0.005                  # hole radius, half-width, thickness [m]
S0 = 60e6                                    # remote tension [Pa]

mesh = Mesh.plate_with_hole_quarter(radius=A, half_width=W, n_radial=40,
                                    n_circumferential=48, grading=8.0)
aluminium = Material(E=70e9, nu=0.33, thickness=T, yield_strength=275e6, name="Al 6061-T6")
model = FEAModel(mesh, aluminium)

model.fix(mesh.nodes_at(x=0.0), "x")
model.fix(mesh.nodes_at(y=0.0), "y")
edge = mesh.nodes_at(x=W)
edge = edge[np.argsort(mesh.nodes[edge, 1])]
model.add_edge_traction(edge, (S0, 0.0))

r = model.solve()
crown = mesh.closest_node(0.0, A)
print(f"stress concentration factor sigma_xx(0, a) / s0 = {r.nodal_stresses[crown, 0] / S0:.3f} "
      "(textbook ~3 for an infinite plate)")
print(f"minimum safety factor = {r.safety_factor.min():.2f}")
plot_results(model, r, filename=os.path.join(OUTPUT, "03_plate_with_hole.png"),
             title="Plate with hole", scale=150)
print("picture saved to output/03_plate_with_hole.png")
