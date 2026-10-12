"""
Example 1: the simplest possible job. A cantilever (diving-board) beam.

    wall |=================================| <- 2000 N down at the tip
         x = 0                          x = 1.0 m

Change the numbers in the SETTINGS block and run the file again.
"""

import os

from _setup import OUTPUT

from engine import FEAModel, Material, Mesh
from visualizer import plot_results

# ----------------------------- SETTINGS (edit me) -----------------------------
LENGTH = 1.0            # m
HEIGHT = 0.1            # m
THICKNESS = 0.02        # m (out-of-plane depth)
E = 210e9               # Pa  (steel 210e9, aluminium 70e9, titanium 114e9)
NU = 0.3                # Poisson's ratio
YIELD = 250e6           # Pa  (used for the safety factor)
FORCE = -2000.0         # N   (negative = downward)
NX, NY = 100, 10        # elements along the length / through the height
# ------------------------------------------------------------------------------

mesh = Mesh.structured_rectangle(length=LENGTH, height=HEIGHT, nx=NX, ny=NY)
steel = Material(E=E, nu=NU, thickness=THICKNESS, plane="stress", yield_strength=YIELD)
model = FEAModel(mesh, steel)

model.fix(mesh.nodes_at(x=0.0), "xy")                       # clamp the left edge
tip = mesh.closest_node(LENGTH, HEIGHT / 2)                 # middle of the right edge
model.add_nodal_force(tip, fy=FORCE)                        # push it

r = model.solve()

I = THICKNESS * HEIGHT**3 / 12
hand_calc = abs(FORCE) * LENGTH**3 / (3 * E * I)
print(f"tip drop (FEA)        : {r.displacements[tip, 1]:+.6e} m")
print(f"tip drop (hand calc)  : {-hand_calc:+.6e} m   (P L^3 / 3EI)")
print(f"max Von Mises stress  : {r.max_von_mises / 1e6:.2f} MPa")
print(f"minimum safety factor : {r.safety_factor.min():.2f}")

plot_results(model, r, filename=os.path.join(OUTPUT, "01_my_beam.png"), title="My beam")
print("picture saved to output/01_my_beam.png")
