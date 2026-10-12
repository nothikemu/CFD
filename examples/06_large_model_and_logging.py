"""
Example 6: a big model (>50 000 unknowns) with the full step-by-step log.

Above FEAModel.ITERATIVE_THRESHOLD free DOFs the engine automatically
switches from the direct solver to preconditioned conjugate gradient.
Turning on logging prints every step of the maths.
"""

import logging
import os
import sys

from _setup import OUTPUT

from engine import FEAModel, Material, Mesh
from visualizer import plot_results

logging.basicConfig(level=logging.INFO, stream=sys.stdout, format="%(message)s")
logging.getLogger("matplotlib").setLevel(logging.WARNING)

mesh = Mesh.structured_rectangle(2.0, 0.1, nx=640, ny=40)        # 52 480 free DOFs
model = FEAModel(
    mesh,
    Material(210e9, 0.3, thickness=0.05),
    iterative_threshold=50_000,   # switch point (lower it to force CG on small models)
    cg_rtol=1e-10,                # CG stopping tolerance
    cg_maxiter=None,              # None = automatic
    direct_fallback=True,         # if CG stalls, finish with the direct solver
)
model.fix(mesh.nodes_at(x=0.0), "xy")
model.add_nodal_force(mesh.closest_node(2.0, 0.05), fy=-5000)
r = model.solve()

print("solver used:", r.diagnostics["solver"], "| iterations:", r.diagnostics["iterations"])
plot_results(model, r, filename=os.path.join(OUTPUT, "06_large_model.png"), title="Large model")
