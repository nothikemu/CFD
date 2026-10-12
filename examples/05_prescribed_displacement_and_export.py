"""
Example 5: move a boundary by a set amount (instead of pushing with a force),
read the reaction force needed, and export every result to CSV files.

A 0.2 x 0.1 m block, left edge held, right edge pulled 0.1 mm to the right.
"""

import os

import numpy as np

from _setup import OUTPUT

from engine import FEAModel, Material, Mesh

mesh = Mesh.structured_rectangle(0.2, 0.1, nx=20, ny=10)
model = FEAModel(mesh, Material(E=70e9, nu=0.33, thickness=0.01))

model.fix(mesh.nodes_at(x=0.0), "x")              # left edge cannot move sideways
model.fix(mesh.closest_node(0.0, 0.0), "y")       # one point stops vertical sliding
right = mesh.nodes_at(x=0.2)
model.fix(right, "x", value=1e-4)                 # right edge moved +0.1 mm

r = model.solve()
pull = r.reactions[2 * right].sum()
print(f"force needed to stretch the block 0.1 mm: {pull:.1f} N")

# ---- export nodal results
table = np.column_stack([
    np.arange(mesh.n_nodes), mesh.nodes, r.displacements,
    r.nodal_stresses, r.von_mises, r.reactions.reshape(-1, 2),
])
header = "node,x,y,ux,uy,sxx,syy,txy,von_mises,Rx,Ry"
path = os.path.join(OUTPUT, "05_nodal_results.csv")
np.savetxt(path, table, delimiter=",", header=header, comments="", fmt="%.10g")

# ---- export Gauss-point results
gp = np.column_stack([
    np.repeat(np.arange(mesh.n_elements), 4), r.gauss_points.reshape(-1, 2),
    r.gauss_strains.reshape(-1, 3), r.gauss_stresses.reshape(-1, 3),
    r.gauss_von_mises.ravel(),
])
path_gp = os.path.join(OUTPUT, "05_gauss_results.csv")
np.savetxt(path_gp, gp, delimiter=",", comments="", fmt="%.10g",
           header="element,x,y,exx,eyy,gxy,sxx,syy,txy,von_mises")

print("diagnostics:", {k: r.diagnostics[k] for k in ("solver", "relative_residual", "n_free")})
print(f"CSV files written: {path}, {path_gp}")
