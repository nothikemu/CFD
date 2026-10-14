[← Documentation home](README.md) · [Main README](../README.md)

# User guide: building a model

## Building your own model, step by step

Every model is the same five steps: **shape → material → supports → loads → solve**.

```python
from engine import Mesh, Material, FEAModel
from visualizer import plot_results

mesh  = Mesh.structured_rectangle(length=1.0, height=0.1, nx=100, ny=10)   # 1. shape
steel = Material(E=210e9, nu=0.3, thickness=0.02, yield_strength=250e6)    # 2. material
model = FEAModel(mesh, steel)
model.fix(mesh.nodes_at(x=0.0), "xy")                                      # 3. supports
model.add_nodal_force(mesh.closest_node(1.0, 0.05), fy=-2000)              # 4. loads
r = model.solve()                                                          # 5. solve
plot_results(model, r, filename="output/my_beam.png")
print(r.max_von_mises, r.safety_factor.min())
```

> **Where to put your script:** in the `CFD` folder (next to `app.py`), so `from engine import ...` works. If you put it in `examples/`, add `from _setup import OUTPUT` as the first import, as the examples do. To run it from anywhere else, first do `import sys; sys.path.insert(0, "/path/to/CFD")`.

### Units: pick one consistent set

The solver has no built-in units. It just needs them to match. Two common sets:

| Quantity | SI set | Millimetre set |
|---|---|---|
| Length | m | mm |
| Force | N | N |
| Stress, `E`, pressure | Pa (N/m²) | MPa (N/mm²) |
| Steel `E` | `210e9` | `210000` |
| Results come out in | m, Pa | mm, MPa |

Mixing sets (for example mm for lengths with Pa for `E`) gives answers that are wrong by factors of 1000 or more.

Typical material values:

| Material | `E` (Pa) | `nu` | Yield (Pa) |
|---|---|---|---|
| Structural steel | 200e9 to 210e9 | 0.30 | 250e6 to 355e6 |
| Stainless steel | 193e9 | 0.29 | 215e6 |
| Aluminium 6061-T6 | 69e9 | 0.33 | 275e6 |
| Titanium Ti-6Al-4V | 114e9 | 0.34 | 880e6 |
| Concrete | 25e9 to 40e9 | 0.20 | (brittle; use with care) |
| Wood (along grain, approximate) | 10e9 to 14e9 | 0.30 | (anisotropic; approximate only) |

### Shape (mesh)

A mesh is a set of **nodes** (points) joined into **elements** (four-sided cells). There are three ways to make one.

**a) Rectangle (most common)**
```python
mesh = Mesh.structured_rectangle(length, height, nx, ny, origin=(x0, y0))
```
- `nx`, `ny`: number of elements along x and y.
- `origin`: bottom-left corner, default `(0, 0)`. Use `origin=(0, -height/2)` to put the beam's centre line on y = 0.

**b) Quarter plate with a circular hole**
```python
mesh = Mesh.plate_with_hole_quarter(radius, half_width, n_radial, n_circumferential, grading=1.0)
```
- Covers the square `[0, half_width] × [0, half_width]` minus a hole of `radius` at the origin.
- `n_circumferential` must be even.
- `grading > 1` packs small elements near the hole, where stress changes fastest (8 is a good value).

**c) Any shape: supply the points yourself**
```python
nodes = [[x0, y0], [x1, y1], ...]          # one row per node; node numbers start at 0
elements = [[a, b, c, d], ...]             # 4 node numbers per element
mesh = Mesh(nodes, elements)
```

Rules for your own mesh (they are checked automatically):

1. **List each element's 4 nodes counter-clockwise** (anticlockwise) as you walk around it, for example bottom-left → bottom-right → top-right → top-left. Clockwise elements are rejected with `MeshTopologyError`.
2. Elements must be **convex**: every inside angle under 180°. Badly bent elements raise `DistortedElementError`.
3. Neighbouring elements must **share** node numbers along their common edge. Don't duplicate points, or the parts won't be connected.
4. Node numbers must exist (0 … number of nodes − 1), and a node can't appear twice in one element.
5. Only quadrilaterals. No triangles; a triangle can't be faked by repeating a node.

See `examples/04_custom_mesh.py` for a complete hand-built mesh. Meshes from other tools (Gmsh, etc.) work too: load their quad node coordinates and connectivity into two arrays and pass them to `Mesh(...)`.

Useful mesh properties:

```python
mesh.n_nodes, mesh.n_elements, mesh.n_dofs   # sizes
mesh.nodes                                   # (n_nodes, 2) coordinates
mesh.elements                                # (n_elements, 4) connectivity
mesh.bounding_box()                          # (xmin, xmax, ymin, ymax)
print(mesh)                                  # Mesh(n_nodes=..., n_elements=..., n_dofs=...)
```

### Picking nodes (for supports and loads)

| Call | Returns | Example |
|---|---|---|
| `mesh.nodes_at(x=value)` | All nodes on the vertical line x = value | left edge: `mesh.nodes_at(x=0.0)` |
| `mesh.nodes_at(y=value)` | All nodes on the horizontal line y = value | bottom edge: `mesh.nodes_at(y=0.0)` |
| `mesh.nodes_at(x=a, y=b)` | The node(s) exactly at (a, b) | corner |
| `mesh.nodes_at(x=..., tol=1e-6)` | Same, with your own matching tolerance | |
| `mesh.closest_node(x, y)` | The single node nearest (x, y) | load point: `mesh.closest_node(1.0, 0.05)` |
| `mesh.nodes_where(rule)` | All nodes where `rule(x, y)` is True (works on arrays) | `mesh.nodes_where(lambda x, y: (x > 0.5) & (y < 0.01))` |

To get the nodes along an edge **in order** (needed for edge loads), sort them:
```python
edge = mesh.nodes_at(x=1.0)
edge = edge[np.argsort(mesh.nodes[edge, 1])]   # sort by y
```

### Material

```python
Material(E, nu, thickness=1.0, plane="stress", yield_strength=None, name="material")
```

| Argument | Meaning | Allowed |
|---|---|---|
| `E` | Young's modulus (stiffness) | > 0 |
| `nu` | Poisson's ratio | between −1 and 0.5 (not 0.5 itself) |
| `thickness` | Out-of-plane thickness. Forces and areas are multiplied by it. For plane strain, use 1.0 to get results "per unit depth". | > 0 |
| `plane` | `"stress"` for thin flat parts (sheet metal, brackets, beams); `"strain"` for long or thick bodies that can't stretch in the third direction (dams, tunnels, retaining walls, long pipes) | `"stress"`, `"strain"` (also `"plane_stress"`, `"plane strain"`, …) |
| `yield_strength` | If given, results include a **safety factor** = yield / Von Mises | > 0 or `None` |
| `name` | Label used in logs | any text |

You can switch the assumption later with `steel.plane = "strain"`. If the model was already assembled, build a new `FEAModel` afterwards.

### Supports (Dirichlet boundary conditions)

```python
model.fix(nodes, dofs="xy", value=0.0)
```

| Want | Call |
|---|---|
| Clamp / weld (no movement) | `model.fix(nodes, "xy")` |
| Roller that can slide **horizontally** (holds y only) | `model.fix(nodes, "y")` |
| Roller that can slide **vertically** (holds x only) | `model.fix(nodes, "x")` |
| Pin (a single point held in x and y) | `model.fix(node, "xy")` |
| **Symmetry line** at x = const | `model.fix(mesh.nodes_at(x=...), "x")` |
| **Symmetry line** at y = const | `model.fix(mesh.nodes_at(y=...), "y")` |
| **Move** nodes by a set amount | `model.fix(nodes, "x", value=0.001)` or `model.fix(node, "xy", value=(dx, dy))` |

`nodes` can be one node number or a list/array. Calling `fix` again on the same DOF overwrites it. `model.clear_constraints()` removes all supports.

**The golden rule:** the supports must stop the part from **sliding in x, sliding in y and spinning**. The minimum is 3 restrained directions that aren't all parallel and don't all pass through one point, for example a pin (x+y) at one point plus a y-roller at another. If you forget, the solver stops with a `SingularSystemError` that says which motion is still free.

### Loads (Neumann boundary conditions)

**Point force on a node**
```python
model.add_nodal_force(node, fx=0.0, fy=0.0)      # units: force (N)
```
Calls add up, so you can call it repeatedly on the same node.

**Pressure / traction along an edge** (force per area, Pa)
```python
model.add_edge_traction(edge_nodes, (tx, ty))
```
- `edge_nodes` is the **ordered** chain of nodes along a boundary (see [Picking nodes](#picking-nodes-for-supports-and-loads)), and each neighbouring pair must be an element edge.
- `(tx, ty)` is the force per unit area in x and y. Downward pressure on a top face is `(0, -p)`, and tension pulling a right face is `(+s, 0)`.
- The edge length × `thickness` is accounted for automatically, and the method returns the size of the total force added.

**Varying pressure** (e.g. water pressure growing with depth): pass a function of position:
```python
rho_g, surface = 9810.0, 5.0
model.add_edge_traction(wall_nodes, lambda x, y: (rho_g * (surface - y), 0.0))
```
The function receives `x, y` and returns `(tx, ty)`. It's integrated exactly for anything up to a cubic variation along each element edge.

`model.clear_loads()` removes all loads. Loads are applied to the model's force vector `model.F` (2 entries per node: x, y).

**Gravity / self-weight** has no dedicated call. Approximate it by spreading each element's weight onto its four nodes:
```python
rho, g = 7850.0, 9.81
area = mesh.signed_areas()
for e, a in zip(mesh.elements, area):
    for n in e:
        model.add_nodal_force(n, fy=-rho * g * a * steel.thickness / 4)
```

### Solve

```python
r = model.solve()
```
This assembles, checks, solves and post-processes in one call, and returns an `FEAResults` object ([Results](results-and-plots.md)). To re-solve after changing loads or supports, just call `model.solve()` again. The stiffness matrix is reused automatically because it only depends on mesh + material.
