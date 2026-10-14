[← Documentation home](README.md) · [Main README](../README.md)

# API reference

## Full API reference

All public names are importable from `engine`:

```python
from engine import (Mesh, Material, FEAModel, FEAResults, Q4Element,
                    von_mises, cantilever_tip_deflection,
                    FEAError, MeshTopologyError, DistortedElementError,
                    SingularSystemError, SolverConvergenceError, ElementIntegrityError)
```

### `Mesh(nodes, elements, validate=True)`
| Member | Description |
|---|---|
| `nodes`, `elements` | Coordinate array (n, 2) float; connectivity (m, 4) int |
| `n_nodes`, `n_elements`, `n_dofs` | Sizes (`n_dofs = 2 · n_nodes`) |
| `element_dofs` | (m, 8) global DOF numbers per element, ordered `[u1, v1, …, u4, v4]` |
| `element_coords` | (m, 4, 2) coordinates per element |
| `node_dofs(n)` *(static)* | `(2n, 2n+1)` |
| `dof_index(n, "x"/"y")` *(static)* | Single DOF number |
| `signed_areas()` | Element areas (positive = counter-clockwise) |
| `corner_cross_products()` | Corner convexity values (all positive = convex) |
| `bounding_box()` | `(xmin, xmax, ymin, ymax)` |
| `characteristic_length()` | Bounding-box diagonal |
| `validate()` | Re-run topology checks |
| `nodes_where(rule)`, `nodes_at(x, y, tol)`, `closest_node(x, y)` | Node selection ([user guide](user-guide.md#picking-nodes-for-supports-and-loads)) |
| `structured_rectangle(length, height, nx, ny, origin=(0,0))` *(classmethod)* | Rectangular grid |
| `plate_with_hole_quarter(radius, half_width, n_radial, n_circumferential, grading=1.0)` *(classmethod)* | Curved graded mesh |

### `Material(E, nu, thickness=1.0, plane="stress", yield_strength=None, name="material")`
| Member | Description |
|---|---|
| `E`, `nu`, `thickness`, `yield_strength`, `name` | Inputs |
| `plane` | `"stress"` or `"strain"`; can be reassigned |
| `D` | 3×3 constitutive matrix |
| `shear_modulus` | G = E / 2(1+ν) |
| `out_of_plane_stress(stress)` | σzz (ν(σxx+σyy) in plane strain, 0 in plane stress) |
| `PLANE_STRESS`, `PLANE_STRAIN` | Constants `"stress"`, `"strain"` |

### `FEAModel(mesh, material, iterative_threshold=50000, cg_rtol=1e-10, cg_maxiter=None, direct_fallback=True)`
| Member | Description |
|---|---|
| `fix(nodes, dofs="xy", value=0.0)` | Supports / prescribed displacement |
| `add_nodal_force(node, fx=0, fy=0)` | Point load (accumulates) |
| `add_edge_traction(edge_nodes, traction)` | Edge traction: constant `(tx, ty)` or `f(x, y) -> (tx, ty)`; returns resultant magnitude |
| `clear_loads()`, `clear_constraints()` | Reset |
| `assemble()` | Build and return global CSR stiffness `K` |
| `solve()` | Full analysis → `FEAResults` |
| `F` | Global force vector |
| `K`, `Ke`, `detJ`, `results` | Assembled matrix, element matrices, Gauss-point Jacobians, last results |
| `constrained_dofs` | Sorted fixed DOF numbers |
| `ITERATIVE_THRESHOLD` | Class default 50,000 |

### `FEAResults`
Fields are listed in [Results](results-and-plots.md#reading-and-exporting-results). Also has `max_displacement` and `max_von_mises`.

`diagnostics` keys: `n_nodes`, `n_elements`, `n_dofs`, `n_free`, `n_constrained`, `nnz`, `solver`, `iterations` (−1 for direct), `relative_residual`, `equilibrium_fx`, `equilibrium_fy`, `equilibrium_mz`, `detJ_min`, `detJ_max`, `solve_time`, `total_time`.

### `Q4Element` (low-level, for study or extension)
| Member | Description |
|---|---|
| `NODE_NATURAL` | Corner natural coordinates (−1,−1), (1,−1), (1,1), (−1,1) |
| `GAUSS_POINTS`, `GAUSS_WEIGHTS`, `GAUSS_ABSCISSA` | 2×2 rule: ±1/√3, weights 1 |
| `shape_functions(xi, eta)` | N₁…N₄ |
| `shape_derivatives(xi, eta)` | 2×4 derivatives w.r.t. ξ, η |
| `jacobian(coords, xi, eta)` | (m, 2, 2) Jacobian |
| `kinematics(coords, xi, eta)` | `(detJ, B)` with the distortion check |
| `stiffness(coords, D, thickness, verify=True)` | `(Ke, detJ)` for a batch of elements; a single element `(4, 2)` works too |
| `verify_stiffness(Ke)` | Symmetry and positive-semi-definiteness check |
| `extrapolation_matrix()` | 4×4 Gauss → corner stress matrix |

### Helpers
- `von_mises(stress, szz=None)`: Von Mises from `[..., (σxx, σyy, τxy)]` arrays.
- `cantilever_tip_deflection(P, L, E, I)`: textbook PL³/3EI.

### Exceptions
All inherit from `FEAError` (a `RuntimeError`), so `except FEAError:` catches every solver problem.
`MeshTopologyError`, `DistortedElementError`, `SingularSystemError`, `SolverConvergenceError`, `ElementIntegrityError`.
Bad numeric inputs to `Material` / `Mesh` factories raise `ValueError`; bad node numbers raise `IndexError`.

### Visualizer functions
`plot_results(...)`, `plot_stress_components(...)` ([plots](results-and-plots.md#pictures-visualizerpy)) plus helpers:
- `auto_scale_factor(mesh, results, target_fraction=0.10)`
- `quad_triangulation(mesh, coords)`: splits quads into triangles for Matplotlib contouring.
