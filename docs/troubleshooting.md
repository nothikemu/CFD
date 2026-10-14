[← Documentation home](README.md) · [Main README](../README.md)

# Troubleshooting

## Troubleshooting (every error message explained)

| You see | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError: No module named 'numpy'` (or scipy/matplotlib) | Packages not installed in this Python | `python3 -m pip install -r requirements.txt` (activate your `.venv` first if you made one) |
| `ModuleNotFoundError: No module named 'engine'` | Your script isn't in the `CFD` folder | Move it next to `app.py`, or add `sys.path.insert(0, "/path/to/CFD")` before the imports |
| `python3: command not found` (Windows) | Windows calls it `python` | Use `python` |
| `MeshTopologyError: ... clockwise or degenerate ...` | Element nodes listed clockwise, or element has zero area | Reverse the order of that element's 4 nodes, e.g. `[a, b, c, d]` → `[a, d, c, b]` |
| `MeshTopologyError: ... outside [0, N]` | Connectivity refers to a node number that doesn't exist | Node numbers start at 0. Check your element list. |
| `MeshTopologyError: ... same node more than once` | A node is repeated inside an element (e.g. faking a triangle) | Use 4 distinct nodes; split triangles into quads |
| `MeshTopologyError: nodes must have shape (n_nodes, 2)` / `elements must have shape (n_elements, 4)` | Input arrays are the wrong shape | Nodes: rows of `[x, y]`. Elements: rows of 4 integers. |
| `WARNING ... non-convex element(s)` | An element has an inside angle ≥ 180° | Reshape those elements. If it's severe, the next error appears. |
| `DistortedElementError: Non-positive or near-zero Jacobian ...` | Element is inverted, collapsed or too distorted | Fix the listed element ids: make them convex and counter-clockwise |
| `SingularSystemError: no Dirichlet constraints defined` | No supports at all | Add `model.fix(...)` |
| `SingularSystemError: ... unrestrained rigid-body mode(s)` | Part can still slide or spin | Add supports in the missing direction. Check that rollers aren't all parallel and don't all pass through one point. |
| `SingularSystemError: ... mesh component N ...` | Your mesh has **separate, unconnected pieces**, and one isn't supported | Make the pieces share nodes, or support each piece |
| `SolverConvergenceError` | Iterative solver didn't converge with `direct_fallback=False` | Use `direct_fallback=True`, raise `cg_maxiter`, or loosen `cg_rtol` |
| `WARNING small LU pivot ratio` | Extremely slender part or huge stiffness contrast; results are probably fine but less precise | Check equilibrium in the log; refine or rescale units |
| `ValueError: Poisson's ratio must satisfy -1 < nu < 0.5` | `nu` of 0.5 or more | Use 0.499 at most (rubber-like materials aren't well suited to Q4) |
| `IndexError: node N out of range` | Node number in `fix`/`add_nodal_force` doesn't exist | Use `mesh.closest_node` / `mesh.nodes_at` to find node numbers |
| `ValueError: zero-length traction segment` | Same node twice in an edge chain | Remove duplicates; sort the edge nodes |
| Edge load looks wrong / too big | Edge nodes not in order along the edge | Sort them: `edge[np.argsort(mesh.nodes[edge, 1])]` (or by x for horizontal edges) |
| Picture looks flat / no visible bending | Movement is tiny compared with the part | Pass `scale=` (e.g. 100) to `plot_results` or `--scale` to `app.py` |
| Picture is wildly distorted | Scale too large | Lower `scale` |
| No window appears on Linux server | No display available | Expected: pictures are saved in `output/` |
| Answers off by 1000× | Units mixed (mm with Pa, etc.) | See [Units](user-guide.md#units-pick-one-consistent-set) |
| Deflection a bit smaller than hand calc | Too few elements through the thickness (locking) | Refine: [Accuracy](accuracy-and-validation.md#getting-accurate-answers) |
