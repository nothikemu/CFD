[← Documentation home](README.md) · [Main README](../README.md)

# Results, export and plots

## Reading and exporting results

`r = model.solve()` gives you:

| Attribute | Shape | Meaning |
|---|---|---|
| `r.displacements` | (nodes, 2) | Movement of each node `[ux, uy]` |
| `r.u` | (2·nodes,) | Same as one long vector `[ux0, uy0, ux1, uy1, ...]` |
| `r.reactions` | (2·nodes,) | Support forces (non-zero only where you used `fix`) |
| `r.external_forces` | (2·nodes,) | The loads you applied |
| `r.nodal_stresses` | (nodes, 3) | Smoothed stress at each node `[σxx, σyy, τxy]` |
| `r.nodal_szz` | (nodes,) | Out-of-plane stress σzz (zero for plane stress) |
| `r.von_mises` | (nodes,) | Von Mises equivalent stress at each node |
| `r.safety_factor` | (nodes,) or `None` | yield_strength / Von Mises (needs `yield_strength`) |
| `r.gauss_points` | (elements, 4, 2) | x, y of the 4 internal sample points per element |
| `r.gauss_strains` | (elements, 4, 3) | Strains `[εxx, εyy, γxy]` at those points |
| `r.gauss_stresses` | (elements, 4, 3) | Stresses `[σxx, σyy, τxy]` at those points (the most accurate stresses) |
| `r.gauss_von_mises` | (elements, 4) | Von Mises at those points |
| `r.element_nodal_stresses` | (elements, 4, 3) | Stress at each element's corners before averaging |
| `r.max_displacement` | number | Largest movement anywhere |
| `r.max_von_mises` | number | Largest nodal Von Mises stress |
| `r.diagnostics` | dict | Solver name, iterations, residual, equilibrium errors, det(J) range, timings, sizes |

Everyday recipes:

```python
node = mesh.closest_node(1.0, 0.05)
ux, uy = r.displacements[node]                      # movement of one point
worst = int(np.argmax(r.von_mises))                 # most stressed node
print(mesh.nodes[worst], r.von_mises[worst])
left = mesh.nodes_at(x=0.0)
print("support force y:", r.reactions[2 * left + 1].sum())   # x-reactions: r.reactions[2 * left]
print("yields?", r.safety_factor.min() < 1)
print(r.diagnostics)
```

**Export to CSV** (opens in Excel / Sheets):

```python
import numpy as np
table = np.column_stack([np.arange(mesh.n_nodes), mesh.nodes, r.displacements,
                         r.nodal_stresses, r.von_mises, r.reactions.reshape(-1, 2)])
np.savetxt("output/results.csv", table, delimiter=",", comments="", fmt="%.10g",
           header="node,x,y,ux,uy,sxx,syy,txy,von_mises,Rx,Ry")
```
A full version (with Gauss-point output too) is in `examples/05_prescribed_displacement_and_export.py`.

## Pictures (`visualizer.py`)

### Main picture: mesh + Von Mises on the deformed shape

```python
from visualizer import plot_results
fig, used_scale = plot_results(model, r, filename="output/my_result.png")
```

| Argument | Default | Meaning |
|---|---|---|
| `model` | (required) | The solved `FEAModel` |
| `results` | `model.results` | Results to draw |
| `scale` | automatic | How many times to exaggerate the movement. Automatic makes the largest movement 10 % of the part size. `scale=1` shows true size. |
| `filename` | `"fea_results.png"` | Where to save (folders are created). `None` = don't save. |
| `title` | `"Q4 Finite Element Analysis"` | Heading text |
| `stress_unit` / `length_unit` | `"Pa"` / `"m"` | Axis/colour-bar labels only (no conversion) |
| `n_levels` | 32 | Number of colour bands |
| `cmap` | `"turbo"` | Any Matplotlib colour map (`"viridis"`, `"jet"`, `"plasma"`, …) |
| `dpi` | 220 | Image resolution |
| `figsize` | automatic | `(width, height)` in inches |
| `show` | `False` | `True` opens a window (desktop) or shows inline (notebooks) |

Left (or top) panel: the undeformed mesh. Red squares are clamped nodes, orange/green triangles are rollers, dark arrows are loads.
Right (or bottom) panel: colour = Von Mises stress, drawn on the exaggerated bent shape. The dashed line is the original outline and the white circle marks the maximum.
Long, slender parts are stacked vertically, and squarish parts are drawn side by side.

### Separate stress components

```python
from visualizer import plot_stress_components
plot_stress_components(model, r, filename="output/components.png", stress_unit="Pa", cmap="RdBu_r", dpi=200)
```
Draws σxx, σyy and τxy with a red (tension) / blue (compression) colour scale.

### Showing on screen and customising

- In a notebook, use `plot_results(..., show=True)` or call `plt.show()`.
- On a desktop, `show=True` opens a window.
- On a Linux machine **without a screen** (server, container), pictures are only saved to files. That happens automatically.
- Both functions return Matplotlib figures, so you can add things before saving: `fig.axes[1].set_xlim(...)`, `fig.savefig("x.pdf")`.
