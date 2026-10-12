# Q4 FEA Structural Solver

A 2D finite element program written from scratch in Python. It works out how flat parts such as beams, brackets and plates **bend and stretch** under load, and **how much stress** they carry. It also tells you whether they'll yield. No commercial or black-box solver is used: every step of the method is plain, readable Python (NumPy, SciPy sparse matrices, Matplotlib).

![example output](https://img.shields.io/badge/validated-Euler--Bernoulli%20%3C2%25-brightgreen) ![tests](https://img.shields.io/badge/tests-24%20passing-brightgreen) ![python](https://img.shields.io/badge/python-3.9%2B-blue)

---

## Table of contents

1. [What it can do](#1-what-it-can-do)
2. [Fastest start: no terminal (Google Colab)](#2-fastest-start-no-terminal-google-colab)
3. [Installing on your own computer](#3-installing-on-your-own-computer)
4. [Running the built-in benchmarks (`app.py`)](#4-running-the-built-in-benchmarks-apppy)
5. [Ready-made examples (`examples/`)](#5-ready-made-examples-examples)
6. [Building your own model, step by step](#6-building-your-own-model-step-by-step)
7. [Reading and exporting results](#7-reading-and-exporting-results)
8. [Pictures (`visualizer.py`)](#8-pictures-visualizerpy)
9. [Configuration: solver, accuracy, logging](#9-configuration-solver-accuracy-logging)
10. [Getting accurate answers](#10-getting-accurate-answers)
11. [Model size and speed](#11-model-size-and-speed)
12. [Validation tests](#12-validation-tests)
13. [Troubleshooting (every error message explained)](#13-troubleshooting-every-error-message-explained)
14. [Full API reference](#14-full-api-reference)
15. [The maths inside](#15-the-maths-inside)
16. [Limitations](#16-limitations)
17. [Project layout](#17-project-layout)

---

## 1. What it can do

| Feature | Details |
|---|---|
| Element | 4-node isoparametric quadrilateral (Q4), full 2×2 Gauss integration |
| Meshes | Built-in rectangle and plate-with-hole generators, or **any mesh you supply** (coordinates + connectivity) |
| Material | Isotropic linear elastic: Young's modulus `E`, Poisson's ratio `nu`, thickness, optional yield strength |
| 2D assumption | **Plane stress** (thin plates, sheet parts) or **plane strain** (thick/long bodies: dams, tunnels, retaining walls) |
| Supports | Fix x, y or both on any nodes: clamps, pins, rollers, symmetry planes. Prescribed (non-zero) displacements too |
| Loads | Point forces on nodes; distributed pressure/traction on edges, constant or varying with position |
| Solver | Sparse direct solver (SuperLU). Above 50,000 unknowns it switches automatically to preconditioned conjugate gradient, with a direct fallback |
| Results | Displacements, reaction forces, strains, stresses (σxx, σyy, τxy, σzz), Von Mises stress, safety factor |
| Checks | Counter-clockwise numbering, distorted elements (Jacobian), missing supports, element symmetry/positivity, global equilibrium |
| Output | Step-by-step maths log, high-resolution PNG contour plots, CSV export |

---

## 2. Fastest start: no terminal (Google Colab)

Colab runs Python in your web browser for free. You don't need to install anything.

1. Go to <https://colab.research.google.com>.
2. Click **File → Open notebook → GitHub**, paste `https://github.com/nothikemu/CFD`, and pick
   `examples/quickstart.ipynb` on branch `claude/sharp-carson-1l4xvs`.
   *(Alternatively, click **New notebook** and paste the cells from the box below.)*
3. Press **▶** on each grey box from top to bottom (or press **Shift + Enter**).
4. In the **SETTINGS** box, change the length, height, material or force, then press **▶** again.

> **Private repository?** Colab will ask you to sign in to GitHub. If the repo is public you won't be asked.

Manual Colab cells:

```python
# Cell 1: download + install
!git clone -b claude/sharp-carson-1l4xvs https://github.com/nothikemu/CFD.git
%cd CFD
!pip -q install -r requirements.txt
```

```python
# Cell 2: check everything works
!python -m pytest test_validation.py -q
```

```python
# Cell 3: run the benchmarks and show the pictures
!python app.py
from IPython.display import Image, display
display(Image("output/cantilever_von_mises.png"))
```

---

## 3. Installing on your own computer

### 3.1 You need Python 3.9 or newer

Check by opening a terminal (**Mac:** Terminal app · **Windows:** PowerShell · **Linux:** any terminal):

```bash
python3 --version        # Mac / Linux
python --version         # Windows
```

If it's missing or older than 3.9, install it from <https://www.python.org/downloads/>. On Windows, tick **"Add python.exe to PATH"** in the installer.

> Throughout this guide, **Windows users type `python` wherever it says `python3`**.

### 3.2 Download the code

**With git:**
```bash
git clone -b claude/sharp-carson-1l4xvs https://github.com/nothikemu/CFD.git
cd CFD
```

**Without git:** on the GitHub page, switch the branch dropdown to `claude/sharp-carson-1l4xvs` and click **Code → Download ZIP**. Unzip it, then open a terminal inside the unzipped folder.

### 3.3 (Recommended) make a private Python environment

This keeps the solver's packages separate from everything else on your computer.

```bash
# Mac / Linux
python3 -m venv .venv
source .venv/bin/activate
```
```powershell
# Windows PowerShell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

If Windows refuses with "running scripts is disabled", run this once, then activate again:
`Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`

Run the activate line again every time you open a new terminal. To leave the environment, type `deactivate`.

### 3.4 Install the packages

```bash
python3 -m pip install -r requirements.txt
```

This installs `numpy`, `scipy`, `matplotlib` and `pytest`. Nothing else is needed.

### 3.5 Check it works

```bash
python3 -m pytest test_validation.py -v
```

You should see **`24 passed`**.

### 3.6 Using an editor instead of a terminal (VS Code)

1. Install [VS Code](https://code.visualstudio.com/) and its **Python** extension.
2. **File → Open Folder…** and choose the `CFD` folder.
3. Pick your interpreter: press `Ctrl/Cmd + Shift + P`, run **Python: Select Interpreter**, and choose `.venv` if you made one.
4. Open any `.py` file (for example `examples/01_my_beam.py`) and click **▶ Run Python File** in the top-right corner.
5. Pictures appear in the `output/` folder in the sidebar. Click one to view it.

PyCharm, Spyder and Thonny work the same way: open the folder, open a file, press Run.

---

## 4. Running the built-in benchmarks (`app.py`)

```bash
python3 app.py
```

This runs two problems and prints every step of the calculation:

1. **Cantilever beam:** a 2 m steel beam clamped at one end with 5 kN down at the tip. The result is compared with the textbook formula δ = PL³/3EI and passes if within 2 %.
2. **Plate with a hole:** an aluminium plate in tension, compared with the textbook stress-concentration factor (≈3).

Pictures go to `output/`:

| File | Shows |
|---|---|
| `output/cantilever_von_mises.png` | Undeformed mesh + Von Mises stress on the bent shape |
| `output/cantilever_components.png` | σxx, σyy and τxy separately |
| `output/plate_hole_von_mises.png` | Plate-with-hole mesh and stress |

### All `app.py` options

| Option | What it does | Example |
|---|---|---|
| `--nx N` | Elements along the beam (default 160) | `python3 app.py --nx 320` |
| `--ny N` | Elements through the beam height (default 8) | `python3 app.py --ny 16` |
| `--plane stress\|strain` | 2D assumption for both benchmarks (default `stress`) | `python3 app.py --plane strain` |
| `--scale X` | Exaggerate the bend X times in pictures (default: automatic) | `python3 app.py --scale 50` |
| `--no-plot` | Skip pictures (faster) | `python3 app.py --no-plot` |
| `--run-tests` | Run the validation suite afterwards | `python3 app.py --run-tests` |
| `--quiet` | Print only warnings/errors | `python3 app.py --quiet` |
| `--help` | List options | `python3 app.py --help` |

Options can be combined, for example `python3 app.py --nx 320 --ny 16 --scale 30 --run-tests`.
The program exits with code `0` if the cantilever is within 2 % (and, with `--run-tests`, all tests pass), and `1` otherwise. That makes it usable in CI.

### Understanding the log

Every solve prints numbered steps:

| Step | Meaning | What "healthy" looks like |
|---|---|---|
| `[STEP 1]` | Mesh size: nodes, elements, degrees of freedom (2 per node) | |
| `[STEP 2]` | The 3×3 material stiffness matrix **D** | |
| `[STEP 3]` | Jacobian determinants at all Gauss points; domain area two ways; eigenvalues of the first element matrix | all det(J) positive; both areas equal; exactly 3 eigenvalues ≈ 0 (the rigid-body motions) |
| `[STEP 4]` | Global sparse matrix: size, non-zeros, memory vs dense, symmetry error | symmetry error tiny compared with the matrix values |
| `[STEP 5]` | Free vs fixed unknowns; rigid-body restraint check | "verified" |
| `[STEP 6]` | Which solver ran, time, residual | relative residual ≲ 1e-8 |
| `[STEP 7]` | Reaction forces and global equilibrium (ΣFx, ΣFy, ΣM) | relative imbalance ≲ 1e-8 |
| `[STEP 8]` | Strain and stress ranges at Gauss points | |
| `[STEP 9]` | Nodal Von Mises; safety factor if a yield strength was given | safety factor > 1 means no yielding |

---

## 5. Ready-made examples (`examples/`)

Each file runs on its own, from any folder or with an editor's Run button. Pictures and CSVs go to `output/`.

| File | Teaches | Run |
|---|---|---|
| `01_my_beam.py` | The simplest job: clamp, push, solve, plot. Edit the **SETTINGS** block. | `python3 examples/01_my_beam.py` |
| `02_simply_supported_pressure.py` | Pin + roller supports, **pressure on an edge**, reactions, comparison with 5qL⁴/384EI | `python3 examples/02_simply_supported_pressure.py` |
| `03_plate_with_hole.py` | **Symmetry** (model ¼ of a part), curved graded mesh, stress concentration, safety factor | `python3 examples/03_plate_with_hole.py` |
| `04_custom_mesh.py` | **Your own mesh** from raw coordinates + connectivity (tapered bracket), **plane strain**, custom node selection, angled load | `python3 examples/04_custom_mesh.py` |
| `05_prescribed_displacement_and_export.py` | **Move a boundary by a set amount**, read the force needed, **export results to CSV** | `python3 examples/05_prescribed_displacement_and_export.py` |
| `06_large_model_and_logging.py` | 52,000-unknown model, **iterative solver**, solver settings, full logging | `python3 examples/06_large_model_and_logging.py` |
| `quickstart.ipynb` | Notebook version for Colab / Jupyter | open in Colab or Jupyter |

The quickest way to start your own job is to copy `01_my_beam.py`, rename it and edit it.

---

## 6. Building your own model, step by step

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

### 6.1 Units: pick one consistent set

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

### 6.2 Shape (mesh)

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

### 6.3 Picking nodes (for supports and loads)

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

### 6.4 Material

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

### 6.5 Supports (Dirichlet boundary conditions)

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

### 6.6 Loads (Neumann boundary conditions)

**Point force on a node**
```python
model.add_nodal_force(node, fx=0.0, fy=0.0)      # units: force (N)
```
Calls add up, so you can call it repeatedly on the same node.

**Pressure / traction along an edge** (force per area, Pa)
```python
model.add_edge_traction(edge_nodes, (tx, ty))
```
- `edge_nodes` is the **ordered** chain of nodes along a boundary (see §6.3), and each neighbouring pair must be an element edge.
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

### 6.7 Solve

```python
r = model.solve()
```
This assembles, checks, solves and post-processes in one call, and returns an `FEAResults` object (§7). To re-solve after changing loads or supports, just call `model.solve()` again. The stiffness matrix is reused automatically because it only depends on mesh + material.

---

## 7. Reading and exporting results

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

---

## 8. Pictures (`visualizer.py`)

### 8.1 Main picture: mesh + Von Mises on the deformed shape

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

### 8.2 Separate stress components

```python
from visualizer import plot_stress_components
plot_stress_components(model, r, filename="output/components.png", stress_unit="Pa", cmap="RdBu_r", dpi=200)
```
Draws σxx, σyy and τxy with a red (tension) / blue (compression) colour scale.

### 8.3 Showing on screen and customising

- In a notebook, use `plot_results(..., show=True)` or call `plt.show()`.
- On a desktop, `show=True` opens a window.
- On a Linux machine **without a screen** (server, container), pictures are only saved to files. That happens automatically.
- Both functions return Matplotlib figures, so you can add things before saving: `fig.axes[1].set_xlim(...)`, `fig.savefig("x.pdf")`.

---

## 9. Configuration: solver, accuracy, logging

### 9.1 Solver settings

```python
model = FEAModel(mesh, material,
                 iterative_threshold=50_000,  # above this many free unknowns: use iterative solver
                 cg_rtol=1e-10,               # iterative solver stopping tolerance (relative)
                 cg_maxiter=None,             # max iterations (None = automatic)
                 direct_fallback=True)        # if iterative solver stalls: finish with direct solver
```

| Setting | When to change it |
|---|---|
| `iterative_threshold` | Lower it (e.g. `10`) to force the iterative solver on small models. Raise it (e.g. `10**9`) to always use the direct solver, which is often fastest in 2D if you have enough memory. |
| `cg_rtol` | Smaller = more accurate, more iterations. `1e-8` to `1e-12` is sensible. |
| `cg_maxiter` | Cap run time on huge models. |
| `direct_fallback` | Set `False` if you'd rather get a `SolverConvergenceError` than wait for a direct solve. |

**How the solver chooses:**
- Up to `iterative_threshold` free unknowns, it uses the SuperLU sparse direct factorisation.
- Above that, it uses conjugate gradient. The matrix is scaled symmetrically by its diagonal first, then an incomplete-LU preconditioner is built with a symmetric fill-reducing ordering. If CG doesn't converge, it falls back to the direct solver.
- `r.diagnostics["solver"]` tells you which one ran.

### 9.2 Logging (the step-by-step maths output)

The engine reports through Python's standard `logging` under the name `engine.fea`. `app.py` turns it on automatically. In your own scripts:

```python
import logging, sys
logging.basicConfig(level=logging.INFO, stream=sys.stdout, format="%(message)s")
logging.getLogger("matplotlib").setLevel(logging.WARNING)   # hide Matplotlib chatter
```
- **Silence everything except problems:** `logging.getLogger("engine.fea").setLevel(logging.WARNING)`.
- **Save the log to a file:** add `filename="run.log"` to `basicConfig` (and remove `stream=`).

### 9.3 Advanced knobs

- `Mesh(nodes, elements, validate=False)` skips the up-front mesh checks. The Jacobian checks during the solve still run.
- `FEAModel.assemble()` builds the global stiffness matrix without solving: `K = model.assemble()` (SciPy CSR). Useful for your own analysis, e.g. eigenvalues.
- `model.K`, `model.Ke` (all element matrices, shape (elements, 8, 8)), `model.detJ` (Jacobian determinants, shape (elements, 4)) are available after assembly.
- `model.constrained_dofs` lists every fixed DOF index (node `n` → x is `2n`, y is `2n+1`).
- `Mesh.node_dofs(n)` → `(2n, 2n+1)`; `Mesh.dof_index(n, "y")` → `2n+1`.

---

## 10. Getting accurate answers

1. **Use enough elements, especially through the thickness of bending parts.** Q4 elements are slightly too stiff in bending ("shear locking"). Error against beam theory for the L/h = 20 cantilever:

   | Mesh (along × through) | Error |
   |---|---|
   | 20 × 1 | ≈ −32 % |
   | 80 × 4 | ≈ −2.9 % |
   | 160 × 8 | ≈ −0.7 % |
   | 320 × 16 | ≈ −0.1 % |

   Rule of thumb: **at least 8 elements through the height of anything that bends**, with elements roughly square.
2. **Check convergence:** solve, double `nx` and `ny`, and solve again. If the answer barely changes, the mesh is fine enough.
3. **Keep elements well shaped:** close to square, angles near 90°, side ratios under about 5:1.
4. **Point loads and sharp inside corners create infinite stress in theory.** Very close to them, the stress keeps rising as you refine. Judge stresses a little away from those spots, or spread the load with `add_edge_traction`. Displacements and reactions are not affected.
5. **Gauss-point stresses** (`r.gauss_stresses`) are the most accurate. Nodal stresses are extrapolated and averaged for smooth pictures, and can overshoot slightly where stress changes rapidly.
6. **Fine-mesh the hot spots** (holes, notches, fillets). `plate_with_hole_quarter(..., grading=8)` does this automatically. For your own meshes, place more nodes there.
7. **Linear only:** results are valid for small movements and stresses below yield. A safety factor below 1 means the real part would yield, and beyond that point this program's answer is no longer accurate.

---

## 11. Model size and speed

Measured on a single CPU core:

| Model | Unknowns | Solver | Time |
|---|---|---|---|
| Cantilever 160 × 8 | 2,900 | direct | ≈ 0.06 s |
| Plate with hole 40 × 48 | 4,000 | direct | ≈ 0.08 s |
| Cantilever 640 × 40 | 52,500 | iterative (6 to 7 CG iterations) | ≈ 2 s |

- The global matrix is always sparse: a 2,900-unknown model uses 0.6 MB instead of 67 MB dense.
- Element integration is vectorised over all elements at once.
- Memory, not speed, is usually the limit. A few million unknowns is feasible on a 16 GB machine with the iterative solver.
- For the fastest results on big meshes, turn logging down to `WARNING`; the per-solve checks are cheap.

---

## 12. Validation tests

```bash
python3 -m pytest test_validation.py -v      # all tests, verbose
python3 -m pytest test_validation.py -q      # short summary
python3 -m pytest test_validation.py -k cantilever -s   # only cantilever tests, show prints
python3 test_validation.py                   # same as -v -s
python3 app.py --run-tests                   # benchmarks + tests
```

| Test | Proves |
|---|---|
| `test_cantilever_tip_deflection_matches_euler_bernoulli` | **Main benchmark:** tip deflection within 2 % of PL³/3EI |
| `test_cantilever_plane_strain_matches_effective_modulus` | Plane strain matches beam theory with E/(1−ν²) |
| `test_cantilever_equilibrium_and_root_bending_stress` | Reactions balance the load and moment; bending stress matches Mc/I |
| `test_mesh_convergence_is_monotonic` | Error shrinks steadily as the mesh is refined |
| `test_shape_functions_*`, `test_shape_derivatives_*` | Shape functions sum to 1, equal 1 at their own node, and derivatives match finite differences |
| `test_element_stiffness_symmetric_psd_with_three_rigid_body_modes` | Element matrix is symmetric, non-negative, has exactly 3 zero-energy (rigid) modes, and integrates area exactly |
| `test_jacobian_of_rectangle_is_constant`, `test_extrapolation_inverts_gauss_interpolation` | Geometry mapping and stress extrapolation are exact |
| `test_patch_test_*` | Constant-stress states are reproduced **exactly** on a distorted mesh (the standard FE correctness test) |
| `test_constitutive_matrices`, `test_von_mises_formula` | Material matrices and Von Mises formula |
| `test_clockwise_connectivity_rejected`, `test_bad_connectivity_rejected`, `test_distorted_element_raises_jacobian_error` | Bad meshes are caught |
| `test_unconstrained_model_is_singular`, `test_slender_beam_not_flagged_singular_*` | Missing supports are caught, and valid slender models are not wrongly rejected |
| `test_iterative_solver_*`, `test_large_model_uses_iterative_solver_above_50k_dofs` | Iterative solver agrees with the direct solver, kicks in above 50k unknowns, and handles non-convergence |
| `test_prescribed_displacement_reaction`, `test_edge_traction_consistent_loads` | Prescribed movements and distributed loads are exact |

---

## 13. Troubleshooting (every error message explained)

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
| Answers off by 1000× | Units mixed (mm with Pa, etc.) | See §6.1 |
| Deflection a bit smaller than hand calc | Too few elements through the thickness (locking) | Refine: §10 |

---

## 14. Full API reference

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
| `nodes_where(rule)`, `nodes_at(x, y, tol)`, `closest_node(x, y)` | Node selection (§6.3) |
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
Fields are listed in §7. Also has `max_displacement` and `max_von_mises`.

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
`plot_results(...)`, `plot_stress_components(...)` (§8) plus helpers:
- `auto_scale_factor(mesh, results, target_fraction=0.10)`
- `quad_triangulation(mesh, coords)`: splits quads into triangles for Matplotlib contouring.

---

## 15. The maths inside

For readers who want to check the method. Everything below is implemented explicitly in `engine/fea_solver.py`.

**Shape functions** (natural coordinates ξ, η ∈ [−1, 1], node i at (ξᵢ, ηᵢ)):
`Nᵢ = ¼ (1 + ξ ξᵢ)(1 + η ηᵢ)`

**Jacobian:** `J = [∂N/∂ξ; ∂N/∂η] · [x y]` (2×4 · 4×2). det J is checked at every Gauss point and must be > 10⁻¹⁰ × (element size)².

**Global derivatives:** `[∂N/∂x; ∂N/∂y] = J⁻¹ [∂N/∂ξ; ∂N/∂η]`

**Strain-displacement matrix B (3×8)**, relating `ε = [εxx, εyy, γxy]ᵀ` to `uₑ = [u₁ v₁ … u₄ v₄]ᵀ`:
```
| ∂N1/∂x   0      ∂N2/∂x   0     ... |
|   0     ∂N1/∂y   0     ∂N2/∂y  ... |
| ∂N1/∂y  ∂N1/∂x  ∂N2/∂y  ∂N2/∂x ... |
```

**Material matrix D:**
- Plane stress: `E/(1−ν²) · [[1, ν, 0], [ν, 1, 0], [0, 0, (1−ν)/2]]`
- Plane strain: `E/((1+ν)(1−2ν)) · [[1−ν, ν, 0], [ν, 1−ν, 0], [0, 0, (1−2ν)/2]]`

**Element stiffness** (2×2 Gauss, ξ, η = ±1/√3, weights 1):
`Kₑ = Σ_g Bᵀ D B · det J · t`. The result is checked to be symmetric and positive semi-definite.

**Assembly:** all 64 entries of every Kₑ are scattered into a SciPy COO matrix (duplicates summed), then converted to CSR.

**Boundary conditions by partitioning** (f = free, c = constrained):
`K_ff u_f = F_f − K_fc u_c`. Reactions are `R = K u − F` at constrained DOFs.

**Stress recovery:**
- At each Gauss point, `ε = B uₑ` and `σ = D ε`.
- Stresses are extrapolated to the corners with `σ_node = E σ_gauss`, where `E[i, k] = N_k(√3 ξᵢ, √3 ηᵢ)`. Values are then averaged over the elements sharing each node.

**Von Mises:**
`σ_vm = √( ½[(σxx−σyy)² + (σyy−σzz)² + (σzz−σxx)²] + 3τxy² )`, with `σzz = ν(σxx+σyy)` in plane strain and 0 in plane stress.

**Consistent edge loads:** `fₐ = t ∫ Nₐ(s) T(s) ds` with 2-point Gauss per edge.

**Support check:**
- For every connected piece of mesh, the 3 rigid-body motions sampled at the constrained DOFs must have rank 3. This is exact, with no round-off guessing.
- It works because fully integrated Q4 elements have no spurious zero-energy modes.

---

## 16. Limitations

- **2D only**: plane stress or plane strain. No 3D, axisymmetric, plates/shells or beams-as-lines.
- **Linear elastic, small deformations**: no plasticity, contact, buckling, large rotations or creep.
- **One isotropic material per model** (no composites/orthotropic wood, no multi-material).
- **Static loads only**: no dynamics, vibration modes or heat.
- **Q4 elements only** (no triangles or 8-node quads), so fine meshes are needed in bending.
- Gravity, thermal loads and springs aren't built in (gravity can be added as nodal forces; see §6.6).

---

## 17. Project layout

```
CFD/
├── engine/
│   ├── __init__.py          public API exports
│   └── fea_solver.py        Mesh, Material, Q4Element, FEAModel, FEAResults, helpers, exceptions
├── examples/
│   ├── _setup.py            makes the examples runnable from any folder
│   ├── 01_my_beam.py
│   ├── 02_simply_supported_pressure.py
│   ├── 03_plate_with_hole.py
│   ├── 04_custom_mesh.py
│   ├── 05_prescribed_displacement_and_export.py
│   ├── 06_large_model_and_logging.py
│   └── quickstart.ipynb     Colab / Jupyter notebook
├── app.py                   benchmark driver with step-by-step log (CLI)
├── visualizer.py            Matplotlib plots
├── test_validation.py       pytest validation suite
├── requirements.txt
└── output/                  generated pictures and CSVs (not committed)
```
