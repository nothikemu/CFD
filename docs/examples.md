[← Documentation home](README.md) · [Main README](../README.md)

# Examples

Six ready-to-run scripts live in [`examples/`](../examples). Each one runs on its own, from any folder or with an editor's **Run** button, and writes its pictures and CSV files to `output/`.

```bash
python3 examples/01_my_beam.py
python3 examples/02_simply_supported_pressure.py
python3 examples/03_plate_with_hole.py
python3 examples/04_custom_mesh.py
python3 examples/05_prescribed_displacement_and_export.py
python3 examples/06_large_model_and_logging.py
```

Run them all at once:

```bash
# Mac / Linux
for f in examples/0*.py; do python3 "$f"; done
```
```powershell
# Windows PowerShell
Get-ChildItem examples\0*.py | ForEach-Object { python $_.FullName }
```

The fastest way to start your own job: **copy the example closest to your problem, rename it and edit the numbers.**

| # | File | Features shown |
|---|---|---|
| 1 | `01_my_beam.py` | Rectangle mesh, clamp, point load, safety factor, plot |
| 2 | `02_simply_supported_pressure.py` | Pin + roller supports, **edge pressure**, reactions, component plot |
| 3 | `03_plate_with_hole.py` | **Symmetry**, curved graded mesh, stress concentration |
| 4 | `04_custom_mesh.py` | **Hand-built mesh**, plane strain, `nodes_where`, angled load |
| 5 | `05_prescribed_displacement_and_export.py` | **Prescribed displacement**, reaction force, **CSV export** |
| 6 | `06_large_model_and_logging.py` | **52 000-unknown** model, iterative solver settings, full log |
| — | `quickstart.ipynb` | Notebook version for Google Colab / Jupyter |

---

## 1. Cantilever beam: `01_my_beam.py`

A bar clamped at one end with a point load at the other. Its **SETTINGS** block at the top holds every number you'd want to change:

```python
LENGTH = 1.0            # m
HEIGHT = 0.1            # m
THICKNESS = 0.02        # m (out-of-plane depth)
E = 210e9               # Pa  (steel 210e9, aluminium 70e9, titanium 114e9)
NU = 0.3                # Poisson's ratio
YIELD = 250e6           # Pa  (used for the safety factor)
FORCE = -2000.0         # N   (negative = downward)
NX, NY = 100, 10        # elements along the length / through the height
```

Expected output:
```
tip drop (FEA)        : -1.906742e-03 m
tip drop (hand calc)  : -1.904762e-03 m   (P L^3 / 3EI)
max Von Mises stress  : 65.02 MPa
minimum safety factor : 3.85
picture saved to output/01_my_beam.png
```

![cantilever](images/tutorial_beam.png)

A line-by-line explanation of the same model is in the [Tutorial](tutorial.md).

---

## 2. Beam on two supports with pressure: `02_simply_supported_pressure.py`

A 3 m beam resting on a pin (left) and a roller (right), with uniform pressure on its top face.

Key lines:
```python
model.fix(left, "xy")                         # pin: holds x and y
model.fix(right, "y")                         # roller: holds y, free to slide in x
top = mesh.nodes_at(y=H / 2)
top = top[np.argsort(mesh.nodes[top, 0])]     # edge nodes must be in order
model.add_edge_traction(top, (0.0, -PRESSURE))
```

Expected output:
```
mid-span deflection FEA  : 3.760388e-03 m
mid-span deflection beam : 3.750000e-03 m  (5 q L^4 / 384 E I)
support reactions (y)    : left 15000.0 N, right 15000.0 N, expected 15000.0 N each
max Von Mises            : 59.52 MPa
```
The mid-span deflection is within 0.3 % of beam theory. Each support carries exactly half the total load, and the small extra FEA deflection is the shear deformation that simple beam theory ignores.

![simply supported](images/simply_supported.png)

---

## 3. Plate with a hole: `03_plate_with_hole.py`

A plate pulled from both sides, with a round hole in the middle. Because the part is symmetric, only **one quarter** is modelled. Symmetry is enforced with rollers on the cut lines:

```python
mesh = Mesh.plate_with_hole_quarter(radius=A, half_width=W, n_radial=40,
                                    n_circumferential=48, grading=8.0)
model.fix(mesh.nodes_at(x=0.0), "x")    # symmetry line x = 0
model.fix(mesh.nodes_at(y=0.0), "y")    # symmetry line y = 0
model.add_edge_traction(edge, (S0, 0.0))
```

Expected output:
```
stress concentration factor sigma_xx(0, a) / s0 = 3.132 (textbook ~3 for an infinite plate)
minimum safety factor = 1.50
```
The hole roughly **triples** the stress at its edge, the classic result. This plate is square and finite, so the value sits slightly above the infinite-plate textbook 3.0.

![plate with hole](images/plate_with_hole.png)

---

## 4. Your own mesh: `04_custom_mesh.py`

A tapered concrete bracket built **from raw coordinates and connectivity**, so you can see how to describe any shape yourself:

```python
nodes = []
for j in range(NY + 1):                       # rows bottom -> top
    for i in range(NX + 1):                   # columns left -> right
        x = L * i / NX
        h = H_ROOT + (H_TIP - H_ROOT) * x / L
        nodes.append([x, -h / 2 + h * j / NY])

elements = []
for j in range(NY):
    for i in range(NX):
        n0 = j * (NX + 1) + i                 # bottom-left corner
        elements.append([n0, n0 + 1, n0 + NX + 2, n0 + NX + 1])   # counter-clockwise

mesh = Mesh(nodes, elements)
concrete = Material(E=30e9, nu=0.2, thickness=1.0, plane="strain")
wall = mesh.nodes_where(lambda x, y: np.isclose(x, 0.0))     # custom node selection
model.add_nodal_force(tip, fx=20e3, fy=-50e3)                # angled load
```

Expected output:
```
Mesh(n_nodes=369, n_elements=320, n_dofs=738)
tip displacement: ux = +3.8517e-06 m, uy = -2.6962e-04 m
max Von Mises   : 2.372 MPa
```

![custom mesh](images/custom_mesh.png)

The rules for writing your own `nodes` / `elements` are in the [User guide → mesh](user-guide.md#shape-mesh).

---

## 5. Prescribed displacement and CSV export: `05_prescribed_displacement_and_export.py`

Instead of applying a force, this **moves** the right edge of a block by 0.1 mm and asks what force that takes:

```python
model.fix(mesh.nodes_at(x=0.0), "x")
model.fix(mesh.closest_node(0.0, 0.0), "y")
model.fix(right, "x", value=1e-4)           # move by +0.1 mm
r = model.solve()
pull = r.reactions[2 * right].sum()         # force needed
```

Expected output:
```
force needed to stretch the block 0.1 mm: 35000.0 N
CSV files written: .../output/05_nodal_results.csv, .../output/05_gauss_results.csv
```
Check: E · A · strain = 70e9 × (0.1 × 0.01) × (0.0001 / 0.2) = 35 000 N exactly.

The CSV files open directly in Excel, LibreOffice or Google Sheets:

| File | Columns |
|---|---|
| `05_nodal_results.csv` | `node, x, y, ux, uy, sxx, syy, txy, von_mises, Rx, Ry` |
| `05_gauss_results.csv` | `element, x, y, exx, eyy, gxy, sxx, syy, txy, von_mises` |

---

## 6. Large model and the full log: `06_large_model_and_logging.py`

A 640 × 40 element beam with **52 480 unknowns**. Above 50 000 the engine automatically switches to the iterative conjugate-gradient solver. All solver settings are shown with their defaults:

```python
model = FEAModel(
    mesh,
    Material(210e9, 0.3, thickness=0.05),
    iterative_threshold=50_000,   # switch point
    cg_rtol=1e-10,                # CG stopping tolerance
    cg_maxiter=None,              # None = automatic
    direct_fallback=True,         # if CG stalls, finish with the direct solver
)
```

With logging switched on, it prints every step of the maths. The end looks like:
```
[STEP 6] Linear solve: Preconditioned conjugate gradient (Jacobi-scaled ILU) (6 iterations), ...
...
solver used: Preconditioned conjugate gradient (Jacobi-scaled ILU) | iterations: 6
```
The whole run takes about 2 seconds. See [Configuration](configuration.md) for every solver option.

---

## Notebook: `quickstart.ipynb`

The same workflow as a notebook. Open it in Google Colab (no installation needed) or in Jupyter. Its cells are:

1. download and install;
2. run the tests;
3. your beam, with an editable SETTINGS box and an inline picture;
4. the built-in benchmarks with the full log;
5. a further example.

How to open it: [Installation → Colab](installation.md#fastest-start-no-terminal-google-colab).
