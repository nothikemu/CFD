[← Documentation home](README.md) · [Main README](../README.md)

# Configuration and performance

## Configuration: solver, accuracy, logging

### Solver settings

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

### Logging (the step-by-step maths output)

The engine reports through Python's standard `logging` under the name `engine.fea`. `app.py` turns it on automatically. In your own scripts:

```python
import logging, sys
logging.basicConfig(level=logging.INFO, stream=sys.stdout, format="%(message)s")
logging.getLogger("matplotlib").setLevel(logging.WARNING)   # hide Matplotlib chatter
```
- **Silence everything except problems:** `logging.getLogger("engine.fea").setLevel(logging.WARNING)`.
- **Save the log to a file:** add `filename="run.log"` to `basicConfig` (and remove `stream=`).

### Advanced knobs

- `Mesh(nodes, elements, validate=False)` skips the up-front mesh checks. The Jacobian checks during the solve still run.
- `FEAModel.assemble()` builds the global stiffness matrix without solving: `K = model.assemble()` (SciPy CSR). Useful for your own analysis, e.g. eigenvalues.
- `model.K`, `model.Ke` (all element matrices, shape (elements, 8, 8)), `model.detJ` (Jacobian determinants, shape (elements, 4)) are available after assembly.
- `model.constrained_dofs` lists every fixed DOF index (node `n` → x is `2n`, y is `2n+1`).
- `Mesh.node_dofs(n)` → `(2n, 2n+1)`; `Mesh.dof_index(n, "y")` → `2n+1`.

## Model size and speed

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
