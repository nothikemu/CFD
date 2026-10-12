# Q4 FEA Structural Solver

A 2D linear-elastic finite element solver written from scratch in Python. It uses NumPy and SciPy sparse matrices, and Matplotlib for plots.

| File | Purpose |
|------|---------|
| `engine/fea_solver.py` | `Mesh` (CCW check, DOF map, mesh generators), `Material` (plane stress or strain `D`), `Q4Element` (shape functions, Jacobian, `B`, 2x2 Gauss stiffness), `FEAModel` (COO→CSR assembly, Dirichlet/Neumann BCs, partitioned solve, reactions, stress recovery, Von Mises) |
| `app.py` | Runs the benchmarks with step-by-step logs: a cantilever compared with Euler-Bernoulli, and a plate with a hole compared with Kirsch/Heywood |
| `visualizer.py` | Draws the undeformed wireframe and a Von Mises contour on the scaled deformed shape |
| `test_validation.py` | pytest validation suite (cantilever within 2%, patch tests, element checks, error guards, solver paths) |

```bash
pip install -r requirements.txt
python app.py                 # benchmarks + figures in ./output
python app.py --run-tests     # ...then the validation suite
python -m pytest test_validation.py -v
```

Solver selection: SuperLU direct factorisation up to 50,000 free DOFs. Above that, conjugate gradient with a Jacobi-scaled ILU preconditioner. If CG fails to converge, the solver falls back to the direct method.
