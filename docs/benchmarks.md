[← Documentation home](README.md) · [Main README](../README.md)

# Built-in benchmarks (`app.py`)

## Running the built-in benchmarks (`app.py`)

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
