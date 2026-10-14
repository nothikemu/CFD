# Documentation

Start with the guides from the top of this list; the reference pages are for looking things up.

## Guides

| Page | Read it to… |
|---|---|
| [Installation](installation.md) | Get it running in Google Colab (no install), on Windows / Mac / Linux, or in VS Code |
| [Tutorial](tutorial.md) | Do a complete first analysis step by step, with the output you should see at each step |
| [Examples](examples.md) | Copy a ready-made scenario: beams, pressure, holes, custom shapes, CSV export, large models |
| [User guide](user-guide.md) | Build any model: units, meshes, choosing nodes, materials, supports, loads |
| [Results, export and plots](results-and-plots.md) | Read every result, export to Excel/CSV, customise pictures |
| [Built-in benchmarks (`app.py`)](benchmarks.md) | Use the command-line program and understand its step-by-step log |

## Reference

| Page | Contents |
|---|---|
| [Configuration and performance](configuration.md) | Solver selection and tolerances, logging, advanced options, model size and speed |
| [Accuracy and validation](accuracy-and-validation.md) | How fine the mesh must be, and what each of the 24 tests proves |
| [Troubleshooting](troubleshooting.md) | Every error and warning message, what causes it and how to fix it |
| [API reference](api-reference.md) | Every class, method, argument and result field |
| [Theory and limitations](theory.md) | The equations implemented, and what the solver cannot do |

## Pictures

The images in [`images/`](images) are produced by [`make_images.py`](make_images.py):

```bash
python3 docs/make_images.py
```
