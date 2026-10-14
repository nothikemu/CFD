[← Documentation home](README.md) · [Main README](../README.md)

# Accuracy and validation

## Getting accurate answers

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

## Validation tests

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
