[← Documentation home](README.md) · [Main README](../README.md)

# Theory and limitations

## The maths inside

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

## Limitations

- **2D only**: plane stress or plane strain. No 3D, axisymmetric, plates/shells or beams-as-lines.
- **Linear elastic, small deformations**: no plasticity, contact, buckling, large rotations or creep.
- **One isotropic material per model** (no composites/orthotropic wood, no multi-material).
- **Static loads only**: no dynamics, vibration modes or heat.
- **Q4 elements only** (no triangles or 8-node quads), so fine meshes are needed in bending.
- Gravity, thermal loads and springs aren't built in (gravity can be added as nodal forces; see [Loads](user-guide.md#loads-neumann-boundary-conditions)).
