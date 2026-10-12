"""
Validation suite for the Q4 FEA engine.

Primary benchmark: an Euler-Bernoulli cantilever (L/h = 20) clamped at x = 0
and loaded by a concentrated downward tip force. The FE tip deflection must
match delta = P L^3 / (3 E I) within 2 %.

Supporting checks: element-level properties (symmetry, rigid-body modes,
exact area integration), the constant-stress patch test on a distorted mesh,
constitutive matrices, Jacobian / orientation guards, singular-system
detection, prescribed displacements, the iterative-solver path, traction
integration and stress extrapolation.

Run with:  python -m pytest test_validation.py -v      (or python test_validation.py)
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

from engine import (
    DistortedElementError,
    FEAModel,
    Material,
    Mesh,
    MeshTopologyError,
    Q4Element,
    SingularSystemError,
    SolverConvergenceError,
    cantilever_tip_deflection,
    von_mises,
)

TOLERANCE = 0.02

# Cantilever benchmark data (steel, rectangular section)
L, H, T = 2.0, 0.1, 0.05
E, NU = 210e9, 0.3
P = 5_000.0
I = T * H**3 / 12.0


def build_cantilever(nx: int = 160, ny: int = 8, plane: str = "stress", **model_kwargs):
    mesh = Mesh.structured_rectangle(L, H, nx, ny, origin=(0.0, -0.5 * H))
    model = FEAModel(mesh, Material(E, NU, thickness=T, plane=plane), **model_kwargs)
    model.fix(mesh.nodes_at(x=0.0), "xy")
    tip = mesh.closest_node(L, 0.0)
    model.add_nodal_force(tip, fy=-P)
    return model, tip


# --------------------------------------------------------------------------- #
# Primary benchmark                                                           #
# --------------------------------------------------------------------------- #
def test_cantilever_tip_deflection_matches_euler_bernoulli():
    model, tip = build_cantilever()
    results = model.solve()
    delta_fe = -results.displacements[tip, 1]
    delta_exact = cantilever_tip_deflection(P, L, E, I)
    error = abs(delta_fe - delta_exact) / delta_exact
    print(
        f"\n  Euler-Bernoulli delta = {delta_exact:.6e} m | FE delta = {delta_fe:.6e} m "
        f"| error = {100 * error:.4f} %"
    )
    assert delta_fe > 0.0, "tip must deflect downward"
    assert error < TOLERANCE, f"tip deflection error {100 * error:.3f} % exceeds 2 %"


def test_cantilever_plane_strain_matches_effective_modulus():
    model, tip = build_cantilever(plane="strain")
    delta_fe = -model.solve().displacements[tip, 1]
    delta_exact = cantilever_tip_deflection(P, L, E / (1.0 - NU**2), I)
    assert abs(delta_fe - delta_exact) / delta_exact < TOLERANCE


def test_cantilever_equilibrium_and_root_bending_stress():
    model, _ = build_cantilever()
    r = model.solve()
    np.testing.assert_allclose(r.reactions[1::2].sum(), P, rtol=1e-9)
    assert abs(r.reactions[0::2].sum()) < 1e-8 * P
    # reaction moment about the root equals P L
    root = model.mesh.nodes_at(x=0.0)
    m_root = -np.sum(model.mesh.nodes[root, 1] * r.reactions[2 * root])
    np.testing.assert_allclose(abs(m_root), P * L, rtol=1e-9)
    assert abs(r.diagnostics["equilibrium_fx"]) < 1e-8 * P
    assert abs(r.diagnostics["equilibrium_fy"]) < 1e-8 * P
    # bending stress one element away from the clamp vs M c / I (Gauss points sit inside)
    x_probe = 0.25 * L
    sigma_beam = P * (L - x_probe) * (0.5 * H) / I
    gx = r.gauss_points[..., 0]
    near = np.abs(gx - x_probe) < 0.5 * L / 160
    sigma_fe = np.abs(r.gauss_stresses[..., 0][near]).max()
    # Gauss points lie at +/- h/(2*8) * (1 - 1/sqrt 3) inside the fibre -> compare to the linear profile
    y_gauss = np.abs(r.gauss_points[..., 1][near]).max()
    sigma_expected = sigma_beam * y_gauss / (0.5 * H)
    assert abs(sigma_fe - sigma_expected) / sigma_expected < TOLERANCE


def test_mesh_convergence_is_monotonic():
    delta_exact = cantilever_tip_deflection(P, L, E, I)
    errors = []
    for nx, ny in [(20, 1), (40, 2), (80, 4), (160, 8)]:
        model, tip = build_cantilever(nx, ny)
        errors.append(abs(-model.solve().displacements[tip, 1] - delta_exact) / delta_exact)
    assert all(e1 > e2 for e1, e2 in zip(errors, errors[1:])), errors


# --------------------------------------------------------------------------- #
# Element-level checks                                                        #
# --------------------------------------------------------------------------- #
DISTORTED_QUAD = np.array([[0.0, 0.0], [2.1, 0.3], [1.7, 1.6], [-0.2, 1.1]])


def test_shape_functions_partition_of_unity_and_kronecker_delta():
    for xi, eta in [(-0.3, 0.7), (0.0, 0.0), (0.9, -0.95)]:
        assert np.isclose(Q4Element.shape_functions(xi, eta).sum(), 1.0)
        assert np.allclose(Q4Element.shape_derivatives(xi, eta).sum(axis=1), 0.0)
    for i, (xi, eta) in enumerate(Q4Element.NODE_NATURAL):
        np.testing.assert_allclose(Q4Element.shape_functions(xi, eta), np.eye(4)[i], atol=1e-15)


def test_shape_derivatives_match_finite_differences():
    xi, eta, h = 0.31, -0.52, 1e-6
    dN = Q4Element.shape_derivatives(xi, eta)
    fd_xi = (Q4Element.shape_functions(xi + h, eta) - Q4Element.shape_functions(xi - h, eta)) / (2 * h)
    fd_eta = (Q4Element.shape_functions(xi, eta + h) - Q4Element.shape_functions(xi, eta - h)) / (2 * h)
    np.testing.assert_allclose(dN[0], fd_xi, atol=1e-9)
    np.testing.assert_allclose(dN[1], fd_eta, atol=1e-9)


def test_element_stiffness_symmetric_psd_with_three_rigid_body_modes():
    D = Material(E, NU).D
    Ke, detJ = Q4Element.stiffness(DISTORTED_QUAD, D, 0.01)
    Ke = Ke[0]
    np.testing.assert_allclose(Ke, Ke.T, rtol=0, atol=1e-12 * np.abs(Ke).max())
    eig = np.linalg.eigvalsh(Ke)
    tol = 1e-10 * eig.max()
    assert np.sum(np.abs(eig) < tol) == 3, eig
    assert np.all(eig > -tol)
    # rigid translations and rotation produce zero nodal forces
    x, y = DISTORTED_QUAD[:, 0], DISTORTED_QUAD[:, 1]
    for mode in (
        np.tile([1.0, 0.0], 4),
        np.tile([0.0, 1.0], 4),
        np.column_stack([-y, x]).ravel(),
    ):
        assert np.linalg.norm(Ke @ mode) < 1e-9 * np.abs(Ke).max()
    # sum_g w_g det(J_g) integrates the element area exactly
    area = Mesh(DISTORTED_QUAD, [[0, 1, 2, 3]]).signed_areas()[0]
    np.testing.assert_allclose(detJ.sum(), area, rtol=1e-13)


def test_jacobian_of_rectangle_is_constant():
    rect = np.array([[0.0, 0.0], [4.0, 0.0], [4.0, 2.0], [0.0, 2.0]])
    for xi, eta in Q4Element.GAUSS_POINTS:
        J = Q4Element.jacobian(rect, xi, eta)[0]
        np.testing.assert_allclose(J, [[2.0, 0.0], [0.0, 1.0]], atol=1e-15)


def test_extrapolation_inverts_gauss_interpolation():
    Eg = Q4Element.extrapolation_matrix()
    A = np.array([Q4Element.shape_functions(xi, eta) for xi, eta in Q4Element.GAUSS_POINTS])
    np.testing.assert_allclose(Eg @ A, np.eye(4), atol=1e-13)


# --------------------------------------------------------------------------- #
# Patch test (completeness / consistency on a distorted mesh)                 #
# --------------------------------------------------------------------------- #
def _distorted_patch():
    nodes = np.array(
        [
            [0.0, 0.0], [2.0, 0.0], [4.0, 0.0],
            [0.0, 1.5], [2.4, 1.7], [4.0, 2.0],
            [0.0, 4.0], [1.6, 4.0], [4.0, 4.0],
        ]
    )
    elements = np.array([[0, 1, 4, 3], [1, 2, 5, 4], [3, 4, 7, 6], [4, 5, 8, 7]])
    return Mesh(nodes, elements)


@pytest.mark.parametrize("plane", ["stress", "strain"])
def test_patch_test_linear_displacement_reproduced_exactly(plane):
    mesh = _distorted_patch()
    mat = Material(1000.0, 0.25, thickness=1.0, plane=plane)
    model = FEAModel(mesh, mat)
    a = np.array([1e-3, 2e-3, -1e-3, 0.5e-3, 1.5e-3, -0.7e-3])

    def exact(xy):
        x, y = xy[:, 0], xy[:, 1]
        return np.column_stack([a[0] + a[1] * x + a[2] * y, a[3] + a[4] * x + a[5] * y])

    boundary = [n for n in range(mesh.n_nodes) if n != 4]
    for n in boundary:
        model.fix(n, "xy", value=tuple(exact(mesh.nodes[[n]])[0]))
    r = model.solve()
    np.testing.assert_allclose(r.displacements, exact(mesh.nodes), atol=1e-14)
    eps = np.array([a[1], a[5], a[2] + a[4]])
    sigma = mat.D @ eps
    np.testing.assert_allclose(r.gauss_stresses.reshape(-1, 3), np.tile(sigma, (16, 1)), rtol=1e-9)
    np.testing.assert_allclose(r.nodal_stresses, np.tile(sigma, (9, 1)), rtol=1e-9)


def test_patch_test_uniform_traction():
    mesh = _distorted_patch()
    mat = Material(1000.0, 0.25, thickness=0.5)
    model = FEAModel(mesh, mat)
    s0 = 7.0
    model.fix(mesh.nodes_at(x=0.0), "x")
    model.fix(0, "y")
    model.add_edge_traction([2, 5, 8], (s0, 0.0))
    r = model.solve()
    np.testing.assert_allclose(r.gauss_stresses[..., 0], s0, rtol=1e-10)
    np.testing.assert_allclose(r.gauss_stresses[..., 1:], 0.0, atol=1e-10 * s0)
    np.testing.assert_allclose(r.displacements[8, 0], s0 / mat.E * 4.0, rtol=1e-10)
    np.testing.assert_allclose(r.von_mises, s0, rtol=1e-10)


# --------------------------------------------------------------------------- #
# Material                                                                    #
# --------------------------------------------------------------------------- #
def test_constitutive_matrices():
    m = Material(200.0, 0.25, plane="stress")
    c = 200.0 / (1 - 0.0625)
    np.testing.assert_allclose(m.D, c * np.array([[1, 0.25, 0], [0.25, 1, 0], [0, 0, 0.375]]))
    m.plane = "plane_strain"
    c = 200.0 / (1.25 * 0.5)
    np.testing.assert_allclose(m.D, c * np.array([[0.75, 0.25, 0], [0.25, 0.75, 0], [0, 0, 0.25]]))
    np.testing.assert_allclose(m.D[2, 2], m.shear_modulus)
    with pytest.raises(ValueError):
        Material(200.0, 0.5)
    with pytest.raises(ValueError):
        Material(-1.0, 0.3)
    with pytest.raises(ValueError):
        Material(200.0, 0.3, plane="axisymmetric")


def test_von_mises_formula():
    s = np.array([[100.0, 0.0, 0.0], [0.0, 0.0, 50.0], [80.0, -40.0, 30.0]])
    expected = [100.0, np.sqrt(3) * 50.0, np.sqrt(80**2 + 40**2 + 80 * 40 + 3 * 30**2)]
    np.testing.assert_allclose(von_mises(s), expected)


# --------------------------------------------------------------------------- #
# Guards and error handling                                                   #
# --------------------------------------------------------------------------- #
def test_clockwise_connectivity_rejected():
    nodes = [[0, 0], [1, 0], [1, 1], [0, 1]]
    with pytest.raises(MeshTopologyError, match="counter-clockwise"):
        Mesh(nodes, [[0, 3, 2, 1]])


def test_bad_connectivity_rejected():
    nodes = [[0, 0], [1, 0], [1, 1], [0, 1]]
    with pytest.raises(MeshTopologyError):
        Mesh(nodes, [[0, 1, 2, 7]])
    with pytest.raises(MeshTopologyError):
        Mesh(nodes, [[0, 1, 1, 3]])


def test_distorted_element_raises_jacobian_error():
    # CCW with positive area, but node 3 is pulled inside -> re-entrant corner
    arrow = np.array([[0.0, 0.0], [2.0, 0.0], [0.2, 0.2], [0.0, 2.0]])
    mesh = Mesh(arrow, [[0, 1, 2, 3]])
    assert mesh.signed_areas()[0] > 0
    model = FEAModel(mesh, Material(E, NU))
    model.fix([0, 1], "xy")
    with pytest.raises(DistortedElementError, match="Jacobian"):
        model.solve()


def test_unconstrained_model_is_singular():
    mesh = Mesh.structured_rectangle(1.0, 1.0, 2, 2)
    model = FEAModel(mesh, Material(E, NU))
    model.add_nodal_force(8, fx=1.0)
    with pytest.raises(SingularSystemError):
        model.solve()
    model.fix(0, "xy")  # rotation about node 0 still free
    with pytest.raises(SingularSystemError):
        model.solve()


# --------------------------------------------------------------------------- #
# Solver paths                                                                #
# --------------------------------------------------------------------------- #
def test_iterative_solver_matches_direct_solver():
    direct, tip = build_cantilever(80, 4)
    u_direct = direct.solve().u
    iterative, _ = build_cantilever(80, 4, iterative_threshold=10, cg_rtol=1e-12)
    r = iterative.solve()
    assert "conjugate gradient" in r.diagnostics["solver"]
    np.testing.assert_allclose(r.u, u_direct, rtol=0, atol=1e-8 * np.abs(u_direct).max())


def test_large_model_uses_iterative_solver_above_50k_dofs():
    model, tip = build_cantilever(640, 40)  # 52 480 free DOFs
    r = model.solve()
    assert r.diagnostics["n_free"] > FEAModel.ITERATIVE_THRESHOLD
    assert "conjugate gradient" in r.diagnostics["solver"]
    delta_exact = cantilever_tip_deflection(P, L, E, I)
    assert abs(-r.displacements[tip, 1] - delta_exact) / delta_exact < TOLERANCE


def test_iterative_solver_non_convergence_handling():
    strict, _ = build_cantilever(80, 4, iterative_threshold=10, cg_rtol=1e-14, cg_maxiter=1,
                                 direct_fallback=False)
    with pytest.raises(SolverConvergenceError):
        strict.solve()
    fallback, _ = build_cantilever(80, 4, iterative_threshold=10, cg_rtol=1e-14, cg_maxiter=1)
    r = fallback.solve()
    assert "fallback" in r.diagnostics["solver"]
    assert r.diagnostics["relative_residual"] < 1e-10


def test_slender_beam_not_flagged_singular_and_partial_supports_detected():
    # L/h = 1000: badly conditioned but properly supported -> must solve
    mesh = Mesh.structured_rectangle(1000.0, 1.0, 1000, 2)
    model = FEAModel(mesh, Material(1.0, 0.3))
    model.fix(mesh.nodes_at(x=0.0), "xy")
    model.add_nodal_force(mesh.closest_node(1000.0, 0.5), fy=-1e-6)
    assert np.all(np.isfinite(model.solve().u))
    # rollers on one horizontal edge only leave x-translation free
    model.clear_constraints()
    model.fix(mesh.nodes_at(y=0.0), "y")
    with pytest.raises(SingularSystemError, match="rigid-body"):
        model.solve()


def test_prescribed_displacement_reaction():
    mesh = Mesh.structured_rectangle(2.0, 1.0, 4, 2)
    mat = Material(1000.0, 0.0, thickness=1.0)
    model = FEAModel(mesh, mat)
    model.fix(mesh.nodes_at(x=0.0), "x")
    model.fix(0, "y")
    right = mesh.nodes_at(x=2.0)
    model.fix(right, "x", value=0.01)
    r = model.solve()
    # uniaxial strain 0.005 with nu = 0 -> total force E * eps * A
    np.testing.assert_allclose(r.reactions[2 * right].sum(), 1000.0 * 0.005 * 1.0, rtol=1e-10)


def test_edge_traction_consistent_loads():
    mesh = Mesh.structured_rectangle(1.0, 1.0, 1, 3)
    model = FEAModel(mesh, Material(E, NU, thickness=2.0))
    edge = mesh.nodes_at(x=1.0)
    edge = edge[np.argsort(mesh.nodes[edge, 1])]
    model.add_edge_traction(edge, lambda x, y: (0.0, 3.0 * y**2))
    # consistent loads preserve every moment of the traction that is linear per segment:
    # resultant t * int_0^1 3 y^2 dy = 2 ;  first moment t * int_0^1 3 y^3 dy = 1.5
    np.testing.assert_allclose(model.F[1::2].sum(), 2.0, rtol=1e-13)
    np.testing.assert_allclose(np.sum(mesh.nodes[:, 1] * model.F[1::2]), 1.5, rtol=1e-13)
    assert np.all(model.F[0::2] == 0.0)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v", "-s"]))
