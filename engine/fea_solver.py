"""
2D linear-elastic finite element solver built on 4-node isoparametric
quadrilateral (Q4) elements.

Pipeline
--------
1.  Mesh            : nodal coordinates + CCW element connectivity, DOF map.
2.  Material        : isotropic linear elasticity, plane stress / plane strain D.
3.  Q4Element       : bilinear shape functions, Jacobian, B operator,
                      2x2 Gauss-Legendre stiffness integration (vectorised over
                      all elements at once).
4.  FEAModel        : sparse COO -> CSR global assembly, Dirichlet / Neumann
                      boundary conditions, partitioned solve (direct or
                      ILU-preconditioned CG), reactions, equilibrium audit,
                      Gauss-point stress recovery, nodal extrapolation and
                      Von Mises evaluation.

Sign / ordering conventions
---------------------------
*   Global DOF numbering : node n owns DOFs (2n, 2n+1) = (u_x, u_y).
*   Element DOF ordering : [u1, v1, u2, v2, u3, v3, u4, v4].
*   Local node ordering  : counter-clockwise, natural coordinates
                           1:(-1,-1)  2:(+1,-1)  3:(+1,+1)  4:(-1,+1).
*   Strain / stress      : Voigt vectors [e_xx, e_yy, gamma_xy] and
                           [s_xx, s_yy, t_xy] (engineering shear strain).
"""

from __future__ import annotations

import logging
import time
import warnings
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import scipy.sparse as sp
import scipy.sparse.csgraph  # noqa: F401  (registers sp.csgraph)
import scipy.sparse.linalg as spla

logger = logging.getLogger("engine.fea")

ArrayLike = Union[Sequence[float], np.ndarray]


# --------------------------------------------------------------------------- #
# Exceptions                                                                  #
# --------------------------------------------------------------------------- #
class FEAError(RuntimeError):
    """Base class for all solver errors."""


class MeshTopologyError(FEAError):
    """Raised for malformed node / connectivity input or clockwise elements."""


class DistortedElementError(FEAError):
    """Raised when the Jacobian determinant is negative or near zero."""


class SingularSystemError(FEAError):
    """Raised when the reduced stiffness matrix is singular (rigid-body modes)."""


class SolverConvergenceError(FEAError):
    """Raised when the iterative solver fails to reach the requested tolerance."""


class ElementIntegrityError(FEAError):
    """Raised when an element stiffness matrix is not symmetric / PSD."""


# --------------------------------------------------------------------------- #
# Material                                                                    #
# --------------------------------------------------------------------------- #
class Material:
    """
    Isotropic, homogeneous, linear-elastic material for 2D continuum analysis.

    Parameters
    ----------
    E : float
        Young's modulus (> 0).
    nu : float
        Poisson's ratio (-1 < nu < 0.5).
    thickness : float
        Out-of-plane thickness t (> 0). For plane strain this is the unit
        depth the forces are referred to.
    plane : {"stress", "strain"}
        Constitutive assumption.
    yield_strength : float, optional
        Used to report Von Mises safety factors.
    """

    PLANE_STRESS = "stress"
    PLANE_STRAIN = "strain"

    def __init__(
        self,
        E: float,
        nu: float,
        thickness: float = 1.0,
        plane: str = PLANE_STRESS,
        yield_strength: Optional[float] = None,
        name: str = "material",
    ) -> None:
        if not np.isfinite(E) or E <= 0.0:
            raise ValueError(f"Young's modulus must be positive and finite, got {E!r}")
        if not (-1.0 < nu < 0.5):
            raise ValueError(f"Poisson's ratio must satisfy -1 < nu < 0.5, got {nu!r}")
        if not np.isfinite(thickness) or thickness <= 0.0:
            raise ValueError(f"Thickness must be positive and finite, got {thickness!r}")
        if yield_strength is not None and yield_strength <= 0.0:
            raise ValueError(f"Yield strength must be positive, got {yield_strength!r}")
        self.E = float(E)
        self.nu = float(nu)
        self.thickness = float(thickness)
        self.yield_strength = None if yield_strength is None else float(yield_strength)
        self.name = name
        self._plane = self._normalise_plane(plane)

    @staticmethod
    def _normalise_plane(plane: str) -> str:
        key = str(plane).strip().lower().replace("_", " ").replace("-", " ")
        if key in ("stress", "plane stress", "ps"):
            return Material.PLANE_STRESS
        if key in ("strain", "plane strain", "pe"):
            return Material.PLANE_STRAIN
        raise ValueError(f"Unknown plane formulation {plane!r}; use 'stress' or 'strain'")

    @property
    def plane(self) -> str:
        return self._plane

    @plane.setter
    def plane(self, value: str) -> None:
        self._plane = self._normalise_plane(value)

    @property
    def shear_modulus(self) -> float:
        return self.E / (2.0 * (1.0 + self.nu))

    @property
    def D(self) -> np.ndarray:
        """3x3 constitutive matrix relating [e_xx, e_yy, g_xy] to [s_xx, s_yy, t_xy]."""
        E, nu = self.E, self.nu
        if self._plane == Material.PLANE_STRESS:
            c = E / (1.0 - nu * nu)
            return c * np.array(
                [
                    [1.0, nu, 0.0],
                    [nu, 1.0, 0.0],
                    [0.0, 0.0, 0.5 * (1.0 - nu)],
                ]
            )
        c = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
        return c * np.array(
            [
                [1.0 - nu, nu, 0.0],
                [nu, 1.0 - nu, 0.0],
                [0.0, 0.0, 0.5 * (1.0 - 2.0 * nu)],
            ]
        )

    def out_of_plane_stress(self, stress: np.ndarray) -> np.ndarray:
        """sigma_zz for a stress array whose last axis is [s_xx, s_yy, t_xy]."""
        if self._plane == Material.PLANE_STRAIN:
            return self.nu * (stress[..., 0] + stress[..., 1])
        return np.zeros(stress.shape[:-1])

    def __repr__(self) -> str:
        return (
            f"Material(name={self.name!r}, E={self.E:.6g}, nu={self.nu:.4g}, "
            f"t={self.thickness:.6g}, plane={self._plane!r})"
        )


# --------------------------------------------------------------------------- #
# Mesh                                                                        #
# --------------------------------------------------------------------------- #
class Mesh:
    """
    Q4 mesh: nodal coordinates (n_nodes, 2) and connectivity (n_elems, 4).

    Every element must be numbered counter-clockwise; this is checked through
    the signed (shoelace) area of each element. Non-convex elements are
    reported; their Gauss-point Jacobians are checked again during integration.
    """

    DOFS_PER_NODE = 2

    def __init__(self, nodes: ArrayLike, elements: ArrayLike, validate: bool = True) -> None:
        nodes = np.array(nodes, dtype=np.float64, copy=True)
        elements = np.array(elements, copy=True)
        if nodes.ndim != 2 or nodes.shape[1] != 2:
            raise MeshTopologyError(f"nodes must have shape (n_nodes, 2), got {nodes.shape}")
        if nodes.shape[0] < 4:
            raise MeshTopologyError("a Q4 mesh needs at least four nodes")
        if not np.all(np.isfinite(nodes)):
            raise MeshTopologyError("nodal coordinates contain NaN or Inf")
        if elements.ndim != 2 or elements.shape[1] != 4:
            raise MeshTopologyError(
                f"elements must have shape (n_elements, 4), got {elements.shape}"
            )
        if elements.shape[0] == 0:
            raise MeshTopologyError("mesh contains no elements")
        if not np.issubdtype(elements.dtype, np.integer):
            if not np.all(np.equal(np.mod(elements, 1), 0)):
                raise MeshTopologyError("element connectivity must contain integer node ids")
        elements = elements.astype(np.int64)

        self.nodes: np.ndarray = nodes
        self.elements: np.ndarray = elements
        if validate:
            self.validate()

    # ---------------------------------------------------------------- sizes
    @property
    def n_nodes(self) -> int:
        return self.nodes.shape[0]

    @property
    def n_elements(self) -> int:
        return self.elements.shape[0]

    @property
    def n_dofs(self) -> int:
        return self.DOFS_PER_NODE * self.n_nodes

    # ------------------------------------------------------------- DOF map
    @staticmethod
    def node_dofs(node: int) -> Tuple[int, int]:
        """Global (x, y) DOF indices of a node."""
        node = int(node)
        return 2 * node, 2 * node + 1

    @staticmethod
    def dof_index(node: int, component: str) -> int:
        component = component.lower()
        if component not in ("x", "y"):
            raise ValueError(f"DOF component must be 'x' or 'y', got {component!r}")
        return 2 * int(node) + (0 if component == "x" else 1)

    @property
    def element_dofs(self) -> np.ndarray:
        """(n_elements, 8) global DOF indices in [u1, v1, ..., u4, v4] order."""
        e = self.elements
        dofs = np.empty((e.shape[0], 8), dtype=np.int64)
        dofs[:, 0::2] = 2 * e
        dofs[:, 1::2] = 2 * e + 1
        return dofs

    @property
    def element_coords(self) -> np.ndarray:
        """(n_elements, 4, 2) nodal coordinates per element."""
        return self.nodes[self.elements]

    # ----------------------------------------------------------- geometry
    def signed_areas(self) -> np.ndarray:
        """Shoelace signed area of each element (positive for CCW numbering)."""
        x = self.element_coords[..., 0]
        y = self.element_coords[..., 1]
        return 0.5 * np.sum(x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y, axis=1)

    def corner_cross_products(self) -> np.ndarray:
        """
        z-component of (x_{i+1} - x_i) x (x_{i-1} - x_i) at every corner.
        All positive <=> element is strictly convex and CCW, which in turn
        guarantees det(J) > 0 everywhere inside the element.
        """
        c = self.element_coords
        fwd = np.roll(c, -1, axis=1) - c
        bwd = np.roll(c, 1, axis=1) - c
        return fwd[..., 0] * bwd[..., 1] - fwd[..., 1] * bwd[..., 0]

    def bounding_box(self) -> Tuple[float, float, float, float]:
        xmin, ymin = self.nodes.min(axis=0)
        xmax, ymax = self.nodes.max(axis=0)
        return float(xmin), float(xmax), float(ymin), float(ymax)

    def characteristic_length(self) -> float:
        xmin, xmax, ymin, ymax = self.bounding_box()
        return float(np.hypot(xmax - xmin, ymax - ymin))

    def validate(self) -> None:
        e = self.elements
        if e.min() < 0 or e.max() >= self.n_nodes:
            raise MeshTopologyError(
                f"connectivity references node ids outside [0, {self.n_nodes - 1}]"
            )
        sorted_rows = np.sort(e, axis=1)
        repeated = np.any(sorted_rows[:, 1:] == sorted_rows[:, :-1], axis=1)
        if np.any(repeated):
            bad = np.flatnonzero(repeated)
            raise MeshTopologyError(
                f"elements {bad[:10].tolist()} reference the same node more than once"
            )
        areas = self.signed_areas()
        scale = self.characteristic_length() ** 2
        clockwise = areas <= 1e-14 * scale
        if np.any(clockwise):
            bad = np.flatnonzero(clockwise)
            raise MeshTopologyError(
                f"{bad.size} element(s) are clockwise or degenerate (signed area <= 0); "
                f"first offenders: {bad[:10].tolist()} with areas {areas[bad[:5]].tolist()}. "
                "Connectivity must be numbered counter-clockwise."
            )
        nonconvex = np.any(self.corner_cross_products() <= 0.0, axis=1)
        if np.any(nonconvex):
            logger.warning(
                "%d non-convex element(s) detected (first: %s); Jacobians will be "
                "checked at every Gauss point.",
                int(nonconvex.sum()),
                np.flatnonzero(nonconvex)[:10].tolist(),
            )
        used = np.zeros(self.n_nodes, dtype=bool)
        used[e.ravel()] = True
        if not np.all(used):
            logger.warning(
                "%d node(s) are not attached to any element; their DOFs will be "
                "constrained automatically at solve time.",
                int((~used).sum()),
            )

    # ------------------------------------------------------- node queries
    def nodes_where(self, predicate: Callable[[np.ndarray, np.ndarray], np.ndarray]) -> np.ndarray:
        """Node ids for which predicate(x, y) is True (vectorised)."""
        mask = np.asarray(predicate(self.nodes[:, 0], self.nodes[:, 1]), dtype=bool)
        return np.flatnonzero(mask)

    def nodes_at(
        self, x: Optional[float] = None, y: Optional[float] = None, tol: Optional[float] = None
    ) -> np.ndarray:
        """Node ids lying on the line x = const and/or y = const."""
        if x is None and y is None:
            raise ValueError("specify x and/or y")
        tol = 1e-9 * self.characteristic_length() if tol is None else tol
        mask = np.ones(self.n_nodes, dtype=bool)
        if x is not None:
            mask &= np.abs(self.nodes[:, 0] - x) <= tol
        if y is not None:
            mask &= np.abs(self.nodes[:, 1] - y) <= tol
        return np.flatnonzero(mask)

    def closest_node(self, x: float, y: float) -> int:
        return int(np.argmin(np.hypot(self.nodes[:, 0] - x, self.nodes[:, 1] - y)))

    # ---------------------------------------------------------- factories
    @classmethod
    def structured_rectangle(
        cls,
        length: float,
        height: float,
        nx: int,
        ny: int,
        origin: Tuple[float, float] = (0.0, 0.0),
    ) -> "Mesh":
        """Uniform nx-by-ny Q4 grid on [x0, x0+length] x [y0, y0+height]."""
        if nx < 1 or ny < 1:
            raise ValueError("nx and ny must be >= 1")
        if length <= 0 or height <= 0:
            raise ValueError("length and height must be positive")
        xs = origin[0] + np.linspace(0.0, length, nx + 1)
        ys = origin[1] + np.linspace(0.0, height, ny + 1)
        X, Y = np.meshgrid(xs, ys)  # row j -> y index, column i -> x index
        nodes = np.column_stack([X.ravel(), Y.ravel()])
        i, j = np.meshgrid(np.arange(nx), np.arange(ny))
        n0 = (j * (nx + 1) + i).ravel()
        elements = np.column_stack([n0, n0 + 1, n0 + nx + 2, n0 + nx + 1])
        return cls(nodes, elements)

    @classmethod
    def plate_with_hole_quarter(
        cls,
        radius: float,
        half_width: float,
        n_radial: int,
        n_circumferential: int,
        grading: float = 1.0,
    ) -> "Mesh":
        """
        Mapped (curvilinear, non-rectangular) Q4 mesh of the quarter domain
        [0, W] x [0, W] minus a circular hole of radius a centred at the origin.

        Rays run from the hole boundary (a cos(th), a sin(th)) to the outer
        square boundary; radial spacing is geometrically graded towards the
        hole by `grading` (> 1 clusters elements at the hole).
        """
        if not (0.0 < radius < half_width):
            raise ValueError("require 0 < radius < half_width")
        if n_circumferential % 2:
            raise ValueError("n_circumferential must be even so a ray hits the 45 deg corner")
        if n_radial < 1:
            raise ValueError("n_radial must be >= 1")
        theta = np.linspace(0.0, 0.5 * np.pi, n_circumferential + 1)
        inner = radius * np.column_stack([np.cos(theta), np.sin(theta)])
        outer = np.empty_like(inner)
        for k, th in enumerate(theta):
            if th <= 0.25 * np.pi + 1e-14:
                outer[k] = (half_width, half_width * np.tan(th))
            else:
                outer[k] = (half_width / np.tan(th), half_width)
        if abs(grading - 1.0) < 1e-12:
            s = np.linspace(0.0, 1.0, n_radial + 1)
        else:
            q = grading ** (1.0 / max(n_radial - 1, 1))
            steps = q ** np.arange(n_radial)
            s = np.concatenate([[0.0], np.cumsum(steps)]) / steps.sum()
        # node id = k * (n_radial + 1) + i, k: angular index, i: radial index
        nodes = (inner[:, None, :] * (1.0 - s)[None, :, None]) + (
            outer[:, None, :] * s[None, :, None]
        )
        nodes = nodes.reshape(-1, 2)
        nr = n_radial + 1
        i, k = np.meshgrid(np.arange(n_radial), np.arange(n_circumferential))
        n0 = (k * nr + i).ravel()
        # radial-outward then angular-positive is right handed => CCW
        elements = np.column_stack([n0, n0 + 1, n0 + nr + 1, n0 + nr])
        return cls(nodes, elements)

    def __repr__(self) -> str:
        return f"Mesh(n_nodes={self.n_nodes}, n_elements={self.n_elements}, n_dofs={self.n_dofs})"


# --------------------------------------------------------------------------- #
# Q4 element formulation                                                      #
# --------------------------------------------------------------------------- #
class Q4Element:
    """
    Bilinear 4-node isoparametric quadrilateral.

    All routines operate on a batch of elements: `coords` has shape
    (n_elements, 4, 2). A single element is simply a batch of size one.
    """

    NODE_NATURAL = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])
    GAUSS_ABSCISSA = 1.0 / np.sqrt(3.0)
    # Gauss points ordered like the element nodes so stresses can be
    # extrapolated with the same bilinear interpolation.
    GAUSS_POINTS = GAUSS_ABSCISSA * NODE_NATURAL
    GAUSS_WEIGHTS = np.ones(4)

    # ---------------------------------------------------- shape functions
    @staticmethod
    def shape_functions(xi: float, eta: float) -> np.ndarray:
        """N_i(xi, eta) = 1/4 (1 + xi xi_i)(1 + eta eta_i), shape (4,)."""
        return 0.25 * np.array(
            [
                (1.0 - xi) * (1.0 - eta),
                (1.0 + xi) * (1.0 - eta),
                (1.0 + xi) * (1.0 + eta),
                (1.0 - xi) * (1.0 + eta),
            ]
        )

    @staticmethod
    def shape_derivatives(xi: float, eta: float) -> np.ndarray:
        """
        Natural derivatives, shape (2, 4):
            row 0 = dN_i/dxi,  row 1 = dN_i/deta.
        """
        return 0.25 * np.array(
            [
                [-(1.0 - eta), (1.0 - eta), (1.0 + eta), -(1.0 + eta)],
                [-(1.0 - xi), -(1.0 + xi), (1.0 + xi), (1.0 - xi)],
            ]
        )

    # -------------------------------------------------------- kinematics
    @classmethod
    def jacobian(cls, coords: np.ndarray, xi: float, eta: float) -> np.ndarray:
        """
        J = dN_nat (2x4) @ X (4x2) = [[dx/dxi, dy/dxi], [dx/deta, dy/deta]],
        returned with shape (n_elements, 2, 2).
        """
        coords = np.asarray(coords, dtype=np.float64)
        if coords.ndim == 2:
            coords = coords[None]
        return np.einsum("an,enc->eac", cls.shape_derivatives(xi, eta), coords)

    @classmethod
    def kinematics(
        cls,
        coords: np.ndarray,
        xi: float,
        eta: float,
        element_ids: Optional[np.ndarray] = None,
        rel_tol: float = 1e-10,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Evaluate det(J) and the 3x8 B matrix at (xi, eta) for every element.

        Raises DistortedElementError if any det(J) <= rel_tol * reference,
        where the reference is a quarter of the element's bounding-box area
        (det(J) of an undistorted rectangle of the same extent).

        Returns
        -------
        detJ : (n_elements,)
        B    : (n_elements, 3, 8)
        """
        coords = np.asarray(coords, dtype=np.float64)
        if coords.ndim == 2:
            coords = coords[None]
        dN_nat = cls.shape_derivatives(xi, eta)
        J = np.einsum("an,enc->eac", dN_nat, coords)
        detJ = J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0]

        extent = coords.max(axis=1) - coords.min(axis=1)
        reference = 0.25 * extent[:, 0] * extent[:, 1]
        reference = np.where(reference > 0.0, reference, 0.25 * np.sum(extent**2, axis=1))
        bad = detJ <= rel_tol * reference
        if np.any(bad):
            idx = np.flatnonzero(bad)
            ids = idx if element_ids is None else np.asarray(element_ids)[idx]
            raise DistortedElementError(
                f"Non-positive or near-zero Jacobian determinant at natural point "
                f"(xi={xi:+.6f}, eta={eta:+.6f}) in {idx.size} element(s); "
                f"element ids {ids[:10].tolist()} have det(J) = "
                f"{detJ[idx[:10]].tolist()}. The element is inverted, collapsed or "
                f"excessively distorted (check CCW numbering and interior angles < 180 deg)."
            )

        inv_det = 1.0 / detJ
        invJ = np.empty_like(J)
        invJ[:, 0, 0] = J[:, 1, 1] * inv_det
        invJ[:, 0, 1] = -J[:, 0, 1] * inv_det
        invJ[:, 1, 0] = -J[:, 1, 0] * inv_det
        invJ[:, 1, 1] = J[:, 0, 0] * inv_det

        # [dN/dx; dN/dy] = J^{-1} [dN/dxi; dN/deta]
        dN_glob = np.einsum("eab,bn->ean", invJ, dN_nat)

        B = np.zeros((coords.shape[0], 3, 8))
        B[:, 0, 0::2] = dN_glob[:, 0, :]
        B[:, 1, 1::2] = dN_glob[:, 1, :]
        B[:, 2, 0::2] = dN_glob[:, 1, :]
        B[:, 2, 1::2] = dN_glob[:, 0, :]
        return detJ, B

    # ------------------------------------------------------- stiffness
    @classmethod
    def stiffness(
        cls,
        coords: np.ndarray,
        D: np.ndarray,
        thickness: float,
        element_ids: Optional[np.ndarray] = None,
        verify: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        K_e = sum_{g=1..4} w_g * B_g^T D B_g * det(J_g) * t   (2x2 Gauss, w_g = 1).

        Returns
        -------
        Ke   : (n_elements, 8, 8)
        detJ : (n_elements, 4) Jacobian determinants at the Gauss points.
        """
        coords = np.asarray(coords, dtype=np.float64)
        if coords.ndim == 2:
            coords = coords[None]
        n_el = coords.shape[0]
        Ke = np.zeros((n_el, 8, 8))
        detJ_all = np.empty((n_el, 4))
        for g, ((xi, eta), w) in enumerate(zip(cls.GAUSS_POINTS, cls.GAUSS_WEIGHTS)):
            detJ, B = cls.kinematics(coords, xi, eta, element_ids=element_ids)
            detJ_all[:, g] = detJ
            Ke += np.einsum("eki,kl,elj->eij", B, D, B) * (w * thickness * detJ)[:, None, None]

        if verify:
            cls.verify_stiffness(Ke, element_ids)
        # remove round-off asymmetry so downstream solvers see an exactly symmetric K
        Ke = 0.5 * (Ke + np.transpose(Ke, (0, 2, 1)))
        return Ke, detJ_all

    @staticmethod
    def verify_stiffness(
        Ke: np.ndarray, element_ids: Optional[np.ndarray] = None, rel_tol: float = 1e-9
    ) -> None:
        """Check symmetry and positive semi-definiteness of every element matrix."""
        scale = np.max(np.abs(Ke), axis=(1, 2))
        asym = np.max(np.abs(Ke - np.transpose(Ke, (0, 2, 1))), axis=(1, 2))
        bad = asym > rel_tol * scale
        if np.any(bad):
            idx = np.flatnonzero(bad)
            ids = idx if element_ids is None else np.asarray(element_ids)[idx]
            raise ElementIntegrityError(
                f"element stiffness not symmetric for elements {ids[:10].tolist()} "
                f"(max |K - K^T| / max|K| = {(asym[idx] / scale[idx]).max():.3e})"
            )
        sym = 0.5 * (Ke + np.transpose(Ke, (0, 2, 1)))
        eig_min = np.linalg.eigvalsh(sym)[:, 0]
        bad = eig_min < -rel_tol * scale
        if np.any(bad):
            idx = np.flatnonzero(bad)
            ids = idx if element_ids is None else np.asarray(element_ids)[idx]
            raise ElementIntegrityError(
                f"element stiffness not positive semi-definite for elements "
                f"{ids[:10].tolist()} (min eigenvalue {eig_min[idx].min():.3e})"
            )

    # ---------------------------------------------------- extrapolation
    @classmethod
    def extrapolation_matrix(cls) -> np.ndarray:
        """
        4x4 matrix E with sigma_node = E @ sigma_gauss.

        The Gauss points form a smaller 'element' with natural coordinates
        r = sqrt(3) xi; the element corners sit at r = +/- sqrt(3) in that
        frame, so E[i, k] = N_k(sqrt(3) xi_i, sqrt(3) eta_i).
        """
        s3 = np.sqrt(3.0)
        return np.array([cls.shape_functions(s3 * xi, s3 * eta) for xi, eta in cls.NODE_NATURAL])


# --------------------------------------------------------------------------- #
# Results container                                                           #
# --------------------------------------------------------------------------- #
@dataclass
class FEAResults:
    displacements: np.ndarray        # (n_nodes, 2)
    u: np.ndarray                    # (n_dofs,)
    reactions: np.ndarray            # (n_dofs,) non-zero only at constrained DOFs
    external_forces: np.ndarray      # (n_dofs,)
    gauss_points: np.ndarray         # (n_elements, 4, 2) physical Gauss point coords
    gauss_strains: np.ndarray        # (n_elements, 4, 3)
    gauss_stresses: np.ndarray       # (n_elements, 4, 3)
    element_nodal_stresses: np.ndarray  # (n_elements, 4, 3) extrapolated, unaveraged
    nodal_stresses: np.ndarray       # (n_nodes, 3) averaged [s_xx, s_yy, t_xy]
    nodal_szz: np.ndarray            # (n_nodes,) out-of-plane stress
    von_mises: np.ndarray            # (n_nodes,)
    gauss_von_mises: np.ndarray      # (n_elements, 4)
    safety_factor: Optional[np.ndarray]
    diagnostics: Dict[str, float] = field(default_factory=dict)

    @property
    def max_displacement(self) -> float:
        return float(np.max(np.linalg.norm(self.displacements, axis=1)))

    @property
    def max_von_mises(self) -> float:
        return float(np.max(self.von_mises))


# --------------------------------------------------------------------------- #
# Model                                                                       #
# --------------------------------------------------------------------------- #
class FEAModel:
    """
    Assembles and solves K u = F for a Q4 mesh with one material.

    Boundary conditions
    -------------------
    * fix(...)                : Dirichlet (prescribed displacement, default 0)
    * add_nodal_force(...)    : Neumann, concentrated nodal force
    * add_edge_traction(...)  : Neumann, distributed traction over boundary
                                edges (consistent nodal loads, 2-point Gauss)
    """

    ITERATIVE_THRESHOLD = 50_000

    def __init__(
        self,
        mesh: Mesh,
        material: Material,
        iterative_threshold: int = ITERATIVE_THRESHOLD,
        cg_rtol: float = 1e-10,
        cg_maxiter: Optional[int] = None,
        direct_fallback: bool = True,
    ) -> None:
        self.mesh = mesh
        self.material = material
        self.iterative_threshold = int(iterative_threshold)
        self.cg_rtol = float(cg_rtol)
        self.cg_maxiter = cg_maxiter
        self.direct_fallback = bool(direct_fallback)

        self.F = np.zeros(mesh.n_dofs)
        self._prescribed: Dict[int, float] = {}

        self.K: Optional[sp.csr_matrix] = None
        self.Ke: Optional[np.ndarray] = None
        self.detJ: Optional[np.ndarray] = None
        self.results: Optional[FEAResults] = None

    # ================================================= boundary conditions
    def fix(
        self,
        nodes: Union[int, Iterable[int]],
        dofs: str = "xy",
        value: Union[float, Tuple[float, float]] = 0.0,
    ) -> None:
        """Prescribe displacement on `dofs` ('x', 'y' or 'xy') of the given nodes."""
        dofs = dofs.lower()
        if not dofs or any(c not in "xy" for c in dofs):
            raise ValueError(f"dofs must be a combination of 'x' and 'y', got {dofs!r}")
        if np.isscalar(value):
            values = {"x": float(value), "y": float(value)}
        else:
            vx, vy = value
            values = {"x": float(vx), "y": float(vy)}
        for n in np.atleast_1d(np.asarray(nodes, dtype=np.int64)):
            self._check_node(n)
            for c in dofs:
                self._prescribed[Mesh.dof_index(n, c)] = values[c]
        self.results = None

    def add_nodal_force(self, node: int, fx: float = 0.0, fy: float = 0.0) -> None:
        """Concentrated force (force units) at a node; accumulates."""
        self._check_node(node)
        ix, iy = Mesh.node_dofs(node)
        self.F[ix] += float(fx)
        self.F[iy] += float(fy)
        self.results = None

    def add_edge_traction(
        self,
        edge_nodes: Sequence[int],
        traction: Union[
            Tuple[float, float],
            Callable[[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]],
        ],
    ) -> float:
        """
        Distributed traction (force / area) on a chain of boundary nodes.

        `edge_nodes` is an ordered polyline of nodes n0-n1-...-nk; every
        segment is a straight element edge. `traction` is either a constant
        (tx, ty) or a callable (x, y) -> (tx, ty) evaluated at the Gauss points
        of each segment. The consistent nodal load of a segment of length L is

            f_a = t * L/2 * sum_g w_g N_a(s_g) T(x(s_g)),   a = 1, 2

        which is exact for tractions up to cubic order along the segment.
        Returns the resultant force vector magnitude added.
        """
        chain = np.asarray(edge_nodes, dtype=np.int64)
        if chain.ndim != 1 or chain.size < 2:
            raise ValueError("edge_nodes must be an ordered sequence of at least two nodes")
        for n in chain:
            self._check_node(n)
        t = self.material.thickness
        gp = np.array([-1.0, 1.0]) / np.sqrt(3.0)
        total = np.zeros(2)
        for a, b in zip(chain[:-1], chain[1:]):
            xa, xb = self.mesh.nodes[a], self.mesh.nodes[b]
            length = float(np.hypot(*(xb - xa)))
            if length <= 0.0:
                raise ValueError(f"zero-length traction segment between nodes {a} and {b}")
            fa = np.zeros(2)
            fb = np.zeros(2)
            for s in gp:
                Na, Nb = 0.5 * (1.0 - s), 0.5 * (1.0 + s)
                xg = Na * xa + Nb * xb
                if callable(traction):
                    tx, ty = traction(np.array(xg[0]), np.array(xg[1]))
                    T = np.array([float(tx), float(ty)])
                else:
                    T = np.asarray(traction, dtype=np.float64)
                jac = 0.5 * length  # ds/dxi, weight 1
                fa += Na * T * jac * t
                fb += Nb * T * jac * t
            self.F[2 * a : 2 * a + 2] += fa
            self.F[2 * b : 2 * b + 2] += fb
            total += fa + fb
        self.results = None
        return float(np.hypot(*total))

    def clear_loads(self) -> None:
        self.F[:] = 0.0
        self.results = None

    def clear_constraints(self) -> None:
        self._prescribed.clear()
        self.results = None

    @property
    def constrained_dofs(self) -> np.ndarray:
        return np.array(sorted(self._prescribed), dtype=np.int64)

    def _check_node(self, node: int) -> None:
        if not (0 <= int(node) < self.mesh.n_nodes):
            raise IndexError(f"node {node} out of range [0, {self.mesh.n_nodes - 1}]")

    # =========================================================== assembly
    def assemble(self) -> sp.csr_matrix:
        """Integrate every element and assemble the global sparse stiffness K."""
        mesh = self.mesh
        t0 = time.perf_counter()
        D = self.material.D
        logger.info("[STEP 2] Constitutive matrix D (plane %s):\n%s",
                    self.material.plane, np.array2string(D, precision=6, max_line_width=120))

        Ke, detJ = Q4Element.stiffness(
            mesh.element_coords, D, self.material.thickness,
            element_ids=np.arange(mesh.n_elements),
        )
        t1 = time.perf_counter()
        area_quadrature = float(np.sum(detJ))  # sum_g w_g detJ_g = element area
        logger.info(
            "[STEP 3] Q4 kinematics + 2x2 Gauss integration on %d elements "
            "(%d Gauss points, xi,eta = +/-1/sqrt(3) = +/-%.10f, w = 1.0)",
            mesh.n_elements, 4 * mesh.n_elements, Q4Element.GAUSS_ABSCISSA,
        )
        logger.info(
            "         det(J): min = %.6e, max = %.6e, all positive = %s",
            detJ.min(), detJ.max(), bool(np.all(detJ > 0)),
        )
        logger.info(
            "         domain area via sum(w*detJ) = %.10e, via shoelace = %.10e",
            area_quadrature, float(np.sum(mesh.signed_areas())),
        )
        eig = np.linalg.eigvalsh(Ke[0])
        logger.info(
            "         Ke[0] eigenvalues (3 rigid-body zeros expected): %s",
            np.array2string(eig, precision=4, max_line_width=140),
        )
        logger.info("         element integration time: %.4f s", t1 - t0)

        dofs = mesh.element_dofs
        rows = np.repeat(dofs, 8, axis=1).ravel()
        cols = np.tile(dofs, (1, 8)).ravel()
        K = sp.coo_matrix(
            (Ke.ravel(), (rows, cols)), shape=(mesh.n_dofs, mesh.n_dofs)
        ).tocsr()
        K.sum_duplicates()
        K.eliminate_zeros()
        t2 = time.perf_counter()

        asym = abs(K - K.T).max() if K.nnz else 0.0
        density = K.nnz / float(mesh.n_dofs) ** 2
        logger.info(
            "[STEP 4] Global sparse assembly (COO -> CSR): %d x %d, nnz = %d "
            "(density %.4e, %.2f MB vs %.2f MB dense), max|K-K^T| = %.3e, %.4f s",
            mesh.n_dofs, mesh.n_dofs, K.nnz, density,
            (K.data.nbytes + K.indices.nbytes + K.indptr.nbytes) / 1e6,
            8.0 * mesh.n_dofs**2 / 1e6, asym, t2 - t1,
        )
        self.K, self.Ke, self.detJ = K, Ke, detJ
        return K

    # ============================================================== solve
    def solve(self) -> FEAResults:
        mesh = self.mesh
        t_start = time.perf_counter()
        logger.info("=" * 78)
        logger.info("FEA SOLVE  |  %r  |  %r", mesh, self.material)
        logger.info("=" * 78)
        logger.info(
            "[STEP 1] Mesh: %d nodes, %d Q4 elements, %d DOFs (2 per node: u_x=2n, u_y=2n+1)",
            mesh.n_nodes, mesh.n_elements, mesh.n_dofs,
        )
        if self.K is None:
            self.assemble()
        K = self.K

        # unattached nodes would leave zero rows -> constrain them
        attached = np.zeros(mesh.n_nodes, dtype=bool)
        attached[mesh.elements.ravel()] = True
        prescribed = dict(self._prescribed)
        for n in np.flatnonzero(~attached):
            prescribed.setdefault(2 * n, 0.0)
            prescribed.setdefault(2 * n + 1, 0.0)

        if not prescribed:
            raise SingularSystemError(
                "no Dirichlet constraints defined: the structure has 3 unrestrained "
                "rigid-body modes and K is singular"
            )
        c = np.array(sorted(prescribed), dtype=np.int64)
        u_c = np.array([prescribed[i] for i in c])
        self._check_rigid_body_restraint(c)
        free_mask = np.ones(mesh.n_dofs, dtype=bool)
        free_mask[c] = False
        f = np.flatnonzero(free_mask)
        n_free = f.size
        if n_free == 0:
            raise SingularSystemError("every DOF is constrained; nothing to solve")

        # ---- partition  [K_ff K_fc; K_cf K_cc] [u_f; u_c] = [F_f; F_c + R_c]
        K_ff = K[f][:, f].tocsc()
        K_fc = K[f][:, c]
        rhs = self.F[f] - K_fc @ u_c
        logger.info(
            "[STEP 5] Partitioning: %d free DOFs, %d constrained DOFs "
            "(%d with non-zero prescribed value); |F_ext| = %.6e",
            n_free, c.size, int(np.count_nonzero(u_c)), np.linalg.norm(self.F),
        )
        logger.info("         K_ff: %d x %d, nnz = %d", n_free, n_free, K_ff.nnz)

        t_solve = time.perf_counter()
        if n_free > self.iterative_threshold:
            u_f, solver_name, iterations = self._solve_iterative(K_ff, rhs)
        else:
            u_f, solver_name, iterations = self._solve_direct(K_ff, rhs)
        t_solve = time.perf_counter() - t_solve

        u = np.zeros(mesh.n_dofs)
        u[f] = u_f
        u[c] = u_c

        residual = K_ff @ u_f - rhs
        rel_res = np.linalg.norm(residual) / max(np.linalg.norm(rhs), np.finfo(float).tiny)
        logger.info(
            "[STEP 6] Linear solve: %s%s, %.4f s, relative residual ||K_ff u_f - b||/||b|| = %.3e",
            solver_name, "" if iterations is None else f" ({iterations} iterations)",
            t_solve, rel_res,
        )

        # ---- reactions and global equilibrium
        Ku = K @ u
        reactions = np.zeros(mesh.n_dofs)
        reactions[c] = Ku[c] - self.F[c]
        total = self.F + reactions
        x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
        sum_fx, sum_fy = total[0::2].sum(), total[1::2].sum()
        moment = np.sum(x * total[1::2] - y * total[0::2])
        load_scale = max(np.abs(self.F).sum(), np.abs(reactions).sum(), np.finfo(float).tiny)
        moment_scale = load_scale * mesh.characteristic_length()
        logger.info(
            "[STEP 7] Reactions R = K u - F at constrained DOFs: "
            "sum R_x = %.6e, sum R_y = %.6e",
            reactions[0::2].sum(), reactions[1::2].sum(),
        )
        logger.info(
            "         Static equilibrium (applied + reactions): "
            "sum F_x = %.3e, sum F_y = %.3e, sum M_z = %.3e  (relative %.3e)",
            sum_fx, sum_fy, moment,
            max(abs(sum_fx), abs(sum_fy)) / load_scale + abs(moment) / moment_scale,
        )
        if not np.all(np.isfinite(u)):
            raise SingularSystemError("solution contains NaN/Inf; system is singular")

        # ---- stress recovery
        results = self._recover(u, reactions)
        results.diagnostics.update(
            {
                "n_nodes": mesh.n_nodes,
                "n_elements": mesh.n_elements,
                "n_dofs": mesh.n_dofs,
                "n_free": n_free,
                "n_constrained": int(c.size),
                "nnz": int(K.nnz),
                "solver": solver_name,
                "iterations": -1 if iterations is None else iterations,
                "relative_residual": float(rel_res),
                "equilibrium_fx": float(sum_fx),
                "equilibrium_fy": float(sum_fy),
                "equilibrium_mz": float(moment),
                "detJ_min": float(self.detJ.min()),
                "detJ_max": float(self.detJ.max()),
                "solve_time": float(t_solve),
                "total_time": float(time.perf_counter() - t_start),
            }
        )
        logger.info(
            "[DONE]   max |u| = %.6e, max Von Mises = %.6e, total wall time %.4f s",
            results.max_displacement, results.max_von_mises, results.diagnostics["total_time"],
        )
        self.results = results
        return results

    def _check_rigid_body_restraint(self, constrained: np.ndarray) -> None:
        """
        Exact singularity test, independent of round-off.

        A fully integrated Q4 mesh has no spurious zero-energy modes, so K_ff is
        singular if and only if some connected group of elements has a
        rigid-body motion (u_x, u_y) = (a - b y, c + b x) that does not move
        any of its constrained DOFs. For each connected component the 3 rigid
        modes, sampled at the constrained DOFs, must therefore have rank 3.
        """
        mesh = self.mesh
        e = mesh.elements
        pairs = np.vstack([e[:, [i, (i + 1) % 4]] for i in range(4)])
        adjacency = sp.coo_matrix(
            (np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(mesh.n_nodes,) * 2
        ).tocsr()
        n_comp, labels = sp.csgraph.connected_components(adjacency, directed=False)

        is_constrained = np.zeros(mesh.n_dofs, dtype=bool)
        is_constrained[constrained] = True
        attached = np.zeros(mesh.n_nodes, dtype=bool)
        attached[e.ravel()] = True
        for comp in range(n_comp):
            comp_nodes = np.flatnonzero((labels == comp) & attached)
            if comp_nodes.size == 0:
                continue
            xy = mesh.nodes[comp_nodes]
            xy = (xy - xy.mean(axis=0)) / max(np.ptp(xy, axis=0).max(), np.finfo(float).tiny)
            rows = []
            for k, n in enumerate(comp_nodes):
                if is_constrained[2 * n]:
                    rows.append([1.0, 0.0, -xy[k, 1]])
                if is_constrained[2 * n + 1]:
                    rows.append([0.0, 1.0, xy[k, 0]])
            rank = np.linalg.matrix_rank(np.array(rows), tol=1e-10) if rows else 0
            if rank < 3:
                names = ["x-translation", "y-translation", "rotation"]
                free_modes = 3 - rank
                raise SingularSystemError(
                    f"K_ff is singular: mesh component {comp} ({comp_nodes.size} nodes, "
                    f"{len(rows)} constrained DOFs) has {free_modes} unrestrained rigid-body "
                    f"mode(s) among {names}. Add Dirichlet constraints that prevent "
                    "translation in x and y and rotation."
                )
        logger.info(
            "         Rigid-body restraint verified on %d connected component(s) "
            "(rank of constrained rigid-mode matrix = 3)", n_comp,
        )

    def _solve_direct(
        self, A: sp.csc_matrix, b: np.ndarray, pivot_tol: float = 1e-11
    ) -> Tuple[np.ndarray, str, None]:
        """
        Sparse LU (SuperLU) factorisation and solve. A rigid-body mode left
        unrestrained shows up as a pivot that is round-off small relative to
        the largest pivot, rather than an exact zero, so pivots are screened
        against `pivot_tol` before the solution is accepted.
        """
        singular_msg = (
            "reduced stiffness matrix K_ff is singular: the constraints do not suppress "
            "all rigid-body modes (2 translations + 1 rotation), or the mesh contains a mechanism"
        )
        with warnings.catch_warnings():
            warnings.simplefilter("error", spla.MatrixRankWarning)
            try:
                lu = spla.splu(A)
            except (spla.MatrixRankWarning, RuntimeError) as exc:
                raise SingularSystemError(f"{singular_msg}. Underlying error: {exc}") from exc
        pivots = np.abs(lu.U.diagonal())
        ratio = pivots.min() / pivots.max()
        logger.info("         LU pivot ratio min|U_ii| / max|U_ii| = %.3e", ratio)
        if not np.isfinite(ratio) or ratio == 0.0:
            raise SingularSystemError(f"{singular_msg} (zero pivot)")
        if ratio < pivot_tol:
            logger.warning(
                "small LU pivot ratio %.3e: K_ff is ill-conditioned (very slender "
                "geometry or extreme stiffness contrast); verify the results", ratio,
            )
        x = lu.solve(b)
        if not np.all(np.isfinite(x)):
            raise SingularSystemError("direct solve produced NaN/Inf: K_ff is singular")
        return np.atleast_1d(x), "SciPy sparse direct (SuperLU)", None

    def _solve_iterative(self, A: sp.csc_matrix, b: np.ndarray) -> Tuple[np.ndarray, str, int]:
        """
        Preconditioned conjugate gradient for large systems.

        1.  Symmetric Jacobi equilibration  A_s = S A S,  S = diag(A)^(-1/2),
            which keeps A_s SPD with a unit diagonal and removes the stiffness
            scale spread between DOFs.
        2.  Incomplete LU of A_s with a symmetric fill-reducing ordering
            (minimum degree on A^T + A) and diagonal pivoting, so the factors
            stay close to an incomplete Cholesky and CG remains effective.
        3.  CG on A_s y = S b, then u = S y.

        If ILU fails the preconditioner degrades to the (already applied)
        Jacobi scaling; if CG does not reach the tolerance the system is handed
        to the sparse direct solver (or SolverConvergenceError is raised when
        direct_fallback=False) so a result is never silently inaccurate.
        """
        n = A.shape[0]
        diag = A.diagonal()
        if np.any(diag <= 0.0):
            raise SingularSystemError("non-positive diagonal in K_ff; system is singular")
        s = 1.0 / np.sqrt(diag)
        S = sp.diags(s)
        A_s = (S @ A @ S).tocsc()
        b_s = s * b

        precond_name = "Jacobi-scaled ILU"
        t0 = time.perf_counter()
        try:
            ilu = spla.spilu(
                A_s, drop_tol=1e-6, fill_factor=30.0,
                permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
            )
            M = spla.LinearOperator((n, n), matvec=ilu.solve, dtype=np.float64)
            logger.info(
                "         ILU preconditioner: nnz(L+U) = %d (%.2fx nnz(K_ff)), %.3f s",
                ilu.nnz, ilu.nnz / A.nnz, time.perf_counter() - t0,
            )
        except RuntimeError as exc:
            logger.warning("ILU factorisation failed (%s); using Jacobi preconditioning only", exc)
            M = None
            precond_name = "Jacobi"

        iterations = 0

        def _count(_xk: np.ndarray) -> None:
            nonlocal iterations
            iterations += 1

        maxiter = self.cg_maxiter if self.cg_maxiter is not None else max(n // 10, 2000)
        y, info = spla.cg(A_s, b_s, rtol=self.cg_rtol, atol=0.0, maxiter=maxiter, M=M,
                          callback=_count)
        if info < 0:
            raise SingularSystemError(f"CG breakdown (info={info}); K_ff is likely singular")
        if info > 0:
            if not self.direct_fallback:
                raise SolverConvergenceError(
                    f"preconditioned CG did not converge in {info} iterations "
                    f"(rtol={self.cg_rtol:g})"
                )
            logger.warning(
                "preconditioned CG did not converge in %d iterations (rtol=%g); "
                "falling back to the sparse direct solver", info, self.cg_rtol,
            )
            x, name, _ = self._solve_direct(A, b)
            return x, f"{name} [fallback after {iterations} CG iterations]", iterations
        return s * y, f"Preconditioned conjugate gradient ({precond_name})", iterations

    # =================================================== post-processing
    def _recover(self, u: np.ndarray, reactions: np.ndarray) -> FEAResults:
        mesh, mat = self.mesh, self.material
        D = mat.D
        coords = mesh.element_coords
        ue = u[mesh.element_dofs]  # (n_el, 8)

        n_el = mesh.n_elements
        strains = np.empty((n_el, 4, 3))
        gauss_xy = np.empty((n_el, 4, 2))
        for g, (xi, eta) in enumerate(Q4Element.GAUSS_POINTS):
            _, B = Q4Element.kinematics(coords, xi, eta, element_ids=np.arange(n_el))
            strains[:, g, :] = np.einsum("eij,ej->ei", B, ue)
            gauss_xy[:, g, :] = Q4Element.shape_functions(xi, eta) @ coords
        stresses = strains @ D.T  # sigma = D eps

        E = Q4Element.extrapolation_matrix()
        elem_nodal = np.einsum("ig,egc->eic", E, stresses)

        sums = np.zeros((mesh.n_nodes, 3))
        counts = np.zeros(mesh.n_nodes)
        np.add.at(sums, mesh.elements.ravel(), elem_nodal.reshape(-1, 3))
        np.add.at(counts, mesh.elements.ravel(), 1.0)
        nodal = sums / np.maximum(counts, 1.0)[:, None]

        szz_nodal = mat.out_of_plane_stress(nodal)
        vm = von_mises(nodal, szz_nodal)
        vm_gauss = von_mises(stresses, mat.out_of_plane_stress(stresses))

        safety = None
        if mat.yield_strength is not None:
            safety = np.where(vm > 0.0, mat.yield_strength / np.maximum(vm, 1e-300), np.inf)

        logger.info(
            "[STEP 8] Strain recovery eps = B u_e at %d Gauss points; sigma = D eps", 4 * n_el
        )
        logger.info(
            "         Gauss-point extremes: s_xx [%.4e, %.4e], s_yy [%.4e, %.4e], t_xy [%.4e, %.4e]",
            stresses[..., 0].min(), stresses[..., 0].max(),
            stresses[..., 1].min(), stresses[..., 1].max(),
            stresses[..., 2].min(), stresses[..., 2].max(),
        )
        logger.info(
            "[STEP 9] Bilinear Gauss->node extrapolation (E = N(sqrt(3) xi_i)), nodal averaging; "
            "Von Mises: nodal max %.6e, Gauss max %.6e",
            vm.max(), vm_gauss.max(),
        )
        if safety is not None:
            logger.info(
                "         Yield strength %.4e -> minimum safety factor %.4f at node %d",
                mat.yield_strength, safety.min(), int(np.argmin(safety)),
            )

        return FEAResults(
            displacements=u.reshape(-1, 2).copy(),
            u=u,
            reactions=reactions,
            external_forces=self.F.copy(),
            gauss_points=gauss_xy,
            gauss_strains=strains,
            gauss_stresses=stresses,
            element_nodal_stresses=elem_nodal,
            nodal_stresses=nodal,
            nodal_szz=szz_nodal,
            von_mises=vm,
            gauss_von_mises=vm_gauss,
            safety_factor=safety,
        )


# --------------------------------------------------------------------------- #
# Helpers                                                                     #
# --------------------------------------------------------------------------- #
def von_mises(stress: np.ndarray, szz: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Von Mises equivalent stress for Voigt arrays [..., (s_xx, s_yy, t_xy)]:

        s_vm = sqrt( 1/2 [(sx-sy)^2 + (sy-sz)^2 + (sz-sx)^2] + 3 t_xy^2 )
    """
    sx, sy, txy = stress[..., 0], stress[..., 1], stress[..., 2]
    sz = np.zeros_like(sx) if szz is None else szz
    return np.sqrt(
        0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2) + 3.0 * txy**2
    )


def cantilever_tip_deflection(P: float, L: float, E: float, I: float) -> float:
    """Euler-Bernoulli tip deflection of an end-loaded cantilever, delta = P L^3 / (3 E I)."""
    return P * L**3 / (3.0 * E * I)
