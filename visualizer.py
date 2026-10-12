"""
Matplotlib post-processing for the Q4 FEA engine.

Produces a two-panel engineering figure:
  left  : undeformed wireframe mesh with boundary-condition / load markers
  right : continuous Von Mises contour on the deformed geometry (scaled),
          with the undeformed outline ghosted underneath.
"""

from __future__ import annotations

import os
import sys
from typing import Optional, Tuple

import matplotlib

# Headless Linux (servers, CI, containers) has no screen: render to files only.
if (
    sys.platform.startswith("linux")
    and os.environ.get("MPLBACKEND") is None
    and not os.environ.get("DISPLAY")
    and not os.environ.get("WAYLAND_DISPLAY")
):
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.collections import PolyCollection

from engine import FEAModel, FEAResults, Mesh


def auto_scale_factor(mesh: Mesh, results: FEAResults, target_fraction: float = 0.10) -> float:
    """Scale so the largest displacement is drawn as `target_fraction` of the model size."""
    xmin, xmax, ymin, ymax = mesh.bounding_box()
    size = max(xmax - xmin, ymax - ymin)
    umax = results.max_displacement
    return 1.0 if umax <= 0.0 else target_fraction * size / umax


def quad_triangulation(mesh: Mesh, coords: np.ndarray) -> mtri.Triangulation:
    """Split every Q4 into two triangles (1-2-3, 1-3-4) for Gouraud contouring."""
    e = mesh.elements
    tris = np.vstack([e[:, [0, 1, 2]], e[:, [0, 2, 3]]])
    return mtri.Triangulation(coords[:, 0], coords[:, 1], tris)


def _wireframe(ax, mesh: Mesh, coords: np.ndarray, **kwargs) -> PolyCollection:
    polys = coords[mesh.elements]
    pc = PolyCollection(polys, **kwargs)
    ax.add_collection(pc)
    return pc


def _outline_edges(mesh: Mesh) -> np.ndarray:
    """Boundary edges = edges owned by exactly one element, shape (n_edges, 2)."""
    e = mesh.elements
    edges = np.vstack([e[:, [0, 1]], e[:, [1, 2]], e[:, [2, 3]], e[:, [3, 0]]])
    key = np.sort(edges, axis=1)
    _, idx, counts = np.unique(key, axis=0, return_index=True, return_counts=True)
    return edges[idx[counts == 1]]


def _set_limits(ax, *coord_sets: np.ndarray, pad: float = 0.06) -> None:
    pts = np.vstack(coord_sets)
    xmin, ymin = pts.min(axis=0)
    xmax, ymax = pts.max(axis=0)
    span = max(xmax - xmin, ymax - ymin)
    ax.set_xlim(xmin - pad * span, xmax + pad * span)
    ax.set_ylim(ymin - pad * span, ymax + pad * span)
    ax.set_aspect("equal")


def plot_results(
    model: FEAModel,
    results: Optional[FEAResults] = None,
    scale: Optional[float] = None,
    filename: Optional[str] = "fea_results.png",
    title: str = "Q4 Finite Element Analysis",
    stress_unit: str = "Pa",
    length_unit: str = "m",
    n_levels: int = 32,
    cmap: str = "turbo",
    dpi: int = 220,
    figsize: Optional[Tuple[float, float]] = None,
    show: bool = False,
) -> Tuple[plt.Figure, float]:
    """
    Render the undeformed mesh and the Von Mises field on the deformed shape.

    Returns the figure and the deformation scale factor actually used.
    """
    mesh = model.mesh
    results = results if results is not None else model.results
    if results is None:
        raise ValueError("model has not been solved; call model.solve() first")
    scale = auto_scale_factor(mesh, results) if scale is None else float(scale)

    undeformed = mesh.nodes
    deformed = mesh.nodes + scale * results.displacements

    xmin, xmax, ymin, ymax = mesh.bounding_box()
    aspect = (xmax - xmin) / max(ymax - ymin, 1e-300)
    wide = aspect > 2.5  # long, slender parts stack the panels vertically
    if figsize is None:
        figsize = (14, 6.5) if wide else (15, 7)
    fig, axes = plt.subplots(
        2 if wide else 1, 1 if wide else 2,
        figsize=figsize,
        constrained_layout=True,
    )
    ax0, ax1 = axes

    # ------------------------------------------------------- undeformed
    _wireframe(ax0, mesh, undeformed, facecolors="#f4f6f8", edgecolors="#33415c", linewidths=0.45)
    constrained = model.constrained_dofs
    if constrained.size:
        fixed_x = np.unique(constrained[constrained % 2 == 0] // 2)
        fixed_y = np.unique(constrained[constrained % 2 == 1] // 2)
        both = np.intersect1d(fixed_x, fixed_y)
        only_x = np.setdiff1d(fixed_x, both)
        only_y = np.setdiff1d(fixed_y, both)
        if both.size:
            ax0.plot(*undeformed[both].T, "s", ms=4, color="#c1121f", label="fixed u_x, u_y")
        if only_x.size:
            ax0.plot(*undeformed[only_x].T, ">", ms=4, color="#f77f00", label="roller (u_x = 0)")
        if only_y.size:
            ax0.plot(*undeformed[only_y].T, "^", ms=4, color="#2a9d8f", label="roller (u_y = 0)")
    F = results.external_forces.reshape(-1, 2)
    loaded = np.flatnonzero(np.linalg.norm(F, axis=1) > 0.0)
    if loaded.size:
        span = max(xmax - xmin, ymax - ymin)
        fmag = np.linalg.norm(F[loaded], axis=1).max()
        arrows = F[loaded] / fmag * 0.08 * span
        ax0.quiver(
            undeformed[loaded, 0] - arrows[:, 0], undeformed[loaded, 1] - arrows[:, 1],
            arrows[:, 0], arrows[:, 1], angles="xy", scale_units="xy", scale=1.0,
            color="#6a040f", width=0.003, label="applied load",
        )
    ax0.set_title(
        f"Undeformed mesh: {mesh.n_nodes} nodes, {mesh.n_elements} Q4 elements, "
        f"{mesh.n_dofs} DOFs", fontsize=11,
    )
    ax0.set_xlabel(f"x [{length_unit}]")
    ax0.set_ylabel(f"y [{length_unit}]")
    ax0.legend(loc="best", fontsize=8, framealpha=0.9)
    _set_limits(ax0, undeformed)

    # ------------------------------------------------- deformed contour
    tri = quad_triangulation(mesh, deformed)
    vm = results.von_mises
    vmax = float(vm.max())
    vmin = float(vm.min())
    if vmax - vmin <= 1e-12 * max(abs(vmax), 1.0):
        vmax = vmin + max(abs(vmin), 1.0) * 1e-6
    levels = np.linspace(vmin, vmax, n_levels + 1)
    cf = ax1.tricontourf(tri, vm, levels=levels, cmap=cmap, extend="neither")
    _wireframe(ax1, mesh, deformed, facecolors="none", edgecolors=(0, 0, 0, 0.18), linewidths=0.3)
    for a, b in _outline_edges(mesh):
        ax1.plot(*undeformed[[a, b]].T, "--", color="#6c757d", lw=0.8)
    cbar = fig.colorbar(cf, ax=ax1, shrink=0.92, pad=0.01)
    cbar.set_label(f"Von Mises stress [{stress_unit}]")
    imax = int(np.argmax(vm))
    ax1.plot(*deformed[imax], "o", mfc="none", mec="white", ms=9, mew=1.5)
    ax1.annotate(
        f"max {vmax:.4g} {stress_unit}", deformed[imax], xytext=(10, 10),
        textcoords="offset points", fontsize=8, color="black",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.5", alpha=0.9),
    )
    ax1.set_title(
        f"Von Mises on deformed shape (scale x{scale:.4g}, max |u| = "
        f"{results.max_displacement:.4e} {length_unit}); dashed = undeformed",
        fontsize=11,
    )
    ax1.set_xlabel(f"x [{length_unit}]")
    ax1.set_ylabel(f"y [{length_unit}]")
    _set_limits(ax1, undeformed, deformed)

    fig.suptitle(
        f"{title}  |  E = {model.material.E:.4g}, nu = {model.material.nu:.3g}, "
        f"t = {model.material.thickness:.4g}, plane {model.material.plane}",
        fontsize=13, fontweight="bold",
    )
    if filename:
        directory = os.path.dirname(os.path.abspath(filename))
        os.makedirs(directory, exist_ok=True)
        fig.savefig(filename, dpi=dpi)
    if show:
        plt.show()
    return fig, scale


def plot_stress_components(
    model: FEAModel,
    results: Optional[FEAResults] = None,
    filename: Optional[str] = "fea_stress_components.png",
    stress_unit: str = "Pa",
    cmap: str = "RdBu_r",
    dpi: int = 200,
) -> plt.Figure:
    """Three diverging-colour panels for s_xx, s_yy, t_xy on the undeformed mesh."""
    mesh = model.mesh
    results = results if results is not None else model.results
    if results is None:
        raise ValueError("model has not been solved; call model.solve() first")
    tri = quad_triangulation(mesh, mesh.nodes)
    xmin, xmax, ymin, ymax = mesh.bounding_box()
    wide = (xmax - xmin) / max(ymax - ymin, 1e-300) > 2.5
    fig, axes = plt.subplots(
        3 if wide else 1, 1 if wide else 3,
        figsize=(14, 9) if wide else (18, 5.5), constrained_layout=True,
    )
    names = [r"$\sigma_{xx}$", r"$\sigma_{yy}$", r"$\tau_{xy}$"]
    for k, (ax, name) in enumerate(zip(axes, names)):
        field = results.nodal_stresses[:, k]
        lim = float(np.abs(field).max()) or 1.0
        cf = ax.tricontourf(tri, field, levels=np.linspace(-lim, lim, 33), cmap=cmap)
        fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.01).set_label(f"{name} [{stress_unit}]")
        ax.set_title(f"{name} (nodal, extrapolated + averaged)")
        _set_limits(ax, mesh.nodes)
    if filename:
        os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
        fig.savefig(filename, dpi=dpi)
    return fig
