#!/usr/bin/env python3
"""
Benchmark driver for the Q4 FEA engine.

Runs:
  1. Euler-Bernoulli cantilever benchmark (end point load), compared against
     delta = P L^3 / (3 E I) and the Timoshenko correction.
  2. Quarter plate with a circular hole under uniaxial tension on a
     curvilinear mapped mesh, compared against the Kirsch / Heywood stress
     concentration factor.

Figures are written to ./output/. Use --run-tests to execute the pytest
validation suite afterwards.
"""

from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys

import numpy as np

from engine import FEAModel, Material, Mesh, cantilever_tip_deflection

logger = logging.getLogger("app")
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")


def configure_logging(verbose: bool = True) -> None:
    logging.basicConfig(
        level=logging.INFO if verbose else logging.WARNING,
        format="%(asctime)s | %(name)-10s | %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stdout,
        force=True,
    )
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)


def banner(text: str) -> None:
    logger.info("#" * 78)
    logger.info("# %s", text)
    logger.info("#" * 78)


# --------------------------------------------------------------------------- #
def run_cantilever(nx: int, ny: int, plane: str, scale, plot: bool) -> float:
    banner("BENCHMARK 1: end-loaded cantilever vs Euler-Bernoulli")
    L, h, t = 2.0, 0.1, 0.05            # m
    E, nu = 210e9, 0.3                  # Pa, -
    P = 5_000.0                         # N, downward
    I = t * h**3 / 12.0

    mesh = Mesh.structured_rectangle(L, h, nx, ny, origin=(0.0, -0.5 * h))
    material = Material(E, nu, thickness=t, plane=plane, yield_strength=250e6, name="S250 steel")
    model = FEAModel(mesh, material)

    root = mesh.nodes_at(x=0.0)
    model.fix(root, "xy")
    tip_node = mesh.closest_node(L, 0.0)
    model.add_nodal_force(tip_node, fy=-P)
    logger.info(
        "Geometry L = %.4g m, h = %.4g m, t = %.4g m (L/h = %.1f); I = t h^3 / 12 = %.6e m^4",
        L, h, t, L / h, I,
    )
    logger.info("Root x = 0: %d nodes fully clamped; tip node %d at (%.3f, %.3f) carries F_y = %.1f N",
                root.size, tip_node, *mesh.nodes[tip_node], -P)

    results = model.solve()

    E_eff = E / (1.0 - nu**2) if material.plane == Material.PLANE_STRAIN else E
    delta_eb = cantilever_tip_deflection(P, L, E_eff, I)
    kappa = 5.0 / 6.0
    G = E / (2.0 * (1.0 + nu))
    delta_timo = delta_eb + P * L / (kappa * G * h * t)
    delta_fe = -results.displacements[tip_node, 1]
    err = (delta_fe - delta_eb) / delta_eb
    sigma_root = P * L * (0.5 * h) / I

    logger.info("-" * 78)
    logger.info("Analytical (Euler-Bernoulli)  delta = P L^3 / (3 E I) = %.6e m", delta_eb)
    logger.info("Analytical (Timoshenko)       delta + P L / (k G A)  = %.6e m", delta_timo)
    logger.info("Finite element tip deflection                         = %.6e m", delta_fe)
    logger.info("Relative error vs Euler-Bernoulli                     = %+.4f %%", 100 * err)
    logger.info("Relative error vs Timoshenko                          = %+.4f %%",
                100 * (delta_fe - delta_timo) / delta_timo)
    logger.info("Root bending stress M c / I = %.4e Pa; FE max |s_xx| at Gauss pts = %.4e Pa",
                sigma_root, np.abs(results.gauss_stresses[..., 0]).max())
    logger.info("Total vertical reaction = %.4f N (applied %.4f N)",
                results.reactions[1::2].sum(), P)
    logger.info("-" * 78)

    if plot:
        from visualizer import plot_results, plot_stress_components

        path = os.path.join(OUTPUT_DIR, "cantilever_von_mises.png")
        _, used = plot_results(model, results, scale=scale, filename=path,
                               title="Cantilever benchmark (tip load)")
        plot_stress_components(model, results,
                               filename=os.path.join(OUTPUT_DIR, "cantilever_components.png"))
        logger.info("Figures written to %s (deformation scale x%.4g)", OUTPUT_DIR, used)
    return err


# --------------------------------------------------------------------------- #
def run_plate_with_hole(plane: str, scale, plot: bool) -> float:
    banner("BENCHMARK 2: plate with circular hole (Kirsch stress concentration)")
    a, W, t = 0.01, 0.10, 0.005         # hole radius, half-width, thickness [m]
    E, nu = 70e9, 0.33                  # aluminium
    s0 = 60e6                           # remote tension [Pa]

    mesh = Mesh.plate_with_hole_quarter(a, W, n_radial=40, n_circumferential=48, grading=8.0)
    material = Material(E, nu, thickness=t, plane=plane, yield_strength=275e6, name="Al 6061-T6")
    model = FEAModel(mesh, material)

    model.fix(mesh.nodes_at(x=0.0), "x")   # symmetry about the y-axis
    model.fix(mesh.nodes_at(y=0.0), "y")   # symmetry about the x-axis
    right = mesh.nodes_at(x=W)
    right = right[np.argsort(mesh.nodes[right, 1])]
    model.add_edge_traction(right, (s0, 0.0))
    logger.info("Quarter model: a = %.3g m, W = %.3g m, remote s_xx = %.3e Pa on x = W edge",
                a, W, s0)

    results = model.solve()

    top = mesh.closest_node(0.0, a)
    kt_fe = results.nodal_stresses[top, 0] / s0
    ratio = 2 * a / (2 * W)
    kt_heywood = (2.0 + (1.0 - ratio) ** 3) / (1.0 - ratio)  # gross-section
    logger.info("-" * 78)
    logger.info("Kirsch (infinite plate) Kt = 3.000; Heywood finite-width gross Kt = %.4f", kt_heywood)
    logger.info("FE nodal s_xx at hole crown (0, a) / s0 = %.4f  (diff vs Heywood %+.2f %%)",
                kt_fe, 100 * (kt_fe - kt_heywood) / kt_heywood)
    logger.info("(Heywood assumes a long strip; this plate is square, so its free top edge "
                "raises Kt by a few percent. Informational, not a pass/fail check.)")
    logger.info("Max Von Mises %.4e Pa; minimum safety factor vs %.0f MPa yield = %.3f",
                results.max_von_mises, material.yield_strength / 1e6, results.safety_factor.min())
    logger.info("-" * 78)

    if plot:
        from visualizer import plot_results

        path = os.path.join(OUTPUT_DIR, "plate_hole_von_mises.png")
        _, used = plot_results(model, results, scale=scale, filename=path,
                               title="Quarter plate with hole (uniaxial tension)")
        logger.info("Figure written to %s (deformation scale x%.4g)", path, used)
    return kt_fe


# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Q4 FEA engine benchmarks")
    parser.add_argument("--nx", type=int, default=160, help="cantilever elements along the length")
    parser.add_argument("--ny", type=int, default=8, help="cantilever elements through the depth")
    parser.add_argument("--plane", choices=["stress", "strain"], default="stress")
    parser.add_argument("--scale", type=float, default=None,
                        help="deformation scale factor for plots (default: automatic)")
    parser.add_argument("--no-plot", action="store_true", help="skip figure generation")
    parser.add_argument("--run-tests", action="store_true", help="run test_validation.py afterwards")
    parser.add_argument("--quiet", action="store_true", help="only warnings and errors")
    args = parser.parse_args(argv)

    configure_logging(not args.quiet)
    err = run_cantilever(args.nx, args.ny, args.plane, args.scale, not args.no_plot)
    run_plate_with_hole(args.plane, args.scale, not args.no_plot)

    status = 0
    if abs(err) > 0.02:
        logger.error("Cantilever error %.3f %% exceeds 2 %% tolerance", 100 * err)
        status = 1
    else:
        logger.info("Cantilever benchmark PASSED (|error| = %.3f %% <= 2 %%)", 100 * abs(err))

    if args.run_tests:
        banner("Running validation suite (pytest test_validation.py)")
        here = os.path.dirname(os.path.abspath(__file__))
        rc = subprocess.call(
            [sys.executable, "-m", "pytest", os.path.join(here, "test_validation.py"), "-v"],
            cwd=here,
        )
        status = status or rc
    return status


if __name__ == "__main__":
    sys.exit(main())
