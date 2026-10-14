"""
Regenerate the pictures used in the documentation (docs/images/*.png).

    python3 docs/make_images.py
"""

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from engine import FEAModel, Material, Mesh  # noqa: E402
from visualizer import plot_results, plot_stress_components  # noqa: E402

IMAGES = os.path.join(ROOT, "docs", "images")
DPI = 110


def save(fig_and_scale, name):
    fig = fig_and_scale[0] if isinstance(fig_and_scale, tuple) else fig_and_scale
    path = os.path.join(IMAGES, name)
    fig.savefig(path, dpi=DPI)
    print("wrote", os.path.relpath(path, ROOT))


def tutorial_beam():
    mesh = Mesh.structured_rectangle(length=1.0, height=0.1, nx=100, ny=10)
    steel = Material(E=210e9, nu=0.3, thickness=0.02, yield_strength=250e6, name="steel")
    model = FEAModel(mesh, steel)
    model.fix(mesh.nodes_at(x=0.0), "xy")
    model.add_nodal_force(mesh.closest_node(1.0, 0.05), fy=-2000)
    r = model.solve()
    save(plot_results(model, r, filename=None, title="Tutorial beam"), "tutorial_beam.png")
    save(plot_stress_components(model, r, filename=None), "tutorial_components.png")


def simply_supported():
    L, H, T = 3.0, 0.15, 0.05
    mesh = Mesh.structured_rectangle(L, H, nx=240, ny=12, origin=(0.0, -H / 2))
    model = FEAModel(mesh, Material(200e9, 0.3, thickness=T))
    model.fix(mesh.closest_node(0.0, 0.0), "xy")
    model.fix(mesh.closest_node(L, 0.0), "y")
    top = mesh.nodes_at(y=H / 2)
    model.add_edge_traction(top[np.argsort(mesh.nodes[top, 0])], (0.0, -2.0e5))
    r = model.solve()
    save(plot_results(model, r, filename=None, title="Simply supported beam, uniform pressure"),
         "simply_supported.png")


def plate_with_hole():
    a, W = 0.01, 0.10
    mesh = Mesh.plate_with_hole_quarter(a, W, 40, 48, grading=8.0)
    model = FEAModel(mesh, Material(70e9, 0.33, thickness=0.005, yield_strength=275e6))
    model.fix(mesh.nodes_at(x=0.0), "x")
    model.fix(mesh.nodes_at(y=0.0), "y")
    edge = mesh.nodes_at(x=W)
    model.add_edge_traction(edge[np.argsort(mesh.nodes[edge, 1])], (60e6, 0.0))
    r = model.solve()
    save(plot_results(model, r, filename=None, title="Plate with hole", scale=150),
         "plate_with_hole.png")


def custom_mesh():
    L, h0, h1, nx, ny = 1.0, 0.4, 0.1, 40, 8
    nodes = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            x = L * i / nx
            h = h0 + (h1 - h0) * x / L
            nodes.append([x, -h / 2 + h * j / ny])
    elements = []
    for j in range(ny):
        for i in range(nx):
            n0 = j * (nx + 1) + i
            elements.append([n0, n0 + 1, n0 + nx + 2, n0 + nx + 1])
    mesh = Mesh(nodes, elements)
    model = FEAModel(mesh, Material(30e9, 0.2, thickness=1.0, plane="strain"))
    model.fix(mesh.nodes_at(x=0.0), "xy")
    model.add_nodal_force(mesh.closest_node(L, 0.0), fx=20e3, fy=-50e3)
    r = model.solve()
    save(plot_results(model, r, filename=None, title="Tapered bracket (plane strain)"),
         "custom_mesh.png")


if __name__ == "__main__":
    os.makedirs(IMAGES, exist_ok=True)
    tutorial_beam()
    simply_supported()
    plate_with_hole()
    custom_mesh()
