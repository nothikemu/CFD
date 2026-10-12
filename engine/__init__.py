"""Standalone 2D Q4 finite element structural solver."""

from .fea_solver import (
    DistortedElementError,
    ElementIntegrityError,
    FEAError,
    FEAModel,
    FEAResults,
    Material,
    Mesh,
    MeshTopologyError,
    Q4Element,
    SingularSystemError,
    SolverConvergenceError,
    cantilever_tip_deflection,
    von_mises,
)

__all__ = [
    "DistortedElementError",
    "ElementIntegrityError",
    "FEAError",
    "FEAModel",
    "FEAResults",
    "Material",
    "Mesh",
    "MeshTopologyError",
    "Q4Element",
    "SingularSystemError",
    "SolverConvergenceError",
    "cantilever_tip_deflection",
    "von_mises",
]

__version__ = "1.0.0"
