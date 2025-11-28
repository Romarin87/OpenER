"""
Workflow package for elementary reaction TS processing.

The modules cover TS optimization, saddle point checks, SOAP-based deduplication,
IRC tracing, minima optimization, and SMILES validation.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .pipeline import TransitionStatePipeline

__all__ = ["TransitionStatePipeline"]


def __getattr__(name: str):
    if name == "TransitionStatePipeline":
        from .pipeline import TransitionStatePipeline as _TSP

        return _TSP
    raise AttributeError(name)
