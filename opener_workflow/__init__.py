"""
Workflow package for elementary reaction TS processing.

The modules cover TS optimization, saddle point checks, SOAP-based deduplication,
IRC tracing, minima optimization, and SMILES validation.
"""

from .pipeline import TransitionStatePipeline

__all__ = ["TransitionStatePipeline"]
