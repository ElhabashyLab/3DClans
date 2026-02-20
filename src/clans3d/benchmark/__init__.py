"""
Benchmark module for Clans-3D structural similarity tools.

This module provides comprehensive benchmarking of the full Clans-3D pipeline,
measuring timing for PDB download, score computation, and CLANS file generation.
"""

from clans3d.benchmark.benchmark import Benchmark
from clans3d.benchmark.benchmark_result import BenchmarkResult

__all__ = ["Benchmark", "BenchmarkResult"]
