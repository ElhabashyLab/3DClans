"""
Benchmark module for 3DClans structural similarity tools.

This module provides comprehensive benchmarking of the full 3DClans pipeline,
measuring timing for structure download, score computation, and CLANS file generation.
"""

from clans3d.benchmark.benchmark import Benchmark
from clans3d.benchmark.benchmark_result import BenchmarkResult

__all__ = ["Benchmark", "BenchmarkResult"]
