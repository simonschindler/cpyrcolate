This is a huge milestone! Moving from a collection of scripts to a properly structured, compiled Python package is a massive step up in your development workflow.

Here is a comprehensive and professional `README.md` that explains exactly what `cpyrcolate` does, how to install it (with and without the dev tools), and how to use the API we just built.

---

# cpyrcolate

**cpyrcolate** is a high-performance Python package for computing bond percolation statistics on general graphs. It relies on a highly optimized C backend (via Cython) implementing the fast Monte Carlo algorithm by Newman and Ziff.

Unlike traditional grid-based percolation solvers, `cpyrcolate` operates directly on generic NumPy edge lists, making it ideal for studying complex spatial networks, geometric graphs (like Delaunay or k-NN graphs), and arbitrary network topologies.

## Features

* **Blazing Fast Backend:** The core Union-Find logic is written in pure C, completely bypassing Python's `for`-loop overhead.
* **NumPy Native:** Avoids the memory and performance bottlenecks of building full `NetworkX` graph objects for every run.
* **Canonical & Microcanonical Ensembles:** Computes both raw edge-addition paths (microcanonical) and expected behaviors at specific occupation probabilities (canonical).
* **Spatial Spanning Clusters:** Built-in geometric boundary detection to calculate the probability of a spanning cluster connecting two sides of a spatial graph.

## Installation

### Standard Installation

If you just want to use the package in another project, you can install it directly. This will only install the core dependencies (`numpy` and `scipy`).

```bash
pip install .
# or using uv:
uv pip install .

```

### Development Installation

If you are modifying the C/Cython code or want to run the tests and visualization scripts, you should install it in "editable" mode with the `dev` optional dependencies (`networkx`, `matplotlib`, `pytest`).

```bash
uv sync --extra dev
# or using standard pip:
pip install -e ".[dev]"

```

*Note: Because this package contains C extensions, you need a functioning C compiler (like `gcc` or `clang`) installed on your system.*

## Quickstart & Usage

The package exposes two primary functions: `compute_percolation_single` for single-run trajectories and `compute_percolation_statistics` for averaged canonical statistics.

### 1. Single Run (Microcanonical)

Use this to track the size of the largest connected component as edges are added one by one to an empty graph.

```python
import numpy as np
from cpyrcolate import compute_percolation_single

# 1. Define your graph as a NumPy edge list (M, 2)
edges = np.array([
    [0, 1], 
    [1, 2], 
    [2, 3], 
    [3, 0]
], dtype=np.int32)

# 2. Run the percolation simulation
stats = compute_percolation_single(edges=edges)

# The output contains the size of the largest cluster after each edge is added
print("Max cluster size history:", stats['max_cluster_size'])

```

### 2. Averaged Statistics (Canonical)

Use this to calculate the expected percolation strength (normalized size of the giant component) at specific bond occupation probabilities $p$.

```python
from cpyrcolate import compute_percolation_statistics

# 1. Define the probabilities you want to evaluate (0 to 1)
probabilities = np.linspace(0, 1, 50)

# 2. Compute statistics averaged over 40 independent runs
# This convolves the microcanonical runs with the Binomial PMF
perc_stats = compute_percolation_statistics(
    edges=edges, 
    ps=probabilities, 
    runs=40
)

print("Percolation Strength:", perc_stats['max_cluster_size'])

```

### 3. Spanning Cluster Detection

If your nodes have spatial coordinates, `cpyrcolate` can automatically detect if a cluster spans from one side of the space to the other.

```python
# Node coordinates (N, 2)
coords = np.random.rand(100, 2)
edges = ... # Your geometric edge list

stats = compute_percolation_statistics(
    edges=edges,
    ps=np.linspace(0, 1, 50),
    spanning_cluster=True,
    coords=coords,
    axis=0,        # 0 for Left-Right spanning, 1 for Top-Bottom
    margin=0.05    # Use the outer 5% of coordinates as boundaries
)

print("Probability of Spanning:", stats['spanning_cluster'])

```

## Repository Structure

* `src/cpyrcolate/`: The core Python package.
* `percolate_core.c` / `.h`: The pure C Newman-Ziff Union-Find implementation.
* `percolate_cy.pyx`: The Cython wrapper handling memory views.
* `core.py`: The Python API handling data prep and canonical math.


* `tests/`: Testing scripts for both the C-backend and the final Python API.
* Use `make` inside the tests directory to quickly compile and verify the pure C code.



## Acknowledgements

The core Union-Find algorithm implemented in C is based on the fast Monte Carlo algorithm proposed by M. E. J. Newman and R. M. Ziff (2001).


