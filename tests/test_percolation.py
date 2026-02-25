import pytest
import numpy as np
import networkx as nx
from cpyrcolate import compute_percolation_single


def get_theoretical_S(c_values):
    """
    Numerically solves the Erdős-Rényi giant component equation: S = 1 - exp(-c*S)
    """
    S_theory = np.zeros_like(c_values)
    for i, c in enumerate(c_values):
        if c <= 1.0:
            S_theory[i] = 0.0
        else:
            # Fixed-point iteration to solve for S
            S = 0.5  # initial guess
            for _ in range(100):
                S = 1.0 - np.exp(-c * S)
            S_theory[i] = S
    return S_theory


def test_erdos_renyi_percolation():
    """
    Tests the microcanonical percolation output against the
    analytical solution for an Erdős-Rényi random graph.
    """
    # Use 10,000 nodes for a good balance of speed and statistical smoothness
    N = 10_000

    # We want to test up to an average degree of c = 4
    # Since c = 2*M / N, the total edges M = c * N / 2
    M = int(4 * N / 2)

    # Generate ER graph with a fixed seed for reproducible tests
    G = nx.gnm_random_graph(N, M, seed=42)
    edges = np.array(list(G.edges()), dtype=np.int32)

    # Run the single-run microcanonical percolation from our C-backend
    stats = compute_percolation_single(edges=edges, spanning_cluster=False)

    # Extract measured percolation strength (normalized max cluster size)
    measured_max_cluster = stats["max_cluster_size"]
    S_measured = measured_max_cluster / N

    # Calculate the average degree (c) at every step of the edge addition
    m_steps = np.arange(len(measured_max_cluster))
    c_values = 2.0 * m_steps / N

    # Get the theoretical curve
    S_theory = get_theoretical_S(c_values)

    # --- Assertions ---

    # 1. Sub-critical regime check (c < 0.8)
    # The giant component should not have formed yet.
    sub_critical_mask = c_values < 0.8
    max_sub_critical_S = np.max(S_measured[sub_critical_mask])
    assert max_sub_critical_S < 0.05, (
        f"Giant component appeared too early! Max S={max_sub_critical_S:.3f}"
    )

    # 2. Super-critical regime check (c > 1.2)
    # The simulation should closely track the theoretical curve.
    # We skip the critical window (0.8 < c < 1.2) due to finite-size effects.
    super_critical_mask = c_values > 1.2

    # Calculate Mean Absolute Error (MAE) between simulation and theory
    mae = np.mean(
        np.abs(S_measured[super_critical_mask] - S_theory[super_critical_mask])
    )

    # We expect the error to be less than 2%
    assert mae < 0.02, (
        f"Measured percolation strength deviates from theory (MAE={mae:.4f})"
    )
