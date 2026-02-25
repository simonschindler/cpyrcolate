import numpy as np
import scipy.stats
from . import percolate_cy


def _prepare_inputs(edges, coords=None, spanning_cluster=False, axis=0, margin=0.05):
    """Internal helper to format edges and boundaries for the C backend."""
    num_main_edges = edges.shape[0]
    num_nodes_main = np.max(edges) + 1 if num_main_edges > 0 else 0
    if coords is not None:
        num_nodes_main = max(num_nodes_main, len(coords))

    if not spanning_cluster:
        edges_u = edges[:, 0].astype(np.int32)
        edges_v = edges[:, 1].astype(np.int32)
        return (
            num_nodes_main,
            num_nodes_main,
            edges_u,
            edges_v,
            num_main_edges,
            0,
            -1,
            -1,
        )

    if coords is None:
        raise ValueError("coords are required for spanning cluster detection.")

    # Find boundary nodes purely in NumPy
    min_val = np.min(coords[:, axis])
    max_val = np.max(coords[:, axis])
    rng = max_val - min_val

    side_0_nodes = np.where(coords[:, axis] <= min_val + margin * rng)[0]
    side_1_nodes = np.where(coords[:, axis] >= max_val - margin * rng)[0]

    aux_0 = num_nodes_main
    aux_1 = num_nodes_main + 1
    total_nodes = num_nodes_main + 2

    # Create static boundary edges
    aux_0_u = np.full(len(side_0_nodes), aux_0, dtype=np.int32)
    aux_0_v = side_0_nodes.astype(np.int32)
    aux_1_u = np.full(len(side_1_nodes), aux_1, dtype=np.int32)
    aux_1_v = side_1_nodes.astype(np.int32)

    # Prepend boundary edges so they are added first
    edges_u = np.concatenate([aux_0_u, aux_1_u, edges[:, 0]]).astype(np.int32)
    edges_v = np.concatenate([aux_0_v, aux_1_v, edges[:, 1]]).astype(np.int32)

    num_bound_edges = len(side_0_nodes) + len(side_1_nodes)

    return (
        num_nodes_main,
        total_nodes,
        edges_u,
        edges_v,
        num_main_edges,
        num_bound_edges,
        aux_0,
        aux_1,
    )


def compute_percolation_single(
    edges, spanning_cluster=False, coords=None, axis=0, margin=0.05
):
    """Run a single microcanonical percolation realization."""
    (
        num_nodes_main,
        total_nodes,
        edges_u,
        edges_v,
        num_main_edges,
        num_bound_edges,
        aux_0,
        aux_1,
    ) = _prepare_inputs(edges, coords, spanning_cluster, axis, margin)

    # Force boundary edges to be added first, then randomly permute the main graph edges
    order_bound = np.arange(num_bound_edges, dtype=np.int32)
    order_main = (
        np.random.permutation(num_main_edges).astype(np.int32) + num_bound_edges
    )
    order = np.concatenate([order_bound, order_main])

    out_max_size, out_span = percolate_cy.run_percolation(
        num_nodes_main, total_nodes, edges_u, edges_v, order, aux_0, aux_1
    )

    # Slice off the initial boundary edge additions
    return {
        "max_cluster_size": out_max_size[num_bound_edges:],
        "spanning_cluster": out_span[num_bound_edges:],
        "N": num_nodes_main,
        "M": num_main_edges,
    }


def compute_percolation_statistics(
    edges, ps, runs=40, spanning_cluster=False, coords=None, axis=0, margin=0.05
):
    """Compute canonical ensemble percolation statistics averaged over multiple runs."""
    (
        num_nodes_main,
        total_nodes,
        edges_u,
        edges_v,
        num_main_edges,
        num_bound_edges,
        aux_0,
        aux_1,
    ) = _prepare_inputs(edges, coords, spanning_cluster, axis, margin)

    all_max_sizes = np.zeros((runs, num_main_edges + 1))
    all_spans = np.zeros((runs, num_main_edges + 1))
    order_bound = np.arange(num_bound_edges, dtype=np.int32)

    # 1. Run Microcanonical Ensembles
    for r in range(runs):
        order_main = (
            np.random.permutation(num_main_edges).astype(np.int32) + num_bound_edges
        )
        order = np.concatenate([order_bound, order_main])

        out_max_size, out_span = percolate_cy.run_percolation(
            num_nodes_main, total_nodes, edges_u, edges_v, order, aux_0, aux_1
        )
        all_max_sizes[r, :] = out_max_size[num_bound_edges:]
        all_spans[r, :] = out_span[num_bound_edges:]

    # Calculate microcanonical averages
    micro_max_size = (
        np.mean(all_max_sizes, axis=0) / num_nodes_main
        if num_nodes_main > 0
        else np.zeros(num_main_edges + 1)
    )

    # 1-sigma standard error for max cluster size
    ddof = 1 if runs > 1 else 0
    stderr_max = (
        (np.std(all_max_sizes, axis=0, ddof=ddof) / np.sqrt(runs)) / num_nodes_main
        if num_nodes_main > 0
        else np.zeros(num_main_edges + 1)
    )

    # Spanning probability using Bayesian posterior mean
    k_span = np.sum(all_spans, axis=0)
    micro_span = (k_span + 1.0) / (runs + 2.0)
    span_ci_low = scipy.stats.beta.ppf(0.15865, k_span + 1, runs - k_span + 1)
    span_ci_high = scipy.stats.beta.ppf(1 - 0.15865, k_span + 1, runs - k_span + 1)

    # 2. Convert to Canonical Ensemble
    ret = {
        "ps": ps,
        "N": num_nodes_main,
        "M": num_main_edges,
        "max_cluster_size": np.zeros(len(ps)),
        "max_cluster_size_ci": np.zeros((len(ps), 2)),
        "spanning_cluster": np.zeros(len(ps)),
        "spanning_cluster_ci": np.zeros((len(ps), 2)),
    }

    n_array = np.arange(num_main_edges + 1)

    # Convolve with the Binomial PMF
    for i, p in enumerate(ps):
        binomials = scipy.stats.binom.pmf(n_array, num_main_edges, p)

        ret["max_cluster_size"][i] = np.sum(binomials * micro_max_size)
        ret["max_cluster_size_ci"][i, 0] = np.sum(
            binomials * (micro_max_size - stderr_max)
        )
        ret["max_cluster_size_ci"][i, 1] = np.sum(
            binomials * (micro_max_size + stderr_max)
        )

        if spanning_cluster:
            ret["spanning_cluster"][i] = np.sum(binomials * micro_span)
            ret["spanning_cluster_ci"][i, 0] = np.sum(binomials * span_ci_low)
            ret["spanning_cluster_ci"][i, 1] = np.sum(binomials * span_ci_high)

    return ret
