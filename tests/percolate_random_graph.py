import numpy as np
import matplotlib.pyplot as plt
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


if __name__ == "__main__":
    print("Generating random graph...")
    N = 50_000

    # We want to go up to an average degree (c) of 4.
    # Since c = 2*M / N, the total edges M = c * N / 2
    M = int(4 * N / 2)

    # Generate an Erdős-Rényi graph using NetworkX
    G = nx.gnm_random_graph(N, M)
    edges = np.array(list(G.edges()), dtype=np.int32)

    print(f"Running C-backend percolation on {N} nodes and {M} edges...")
    # Run the microcanonical percolation
    stats = compute_percolation_single(edges=edges, spanning_cluster=False)

    # Extract the measured percolation strength
    # Normalize the max cluster size by N to get the fraction S
    measured_max_cluster = stats["max_cluster_size"]
    S_measured = measured_max_cluster / N

    # Calculate the average degree (c) at every step of the edge addition
    m_steps = np.arange(M + 1)
    c_values = 2.0 * m_steps / N

    # Calculate the theoretical expectation
    S_theory = get_theoretical_S(c_values)

    # Plotting
    print("Plotting results...")
    plt.figure(figsize=(10, 6))

    # Plot the theoretical curve
    plt.plot(c_values, S_theory, "k--", linewidth=2, label="Theoretical ER Model")

    # Plot our measured simulation
    plt.plot(
        c_values,
        S_measured,
        "r-",
        alpha=0.8,
        linewidth=2,
        label="Measured (cpyrcolate C-Backend)",
    )

    # Formatting
    plt.axvline(x=1.0, color="gray", linestyle=":", label="Critical Threshold (c=1)")
    plt.title("Percolation on an Erdős–Rényi Random Graph\nTheory vs. Simulation")
    plt.xlabel("Average Degree ($c$)")
    plt.ylabel("Fraction of nodes in Largest Component ($S$)")
    plt.legend(loc="lower right")
    plt.grid(True, linestyle="--", alpha=0.6)

    plt.tight_layout()
    # plt.show()
    # Save the plot for GitHub Actions instead of displaying it interactively
    plt.savefig("er_percolation_plot.png", dpi=300, bbox_inches="tight")
    print("Plot successfully saved to er_percolation_plot.png")
